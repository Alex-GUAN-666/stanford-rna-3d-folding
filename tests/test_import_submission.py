"""Data integrity checks for the new restored-submission import bridge."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from examples.import_submission import import_submission
from rna_folding.io import SUBMISSION_COLUMNS, load_candidate_manifest, write_submission
from rna_folding.selection import Candidate
from rna_folding.sequence import RNASequence


class ImportSubmissionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.sequence_path = self.root / "test_sequences.csv"
        self.sequence_path.write_text(
            "target_id,sequence,description\ntarget_with_underscores,ACGU,first\nother,GAC,second\n",
            encoding="utf-8",
        )
        self.sequences = [RNASequence("target_with_underscores", "ACGU"), RNASequence("other", "GAC")]
        self.submission = self.root / "submission.csv"
        self.destination = self.root / "imported"
        self.predictions = {}
        for index, sequence in enumerate(self.sequences):
            xyz = np.arange(len(sequence.sequence) * 3, dtype=float).reshape(-1, 3) / 8 - 7.5 + index * 100
            self.predictions[sequence.target_id] = [
                Candidate(f"{index}_{slot}", sequence.target_id, "drfold2", xyz + slot * 0.125)
                for slot in range(5)
            ]
        write_submission(self.submission, self.sequences, self.predictions)

    def mutate_rows(self, mutate):
        with self.submission.open(newline="") as stream:
            rows = list(csv.reader(stream))
        mutate(rows)
        with self.submission.open("w", newline="") as stream:
            csv.writer(stream).writerows(rows)

    def run_import(self):
        return import_submission(self.sequence_path, self.submission, "DRfold2", self.destination, "backup")

    def assert_preserved_predictions(self):
        candidates = load_candidate_manifest(self.destination / "candidates.json", self.sequences)
        expected = [candidate for sequence in self.sequences for candidate in self.predictions[sequence.target_id]]
        self.assertEqual(len(candidates), 10)
        for loaded, original in zip(candidates, expected):
            self.assertEqual(loaded.target_id, original.target_id)
            self.assertEqual(loaded.model, "drfold2")
            np.testing.assert_array_equal(loaded.coordinates, original.coordinates)
            self.assertIsNone(loaded.energy)
            self.assertIsNone(loaded.energy_group)
        return candidates

    def test_roundtrip_preserves_values_and_records_source_identity(self):
        result = self.run_import()
        self.assertEqual(result["candidates"], 10)
        candidates = self.assert_preserved_predictions()
        self.assertEqual(candidates[0].candidate_id, "backup__target_with_underscores__slot_01")
        manifest = json.loads((self.destination / "candidates.json").read_text())
        self.assertEqual(set(manifest), {"schema_version", "candidates"})
        self.assertTrue(all(set(entry) == {"candidate_id", "target_id", "model", "path"}
                            for entry in manifest["candidates"]))
        provenance = json.loads((self.destination / "import_provenance.json").read_text())
        for name, path, shape in [("submission", self.submission, [7, 18]),
                                  ("sequences", self.sequence_path, [2, 3])]:
            source = provenance["source_files"][name]
            self.assertEqual(source["filename"], path.name)
            self.assertEqual(source["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual(source["bytes"], path.stat().st_size)
            self.assertEqual(source["shape"], shape)
        self.assertFalse(provenance["model_label"]["identity_verified"])
        self.assertFalse(provenance["coordinates"]["geometry_quality_verified"])
        self.assertEqual(provenance["coordinates"]["transforms_applied"], [])
        self.assertNotIn(str(self.root), json.dumps(provenance))
        self.assertEqual(list(self.destination.glob("*.csv")), [])

    def test_reordered_rows_retain_target_residue_and_slot_correspondence(self):
        self.mutate_rows(lambda rows: rows.__setitem__(slice(1, None), list(reversed(rows[1:]))))
        self.run_import()
        self.assert_preserved_predictions()

    def test_invalid_schema_is_rejected_before_creating_output(self):
        self.mutate_rows(lambda rows: rows[0].__setitem__(3, "x_6"))
        with self.assertRaisesRegex(ValueError, "18-column schema"):
            self.run_import()
        self.assertFalse(self.destination.exists())

    def test_missing_residue_is_rejected_before_creating_output(self):
        self.mutate_rows(lambda rows: rows.pop())
        with self.assertRaisesRegex(ValueError, "missing"):
            self.run_import()
        self.assertFalse(self.destination.exists())

    def test_all_zero_legacy_fallback_is_rejected_for_any_target(self):
        # The first target is valid: discovering a failure later must still
        # leave no partially written first-target arrays behind.
        def zero_last_target(rows):
            for row in rows[1:]:
                if row[0].startswith("other_"):
                    row[3:] = ["0"] * 15
        self.mutate_rows(zero_last_target)
        with self.assertRaisesRegex(ValueError, "known legacy failure fallback"):
            self.run_import()
        self.assertFalse(self.destination.exists())

    def test_one_zero_slot_is_preserved_without_inventing_replacements(self):
        def zero_first_slot(rows):
            for row in rows[1:]:
                row[3:6] = ["0"] * 3
        self.mutate_rows(zero_first_slot)
        self.run_import()
        candidates = load_candidate_manifest(self.destination / "candidates.json", self.sequences)
        np.testing.assert_array_equal(candidates[0].coordinates, np.zeros((4, 3)))
        np.testing.assert_array_equal(candidates[5].coordinates, np.zeros((3, 3)))

    def test_nonempty_output_remains_unchanged(self):
        self.destination.mkdir()
        protected = self.destination / "keep.txt"
        protected.write_bytes(b"user data\x00\xff")
        with self.assertRaisesRegex(ValueError, "new or empty"):
            self.run_import()
        self.assertEqual(protected.read_bytes(), b"user data\x00\xff")
        self.assertEqual(list(self.destination.iterdir()), [protected])

    def test_existing_empty_output_is_supported(self):
        self.destination.mkdir()
        self.run_import()
        self.assert_preserved_predictions()

    def test_empty_model_and_unsafe_prefix_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "nonempty"):
            import_submission(self.sequence_path, self.submission, "  ", self.destination)
        with self.assertRaisesRegex(ValueError, "prefix"):
            import_submission(self.sequence_path, self.submission, "drfold2", self.destination, "../escape")
        self.assertFalse(self.destination.exists())

    def test_cli_import_produces_loadable_bundle(self):
        command = [sys.executable, "-m", "examples.import_submission",
                   "--sequences", str(self.sequence_path), "--submission", str(self.submission),
                   "--model", "drfold2", "--output-dir", str(self.destination), "--prefix", "backup"]
        completed = subprocess.run(command, capture_output=True, text=True, check=False,
                                   cwd=Path(__file__).resolve().parents[1])
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(completed.stdout)["targets"], 2)
        self.assert_preserved_predictions()


if __name__ == "__main__":
    unittest.main()
