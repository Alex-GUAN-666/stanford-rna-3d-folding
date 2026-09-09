"""Corrupt-input and end-to-end checks for the new reconstruction."""

import contextlib
import csv
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from rna_folding.cli import main
from rna_folding.io import (SUBMISSION_COLUMNS, load_candidate_manifest, read_c1_pdb,
                            validate_submission, write_c1_pdb, write_json,
                            write_submission)
from rna_folding.metrics import evaluate_usalign, parse_reference_tm_score
from rna_folding.pipeline import assemble, build_plan, run_demo
from rna_folding.selection import Candidate
from rna_folding.sequence import RNASequence, load_sequences


class IOTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.sequence = RNASequence("target_with_underscores", "ACGU")
        self.xyz = np.array([[0., 0., 0.], [1., 2., 3.], [2., 0., 4.], [4., 3., 1.]])
        self.candidates = [Candidate(f"candidate_{i}", self.sequence.target_id,
                                     "trrosettarna", self.xyz + i)
                           for i in range(5)]
        self.submission = self.root / "submission.csv"
        write_submission(self.submission, [self.sequence],
                         {self.sequence.target_id: self.candidates})
        np.save(self.root / "coordinates.npy", self.xyz)
        self.manifest = self.root / "candidates.json"
        self.manifest_value = {"schema_version": 1, "candidates": [{
            "candidate_id": "candidate_0", "target_id": self.sequence.target_id,
            "model": "trrosettarna", "path": "coordinates.npy",
        }]}
        write_json(self.manifest, self.manifest_value)

    def mutate_rows(self, mutate):
        with self.submission.open(newline="") as stream:
            rows = list(csv.reader(stream))
        mutate(rows)
        with self.submission.open("w", newline="") as stream:
            csv.writer(stream).writerows(rows)

    def test_roundtrip_preserves_candidate_slots(self):
        result = validate_submission(self.submission, [self.sequence])
        self.assertEqual(result["residues"], 4)
        with self.submission.open(newline="") as stream:
            reader = csv.DictReader(stream)
            self.assertEqual(reader.fieldnames, SUBMISSION_COLUMNS)
            rows = list(reader)
        for residue, row in enumerate(rows):
            for slot, candidate in enumerate(self.candidates, 1):
                np.testing.assert_allclose([float(row[f"{axis}_{slot}"]) for axis in "xyz"],
                                           candidate.coordinates[residue])

    def test_duplicate_and_missing_rows_are_rejected(self):
        self.mutate_rows(lambda rows: rows.append(rows[1]))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            validate_submission(self.submission, [self.sequence])
        self.mutate_rows(lambda rows: rows.__delitem__(slice(-2, None)))
        with self.assertRaisesRegex(ValueError, "missing"):
            validate_submission(self.submission, [self.sequence])

    def test_sequence_metadata_must_match(self):
        self.mutate_rows(lambda rows: rows[1].__setitem__(1, "U"))
        with self.assertRaisesRegex(ValueError, "mismatch"):
            validate_submission(self.submission, [self.sequence])

    def test_residue_index_must_be_exact(self):
        self.mutate_rows(lambda rows: rows[1].__setitem__(2, "1.0"))
        with self.assertRaisesRegex(ValueError, "mismatch"):
            validate_submission(self.submission, [self.sequence])

    def test_nonfinite_and_sentinel_coordinates_rejected(self):
        for bad in ["nan", "inf", "-inf", "1e18"]:
            self.mutate_rows(lambda rows: rows[1].__setitem__(3, bad))
            with self.subTest(value=bad), self.assertRaises(ValueError):
                validate_submission(self.submission, [self.sequence])

    def test_wrong_column_count_rejected(self):
        self.mutate_rows(lambda rows: rows[1].append("extra"))
        with self.assertRaisesRegex(ValueError, "number of columns"):
            validate_submission(self.submission, [self.sequence])

    def test_wrong_schema_rejected(self):
        self.mutate_rows(lambda rows: rows[0].__setitem__(3, "x"))
        with self.assertRaisesRegex(ValueError, "schema"):
            validate_submission(self.submission, [self.sequence])

    def test_no_padding_or_repeated_candidate_id(self):
        for bad in [self.candidates[:4], self.candidates[:4] + self.candidates[:1]]:
            with self.subTest(count=len(bad)), self.assertRaisesRegex(ValueError, "five distinct"):
                write_submission(self.submission, [self.sequence], {self.sequence.target_id: bad})

    def test_manifest_load_and_array_length(self):
        loaded = load_candidate_manifest(self.manifest, [self.sequence])
        np.testing.assert_equal(loaded[0].coordinates, self.xyz)
        np.save(self.root / "coordinates.npy", self.xyz[:-1])
        with self.assertRaisesRegex(ValueError, "length"):
            load_candidate_manifest(self.manifest, [self.sequence])

    def test_manifest_duplicate_id_and_unknown_target(self):
        self.manifest_value["candidates"].append(dict(self.manifest_value["candidates"][0]))
        write_json(self.manifest, self.manifest_value)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            load_candidate_manifest(self.manifest, [self.sequence])
        self.manifest_value["candidates"] = self.manifest_value["candidates"][:1]
        self.manifest_value["candidates"][0]["target_id"] = "unknown"
        write_json(self.manifest, self.manifest_value)
        with self.assertRaisesRegex(ValueError, "Unexpected target"):
            load_candidate_manifest(self.manifest, [self.sequence])

    def test_manifest_rejects_path_escape_and_symlink(self):
        for candidate_path in ["../outside.npy", str(self.root / "coordinates.npy")]:
            self.manifest_value["candidates"][0]["path"] = candidate_path
            write_json(self.manifest, self.manifest_value)
            with self.subTest(path=candidate_path), self.assertRaises(ValueError):
                load_candidate_manifest(self.manifest, [self.sequence])
        nested = self.root / "nested"
        nested.mkdir()
        (nested / "link.npy").symlink_to(self.root / "coordinates.npy")
        self.manifest_value["candidates"][0]["path"] = "link.npy"
        write_json(nested / "candidates.json", self.manifest_value)
        with self.assertRaisesRegex(ValueError, "escapes"):
            load_candidate_manifest(nested / "candidates.json", [self.sequence])

    def test_manifest_rejects_pickle_and_nonreal_array(self):
        for array in [self.xyz.astype(object), self.xyz.astype(complex), self.xyz.astype(str)]:
            np.save(self.root / "coordinates.npy", array)
            with self.subTest(dtype=str(array.dtype)), self.assertRaises(ValueError):
                load_candidate_manifest(self.manifest, [self.sequence])

    def test_manifest_rejects_metadata_typos_and_boolean_energy(self):
        entry = self.manifest_value["candidates"][0]
        entry["energy"] = True
        write_json(self.manifest, self.manifest_value)
        with self.assertRaisesRegex(ValueError, "energy"):
            load_candidate_manifest(self.manifest, [self.sequence])
        del entry["energy"]
        entry["energi"] = 1
        write_json(self.manifest, self.manifest_value)
        with self.assertRaisesRegex(ValueError, "unknown"):
            load_candidate_manifest(self.manifest, [self.sequence])

    def test_pdb_roundtrip_and_alternate_priority(self):
        pdb = self.root / "trace.pdb"
        write_c1_pdb(pdb, self.sequence.sequence, self.xyz)
        sequence, xyz = read_c1_pdb(pdb)
        self.assertEqual(sequence, self.sequence.sequence)
        np.testing.assert_allclose(xyz, self.xyz, atol=0.00051)
        lines = pdb.read_text().splitlines(keepends=True)
        alternate = lines[0][:16] + "A" + lines[0][17:30] + f"{99.:8.3f}" + lines[0][38:]
        pdb.write_text(alternate + "".join(lines))
        _, xyz = read_c1_pdb(pdb)
        self.assertEqual(xyz[0, 0], 0.0)

    def test_pdb_rejects_multiple_chains_and_duplicate_atoms(self):
        pdb = self.root / "trace.pdb"
        write_c1_pdb(pdb, self.sequence.sequence, self.xyz)
        lines = pdb.read_text().splitlines(keepends=True)
        pdb.write_text(lines[0] + "".join(lines))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            read_c1_pdb(pdb)
        lines[1] = lines[1][:21] + "B" + lines[1][22:]
        pdb.write_text("".join(lines))
        with self.assertRaisesRegex(ValueError, "single-chain"):
            read_c1_pdb(pdb)

    def test_reference_normalization_parser_and_ambiguity(self):
        for label in ["Chain", "Structure"]:
            output = (f"TM-score= 0.70000 (normalized by length of {label}_1: L=100)\n"
                      f"TM-score= 0.42000 (normalized by length of {label}_2: L=200)\n")
            self.assertEqual(parse_reference_tm_score(output), 0.42)
            with self.assertRaisesRegex(ValueError, "Expected one"):
                parse_reference_tm_score(output + output)
        with self.assertRaises(ValueError):
            parse_reference_tm_score("TM-score= 1.1 (normalized by length of Chain_2: L=4)")
        with self.assertRaises(ValueError):
            parse_reference_tm_score("TM-score= 0.7 (normalized by length of Chain_1: L=4)")

    def test_usalign_command_and_failure(self):
        pdb = self.root / "trace.pdb"
        write_c1_pdb(pdb, self.sequence.sequence, self.xyz)
        output = "TM-score= 1.00000 (normalized by length of Structure_2: L=4)\n"
        completed = subprocess.CompletedProcess([], 0, stdout=output, stderr="")
        with patch("rna_folding.metrics.shutil.which", return_value="/usr/bin/USalign"), \
                patch("rna_folding.metrics.subprocess.run", return_value=completed) as run:
            result = evaluate_usalign(pdb, pdb, "USalign")
            self.assertEqual(result["tm_score_reference"], 1.0)
            self.assertNotIn("-TMscore", run.call_args.args[0])
            self.assertEqual(run.call_args.args[0][-4:], ["-mol", "RNA", "-atom", " C1'"])
            self.assertNotIn("shell", run.call_args.kwargs)
            evaluate_usalign(pdb, pdb, "USalign", by_residue=True)
            self.assertEqual(run.call_args.args[0][-2:], ["-TMscore", "1"])
            run.return_value = subprocess.CompletedProcess([], 1, stdout="", stderr="bad input")
            with self.assertRaisesRegex(RuntimeError, "bad input"):
                evaluate_usalign(pdb, pdb, "USalign")

    def test_cli_reports_bad_input_without_success(self):
        with contextlib.redirect_stderr(io.StringIO()) as err:
            result = main(["validate", "--sequences", str(self.root / "absent.csv"),
                           "--submission", str(self.submission)])
        self.assertEqual(result, 2)
        self.assertIn("error:", err.getvalue())

    def test_only_long_plan_targets_are_segmented(self):
        plan = build_plan([RNASequence("medium", "A" * 400), RNASequence("long", "A" * 520)])
        self.assertEqual(plan["targets"][0]["segments"], [{"start": 0, "end": 400}])
        self.assertEqual(len(plan["targets"][1]["segments"]), 2)

    def test_cli_and_assemble_prevent_input_output_collisions(self):
        sequence_path = self.root / "sequences.csv"
        sequence_path.write_text(f"target_id,sequence\n{self.sequence.target_id},ACGU\n")
        original_sequence = sequence_path.read_bytes()
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(["plan", "--sequences", str(sequence_path),
                                   "--output", str(sequence_path)]), 2)
        self.assertEqual(sequence_path.read_bytes(), original_sequence)
        entry = self.manifest_value["candidates"][0]
        self.manifest_value["candidates"] = [dict(entry, candidate_id=f"c{i}") for i in range(5)]
        write_json(self.manifest, self.manifest_value)
        for protected in [sequence_path, self.manifest, self.root / "coordinates.npy"]:
            original = protected.read_bytes()
            with self.subTest(path=protected), self.assertRaisesRegex(ValueError, "overwrite"):
                assemble(sequence_path, self.manifest, protected)
            with self.subTest(audit=protected), self.assertRaisesRegex(ValueError, "overwrite"):
                assemble(sequence_path, self.manifest, self.submission, protected)
            self.assertEqual(protected.read_bytes(), original)
        with self.assertRaisesRegex(ValueError, "overwrite"):
            assemble(sequence_path, self.manifest, self.submission, self.submission)
        alias = self.root / "alias.csv"
        alias.symlink_to(sequence_path)
        with self.assertRaisesRegex(ValueError, "overwrite"):
            assemble(sequence_path, self.manifest, alias)

    def test_demo_never_overwrites_existing_inputs(self):
        sequence_path = self.root / "sequences.csv"
        sequence_path.write_text("target_id,sequence\nreal_target,ACGU\n")
        original = sequence_path.read_bytes()
        with self.assertRaisesRegex(ValueError, "new or empty"):
            run_demo(self.root)
        self.assertEqual(sequence_path.read_bytes(), original)
        self.assertFalse((self.root / "coordinates").exists())

    def test_synthetic_demo_exercises_three_routes_and_stitching(self):
        directory = self.root / "demo"
        summary = run_demo(directory, seed=2025)
        self.assertTrue(summary["synthetic"])
        self.assertEqual(summary["validation"]["residues"], 752)
        self.assertEqual(summary["segmented_candidates"], 7)
        self.assertLess(summary["max_absolute_stitch_roundtrip_error"], 1e-8)
        sequences = load_sequences(directory / "sequences.csv")
        candidates = load_candidate_manifest(directory / "candidates.json", sequences)
        self.assertEqual(len(candidates), 21)
        audit = json.loads((directory / "audit.json").read_text())
        self.assertEqual({target["preferred_model"] for target in audit["targets"]},
                         {"trrosettarna", "drfold2", "protenix"})
        for target in audit["targets"]:
            chosen = target["selected"]
            self.assertEqual(len({c["candidate_id"] for c in chosen}), 5)
            self.assertEqual(len({c["model"] for c in chosen}), 3)
        first_csv = (directory / "submission.csv").read_bytes()
        second_directory = self.root / "demo-repeat"
        run_demo(second_directory, seed=2025)
        self.assertEqual(first_csv, (second_directory / "submission.csv").read_bytes())


if __name__ == "__main__":
    unittest.main()
