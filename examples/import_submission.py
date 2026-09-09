"""Import a restored notebook's submission CSV into the candidate toolchain.

This bridge is a new implementation written after the notebook backups were
restored. It does not run or verify a folding model. It retains the five source
coordinate slots, without alignment, perturbation, ranking, or invented energy.
Run from the repository root with ``python -m examples.import_submission``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile

import numpy as np

from rna_folding.io import SUBMISSION_COLUMNS, validate_submission, write_json
from rna_folding.sequence import load_sequences


def _source_record(path: Path) -> dict[str, object]:
    """Record source identity without embedding private local directories."""
    contents = path.read_bytes()
    return {
        "filename": path.name,
        "sha256": hashlib.sha256(contents).hexdigest(),
        "bytes": len(contents),
    }


def _check_destination(destination: Path) -> None:
    if destination.is_symlink() or (
        destination.exists()
        and (not destination.is_dir() or any(destination.iterdir()))
    ):
        raise ValueError("Output directory must be new or empty; existing files are never overwritten")


def import_submission(
    sequences_path: str | Path,
    submission_path: str | Path,
    model: str,
    output_dir: str | Path,
    prefix: str = "imported",
) -> dict[str, object]:
    """Validate an existing CSV and atomically publish a candidate bundle.

    The model label is supplied by the caller, stripped and lowercased; it is
    not inferred from coordinates and is not evidence of model provenance.
    Values are parsed as float64 and saved unchanged from that representation.
    CSV textual formatting is not preserved. Repeated coordinate slots are
    retained; this function makes no claim of conformational diversity.
    """
    if not isinstance(model, str) or not model.strip():
        raise ValueError("model must be a nonempty user-supplied label")
    model_label = model.strip().lower()
    if not isinstance(prefix, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", prefix) or prefix in {".", ".."}:
        raise ValueError("prefix must use letters, digits, underscores, dots or hyphens")
    sequences_path = Path(sequences_path)
    submission_path = Path(submission_path)
    destination = Path(output_dir)
    _check_destination(destination)

    sources = {
        "sequences": _source_record(sequences_path),
        "submission": _source_record(submission_path),
    }
    sequences = load_sequences(sequences_path)
    validation = validate_submission(submission_path, sequences)
    with sequences_path.open(encoding="utf-8-sig", newline="") as stream:
        columns = next(csv.reader(stream))
    sources["sequences"].update({"shape": [len(sequences), len(columns)], "columns": columns})
    sources["submission"].update({
        "shape": [validation["residues"], len(SUBMISSION_COLUMNS)],
        "columns": SUBMISSION_COLUMNS,
    })

    # The existing validator permits reordered rows, so residue identity must
    # come from the full ID rather than row position or splitting on underscores.
    with submission_path.open(encoding="utf-8", newline="") as stream:
        rows = {row["ID"]: row for row in csv.DictReader(stream)}
    entries, arrays = [], []
    for target_index, sequence in enumerate(sequences):
        xyz = np.array([
            [float(rows[f"{sequence.target_id}_{resid}"][column])
             for column in SUBMISSION_COLUMNS[3:]]
            for resid in range(1, len(sequence.sequence) + 1)
        ], dtype=np.float64).reshape(len(sequence.sequence), 5, 3)
        if np.all(xyz == 0):
            raise ValueError(
                f"Target {sequence.target_id}: all-zero coordinates across all five "
                "candidates are a known legacy failure fallback and cannot be imported. "
                "This check does not establish general geometry validity."
            )
        for slot in range(1, 6):
            relative_path = f"coordinates/target_{target_index:05d}_slot_{slot:02d}.npy"
            entries.append({
                "candidate_id": f"{prefix}__{sequence.target_id}__slot_{slot:02d}",
                "target_id": sequence.target_id,
                "model": model_label,
                "path": relative_path,
            })
            arrays.append((relative_path, xyz[:, slot - 1, :]))

    # Reject changing input files before writing any output. Hashes refer to the
    # original CSV bytes, not to a normalized or reconstructed copy.
    for key, path in (("sequences", sequences_path), ("submission", submission_path)):
        if _source_record(path)["sha256"] != sources[key]["sha256"]:
            raise ValueError(f"The {key} CSV changed during import; retry with stable inputs")

    provenance = {
        "schema_version": 1,
        "artifact_type": "submission_import",
        "implementation": "New bridge implemented after restoration of the historical notebook backups",
        "source_files": sources,
        "validation": validation,
        "model_label": {
            "provided": model,
            "stored": model_label,
            "identity_verified": False,
            "meaning": "User-supplied label; model identity is not verified by this import",
        },
        "coordinates": {
            "dtype": "float64",
            "shape_per_candidate": "[target_sequence_length, 3]",
            "source_slot_order": [1, 2, 3, 4, 5],
            "numeric_values_preserved": True,
            "transforms_applied": [],
            "energy_metadata_available": False,
            "geometry_quality_verified": False,
        },
        "candidate_count": len(entries),
        "original_csv_copied": False,
    }
    # Stage a complete bundle beside the destination. Remove only an existing
    # empty destination before renaming, since Windows cannot rename over it.
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".rna-import-", dir=destination.parent) as temporary:
        staging = Path(temporary) / "bundle"
        (staging / "coordinates").mkdir(parents=True)
        for relative_path, coordinates in arrays:
            np.save(staging / relative_path, coordinates, allow_pickle=False)
        write_json(staging / "candidates.json", {"schema_version": 1, "candidates": entries})
        write_json(staging / "import_provenance.json", provenance)
        _check_destination(destination)
        if destination.exists():
            # rmdir fails if the directory gained files after the check; it
            # also cannot follow a newly substituted directory symlink.
            destination.rmdir()
        os.rename(staging, destination)

    return {
        "targets": len(sequences),
        "candidates": len(entries),
        "model_label": model_label,
        "manifest": str(destination / "candidates.json"),
        "provenance": str(destination / "import_provenance.json"),
        "scope": "Existing CSV import only; no model execution or score verification",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", type=Path, required=True)
    parser.add_argument("--submission", type=Path, required=True)
    parser.add_argument("--model", required=True,
                        help="User-supplied label, e.g. drfold2; stripped/lowercased, not verified")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="New or empty directory for .npy files, manifest and provenance")
    parser.add_argument("--prefix", default="imported", help="Prefix for candidate IDs")
    args = parser.parse_args(argv)
    try:
        result = import_submission(args.sequences, args.submission, args.model,
                                   args.output_dir, args.prefix)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
