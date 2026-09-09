"""Strict candidate, submission and C1' PDB I/O (new reconstruction code)."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Mapping, Sequence

import numpy as np

from .geometry import validate_coordinates
from .selection import Candidate
from .sequence import RNASequence

SUBMISSION_COLUMNS = ["ID", "resname", "resid"] + [
    f"{axis}_{i}" for i in range(1, 6) for axis in "xyz"
]


def ensure_distinct_outputs(outputs: Sequence[str | Path],
                            inputs: Sequence[str | Path]) -> None:
    """Reject output aliases of inputs or one another before any writes."""
    resolved_inputs = [Path(path).resolve() for path in inputs]
    resolved_outputs: list[Path] = []
    for output in outputs:
        path = Path(output).resolve()
        for existing in resolved_inputs + resolved_outputs:
            aliases = path == existing
            if not aliases and path.exists() and existing.exists():
                aliases = path.samefile(existing)
            if aliases:
                raise ValueError(f"Output path would overwrite an input or another output: {output}")
        resolved_outputs.append(path)


def write_json(path: str | Path, value: object) -> None:
    """Replace a JSON artifact atomically; non-finite JSON numbers are forbidden."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    content = json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                         delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(content)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _sequence_map(sequences: Sequence[RNASequence]) -> dict[str, RNASequence]:
    result = {s.target_id: s for s in sequences}
    if not sequences or len(result) != len(sequences):
        raise ValueError("Sequences must be nonempty and have unique target IDs")
    return result


def load_candidate_manifest(path: str | Path,
                            sequences: Sequence[RNASequence]) -> list[Candidate]:
    """Load numeric arrays from safe paths relative to a versioned JSON manifest.

    Absolute paths, paths escaping the manifest directory (including symlinks),
    duplicate IDs, unexpected targets and unknown fields are rejected.
    """
    path = Path(path).resolve()
    targets = _sequence_map(sequences)
    with path.open(encoding="utf-8") as stream:
        data = json.load(stream)
    if not isinstance(data, dict) or set(data) != {"schema_version", "candidates"}:
        raise ValueError("Manifest must contain only schema_version and candidates")
    if type(data["schema_version"]) is not int or data["schema_version"] != 1:
        raise ValueError("Only manifest schema_version 1 is supported")
    entries = data["candidates"]
    if not isinstance(entries, list) or not entries:
        raise ValueError("Manifest candidates must be a nonempty list")
    required = {"candidate_id", "target_id", "model", "path"}
    allowed = required | {"energy", "energy_group"}
    ids: set[str] = set()
    observed_targets: set[str] = set()
    result = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict) or not required <= set(entry) <= allowed:
            raise ValueError(f"Candidate {index}: missing or unknown manifest fields")
        for key in required:
            if not isinstance(entry[key], str) or not entry[key].strip():
                raise ValueError(f"Candidate {index}: {key} must be a nonempty string")
        candidate_id = entry["candidate_id"]
        if candidate_id in ids:
            raise ValueError(f"Duplicate candidate_id: {candidate_id}")
        ids.add(candidate_id)
        target_id = entry["target_id"]
        if target_id not in targets:
            raise ValueError(f"Unexpected target_id: {target_id}")
        observed_targets.add(target_id)
        relative = Path(entry["path"])
        if relative.is_absolute():
            raise ValueError("Candidate paths must be relative to the manifest")
        array_path = (path.parent / relative).resolve()
        if not array_path.is_relative_to(path.parent):
            raise ValueError(f"Candidate path escapes manifest directory: {relative}")
        if array_path.suffix.lower() != ".npy":
            raise ValueError("Candidate coordinates must use the .npy format")
        energy = entry.get("energy")
        if energy is not None and (isinstance(energy, bool) or
                                   not isinstance(energy, (int, float)) or
                                   not math.isfinite(energy)):
            raise ValueError(f"Candidate {candidate_id}: energy must be finite or null")
        group = entry.get("energy_group")
        if group is not None and (not isinstance(group, str) or not group.strip()):
            raise ValueError(f"Candidate {candidate_id}: energy_group must be nonempty or null")
        coordinates = np.load(array_path, allow_pickle=False)
        if coordinates.dtype.kind not in "fiu":
            raise ValueError("Candidate arrays must contain real numeric coordinates")
        coordinates = validate_coordinates(coordinates, len(targets[target_id].sequence))
        result.append(Candidate(candidate_id, target_id, entry["model"],
                                coordinates, energy=energy, energy_group=group))
    if observed_targets != set(targets):
        raise ValueError(f"Missing candidate targets: {sorted(set(targets) - observed_targets)}")
    return result


def write_submission(path: str | Path, sequences: Sequence[RNASequence],
                     selections: Mapping[str, Sequence[Candidate]]) -> None:
    """Write five explicitly selected candidates per nucleotide, never padding."""
    targets = _sequence_map(sequences)
    if set(selections) != set(targets):
        raise ValueError("Selection target IDs must exactly match the sequences")
    checked = {}
    for target_id, seq in targets.items():
        chosen = list(selections[target_id])
        if len(chosen) != 5 or len({c.candidate_id for c in chosen}) != 5:
            raise ValueError(f"Target {target_id} requires five distinct candidate IDs")
        if any(c.target_id != target_id for c in chosen):
            raise ValueError(f"Candidate target mismatch for {target_id}")
        checked[target_id] = [validate_coordinates(c.coordinates, len(seq.sequence))
                              for c in chosen]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="",
                                         dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            writer = csv.writer(stream)
            writer.writerow(SUBMISSION_COLUMNS)
            for seq in sequences:
                arrays = checked[seq.target_id]
                for offset, resname in enumerate(seq.sequence):
                    coordinates = [format(float(value), ".10g")
                                   for xyz in arrays for value in xyz[offset]]
                    writer.writerow([f"{seq.target_id}_{offset + 1}", resname,
                                     offset + 1, *coordinates])
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def validate_submission(path: str | Path,
                        sequences: Sequence[RNASequence]) -> dict[str, object]:
    """Check the exact 18-column schema, sequence metadata and finite coordinates.

    Input row order may vary. Each expected ID must occur exactly once.
    """
    _sequence_map(sequences)
    expected = {f"{seq.target_id}_{offset + 1}": (base, offset + 1)
                for seq in sequences for offset, base in enumerate(seq.sequence)}
    seen: set[str] = set()
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != SUBMISSION_COLUMNS:
            raise ValueError("Submission must have the exact 18-column schema and order")
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Line {line}: wrong number of columns")
            identifier = row["ID"]
            if identifier not in expected:
                raise ValueError(f"Line {line}: unexpected residue ID {identifier}")
            if identifier in seen:
                raise ValueError(f"Line {line}: duplicate residue ID {identifier}")
            seen.add(identifier)
            base, resid = expected[identifier]
            if row["resname"] != base or row["resid"] != str(resid):
                raise ValueError(f"Line {line}: sequence or residue index mismatch")
            try:
                xyz = np.array([float(row[key]) for key in SUBMISSION_COLUMNS[3:]])
            except ValueError as exc:
                raise ValueError(f"Line {line}: nonnumeric coordinates") from exc
            validate_coordinates(xyz.reshape(5, 3), 5)
    missing = set(expected) - seen
    if missing:
        raise ValueError(f"Submission is missing {len(missing)} residue rows")
    return {"valid": True, "targets": len(sequences), "residues": len(seen),
            "candidates_per_target": 5, "columns": len(SUBMISSION_COLUMNS)}


def write_c1_pdb(path: str | Path, sequence: str, coordinates: np.ndarray) -> None:
    """Export one C1' atom per RNA residue in single-chain, single-model PDB."""
    if not sequence or set(sequence) - set("ACGU"):
        raise ValueError("PDB sequence must be nonempty uppercase A/C/G/U")
    xyz = validate_coordinates(coordinates, len(sequence))
    if len(sequence) > 9999:
        raise ValueError("PDB output is limited to 9999 residues; use another format")
    if any(len(f"{float(value):8.3f}") != 8 for value in xyz.flat):
        raise ValueError("Coordinates exceed PDB fixed-width coordinate limits")
    lines = []
    for i, (base, (x, y, z)) in enumerate(zip(sequence, xyz), 1):
        lines.append(f"ATOM  {i:5d}  C1' {base:>3s} A{i:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n")
    lines.extend(["TER\n", "END\n"])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="ascii")


def read_c1_pdb(path: str | Path) -> tuple[str, np.ndarray]:
    """Read single-chain PDB C1' traces with blank/A alternate locations.

    Blank alternate locations take priority over A. Other alternates are
    ignored; unsupported residues, multiple chains/models and duplicate atom
    records within the same alternate location are rejected. Residue numbering
    may be nonconsecutive; sequence and coordinates follow file order. Original
    residue identifiers are not returned. A read/write roundtrip renumbers
    residues from 1, so retain original PDBs for residue-index evaluation.
    """
    residues: dict[tuple[str, str, str], dict[str, tuple[str, list[float]]]] = {}
    models = 0
    with Path(path).open(encoding="ascii") as stream:
        for line in stream:
            record = line[:6].strip()
            if record == "MODEL":
                models += 1
                if models > 1:
                    raise ValueError("Only single-model PDB files are supported")
            if record not in {"ATOM", "HETATM"} or line[12:16].strip() != "C1'":
                continue
            alternate = line[16:17]
            if alternate not in {" ", "A"}:
                continue
            base = line[17:20].strip()
            if base not in "ACGU" or len(base) != 1:
                raise ValueError(f"Unsupported RNA PDB residue: {base}")
            key = (line[21:22], line[22:26], line[26:27])
            variants = residues.setdefault(key, {})
            if alternate in variants:
                raise ValueError(f"Duplicate C1' atom for residue {key}")
            if variants and any(item[0] != base for item in variants.values()):
                raise ValueError(f"Inconsistent residue names across alternate locations: {key}")
            try:
                xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError as exc:
                raise ValueError("Invalid PDB coordinate fields") from exc
            variants[alternate] = (base, xyz)
    if not residues:
        raise ValueError("PDB contains no supported C1' atoms")
    if len({key[0] for key in residues}) != 1:
        raise ValueError("Only single-chain PDB files are supported")
    chosen = [variants.get(" ", variants.get("A")) for variants in residues.values()]
    sequence = "".join(item[0] for item in chosen)
    xyz = validate_coordinates(np.array([item[1] for item in chosen], dtype=float))
    return sequence, xyz
