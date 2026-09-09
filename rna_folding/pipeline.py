"""Executable orchestration reconstructed after the competition.

Real model weights and archived inference notebooks are external prerequisites;
the demo uses mathematical curves solely to exercise the software contract.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .geometry import stitch_segments
from .io import (ensure_distinct_outputs, load_candidate_manifest, validate_submission, write_c1_pdb,
                 write_json, write_submission)
from .selection import select_candidates
from .sequence import RNASequence, load_sequences, route_model, segment_ranges

RECONSTRUCTION_NOTICE = (
    "New post-competition reconstruction; not archived competition code or a "
    "reproduction of the historical leaderboard score."
)


def build_plan(sequences: list[RNASequence]) -> dict[str, object]:
    """Record routing and half-open segment ranges without running any models."""
    return {
        "schema_version": 1,
        "provenance": RECONSTRUCTION_NOTICE,
        "routing_status": "Report-informed fixed default policy; not a validated optimum.",
        "segment_convention": "zero-based [start, end), end exclusive",
        "targets": [{
            "target_id": seq.target_id,
            "length": len(seq.sequence),
            "preferred_model": route_model(len(seq.sequence)),
            "segments": [{"start": start, "end": end}
                         for start, end in (segment_ranges(len(seq.sequence))
                                            if len(seq.sequence) > 480
                                            else [(0, len(seq.sequence))])],
        } for seq in sequences],
    }


def assemble(sequences_path: str | Path, manifest_path: str | Path,
             output: str | Path, audit_path: str | Path | None = None) -> dict[str, object]:
    """Select five real input candidates per target and export a validated CSV."""
    sequences = load_sequences(sequences_path)
    candidates = load_candidate_manifest(manifest_path, sequences)
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    input_paths = [sequences_path, manifest_path] + [
        manifest_path.parent / entry["path"] for entry in manifest["candidates"]
    ]
    ensure_distinct_outputs([output] + ([audit_path] if audit_path is not None else []), input_paths)
    selected = {seq.target_id: select_candidates(
        [candidate for candidate in candidates if candidate.target_id == seq.target_id],
        seq, count=5
    ) for seq in sequences}
    write_submission(output, sequences, selected)
    validation = validate_submission(output, sequences)
    audit = {
        "schema_version": 1,
        "provenance": RECONSTRUCTION_NOTICE,
        "sequences_file": Path(sequences_path).name,
        "candidate_manifest": Path(manifest_path).name,
        "submission_file": Path(output).name,
        "selection_policy": (
            "Route-prioritized round robin across models. Energy ordering is limited "
            "to candidates from the same model and declared energy_group. "
            "No coordinate averaging or candidate padding."
        ),
        "validation": validation,
        "targets": [{
            "target_id": seq.target_id,
            "length": len(seq.sequence),
            "preferred_model": route_model(len(seq.sequence)),
            "candidate_count": sum(c.target_id == seq.target_id for c in candidates),
            "selected": [{
                "submission_slot": slot,
                "candidate_id": candidate.candidate_id,
                "model": candidate.model,
                "energy": candidate.energy,
                "energy_group": candidate.energy_group,
            } for slot, candidate in enumerate(selected[seq.target_id], 1)],
        } for seq in sequences],
    }
    if audit_path is not None:
        write_json(audit_path, audit)
    return audit


def run_demo(output_dir: str | Path, seed: int = 2025) -> dict[str, object]:
    """Generate deterministic, explicitly synthetic fixtures and run the pipeline."""
    output_dir = Path(output_dir)
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise ValueError("Demo output directory must be new or empty; choose a fresh directory")
    coordinate_dir = output_dir / "coordinates"
    coordinate_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    sequences = [RNASequence(name, ("ACGU" * ((length + 3) // 4))[:length])
                 for name, length in [("synthetic_short", 72),
                                      ("synthetic_medium", 160),
                                      ("synthetic_long", 520)]]
    sequences_path = output_dir / "sequences.csv"
    with sequences_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["target_id", "sequence"])
        writer.writerows((seq.target_id, seq.sequence) for seq in sequences)
    manifest = {"schema_version": 1, "candidates": []}
    max_stitch_error = 0.0
    segmented_count = 0
    pdb_fixture = None
    for seq in sequences:
        length = len(seq.sequence)
        for candidate_index in range(7):
            t = np.arange(length, dtype=float)
            angle = t * (0.52 + candidate_index * 0.009)
            radius = 9.0 + candidate_index * 0.2
            xyz = np.column_stack((radius * np.cos(angle), radius * np.sin(angle),
                                   t * 1.8 + np.sin(t * 0.13) * candidate_index * 0.1))
            xyz += rng.normal(0.0, 0.01, xyz.shape)
            if length > 480:
                segments = []
                for segment_index, (start, end) in enumerate(segment_ranges(length)):
                    piece = xyz[start:end].copy()
                    if segment_index:
                        theta = 0.31 * segment_index
                        rotation = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                                             [np.sin(theta), np.cos(theta), 0.0],
                                             [0.0, 0.0, 1.0]])
                        piece = piece @ rotation + [17.0, -11.0, 8.0]
                    segments.append((start, piece))
                stitched = stitch_segments(segments, length)
                max_stitch_error = max(max_stitch_error, float(np.max(np.abs(stitched - xyz))))
                xyz = stitched
                segmented_count += 1
            identifier = f"{seq.target_id}_{candidate_index + 1}"
            relative = f"coordinates/{identifier}.npy"
            np.save(output_dir / relative, xyz, allow_pickle=False)
            # Model names label routing branches only; every array is synthetic.
            preferred = route_model(length)
            models = [preferred] + [model for model in ("trrosettarna", "drfold2", "protenix")
                                    if model != preferred]
            manifest["candidates"].append({
                "candidate_id": identifier, "target_id": seq.target_id,
                "model": models[candidate_index % len(models)], "path": relative,
                "energy": float(candidate_index),
                "energy_group": f"synthetic-arbitrary-energy-{seq.target_id}",
            })
            if pdb_fixture is None:
                pdb_fixture = (seq.sequence, xyz.copy())
    write_json(output_dir / "candidates.json", manifest)
    write_json(output_dir / "plan.json", build_plan(sequences))
    audit = assemble(sequences_path, output_dir / "candidates.json",
                     output_dir / "submission.csv", output_dir / "audit.json")
    write_c1_pdb(output_dir / "synthetic_trace.pdb", *pdb_fixture)
    summary = {
        "provenance": RECONSTRUCTION_NOTICE,
        "synthetic": True,
        "description": (
            "Mathematical curves, not RNA model predictions. Model names identify "
            "routing branches only; energies are arbitrary software fixtures."
        ),
        "seed": seed,
        "target_lengths": [len(seq.sequence) for seq in sequences],
        "input_candidates": len(manifest["candidates"]),
        "segmented_candidates": segmented_count,
        "max_absolute_stitch_roundtrip_error": max_stitch_error,
        "validation": audit["validation"],
        "artifacts": ["sequences.csv", "candidates.json", "plan.json", "submission.csv",
                      "audit.json", "synthetic_trace.pdb", "demo_summary.json"],
        "interpretation": "Pipeline smoke test only; no biological accuracy or TM-score measured.",
    }
    write_json(output_dir / "demo_summary.json", summary)
    return summary
