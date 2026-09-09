"""Export validated sequences to FASTA and the upstream Protenix job schema.

This creates input files only. It does not install models, fetch weights, run
neural inference, or reproduce historical training. Native output structures
must be mapped back to one C1' coordinate per original residue before assembly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from rna_folding.sequence import load_sequences, route_model, segment_ranges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sequences = load_sequences(args.sequences)
    destination = args.output_dir
    fasta_dir = destination / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    jobs, manifest = [], []
    planned = []
    for record in sequences:
        windows = segment_ranges(len(record.sequence)) if len(record.sequence) > 480 else [(0, len(record.sequence))]
        for index, (start, end) in enumerate(windows):
            name = f"{record.target_id}__window_{index:03d}"
            subsequence = record.sequence[start:end]
            planned.append((fasta_dir / f"{name}.fasta", f">{name}\n{subsequence}\n"))
            jobs.append({"name": name, "sequences": [{"rnaSequence": {"sequence": subsequence, "count": 1}}]})
            manifest.append({"name": name, "target_id": record.target_id, "start": start,
                             "end": end, "indexing": "0-based half-open", "preferred_model": route_model(len(record.sequence))})
    planned.extend([
        (destination / "protenix_jobs.json", json.dumps(jobs, indent=2) + "\n"),
        (destination / "window_manifest.json", json.dumps({"schema_version": 1, "jobs": manifest}, indent=2) + "\n"),
    ])
    existing = [str(path) for path, _ in planned if path.exists()]
    if existing:
        parser.error("output files already exist; choose an unused directory: " + existing[0])
    for path, content in planned:
        path.write_text(content, encoding="utf-8")
    print(json.dumps({"targets": len(sequences), "jobs": len(jobs), "output_dir": str(destination),
                      "scope": "input export only; no folding model was executed"}, indent=2))


if __name__ == "__main__":
    main()
