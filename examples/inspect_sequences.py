"""Summarize a supplied sequence CSV without opening structure labels.

From the repository root::

    python -m examples.inspect_sequences --sequences data/test_sequences.csv \
        --output outputs/sequence_summary.json

The output describes only the input file. Routing counts use the reconstructed
workflow's documented length rules, not benchmark-derived model performance.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from statistics import median
import sys
from typing import Sequence

# Permit both `python -m examples.inspect_sequences` and direct script usage.
if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rna_folding.sequence import RNASequence, load_sequences, route_model


def summarize_sequences(records: Sequence[RNASequence]) -> dict[str, object]:
    """Calculate unweighted target lengths and residue-weighted composition."""
    if not records:
        raise ValueError("At least one sequence is required")
    lengths = [len(record.sequence) for record in records]
    total = sum(lengths)
    bases = Counter(base for record in records for base in record.sequence)
    routes = Counter(route_model(length) for length in lengths)
    return {
        "n_targets": len(records),
        "total_residues": total,
        "length_nt": {
            "min": min(lengths),
            "max": max(lengths),
            "median": median(lengths),
        },
        "route_counts": {
            name: routes[name] for name in ("trrosettarna", "drfold2", "protenix")
        },
        "route_boundaries_nt": {
            "trrosettarna": "1-100",
            "drfold2": "101-480",
            "protenix": ">480",
        },
        "base_counts": {base: bases[base] for base in "ACGU"},
        "base_proportions": {base: bases[base] / total for base in "ACGU"},
        "composition_denominator": "all residues in the supplied sequence CSV",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", required=True, type=Path,
                        help="CSV containing target_id and sequence columns")
    parser.add_argument("--output", required=True, type=Path,
                        help="Destination JSON file")
    args = parser.parse_args(argv)
    try:
        if args.sequences.resolve() == args.output.resolve():
            raise ValueError("Output must differ from the source sequence file")
        summary = summarize_sequences(load_sequences(args.sequences))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Sequence inspection failed: {exc}\n")
    print(f"Wrote {args.output} ({summary['n_targets']} targets, "
          f"{summary['total_residues']} residues)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
