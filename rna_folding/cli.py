"""CLI for the post-competition reconstruction and synthetic demonstration."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .io import ensure_distinct_outputs, validate_submission, write_json
from .metrics import evaluate_usalign
from .pipeline import assemble, build_plan, run_demo
from .sequence import load_sequences


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        description="Reconstructed RNA workflow; demo data are synthetic, not trained model outputs."
    )
    root.add_argument("--version", action="version", version=__version__)
    commands = root.add_subparsers(dest="command", required=True)
    plan = commands.add_parser("plan", help="Record model routing and segmentation without inference")
    plan.add_argument("--sequences", required=True, help="CSV with target_id and sequence columns")
    plan.add_argument("--output", required=True, help="Output plan JSON")
    export = commands.add_parser("assemble", help="Select five supplied candidates and export CSV")
    export.add_argument("--sequences", required=True)
    export.add_argument("--manifest", required=True, help="Version 1 candidate JSON manifest")
    export.add_argument("--output", required=True, help="Output submission CSV")
    export.add_argument("--audit", help="Optional per-target selection provenance JSON")
    validate = commands.add_parser("validate", help="Strictly validate a submission against sequences")
    validate.add_argument("--sequences", required=True)
    validate.add_argument("--submission", required=True)
    demo = commands.add_parser("demo", help="Run a deterministic synthetic pipeline smoke test")
    demo.add_argument("--output-dir", required=True)
    demo.add_argument("--seed", type=int, default=2025)
    evaluate = commands.add_parser("evaluate", help="Run optional USalign on a pair of C1' PDB traces")
    evaluate.add_argument("--predicted", required=True, help="Predicted single-chain PDB")
    evaluate.add_argument("--reference", required=True, help="Reference single-chain PDB")
    evaluate.add_argument("--usalign", required=True, help="Installed USalign executable name/path")
    evaluate.add_argument("--by-residue", action="store_true",
                          help="Force matching residue indices with -TMscore 1")
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "plan":
            ensure_distinct_outputs([args.output], [args.sequences])
            result = build_plan(load_sequences(args.sequences))
            write_json(args.output, result)
        elif args.command == "assemble":
            result = assemble(args.sequences, args.manifest, args.output, args.audit)
        elif args.command == "validate":
            result = validate_submission(args.submission, load_sequences(args.sequences))
        elif args.command == "demo":
            result = run_demo(args.output_dir, args.seed)
        else:
            result = evaluate_usalign(args.predicted, args.reference, args.usalign,
                                      by_residue=args.by_residue)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))
    return 0
