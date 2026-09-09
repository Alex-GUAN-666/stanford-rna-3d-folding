"""Verify the public CPU workflow in the active environment.

Install the package first: python -m pip install .
Then: python scripts/verify_repo.py --output outputs/verification.json
Optional: install requirements-notebooks.txt and add --notebook.
Historical GPU references and the real USalign binary are not executed.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import sysconfig
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/verification.json")
    parser.add_argument("--notebook", action="store_true")
    args = parser.parse_args()
    report = {
        "schema_version": 1,
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "platform": platform.system(),
        "isolated_venv": sys.prefix != sys.base_prefix,
        "scope": "CPU software behavior on synthetic inputs; not a historical score reproduction",
        "checks": [],
        "passed": False,
    }
    # Bind the record to the tested source bytes, without hashing this report.
    paths = [ROOT / "pyproject.toml", ROOT / "requirements.txt", Path(__file__).resolve()]
    for folder in ["rna_folding", "examples", "tests"]:
        paths.extend(sorted((ROOT / folder).glob("*.py")))
    paths += [ROOT / "notebooks/01_workflow_walkthrough.ipynb"]
    report["source_sha256"] = {str(p.relative_to(ROOT)).replace(os.sep, "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        import numpy as np
        report["numpy"] = np.__version__
        report["package_version"] = importlib.metadata.version("rna-folding-toolkit")
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env["PYTHONNOUSERSITE"] = "1"
        with tempfile.TemporaryDirectory(prefix="rna-reader-") as tmp:
            work = Path(tmp)

            def clean(text: str) -> str:
                for old, new in [(str(ROOT), "{checkout}"), (str(work), "{work}"), (sys.executable, "{python}")]:
                    text = text.replace(old, new)
                return text

            def run(name: str, command: list[str], cwd: Path = ROOT) -> str:
                started = time.monotonic()
                result = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True, timeout=180)
                report["checks"].append({"name": name, "command": [clean(s) for s in command],
                                          "returncode": result.returncode, "seconds": round(time.monotonic() - started, 3),
                                          "stdout": clean(result.stdout), "stderr": clean(result.stderr)})
                print(f"{name}: {'PASS' if result.returncode == 0 else 'FAIL'}", flush=True)
                if result.returncode:
                    raise RuntimeError(f"{name} failed; inspect the recorded stdout/stderr")
                return result.stdout

            run("unit_and_integration_tests", [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"])
            # A non-editable installation must also work away from the checkout.
            installed = json.loads(run("installed_package_outside_checkout", [sys.executable, "-c",
                "import json,rna_folding,hashlib; from pathlib import Path; "
                "package=Path(rna_folding.__file__).parent; "
                "print(json.dumps({'version':rna_folding.__version__,'module_path':rna_folding.__file__,"
                "'source_sha256':{'rna_folding/'+p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in package.glob('*.py')}}))"], work))
            report["installed_package_in_site_packages"] = "site-packages" in installed["module_path"].replace("\\", "/")
            if not report["installed_package_in_site_packages"]:
                raise RuntimeError("Install with 'python -m pip install .' (non-editable) to verify the distributed package")
            expected_package_hashes = {name: value for name, value in report["source_sha256"].items() if name.startswith("rna_folding/")}
            report["installed_package_matches_checkout"] = installed["source_sha256"] == expected_package_hashes
            if not report["installed_package_matches_checkout"]:
                raise RuntimeError("Installed package differs from this checkout; rerun 'python -m pip install .'")
            executable = Path(sysconfig.get_path("scripts")) / ("rna-folding.exe" if os.name == "nt" else "rna-folding")
            run("installed_console_entrypoint", [str(executable), "--version"], work)
            demo = work / "demo"
            summary = json.loads(run("installed_console_demo", [str(executable), "demo", "--output-dir", str(demo)], work))
            sequences, submission = demo / "sequences.csv", demo / "submission.csv"
            cli = [sys.executable, "-m", "rna_folding"]
            run("submission_validation", cli + ["validate", "--sequences", str(sequences), "--submission", str(submission)])
            run("routing_plan", cli + ["plan", "--sequences", str(sequences), "--output", str(work / "plan.json")])
            run("sequence_inspection", [sys.executable, "-m", "examples.inspect_sequences", "--sequences", str(sequences), "--output", str(work / "sequence_summary.json")])
            run("model_input_export", [sys.executable, "-m", "examples.export_model_inputs", "--sequences", str(sequences), "--output-dir", str(work / "model_inputs")])
            bundle = work / "imported"
            run("submission_import", [sys.executable, "-m", "examples.import_submission", "--sequences", str(sequences), "--submission", str(submission), "--model", "synthetic_import", "--output-dir", str(bundle)])
            rebuilt = work / "rebuilt.csv"
            run("reassemble_imported_candidates", cli + ["assemble", "--sequences", str(sequences), "--manifest", str(bundle / "candidates.json"), "--output", str(rebuilt), "--audit", str(work / "rebuilt_audit.json")])
            run("rebuilt_submission_validation", cli + ["validate", "--sequences", str(sequences), "--submission", str(rebuilt)])

            def rows(path: Path):
                with path.open(newline="", encoding="utf-8") as stream:
                    return {row["ID"]: row for row in csv.DictReader(stream)}

            before, after = rows(submission), rows(rebuilt)
            if before.keys() != after.keys():
                raise RuntimeError("Import/assembly roundtrip changed residue IDs")
            columns = [f"{axis}_{slot}" for slot in range(1, 6) for axis in "xyz"]
            for key in before:
                if before[key]["resname"] != after[key]["resname"] or before[key]["resid"] != after[key]["resid"]:
                    raise RuntimeError("Roundtrip changed residue metadata")
                np.testing.assert_allclose([float(before[key][c]) for c in columns],
                                           [float(after[key][c]) for c in columns], rtol=0, atol=1e-6)
            expected = {"valid": True, "targets": 3, "residues": 752, "candidates_per_target": 5, "columns": 18}
            if summary["validation"] != expected:
                raise RuntimeError(f"Unexpected demo result: {summary['validation']}")
            report["demo_validation"] = summary["validation"]
            report["import_assembly_roundtrip"] = {"passed": True, "absolute_tolerance": 1e-6, "residue_ids_and_metadata_preserved": True}
            if args.notebook:
                import nbformat
                from nbclient import NotebookClient
                from jupyter_client import KernelManager
                notebook_files = sorted((ROOT / "notebooks").glob("*.ipynb")) + sorted((ROOT / "historical/reference").glob("*.ipynb"))
                for source_notebook in notebook_files:
                    nbformat.validate(nbformat.read(source_notebook, as_version=4))
                report["notebook_schema_validation"] = {
                    "passed": True,
                    "files": [str(p.relative_to(ROOT)).replace(os.sep, "/") for p in notebook_files],
                    "scope": "JSON schema validation only for historical references; no GPU cells executed",
                }
                notebook_path = ROOT / "notebooks/01_workflow_walkthrough.ipynb"
                nb = nbformat.read(notebook_path, as_version=4)
                km = KernelManager(kernel_name="python3")
                km.kernel_spec.argv = [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"]
                started = time.monotonic()
                client = NotebookClient(nb, km=km, timeout=180, resources={"metadata": {"path": str(notebook_path.parent)}})
                executed = client.execute(cleanup_kc=True)
                code_cells = [cell for cell in executed.cells if cell.cell_type == "code"]
                if len(code_cells) != 7 or any(cell.execution_count is None for cell in code_cells):
                    raise RuntimeError("Expected all seven walkthrough code cells to execute")
                executed_path = output.parent / "executed_walkthrough.ipynb"
                nbformat.write(executed, executed_path)
                report["notebook"] = {"passed": True, "code_cells_executed": len(code_cells),
                                      "kernel": "active Python environment", "seconds": round(time.monotonic() - started, 3),
                                      "artifact": executed_path.name}
                print("walkthrough_notebook: PASS", flush=True)
            else:
                report["notebook"] = {"executed": False, "how_to_run": "Install requirements-notebooks.txt and add --notebook"}
        report["passed"] = True
    except Exception as exc:
        report["error"] = f"{type(exc).__name__}: {exc}"
        print(report["error"], file=sys.stderr)
    report["not_executed"] = ["Historical GPU notebooks/scripts", "Training or fine-tuning", "Real USalign binary", "Historical leaderboard score reproduction"]
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Verification report: {output}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
