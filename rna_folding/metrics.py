"""Optional external USalign integration, not an official Kaggle evaluator."""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess

from .io import read_c1_pdb

_REFERENCE_TM = re.compile(
    r"^\s*TM-score\s*=\s*([0-9.eE+-]+)\s*\(normalized by length of "
    r"(?:Chain|Structure)_2\b", re.MULTILINE
)


def parse_reference_tm_score(output: str) -> float:
    """Extract exactly the score normalized by the second (reference) input."""
    matches = _REFERENCE_TM.findall(output)
    if len(matches) != 1:
        raise ValueError("Expected one reference-normalized Chain_2/Structure_2 TM-score")
    score = float(matches[0])
    if not 0.0 <= score <= 1.0:
        raise ValueError("Reference-normalized TM-score must lie in [0, 1]")
    return score


def evaluate_usalign(predicted: str | Path, reference: str | Path,
                     usalign: str | Path, *, by_residue: bool = False,
                     timeout: float = 120.0) -> dict[str, object]:
    """Evaluate single-chain C1' PDB traces with a caller-installed USalign.

    Default USalign alignment can search for structural correspondence.
    ``by_residue`` forces matching residue indices via ``-TMscore 1`` and
    requires matching sequences and residue keys in both input files.
    """
    predicted = Path(predicted).resolve(strict=True)
    reference = Path(reference).resolve(strict=True)
    pred_sequence, pred_xyz = read_c1_pdb(predicted)
    ref_sequence, ref_xyz = read_c1_pdb(reference)
    if len(pred_xyz) < 3 or len(ref_xyz) < 3:
        raise ValueError("USalign inputs require at least three C1' atoms")
    if by_residue:
        if pred_sequence != ref_sequence or _residue_keys(predicted) != _residue_keys(reference):
            raise ValueError("--by-residue requires identical sequences and PDB residue keys")
    executable = shutil.which(str(usalign))
    if executable is None:
        raise ValueError(f"USalign executable not found or not executable: {usalign}")
    command = [executable, str(predicted), str(reference), "-mol", "RNA", "-atom", " C1'"]
    if by_residue:
        command.extend(["-TMscore", "1"])
    try:
        completed = subprocess.run(command, capture_output=True, text=True,
                                   check=False, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"USalign timed out after {timeout:g} seconds") from exc
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()[-2000:]
        raise RuntimeError(f"USalign failed with exit code {completed.returncode}: {detail}")
    return {
        "metric": "USalign reference-normalized TM-score",
        "tm_score_reference": parse_reference_tm_score(completed.stdout),
        "normalization": "length of second input (reference)",
        "alignment_mode": "residue-index correspondence" if by_residue else "structural alignment",
        "predicted_c1_atoms": len(pred_xyz),
        "reference_c1_atoms": len(ref_xyz),
        "command": command,
        "scope": "Local single-pair helper; not the official Kaggle aggregation/evaluator.",
        "usalign_stdout": completed.stdout,
    }


def _residue_keys(path: Path) -> list[tuple[str, str]]:
    keys = []
    seen = set()
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if (line[:6].strip() in {"ATOM", "HETATM"}
                    and line[12:16].strip() == "C1'"
                    and line[16:17] in {" ", "A"}):
                key = (line[22:26], line[26:27])
                if key not in seen:
                    seen.add(key)
                    keys.append(key)
    return keys
