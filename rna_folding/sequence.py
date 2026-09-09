"""Validated sequence input and explicit, report-inspired routing rules.

These explicit rules are a readable starting point for an inference workflow, not a
claim that any one predictor is best for every sequence in its length bucket.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from numbers import Integral
from pathlib import Path
import re


@dataclass(frozen=True)
class RNASequence:
    """A target identifier and a nonempty, uppercase canonical RNA sequence."""

    target_id: str
    sequence: str

    def __post_init__(self) -> None:
        if not isinstance(self.target_id, str) or not re.fullmatch(
            r"[A-Za-z0-9_.-]+", self.target_id
        ):
            raise ValueError("target_id must contain only A-Z, a-z, 0-9, _, ., or -")
        if self.target_id in {".", ".."}:
            raise ValueError("target_id cannot be . or ..")
        if not isinstance(self.sequence, str) or not re.fullmatch(
            r"[ACGU]+", self.sequence
        ):
            raise ValueError("sequence must be nonempty uppercase canonical RNA (ACGU)")


def load_sequences(path: str | Path) -> list[RNASequence]:
    """Read a CSV with target_id and sequence columns; extra columns are allowed.

    Values are not silently uppercased, trimmed, or converted from DNA. A bad
    input should be corrected deliberately, with the source data preserved.
    """

    records: list[RNASequence] = []
    seen: set[str] = set()
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if fields is None or not {"target_id", "sequence"}.issubset(fields):
            raise ValueError("sequence CSV requires target_id and sequence columns")
        if len(fields) != len(set(fields)):
            raise ValueError("sequence CSV has duplicate column names")
        for row in reader:
            try:
                record = RNASequence(row["target_id"], row["sequence"])
            except ValueError as exc:
                raise ValueError(f"sequence CSV line {reader.line_num}: {exc}") from exc
            if record.target_id in seen:
                raise ValueError(f"duplicate target_id: {record.target_id}")
            seen.add(record.target_id)
            records.append(record)
    if not records:
        raise ValueError("sequence CSV contains no targets")
    return records


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def route_model(length: int) -> str:
    """Return the documented <=100, 101-480, and >480 length-bucket model."""

    length = _positive_integer(length, "length")
    if length <= 100:
        return "trrosettarna"
    if length <= 480:
        return "drfold2"
    return "protenix"


def segment_ranges(
    length: int, window: int = 300, overlap: int = 50
) -> list[tuple[int, int]]:
    """Return 0-based half-open windows covering every residue.

    The last window is shifted left when needed to retain a full window instead
    of a short tail. Its overlap may therefore exceed ``overlap``. A sequence
    shorter than ``window`` is represented by one shorter window. An overlap of
    zero is allowed for indexing, although rigid stitching requires at least
    three non-collinear shared points.
    """

    length = _positive_integer(length, "length")
    window = _positive_integer(window, "window")
    if isinstance(overlap, bool) or not isinstance(overlap, Integral):
        raise ValueError("overlap must be an integer")
    overlap = int(overlap)
    if overlap < 0 or overlap >= window:
        raise ValueError("overlap must satisfy 0 <= overlap < window")
    if length <= window:
        return [(0, length)]
    ranges = [(0, window)]
    stride = window - overlap
    while ranges[-1][1] < length:
        start = min(ranges[-1][0] + stride, length - window)
        ranges.append((start, start + window))
    return ranges
