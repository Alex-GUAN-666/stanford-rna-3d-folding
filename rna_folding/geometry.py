"""Coordinate validation and rigid, index-matched segment assembly.

This module handles geometry only. It does not perform folding, relaxation,
energy calculation, sequence alignment, or biological structure validation.
"""

from __future__ import annotations

from numbers import Integral

import numpy as np


def validate_coordinates(xyz: object, length: int | None = None) -> np.ndarray:
    """Return a float64 copy of finite (L, 3) coordinates.

    Coordinates with absolute components >= 1e6 are rejected as likely missing
    value sentinels. This sanity bound is not a physical plausibility check.
    """

    try:
        if np.iscomplexobj(xyz):
            raise ValueError("coordinates must be real-valued, not complex")
        coordinates = np.array(xyz, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise ValueError("coordinates must be numeric with shape (L, 3)") from exc
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or len(coordinates) == 0:
        raise ValueError("coordinates must have nonempty shape (L, 3)")
    if length is not None:
        if isinstance(length, bool) or not isinstance(length, Integral) or length <= 0:
            raise ValueError("length must be a positive integer")
        if len(coordinates) != length:
            raise ValueError(f"coordinate length {len(coordinates)} != sequence length {length}")
    if not np.isfinite(coordinates).all():
        raise ValueError("coordinates must be finite")
    if np.any(np.abs(coordinates) >= 1e6):
        raise ValueError("coordinate magnitude >= 1e6 suggests a missing-value sentinel")
    return coordinates


def _rigid_transform(moving: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if moving.shape != reference.shape:
        raise ValueError("moving and reference coordinates must have the same shape")
    if len(moving) < 3:
        raise ValueError("rigid alignment requires at least 3 overlapping points")
    moving_center = moving.mean(axis=0)
    reference_center = reference.mean(axis=0)
    centered_moving = moving - moving_center
    centered_reference = reference - reference_center
    if np.linalg.matrix_rank(centered_moving) < 2 or np.linalg.matrix_rank(centered_reference) < 2:
        raise ValueError("rigid alignment requires non-collinear overlapping points")
    u, _, vt = np.linalg.svd(centered_moving.T @ centered_reference)
    correction = np.eye(3)
    correction[-1, -1] = 1.0 if np.linalg.det(u @ vt) >= 0 else -1.0
    rotation = u @ correction @ vt
    translation = reference_center - moving_center @ rotation
    return rotation, translation


def kabsch_align(moving: object, reference: object) -> np.ndarray:
    """Align corresponding points using a proper rotation (no reflection).

    At least three non-collinear points are required in each input. Residue
    correspondence is supplied by the caller; this is not a structural aligner
    such as US-align and does not optimize a TM-score objective.
    """

    moving_array = validate_coordinates(moving)
    reference_array = validate_coordinates(reference)
    rotation, translation = _rigid_transform(moving_array, reference_array)
    return validate_coordinates(moving_array @ rotation + translation)


def stitch_segments(segments: list[tuple[int, np.ndarray]], length: int) -> np.ndarray:
    """Align overlapping segments and average their aligned shared residues.

    Each pair is ``(start, coordinates)`` with a 0-based start. Segments are
    sorted by start, aligned to the assembled overlap, then blended with equal
    per-segment weights. The first segment defines the global coordinate frame.
    Missing coverage, invalid bounds, or degenerate overlaps raise ValueError.
    Averaging can distort local geometry; no bond or steric constraints apply.
    """

    if isinstance(length, bool) or not isinstance(length, Integral) or length <= 0:
        raise ValueError("length must be a positive integer")
    if not segments:
        raise ValueError("at least one segment is required")
    checked: list[tuple[int, np.ndarray]] = []
    for item in segments:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise ValueError("each segment must be a (start, coordinates) pair")
        start, xyz = item
        if isinstance(start, bool) or not isinstance(start, Integral) or start < 0:
            raise ValueError("segment start must be a nonnegative integer")
        coordinates = validate_coordinates(xyz)
        if start + len(coordinates) > length:
            raise ValueError("segment extends beyond the target length")
        checked.append((int(start), coordinates))
    checked.sort(key=lambda pair: pair[0])
    if checked[0][0] != 0:
        raise ValueError("segment coverage must start at residue 0")

    assembled = np.zeros((length, 3), dtype=np.float64)
    counts = np.zeros(length, dtype=np.int64)
    covered_end = 0
    for index, (start, coordinates) in enumerate(checked):
        end = start + len(coordinates)
        if start > covered_end:
            raise ValueError(f"gap in segment coverage before residue {start}")
        shared = counts[start:end] > 0
        aligned = coordinates
        if index:
            if shared.sum() < 3:
                raise ValueError("segment alignment requires at least 3 overlapping points")
            rotation, translation = _rigid_transform(
                coordinates[shared], assembled[start:end][shared]
            )
            aligned = coordinates @ rotation + translation
        old_counts = counts[start:end].copy()
        assembled[start:end] = (
            assembled[start:end] * old_counts[:, None] + aligned
        ) / (old_counts[:, None] + 1)
        counts[start:end] += 1
        covered_end = max(covered_end, end)
    if covered_end != length or np.any(counts == 0):
        raise ValueError("segments do not cover the complete target")
    return validate_coordinates(assembled, length=length)
