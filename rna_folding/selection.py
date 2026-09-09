"""Transparent model-diverse selection with explicitly scoped energy ranking."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np

from .geometry import validate_coordinates
from .sequence import RNASequence, route_model


@dataclass(frozen=True)
class Candidate:
    """One externally supplied structure; energy is optional metadata.

    ``energy_group`` identifies a genuinely comparable scoring protocol within
    one model. The caller is responsible for assigning it correctly. This
    package does not calculate energies or calibrate scores between models.
    """

    candidate_id: str
    target_id: str
    model: str
    coordinates: np.ndarray
    energy: float | None = None
    energy_group: str | None = None


def _order_within_groups(candidates: list[Candidate]) -> list[Candidate]:
    """Sort finite scored entries only within their existing comparable slots."""

    ordered = list(candidates)
    group_positions: dict[str, list[int]] = {}
    for index, candidate in enumerate(candidates):
        if candidate.energy is not None and candidate.energy_group is not None:
            group_positions.setdefault(candidate.energy_group, []).append(index)
    for positions in group_positions.values():
        ranked = sorted((candidates[index] for index in positions), key=lambda item: item.energy)
        for position, candidate in zip(positions, ranked):
            ordered[position] = candidate
    return ordered


def select_candidates(
    candidates: list[Candidate], sequence: RNASequence, count: int = 5
) -> list[Candidate]:
    """Select exactly ``count`` candidates in deterministic model round robin.

    The routed model is visited first, then the other documented model names,
    then additional model names alphabetically. Within each model, input order
    is retained except that scored entries sharing a nonempty energy_group are
    sorted ascending among their own slots. Unscored entries, different groups,
    and missing groups do not compete by raw energy. Ties retain input order.

    Every candidate is validated, including candidates that would not make the
    final selection. IDs must be distinct; structures are not deduplicated by
    shape, and model diversity does not guarantee conformational diversity.
    """

    if not isinstance(sequence, RNASequence):
        raise ValueError("sequence must be an RNASequence")
    if isinstance(count, bool) or not isinstance(count, Integral) or count <= 0:
        raise ValueError("count must be a positive integer")
    if len(candidates) < count:
        raise ValueError(f"need at least {count} candidates, received {len(candidates)}")
    by_model: dict[str, list[Candidate]] = {}
    seen: set[str] = set()
    for candidate in candidates:
        if not isinstance(candidate, Candidate):
            raise ValueError("all entries must be Candidate instances")
        if not isinstance(candidate.candidate_id, str) or not candidate.candidate_id.strip():
            raise ValueError("candidate_id must be a nonempty string")
        if candidate.candidate_id in seen:
            raise ValueError(f"duplicate candidate_id: {candidate.candidate_id}")
        seen.add(candidate.candidate_id)
        if candidate.target_id != sequence.target_id:
            raise ValueError(f"candidate {candidate.candidate_id} belongs to a different target")
        if not isinstance(candidate.model, str) or not candidate.model.strip():
            raise ValueError("model must be a nonempty string")
        if candidate.energy is not None:
            if isinstance(candidate.energy, bool) or not isinstance(candidate.energy, Real):
                raise ValueError("energy must be a finite number or None")
            if not np.isfinite(candidate.energy):
                raise ValueError("energy must be finite")
        if candidate.energy_group is not None and (
            not isinstance(candidate.energy_group, str) or not candidate.energy_group.strip()
        ):
            raise ValueError("energy_group must be a nonempty string or None")
        validate_coordinates(candidate.coordinates, length=len(sequence.sequence))
        by_model.setdefault(candidate.model, []).append(candidate)

    preferred = route_model(len(sequence.sequence))
    known = ["trrosettarna", "drfold2", "protenix"]
    priorities = [preferred] + [model for model in known if model != preferred]
    priorities += sorted(set(by_model) - set(known))
    model_order = [model for model in priorities if model in by_model]
    queues = {model: _order_within_groups(by_model[model]) for model in model_order}
    selected: list[Candidate] = []
    round_index = 0
    while len(selected) < count:
        for model in model_order:
            if round_index < len(queues[model]):
                selected.append(queues[model][round_index])
                if len(selected) == count:
                    return selected
        round_index += 1
    return selected
