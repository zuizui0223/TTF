from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .cached_chunked_transfer import CachedChunkedTransferGeometry
from .genetic_empirical_score import EmpiricalSelfScore
from .genetic_self_detectability import GeneticSelfDetectabilityDesign
from .phylogatr_compact_execution import (
    PhylogatrCompactTTFDesign,
    score_phylogatr_compact_world_batch,
)
from .phylogatr_compact_self_detectability import (
    score_phylogatr_compact_self_world_batch,
)


@dataclass(frozen=True)
class CompactEmpiricalPrimaryScore:
    statistic: float
    training_strength: float
    species_scores: Mapping[str, float]


@dataclass(frozen=True)
class _EmpiricalWorld:
    genetic_distance: Mapping[str, np.ndarray]


def _validated_mapping(
    design: PhylogatrCompactTTFDesign,
    genetic_distance: Mapping[str, np.ndarray],
) -> dict[str, np.ndarray]:
    used = design.train_species + design.eval_species
    keys = set(map(str, genetic_distance.keys()))
    missing = set(used) - keys
    extra = keys - set(used)
    if missing:
        raise ValueError(f"genetic distance mapping is missing species: {sorted(missing)}")
    if extra:
        raise ValueError(f"genetic distance mapping has unexpected species: {sorted(extra)}")

    out: dict[str, np.ndarray] = {}
    for name in used:
        values = np.asarray(genetic_distance[name], dtype=float)
        expected = design.template_edges[name].n_edges
        if values.ndim != 1 or len(values) != expected:
            raise ValueError(
                f"genetic distance vector length drift for {name}: "
                f"{values.shape!r}, expected {expected}"
            )
        if not np.isfinite(values).all() or np.any(values < 0.0):
            raise ValueError(
                f"genetic distances must be finite and non-negative for {name}"
            )
        out[name] = values
    return out


def score_phylogatr_compact_distance_mapping(
    design: PhylogatrCompactTTFDesign,
    genetic_distance: Mapping[str, np.ndarray],
    cached_transfer: CachedChunkedTransferGeometry,
) -> CompactEmpiricalPrimaryScore:
    """Score the one frozen empirical mapping via the qualified compact batch path."""
    validated = _validated_mapping(design, genetic_distance)
    scored = score_phylogatr_compact_world_batch(
        design,
        [_EmpiricalWorld(validated)],
        cached_transfer,
    )
    if int(scored.n_eval_species[0]) != len(design.eval_species):
        raise RuntimeError("non-finite empirical held-out species count")
    return CompactEmpiricalPrimaryScore(
        statistic=float(scored.statistics[0]),
        training_strength=float(scored.training_strengths[0]),
        species_scores={
            name: float(scored.species_scores[name][0])
            for name in design.eval_species
        },
    )


def score_phylogatr_compact_self_mapping(
    design: PhylogatrCompactTTFDesign,
    self_design: GeneticSelfDetectabilityDesign,
    genetic_distance: Mapping[str, np.ndarray],
) -> EmpiricalSelfScore:
    """Score the frozen empirical self diagnostic via the qualified compact path."""
    validated = _validated_mapping(design, genetic_distance)
    scored = score_phylogatr_compact_self_world_batch(
        design,
        self_design,
        [_EmpiricalWorld(validated)],
    )
    if int(scored.n_species[0]) != len(design.eval_species):
        raise RuntimeError("non-finite empirical self species count")
    return EmpiricalSelfScore(
        statistic=float(scored.statistics[0]),
        species_scores={
            name: float(scored.species_scores[name][0])
            for name in self_design.species
        },
    )


__all__ = [
    "CompactEmpiricalPrimaryScore",
    "score_phylogatr_compact_distance_mapping",
    "score_phylogatr_compact_self_mapping",
]
