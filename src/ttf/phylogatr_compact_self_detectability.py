from __future__ import annotations

from typing import Sequence

import numpy as np

from .core import spearman_rho
from .genetic_self_detectability import (
    GeneticSelfBatchScore,
    GeneticSelfDetectabilityDesign,
)
from .genetic_simulate import GeneticSyntheticWorld
from .phylogatr_compact_execution import PhylogatrCompactTTFDesign
from .phylogatr_compact_ibd import crossfit_ibd_residuals_compact


def score_phylogatr_compact_self_world_batch(
    design: PhylogatrCompactTTFDesign,
    self_design: GeneticSelfDetectabilityDesign,
    worlds: Sequence[GeneticSyntheticWorld],
) -> GeneticSelfBatchScore:
    """Score frozen self-detectability with exact compact endpoint-safe IBD.

    This mirrors ``score_genetic_self_world_batch`` exactly.  Only the storage and
    execution of endpoint-safe rank-IBD residualization changes; the self-field
    projection, prior, Spearman statistic, species panel, and synthetic worlds are
    unchanged.
    """
    batch = tuple(worlds)
    if not batch:
        raise ValueError("at least one world is required")
    width = len(batch)
    species_scores = {
        name: np.empty(width, dtype=float) for name in self_design.species
    }

    for name in self_design.species:
        response = np.empty((design.template_edges[name].n_edges, width), dtype=float)
        for column, world in enumerate(batch):
            if name not in world.genetic_distance:
                raise ValueError(f"world {column} is missing species {name}")
            result = crossfit_ibd_residuals_compact(
                world.genetic_distance[name],
                design.ibd_designs[name],
            )
            response[:, column] = result.residual_turnover
        prepared = self_design.by_species[name]
        predicted = prepared.prior_offset[:, None] + prepared.projection @ response
        species_scores[name][:] = np.asarray(
            [spearman_rho(predicted[:, j], response[:, j]) for j in range(width)],
            dtype=float,
        )

    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    for values in species_scores.values():
        finite = np.isfinite(values)
        sums[finite] += values[finite]
        counts[finite] += 1
    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite within-species detectability score")
    return GeneticSelfBatchScore(
        statistics=sums / counts,
        species_scores={name: values.copy() for name, values in species_scores.items()},
        n_species=counts,
    )


__all__ = ["score_phylogatr_compact_self_world_batch"]
