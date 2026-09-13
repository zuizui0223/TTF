from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .cached_chunked_transfer import (
    CachedChunkedTransferGeometry,
    prepare_cached_chunked_transfer,
    score_cached_chunked_batch,
)
from .genetic_gate import GeneticTTFDesign
from .genetic_ibd import crossfit_ibd_residuals_prepared
from .genetic_simulate import GeneticSyntheticWorld
from .geometry_control import length_orthogonalized_turnover
from .private_strength import training_private_strength_from_indices


@dataclass(frozen=True)
class GeneticBatchScore:
    statistics: np.ndarray
    training_strengths: np.ndarray
    species_scores: Mapping[str, np.ndarray]
    n_eval_species: np.ndarray


def prepare_genetic_cached_transfer(
    design: GeneticTTFDesign,
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> CachedChunkedTransferGeometry:
    return prepare_cached_chunked_transfer(
        design.prepared,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )


def score_genetic_world_batch(
    design: GeneticTTFDesign,
    worlds: Sequence[GeneticSyntheticWorld],
    cached_transfer: CachedChunkedTransferGeometry,
) -> GeneticBatchScore:
    """Score multiple synthetic worlds together without changing the estimand."""
    batch = tuple(worlds)
    if not batch:
        raise ValueError("at least one world is required")
    width = len(batch)
    used = design.train_species + design.eval_species
    for index, world in enumerate(batch):
        missing = set(used) - set(world.genetic_distance)
        if missing:
            raise ValueError(f"world {index} is missing species: {sorted(missing)}")

    train_response = {
        name: np.empty((design.template_edges[name].n_edges, width), dtype=float)
        for name in design.train_species
    }
    eval_response = {
        name: np.empty((design.template_edges[name].n_edges, width), dtype=float)
        for name in design.eval_species
    }
    strengths = np.empty(width, dtype=float)

    for column, world in enumerate(batch):
        train_column: dict[str, np.ndarray] = {}
        for name in used:
            result = crossfit_ibd_residuals_prepared(
                world.genetic_distance[name],
                design.ibd_designs[name],
            )
            if name in design.train_species:
                response = length_orthogonalized_turnover(
                    result.residual_turnover,
                    design.template_edges[name].length,
                )
                train_response[name][:, column] = response
                train_column[name] = response
            else:
                eval_response[name][:, column] = result.residual_turnover
        strength, _ = training_private_strength_from_indices(
            train_column,
            design.strength_indices,
            design.train_species,
        )
        strengths[column] = float(strength)

    scored = score_cached_chunked_batch(
        cached_transfer,
        train_response,
        eval_response,
    )
    return GeneticBatchScore(
        statistics=np.asarray(scored.statistics, dtype=float).copy(),
        training_strengths=strengths,
        species_scores={
            name: np.asarray(scored.species_scores[name], dtype=float).copy()
            for name in design.eval_species
        },
        n_eval_species=np.asarray(scored.n_eval_species, dtype=np.int64).copy(),
    )


__all__ = [
    "GeneticBatchScore",
    "prepare_genetic_cached_transfer",
    "score_genetic_world_batch",
]
