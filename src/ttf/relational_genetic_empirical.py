from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .core import SpeciesEdges, spearman_rho
from .geometry_control import length_orthogonalized_turnover
from .phylogatr_compact_ibd import (
    crossfit_ibd_residuals_compact,
    prepare_compact_crossfit_ibd_design,
)
from .transfer import fit_boundary_model


@dataclass(frozen=True)
class RelationalSpeciesResponse:
    species: str
    train_turnover: np.ndarray
    eval_turnover: np.ndarray


def post_ibd_responses(
    edges: SpeciesEdges,
    genetic_distance: np.ndarray,
    *,
    min_training_edges: int = 5,
) -> RelationalSpeciesResponse:
    distance = np.asarray(genetic_distance, dtype=float)
    if distance.ndim != 1 or len(distance) != edges.n_edges:
        raise ValueError("genetic distance length must match frozen species edges")
    design = prepare_compact_crossfit_ibd_design(
        edges.length,
        edges.nodes,
        min_training_edges=int(min_training_edges),
    )
    residual = crossfit_ibd_residuals_compact(distance, design).residual_turnover
    train = length_orthogonalized_turnover(residual, edges.length)
    return RelationalSpeciesResponse(
        species=edges.species,
        train_turnover=np.asarray(train, dtype=float),
        eval_turnover=np.asarray(residual, dtype=float),
    )


def source_only_transfer_scores(
    edge_map: Mapping[str, SpeciesEdges],
    response_map: Mapping[str, RelationalSpeciesResponse],
    pairs: list[tuple[str, str]],
    *,
    bandwidth_km: float = 500.0,
    prior_strength: float = 0.25,
    prior_mean: float = 0.0,
    segment_points: int = 5,
) -> np.ndarray:
    if not pairs:
        raise ValueError("at least one source-target pair is required")
    sources = sorted({source for source, _ in pairs})
    targets = {target for _, target in pairs}
    if set(sources) & targets:
        raise ValueError("relational source and target roles must remain species-disjoint")
    missing = (set(sources) | targets) - set(edge_map)
    if missing:
        raise ValueError(f"missing frozen edge geometry for species: {sorted(missing)}")
    missing_response = (set(sources) | targets) - set(response_map)
    if missing_response:
        raise ValueError(f"missing post-IBD response for species: {sorted(missing_response)}")

    models = {}
    for source in sources:
        edges = edge_map[source]
        response = response_map[source]
        train_edges = SpeciesEdges(
            species=source,
            nodes=edges.nodes,
            start=edges.start,
            end=edges.end,
            midpoint=edges.midpoint,
            length=edges.length,
            turnover=response.train_turnover,
        )
        models[source] = fit_boundary_model(
            [train_edges],
            bandwidth=float(bandwidth_km),
            prior_strength=float(prior_strength),
            prior_mean=float(prior_mean),
        )

    out = np.empty(len(pairs), dtype=float)
    for index, (source, target) in enumerate(pairs):
        target_edges = edge_map[target]
        predicted = models[source].predict_segments(
            target_edges.start,
            target_edges.end,
            n_points=int(segment_points),
        )
        score = spearman_rho(predicted, response_map[target].eval_turnover)
        if not np.isfinite(score):
            raise RuntimeError(f"non-finite source-target transfer score for {source} -> {target}")
        out[index] = float(score)
    return out


__all__ = [
    "RelationalSpeciesResponse",
    "post_ibd_responses",
    "source_only_transfer_scores",
]
