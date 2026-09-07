from __future__ import annotations

from dataclasses import replace

import numpy as np

from .core import SpeciesEdges, rank01


def residualize_edge_length(
    edges: SpeciesEdges,
    *,
    degree: int = 1,
) -> SpeciesEdges:
    """Remove species-local edge-length dependence before TTF transfer scoring.

    The nuisance response is the already rank-standardized turnover and the
    predictor is within-species edge-length rank.  Each edge is predicted from
    a polynomial fitted to *all other edges* (leave-one-edge-out cross-fitting),
    so an edge never trains its own nuisance expectation.  Residuals are then
    re-ranked within species to retain TTF's scale invariance.

    This is deliberately a geographic-distance nuisance layer, not a boundary
    model: it has no access to absolute coordinates, other species, known
    barriers, or empirical benchmark outcomes.
    """
    if degree < 1:
        raise ValueError("degree must be >= 1")
    n = edges.n_edges
    if n <= degree + 2:
        raise ValueError("not enough edges for leave-one-out IBD residualization")

    x = rank01(np.asarray(edges.length, dtype=float)) - 0.5
    y = np.asarray(edges.turnover, dtype=float)
    if y.shape != (n,) or not np.isfinite(y).all():
        raise ValueError("edge turnover must be a finite vector")

    design = np.column_stack([np.ones(n)] + [x ** power for power in range(1, degree + 1)])
    prediction = np.empty(n, dtype=float)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        keep[:] = True
        keep[i] = False
        beta, *_ = np.linalg.lstsq(design[keep], y[keep], rcond=None)
        prediction[i] = float(design[i] @ beta)

    residual = y - prediction
    if not np.isfinite(residual).all():
        raise RuntimeError("IBD nuisance fit produced non-finite residuals")
    return replace(edges, turnover=rank01(residual))


def residualize_edge_sets(
    edge_sets: list[SpeciesEdges] | tuple[SpeciesEdges, ...],
    *,
    degree: int = 1,
) -> tuple[SpeciesEdges, ...]:
    return tuple(residualize_edge_length(edges, degree=degree) for edges in edge_sets)
