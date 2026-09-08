from __future__ import annotations

from collections.abc import Mapping

import numpy as np


def _pearson_zero_if_constant(x: np.ndarray, y: np.ndarray) -> float:
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.ndim != 1 or yy.shape != xx.shape or len(xx) < 3:
        raise ValueError("coherence vectors must have equal 1D shape with >=3 values")
    if not np.isfinite(xx).all() or not np.isfinite(yy).all():
        raise ValueError("coherence vectors must be finite")
    dx = xx - float(xx.mean())
    dy = yy - float(yy.mean())
    den = float(np.sqrt(np.dot(dx, dx) * np.dot(dy, dy)))
    if den <= np.finfo(float).eps:
        return 0.0
    return float(np.dot(dx, dy) / den)


def edge_midpoint_neighbor_coherence(
    response: np.ndarray,
    midpoint: np.ndarray,
    *,
    k: int = 4,
) -> float:
    """Scale-free local spatial coherence of one species' edge response.

    The response is already a within-species rank-derived quantity in TTF.  For
    every edge, this statistic averages the responses of its nearest edge
    midpoints and correlates that neighbour average with the focal response.
    Boundary orientation and absolute field location are never supplied.  The
    statistic therefore measures within-species spatial organization, not
    cross-species alignment.
    """
    y = np.asarray(response, dtype=float)
    x = np.asarray(midpoint, dtype=float)
    if y.ndim != 1 or x.ndim != 2 or len(y) != len(x) or len(y) < 3:
        raise ValueError("response/midpoint shapes are incompatible")
    if not np.isfinite(y).all() or not np.isfinite(x).all():
        raise ValueError("response and midpoint must be finite")
    if k < 1:
        raise ValueError("k must be positive")
    use_k = min(int(k), len(y) - 1)
    delta = x[:, None, :] - x[None, :, :]
    distance2 = np.sum(delta * delta, axis=2)
    np.fill_diagonal(distance2, np.inf)
    order = np.argsort(distance2, axis=1, kind="stable")[:, :use_k]
    neighbour_mean = y[order].mean(axis=1)
    return _pearson_zero_if_constant(y, neighbour_mean)


def training_private_strength(
    responses: Mapping[str, np.ndarray],
    midpoints: Mapping[str, np.ndarray],
    train_species: tuple[str, ...] | list[str],
    *,
    k: int = 4,
) -> tuple[float, dict[str, float]]:
    """Equal-species mean nuisance-strength proxy computed from training only."""
    names = tuple(map(str, train_species))
    if not names:
        raise ValueError("at least one training species is required")
    per_species: dict[str, float] = {}
    for name in names:
        if name not in responses or name not in midpoints:
            raise ValueError(f"missing training strength input for {name}")
        per_species[name] = edge_midpoint_neighbor_coherence(
            responses[name], midpoints[name], k=k
        )
    return float(np.mean(list(per_species.values()))), per_species
