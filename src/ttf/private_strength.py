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


def edge_midpoint_neighbor_indices(
    midpoint: np.ndarray,
    *,
    k: int = 4,
) -> np.ndarray:
    """Precompute nearest edge-midpoint neighbours for a fixed species geometry."""
    x = np.asarray(midpoint, dtype=float)
    if x.ndim != 2 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("midpoint must be a finite m x d array with m>=3")
    if k < 1:
        raise ValueError("k must be positive")
    use_k = min(int(k), len(x) - 1)
    delta = x[:, None, :] - x[None, :, :]
    distance2 = np.sum(delta * delta, axis=2)
    np.fill_diagonal(distance2, np.inf)
    return np.argsort(distance2, axis=1, kind="stable")[:, :use_k]


def neighbor_coherence_from_indices(
    response: np.ndarray,
    neighbour_indices: np.ndarray,
) -> float:
    """Compute local response coherence on a precomputed neighbour design."""
    y = np.asarray(response, dtype=float)
    index = np.asarray(neighbour_indices, dtype=np.int64)
    if y.ndim != 1 or len(y) < 3 or not np.isfinite(y).all():
        raise ValueError("response must be finite 1D with >=3 values")
    if index.ndim != 2 or index.shape[0] != len(y) or index.shape[1] < 1:
        raise ValueError("neighbour_indices shape is incompatible with response")
    if np.any(index < 0) or np.any(index >= len(y)):
        raise ValueError("neighbour index out of range")
    neighbour_mean = y[index].mean(axis=1)
    return _pearson_zero_if_constant(y, neighbour_mean)


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
    index = edge_midpoint_neighbor_indices(midpoint, k=k)
    return neighbor_coherence_from_indices(response, index)


def training_private_strength_from_indices(
    responses: Mapping[str, np.ndarray],
    neighbour_indices: Mapping[str, np.ndarray],
    train_species: tuple[str, ...] | list[str],
) -> tuple[float, dict[str, float]]:
    """Exact fast path for equal-species training coherence on frozen geometry."""
    names = tuple(map(str, train_species))
    if not names:
        raise ValueError("at least one training species is required")
    per_species: dict[str, float] = {}
    for name in names:
        if name not in responses or name not in neighbour_indices:
            raise ValueError(f"missing training strength input for {name}")
        per_species[name] = neighbor_coherence_from_indices(
            responses[name], neighbour_indices[name]
        )
    return float(np.mean(list(per_species.values()))), per_species


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
    index = {
        name: edge_midpoint_neighbor_indices(midpoints[name], k=k)
        for name in names
    }
    return training_private_strength_from_indices(responses, index, names)
