from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import average_ranks


@dataclass(frozen=True)
class IndependentRankFactorization:
    """Exact expectation factorization for independent random rank vectors."""

    mean_left_direction: np.ndarray
    mean_right_direction: np.ndarray
    expected_spearman: float


def standardized_rank_direction(values: np.ndarray) -> np.ndarray:
    """Return the unit centered-rank direction used by Spearman correlation.

    Constant vectors map to the zero vector, matching ``spearman_rho`` in
    ``ttf.core``.  Otherwise this is exactly centered average ranks divided by
    their Euclidean norm.
    """

    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("values must be a finite 1D vector with length >= 3")
    ranks = average_ranks(x)
    centered = ranks - float(np.mean(ranks))
    norm = float(np.linalg.norm(centered))
    if norm <= np.finfo(float).eps:
        return np.zeros_like(centered)
    return centered / norm


def spearman_rank_inner_product(x: np.ndarray, y: np.ndarray) -> float:
    """Compute Spearman correlation as an inner product of rank directions."""

    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.ndim != 1 or yy.ndim != 1 or xx.shape != yy.shape:
        raise ValueError("x and y must be equal-length 1D vectors")
    if len(xx) < 3:
        raise ValueError("rank factorization requires at least three entries")
    return float(
        np.dot(
            standardized_rank_direction(xx),
            standardized_rank_direction(yy),
        )
    )


def _probability_weights(n: int, weights: np.ndarray | None) -> np.ndarray:
    if n < 1:
        raise ValueError("support must contain at least one draw")
    if weights is None:
        return np.full(n, 1.0 / n, dtype=float)
    w = np.asarray(weights, dtype=float)
    if w.shape != (n,) or not np.isfinite(w).all() or np.any(w < 0.0):
        raise ValueError("weights must be finite, non-negative, and match support")
    total = float(np.sum(w))
    if total <= 0.0:
        raise ValueError("weights must have positive total mass")
    return w / total


def independent_rank_factorization(
    left_support: np.ndarray,
    right_support: np.ndarray,
    *,
    left_weights: np.ndarray | None = None,
    right_weights: np.ndarray | None = None,
) -> IndependentRankFactorization:
    """Factor E[Spearman(X,Y)] for independent finite-support random vectors.

    Let ``u(v)`` be the standardized centered-rank direction.  Since
    ``Spearman(X,Y) = u(X)^T u(Y)``, conditional independence gives exactly

        E[Spearman(X,Y)] = E[u(X)]^T E[u(Y)].

    The supports here may represent arbitrary finite distributions, including
    tied or constant vectors.  This is an algebraic identity, not a model fit
    and not a correction to the TTF estimator.
    """

    left = np.asarray(left_support, dtype=float)
    right = np.asarray(right_support, dtype=float)
    if left.ndim != 2 or right.ndim != 2 or left.shape[1] != right.shape[1]:
        raise ValueError("supports must be draws x entries with common entry count")
    if left.shape[1] < 3 or not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("supports must be finite with at least three entries per draw")

    lw = _probability_weights(len(left), left_weights)
    rw = _probability_weights(len(right), right_weights)
    left_dirs = np.vstack([standardized_rank_direction(row) for row in left])
    right_dirs = np.vstack([standardized_rank_direction(row) for row in right])
    mean_left = np.sum(left_dirs * lw[:, None], axis=0)
    mean_right = np.sum(right_dirs * rw[:, None], axis=0)
    return IndependentRankFactorization(
        mean_left_direction=mean_left,
        mean_right_direction=mean_right,
        expected_spearman=float(np.dot(mean_left, mean_right)),
    )


__all__ = [
    "IndependentRankFactorization",
    "standardized_rank_direction",
    "spearman_rank_inner_product",
    "independent_rank_factorization",
]
