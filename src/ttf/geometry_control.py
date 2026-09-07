from __future__ import annotations

import numpy as np

from .core import average_ranks, rank01


def _rank_center(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("rank-control inputs must be finite 1D arrays with >=3 values")
    ranks = average_ranks(x)
    return ranks - float(ranks.mean())


def _residualize_centered(y: np.ndarray, nuisance: np.ndarray) -> np.ndarray:
    yy = np.asarray(y, dtype=float)
    zz = np.asarray(nuisance, dtype=float)
    if yy.shape != zz.shape or yy.ndim != 1:
        raise ValueError("response and nuisance must be equal-length 1D arrays")
    yc = yy - float(yy.mean())
    zc = zz - float(zz.mean())
    den = float(np.dot(zc, zc))
    if den <= np.finfo(float).eps:
        return yc
    beta = float(np.dot(zc, yc) / den)
    return yc - beta * zc


def length_orthogonalized_turnover(
    turnover: np.ndarray,
    edge_length: np.ndarray,
) -> np.ndarray:
    """Remove the monotone first-order edge-length component within one species.

    TTF turnover is already a within-species rank score, but this function ranks
    it again defensively so the definition is invariant to any monotone response
    rescaling. Edge length is replaced by its within-species rank. The centered
    turnover rank is linearly projected off centered length rank.

    The residual is rescaled to the population SD of an equally spaced rank01
    vector of the same size. This keeps the marginal scale comparable across
    species, preserving the original equal-species field weighting while making
    the nuisance mean exactly zero. If turnover is perfectly explained by edge
    length, the returned vector is identically zero.
    """
    y = np.asarray(turnover, dtype=float)
    length = np.asarray(edge_length, dtype=float)
    if y.ndim != 1 or length.shape != y.shape or len(y) < 3:
        raise ValueError("turnover and edge_length must be equal 1D arrays with >=3 edges")
    if not np.isfinite(y).all() or not np.isfinite(length).all():
        raise ValueError("turnover and edge_length must be finite")
    if np.any(length < 0):
        raise ValueError("edge lengths cannot be negative")

    y_rank = rank01(y)
    length_rank = rank01(length)
    residual = _residualize_centered(y_rank, length_rank)
    residual_sd = float(np.std(residual, ddof=0))
    if residual_sd <= np.sqrt(np.finfo(float).eps):
        return np.zeros_like(residual)

    n = len(residual)
    target_sd = float(np.sqrt((n * n - 1.0) / (12.0 * n * n)))
    return residual * (target_sd / residual_sd)


def partial_spearman_rho(
    prediction: np.ndarray,
    turnover: np.ndarray,
    edge_length: np.ndarray,
) -> float:
    """Spearman association after linearly partialling edge-length rank.

    Prediction and turnover are independently converted to ranks, then each is
    residualized against the same within-species edge-length rank. The returned
    value is the Pearson correlation of those residual rank vectors, i.e. a
    first-order partial Spearman correlation controlling edge length.
    """
    x = np.asarray(prediction, dtype=float)
    y = np.asarray(turnover, dtype=float)
    z = np.asarray(edge_length, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or z.shape != x.shape:
        raise ValueError("prediction, turnover and edge_length must have equal 1D shape")
    if len(x) < 3:
        return float("nan")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isfinite(z).all():
        raise ValueError("partial Spearman inputs must be finite")

    rx = _rank_center(x)
    ry = _rank_center(y)
    rz = _rank_center(z)
    ex = _residualize_centered(rx, rz)
    ey = _residualize_centered(ry, rz)
    den = float(np.sqrt(np.dot(ex, ex) * np.dot(ey, ey)))
    if den <= np.finfo(float).eps:
        return 0.0
    return float(np.dot(ex, ey) / den)


def spearman_length_association(
    turnover: np.ndarray,
    edge_length: np.ndarray,
) -> float:
    """Within-species Spearman association used only as a geometry diagnostic."""
    y = _rank_center(np.asarray(turnover, dtype=float))
    z = _rank_center(np.asarray(edge_length, dtype=float))
    den = float(np.sqrt(np.dot(y, y) * np.dot(z, z)))
    if den <= np.finfo(float).eps:
        return 0.0
    return float(np.dot(y, z) / den)
