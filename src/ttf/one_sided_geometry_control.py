from __future__ import annotations

from typing import Sequence

import numpy as np

from .core import average_ranks


def centered_rank(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("values must be finite 1D with at least three entries")
    r = average_ranks(x)
    return r - float(r.mean())


def rank_residual(values: np.ndarray, nuisances: Sequence[np.ndarray]) -> np.ndarray:
    """Residualize a rank vector against a frozen nuisance-rank design."""
    r = centered_rank(values)
    if not nuisances:
        return r
    columns = []
    for nuisance in nuisances:
        z = np.asarray(nuisance, dtype=float)
        if z.shape != r.shape or not np.isfinite(z).all():
            raise ValueError("nuisance shape/finite check failed")
        columns.append(centered_rank(z))
    design = np.column_stack(columns)
    keep = np.linalg.norm(design, axis=0) > np.sqrt(np.finfo(float).eps)
    if np.any(keep):
        q = design[:, keep]
        coef, *_ = np.linalg.lstsq(q, r, rcond=None)
        r = r - q @ coef
    r -= float(r.mean())
    return r


def residual_rank_correlation(
    prediction: np.ndarray,
    target: np.ndarray,
    *,
    prediction_nuisances: Sequence[np.ndarray] = (),
    target_nuisances: Sequence[np.ndarray] = (),
) -> float:
    """Correlate prediction/target ranks after side-specific nuisance projection."""
    x = np.asarray(prediction, dtype=float)
    y = np.asarray(target, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 3:
        return float("nan")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("prediction and target must be finite")
    rx = rank_residual(x, prediction_nuisances)
    ry = rank_residual(y, target_nuisances)
    den = float(np.sqrt(np.dot(rx, rx) * np.dot(ry, ry)))
    if den <= np.finfo(float).eps:
        return 0.0
    return float(np.dot(rx, ry) / den)
