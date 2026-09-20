from __future__ import annotations

import hashlib
from typing import Mapping

import numpy as np

from .private_null_inference import envelope_upper_pvalues


REFERENCE_TAG = "lepidoptera-trait-gradient-v02-reference"
EVALUATION_TAG = "lepidoptera-trait-gradient-v02-evaluation"
ALLOWED_CELLS = ("private", "geometry_confounded_trap", "trait_gradient_positive")


def frozen_v02_seed(master_seed: int, tag: str, cell: str, replicate: int) -> int:
    if tag not in {REFERENCE_TAG, EVALUATION_TAG}:
        raise ValueError("unauthorized v0.2 seed namespace")
    if cell not in ALLOWED_CELLS:
        raise ValueError("unknown v0.2 cell")
    if replicate < 0:
        raise ValueError("replicate must be non-negative")
    payload = f"{int(master_seed)}|{tag}|{cell}|{int(replicate)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def target_fe_geometry_residual(
    trait_similarity: np.ndarray,
    geometry_covariates: np.ndarray,
    target_index: np.ndarray,
) -> np.ndarray:
    y = np.asarray(trait_similarity, dtype=float)
    x = np.asarray(geometry_covariates, dtype=float)
    target = np.asarray(target_index, dtype=np.int64)
    if y.ndim != 1 or x.ndim != 2 or len(y) != len(x) or len(target) != len(y):
        raise ValueError("pair arrays drift")
    if not np.isfinite(y).all() or not np.isfinite(x).all():
        raise ValueError("non-finite pair design")
    yc = y.copy()
    xc = x.copy()
    for label in np.unique(target):
        idx = np.flatnonzero(target == label)
        yc[idx] -= float(np.mean(yc[idx]))
        xc[idx] -= np.mean(xc[idx], axis=0)
    scale = np.std(xc, axis=0, ddof=0)
    scale[scale <= np.sqrt(np.finfo(float).eps)] = 1.0
    z = xc / scale
    beta = np.linalg.lstsq(z, yc, rcond=None)[0]
    return yc - z @ beta


def target_equal_mean_correlation(
    pair_transfer: np.ndarray,
    residual_trait_similarity: np.ndarray,
    target_index: np.ndarray,
) -> tuple[float, dict[int, float]]:
    y = np.asarray(pair_transfer, dtype=float)
    x = np.asarray(residual_trait_similarity, dtype=float)
    target = np.asarray(target_index, dtype=np.int64)
    if y.ndim != 1 or x.ndim != 1 or len(y) != len(x) or len(target) != len(y):
        raise ValueError("pair arrays drift")
    per_target: dict[int, float] = {}
    for label in np.unique(target):
        idx = np.flatnonzero(target == label)
        xx = x[idx] - float(np.mean(x[idx]))
        yy = y[idx] - float(np.mean(y[idx]))
        denom = float(np.sqrt(np.dot(xx, xx) * np.dot(yy, yy)))
        if len(idx) >= 3 and denom > np.finfo(float).tiny:
            per_target[int(label)] = float(np.dot(xx, yy) / denom)
    if len(per_target) < 6:
        raise ValueError("fewer than six target correlations are estimable")
    return float(np.mean(list(per_target.values()))), per_target


def wilson_interval(rejections: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials < 1 or rejections < 0 or rejections > trials:
        raise ValueError("invalid binomial counts")
    p = float(rejections) / float(trials)
    denom = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denom
    half = z * np.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denom
    return float(center - half), float(center + half)


def envelope_pvalues(
    statistics: np.ndarray,
    private_reference: np.ndarray,
    geometry_reference: np.ndarray,
) -> np.ndarray:
    p, _ = envelope_upper_pvalues(
        np.asarray(statistics, dtype=float),
        {
            "private": np.asarray(private_reference, dtype=float),
            "geometry_confounded_trap": np.asarray(geometry_reference, dtype=float),
        },
    )
    return np.asarray(p, dtype=float)


__all__ = [
    "ALLOWED_CELLS",
    "EVALUATION_TAG",
    "REFERENCE_TAG",
    "envelope_pvalues",
    "frozen_v02_seed",
    "target_equal_mean_correlation",
    "target_fe_geometry_residual",
    "wilson_interval",
]
