from __future__ import annotations

from typing import Sequence

import numpy as np

from .core import average_ranks, rank01
from .geometry import SpeciesGeometry


def pooled_affine_frame(
    geometries: Sequence[SpeciesGeometry],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the exact pooled affine frame used by the Gate-I boundary simulator."""
    data = tuple(geometries)
    if not data:
        raise ValueError("at least one geometry is required")
    pooled = np.vstack([item.coordinates for item in data])
    center = np.median(pooled, axis=0)
    raw_scale = np.std(pooled, axis=0)
    active = raw_scale > np.sqrt(np.finfo(float).eps)
    if not np.any(active):
        raise ValueError("sampling geometry has no varying coordinate dimension")
    scale = raw_scale.copy()
    scale[~active] = 1.0
    return center, scale, active


def random_private_transition_propensity(
    coordinates: np.ndarray,
    edge_nodes: np.ndarray,
    *,
    center: np.ndarray,
    scale: np.ndarray,
    active: np.ndarray,
    transition_width: float = 0.2,
    n_directions: int = 2048,
    seed: int = 20260908,
) -> np.ndarray:
    """Geometry-only expected edge contrast under random private median hyperplanes.

    This is a diagnostic nuisance, not an empirical trait model. It integrates the
    latent absolute tanh contrast over uniformly random hyperplane normals while
    forcing every hyperplane through that species' median projection, exactly as
    the current fixed-geometry Gate-I private-boundary generator does. Amplitude
    is omitted because it is a positive scalar and therefore irrelevant to the
    rank nuisance in the noise-free latent contrast.
    """
    x = np.asarray(coordinates, dtype=float)
    nodes = np.asarray(edge_nodes, dtype=np.int64)
    c = np.asarray(center, dtype=float)
    s = np.asarray(scale, dtype=float)
    use = np.asarray(active, dtype=bool)
    if x.ndim != 2 or nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ValueError("coordinates must be nxd and edge_nodes must be mx2")
    if c.shape != (x.shape[1],) or s.shape != c.shape or use.shape != c.shape:
        raise ValueError("affine-frame dimensionality mismatch")
    if transition_width <= 0 or n_directions < 32:
        raise ValueError("transition_width must be positive and n_directions >= 32")
    if np.any(nodes < 0) or np.any(nodes >= len(x)):
        raise ValueError("edge node index out of range")

    z = (x - c) / s
    rng = np.random.default_rng(int(seed))
    draws = rng.normal(size=(int(n_directions), int(np.count_nonzero(use))))
    norms = np.linalg.norm(draws, axis=1)
    bad = norms <= np.finfo(float).tiny
    while np.any(bad):
        draws[bad] = rng.normal(size=(int(np.count_nonzero(bad)), draws.shape[1]))
        norms = np.linalg.norm(draws, axis=1)
        bad = norms <= np.finfo(float).tiny
    draws /= norms[:, None]
    normals = np.zeros((int(n_directions), x.shape[1]), dtype=float)
    normals[:, use] = draws

    projection = z @ normals.T
    offsets = np.median(projection, axis=0)
    latent = np.tanh((projection - offsets[None, :]) / float(transition_width))
    contrast = np.abs(latent[nodes[:, 0], :] - latent[nodes[:, 1], :])
    propensity = contrast.mean(axis=1)
    if not np.isfinite(propensity).all():
        raise RuntimeError("private-boundary propensity became non-finite")
    return propensity


def midpoint_radial_distance(
    coordinates: np.ndarray,
    edge_nodes: np.ndarray,
    *,
    center: np.ndarray,
    scale: np.ndarray,
) -> np.ndarray:
    """Species-centered standardized radial distance of each edge midpoint."""
    x = (np.asarray(coordinates, dtype=float) - np.asarray(center, dtype=float)) / np.asarray(scale, dtype=float)
    nodes = np.asarray(edge_nodes, dtype=np.int64)
    midpoint = 0.5 * (x[nodes[:, 0]] + x[nodes[:, 1]])
    species_center = np.median(x, axis=0)
    return np.linalg.norm(midpoint - species_center[None, :], axis=1)


def _centered_ranks(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("rank inputs must be finite 1D vectors with >=3 values")
    r = average_ranks(x)
    return r - float(r.mean())


def rank_residualized_turnover(
    turnover: np.ndarray,
    nuisances: Sequence[np.ndarray],
) -> np.ndarray:
    """Residualize turnover rank against a frozen geometry-only nuisance basis."""
    y = np.asarray(turnover, dtype=float)
    if y.ndim != 1 or len(y) < 3 or not np.isfinite(y).all():
        raise ValueError("turnover must be a finite 1D vector with >=3 edges")
    if not nuisances:
        raise ValueError("at least one nuisance is required")
    yc = rank01(y) - 0.5
    columns = []
    for nuisance in nuisances:
        z = np.asarray(nuisance, dtype=float)
        if z.shape != y.shape or not np.isfinite(z).all():
            raise ValueError("nuisance shape/finite check failed")
        columns.append(_centered_ranks(z))
    design = np.column_stack(columns)
    keep = np.linalg.norm(design, axis=0) > np.sqrt(np.finfo(float).eps)
    if not np.any(keep):
        residual = yc
    else:
        coef, *_ = np.linalg.lstsq(design[:, keep], yc, rcond=None)
        residual = yc - design[:, keep] @ coef
    residual -= float(residual.mean())
    sd = float(np.std(residual, ddof=0))
    if sd <= np.sqrt(np.finfo(float).eps):
        return np.zeros_like(residual)
    n = len(residual)
    target_sd = float(np.sqrt((n * n - 1.0) / (12.0 * n * n)))
    return residual * (target_sd / sd)


def partial_spearman_controls(
    prediction: np.ndarray,
    turnover: np.ndarray,
    nuisances: Sequence[np.ndarray],
) -> float:
    """Spearman correlation after projecting both rank vectors off nuisances."""
    x = np.asarray(prediction, dtype=float)
    y = np.asarray(turnover, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 3:
        return float("nan")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("prediction and turnover must be finite")
    rx = _centered_ranks(x)
    ry = _centered_ranks(y)
    columns = []
    for nuisance in nuisances:
        z = np.asarray(nuisance, dtype=float)
        if z.shape != x.shape or not np.isfinite(z).all():
            raise ValueError("nuisance shape/finite check failed")
        columns.append(_centered_ranks(z))
    if columns:
        design = np.column_stack(columns)
        keep = np.linalg.norm(design, axis=0) > np.sqrt(np.finfo(float).eps)
        if np.any(keep):
            q = design[:, keep]
            bx, *_ = np.linalg.lstsq(q, rx, rcond=None)
            by, *_ = np.linalg.lstsq(q, ry, rcond=None)
            rx = rx - q @ bx
            ry = ry - q @ by
    den = float(np.sqrt(np.dot(rx, rx) * np.dot(ry, ry)))
    if den <= np.finfo(float).eps:
        return 0.0
    return float(np.dot(rx, ry) / den)
