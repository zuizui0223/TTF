from __future__ import annotations

import numpy as np


def standardize_from_reference(
    reference: np.ndarray, *arrays: np.ndarray
) -> tuple[np.ndarray, ...]:
    """Standardize climate arrays using a fixed response-blind reference cloud."""
    ref = np.asarray(reference, dtype=float)
    if ref.ndim != 2 or ref.shape[0] < 2:
        raise ValueError("reference must be a 2D array with at least two rows")
    if np.any(~np.isfinite(ref)):
        raise ValueError("reference must be finite")
    mean = np.mean(ref, axis=0)
    sd = np.std(ref, axis=0, ddof=0)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("reference climate dimensions must have positive finite SD")

    out: list[np.ndarray] = []
    for array in arrays:
        x = np.asarray(array, dtype=float)
        if x.shape[-1] != ref.shape[1]:
            raise ValueError("climate dimension mismatch")
        if np.any(~np.isfinite(x)):
            raise ValueError("climate arrays must be finite")
        out.append((x - mean) / sd)
    return tuple(out)


def nearest_cloud_distance(points: np.ndarray, cloud: np.ndarray) -> np.ndarray:
    """Euclidean distance from each point to its nearest climate analogue."""
    p = np.asarray(points, dtype=float)
    q = np.asarray(cloud, dtype=float)
    if p.ndim != 2 or q.ndim != 2 or p.shape[1] != q.shape[1] or len(q) == 0:
        raise ValueError("points/cloud must be nonempty 2D arrays with equal dimensions")
    if np.any(~np.isfinite(p)) or np.any(~np.isfinite(q)):
        raise ValueError("points/cloud must be finite")
    d2 = np.sum((p[:, None, :] - q[None, :, :]) ** 2, axis=2)
    return np.sqrt(np.min(d2, axis=1))


def historical_analog_memory(
    current_edge_climate: np.ndarray,
    lgm_edge_climate: np.ndarray,
    host_current_cloud: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return M_host, D_current and D_LGM for every frozen genetic edge.

    Edge climate arrays have shape ``(n_edges, n_positions, n_climate_dims)``.
    Positive memory means that the edge geography is farther from the present
    host-resource climate cloud at the LGM than at 0 BP.
    """
    current = np.asarray(current_edge_climate, dtype=float)
    lgm = np.asarray(lgm_edge_climate, dtype=float)
    host = np.asarray(host_current_cloud, dtype=float)
    if current.ndim != 3 or lgm.shape != current.shape:
        raise ValueError(
            "edge climate arrays must have shape (edges, positions, climate_dims)"
        )
    if host.ndim != 2 or host.shape[1] != current.shape[2] or len(host) == 0:
        raise ValueError("host cloud dimension mismatch")
    if np.any(~np.isfinite(current)) or np.any(~np.isfinite(lgm)):
        raise ValueError("edge climate arrays must be finite")

    n_edges, n_positions, n_dims = current.shape
    d_current = nearest_cloud_distance(
        current.reshape(-1, n_dims), host
    ).reshape(n_edges, n_positions).mean(axis=1)
    d_lgm = nearest_cloud_distance(
        lgm.reshape(-1, n_dims), host
    ).reshape(n_edges, n_positions).mean(axis=1)
    return d_lgm - d_current, d_current, d_lgm


def zscore(x: np.ndarray) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    if values.ndim != 1 or np.any(~np.isfinite(values)):
        raise ValueError("vector must be one-dimensional and finite")
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("vector has no finite variation")
    return (values - float(np.mean(values))) / sd


def residualized_unique_fraction(primary: np.ndarray, controls: np.ndarray) -> float:
    """Fraction of standardized primary variance left after linear controls."""
    y = zscore(primary)
    x = np.asarray(controls, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    if x.ndim != 2 or x.shape[0] != len(y) or np.any(~np.isfinite(x)):
        raise ValueError("control matrix is invalid")
    design = np.column_stack([np.ones(len(y)), x])
    coef = np.linalg.lstsq(design, y, rcond=None)[0]
    resid = y - design @ coef
    return float(np.var(resid, ddof=0) / np.var(y, ddof=0))


def partial_host_memory_beta(
    turnover: np.ndarray,
    host_memory: np.ndarray,
    self_memory: np.ndarray,
    current_host_distance: np.ndarray,
) -> float:
    """Species-specific standardized partial coefficient for M_host."""
    y = zscore(turnover)
    host = zscore(host_memory)
    self_hist = zscore(self_memory)
    current = zscore(current_host_distance)
    design = np.column_stack([np.ones(len(y)), host, self_hist, current])
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("predictor matrix is rank deficient")
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    return float(beta[1])


def equal_species_mean(betas) -> float:
    """Equal-species aggregation; edge-rich species receive no extra weight."""
    values = np.asarray(list(betas), dtype=float)
    if values.ndim != 1 or len(values) == 0 or np.any(~np.isfinite(values)):
        raise ValueError("betas must be a nonempty finite vector")
    return float(np.mean(values))


__all__ = [
    "standardize_from_reference",
    "nearest_cloud_distance",
    "historical_analog_memory",
    "zscore",
    "residualized_unique_fraction",
    "partial_host_memory_beta",
    "equal_species_mean",
]