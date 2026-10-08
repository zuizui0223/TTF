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
    # The estimand is exact Euclidean nearest analogue distance.  A KD tree
    # avoids materializing the full n_edge_points x n_host_points tensor,
    # which can exceed runner memory for the frozen 642-species panel.
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        # Bounded-memory exact fallback for small/non-SciPy installations.
        out = np.empty(len(p), dtype=float)
        # Bound intermediate squared distance arrays to <= 200000 elements.
        block = max(1, 200000 // max(1, len(q)))
        for start in range(0, len(p), block):
            end = min(start + block, len(p))
            d2 = np.sum((p[start:end, None, :] - q[None, :, :]) ** 2, axis=2)
            out[start:end] = np.sqrt(np.min(d2, axis=1))
        return out
    tree = cKDTree(q)
    distances, _ = tree.query(p, k=1, workers=1)
    return np.asarray(distances, dtype=float)


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


def great_circle_quadrature_latlon(
    start_latlon: np.ndarray,
    end_latlon: np.ndarray,
    *,
    segment_points: int = 5,
) -> np.ndarray:
    """Interior great-circle quadrature at inherited TTF midpoint fractions.

    Input latitude/longitude arrays have shape ``(n_edges, 2)`` in degrees.
    Output has shape ``(n_edges, segment_points, 2)``.  Fractions are
    ``(i + 0.5) / segment_points``, matching the inherited TTF operator.
    """
    if segment_points < 1:
        raise ValueError("segment_points must be positive")
    a = np.asarray(start_latlon, dtype=float)
    b = np.asarray(end_latlon, dtype=float)
    if a.ndim != 2 or b.shape != a.shape or a.shape[1] != 2:
        raise ValueError("start/end must have shape (n_edges, 2)")
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("start/end must be finite")

    def unit(x: np.ndarray) -> np.ndarray:
        lat = np.deg2rad(x[:, 0])
        lon = np.deg2rad(x[:, 1])
        clat = np.cos(lat)
        return np.column_stack([clat * np.cos(lon), clat * np.sin(lon), np.sin(lat)])

    u = unit(a)
    v = unit(b)
    dot = np.clip(np.sum(u * v, axis=1), -1.0, 1.0)
    omega = np.arccos(dot)
    if np.any(np.isclose(dot, -1.0, atol=1e-12)):
        raise ValueError("antipodal edge has no unique great-circle interpolation")
    t = (np.arange(segment_points, dtype=float) + 0.5) / segment_points
    out = np.empty((len(a), segment_points, 3), dtype=float)
    for j, frac in enumerate(t):
        small = omega < 1e-12
        point = np.empty_like(u)
        if np.any(~small):
            om = omega[~small]
            denom = np.sin(om)
            point[~small] = (
                (np.sin((1.0 - frac) * om) / denom)[:, None] * u[~small]
                + (np.sin(frac * om) / denom)[:, None] * v[~small]
            )
        if np.any(small):
            point[small] = (1.0 - frac) * u[small] + frac * v[small]
        point /= np.linalg.norm(point, axis=1)[:, None]
        out[:, j, :] = point
    lat = np.rad2deg(np.arcsin(np.clip(out[:, :, 2], -1.0, 1.0)))
    lon = np.rad2deg(np.arctan2(out[:, :, 1], out[:, :, 0]))
    return np.stack([lat, lon], axis=2)


__all__.append("great_circle_quadrature_latlon")


def predictor_design_diagnostics(
    host_memory: np.ndarray,
    self_memory: np.ndarray,
    current_host_distance: np.ndarray,
) -> tuple[float, float]:
    """Return host-memory unique fraction and standardized predictor condition number."""
    host = zscore(host_memory)
    self_hist = zscore(self_memory)
    current = zscore(current_host_distance)
    x = np.column_stack([host, self_hist, current])
    unique = residualized_unique_fraction(host, np.column_stack([self_hist, current]))
    condition = float(np.linalg.cond(x))
    if not np.isfinite(condition):
        raise ValueError("non-finite predictor condition number")
    return unique, condition


def studentized_species_mean(betas) -> tuple[float, float, float]:
    """Equal-species mean, its across-species SE, and studentized statistic."""
    x = np.asarray(list(betas), dtype=float)
    if x.ndim != 1 or len(x) < 2 or np.any(~np.isfinite(x)):
        raise ValueError("at least two finite species betas are required")
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=1))
    if sd <= 0 or not np.isfinite(sd):
        raise ValueError("species betas have no finite between-species variation")
    se = sd / np.sqrt(len(x))
    return mean, float(se), float(mean / se)


__all__.extend(["predictor_design_diagnostics", "studentized_species_mean"])