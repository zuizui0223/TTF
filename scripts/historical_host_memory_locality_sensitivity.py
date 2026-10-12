"""Exploratory leave-one-locality sensitivity of already opened HHM edge data.

This module accepts one species' arrays in memory. It reads/writes no files,
opens no sequences, computes no p-values and cannot authorize a new panel.
The caller must verify source identities and the original 259-species domain.

The unscaled-response slope is in p-distance units per ORIGINAL full-sample
predictor SD. Predictor centering/scales remain fixed across deletions. IBD
fits and response SD are recomputed on each retained graph. Thus the optional
standardized_slope divides only Y by its deletion-specific SD while retaining
original X units; it is not a new fully restandardized subset estimand.
"""
from __future__ import annotations

import numpy as np


class NotEstimable(RuntimeError):
    """An identified geometry/model limitation, retained as a result status."""


def _validate_edges(nodes: np.ndarray) -> np.ndarray:
    nodes = np.asarray(nodes)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or len(nodes) == 0:
        raise ValueError("edge_nodes must be a nonempty m-by-2 array")
    if not np.issubdtype(nodes.dtype, np.integer) or np.any(nodes < 0):
        raise ValueError("locality IDs must be nonnegative integers")
    if np.any(nodes[:, 0] == nodes[:, 1]):
        raise ValueError("self edges cannot be repaired by this diagnostic")
    if len(np.unique(np.sort(nodes, axis=1), axis=0)) != len(nodes):
        raise ValueError("duplicate undirected edges")
    return nodes.astype(np.int64, copy=False)


def endpoint_disjoint_residual(
    edge_nodes: np.ndarray, chord_distance_km: np.ndarray, p_distance: np.ndarray
) -> np.ndarray:
    """Re-fit each edge on edges sharing neither endpoint, in linear time.

    Every fit needs at least five training edges and nonsingular log-distance
    variation. One invalid fit makes this graph variant NOT_ESTIMABLE; no
    invalid edge is discarded or assigned a fabricated residual.
    """
    nodes = _validate_edges(edge_nodes)
    distance = np.asarray(chord_distance_km, dtype=float)
    y = np.asarray(p_distance, dtype=float)
    m = len(nodes)
    if distance.shape != (m,) or y.shape != (m,):
        raise ValueError("edge, distance and response rows are misaligned")
    if not np.isfinite(distance).all() or not np.isfinite(y).all() or np.any(distance < 0):
        raise ValueError("nonfinite or invalid input")
    if np.any((y < 0) | (y > 1)):
        raise ValueError("input must be uncorrected edge p-distance in [0, 1]")

    x = np.log1p(distance)
    scale = float(x.std(ddof=0))
    if scale <= 1e-12:
        raise NotEstimable("INVARIANT_GEOGRAPHIC_DISTANCE")
    x = (x - x.mean()) / scale
    _, inverse = np.unique(nodes, return_inverse=True)
    compact = inverse.reshape(nodes.shape)
    left, right = compact[:, 0], compact[:, 1]
    ids = np.r_[left, right]

    def node_sum(v: np.ndarray) -> np.ndarray:
        return np.bincount(ids, weights=np.r_[v, v])

    degree = np.bincount(ids)
    count = m - degree[left] - degree[right] + 1
    if np.any(count < 5):
        raise NotEstimable("INSUFFICIENT_ENDPOINT_DISJOINT_EDGES")

    def disjoint_sum(v: np.ndarray) -> np.ndarray:
        local = node_sum(v)
        return v.sum() - local[left] - local[right] + v

    sx = disjoint_sum(x)
    sxx = disjoint_sum(x * x)
    sy = disjoint_sum(y)
    sxy = disjoint_sum(x * y)
    determinant = count * sxx - sx * sx
    tolerance = 1e-12 * np.maximum(np.abs(count * sxx), 1.0)
    if np.any(determinant <= tolerance):
        raise NotEstimable("SINGULAR_ENDPOINT_DISJOINT_FIT")
    slope = (count * sxy - sx * sy) / determinant
    intercept = (sy - slope * sx) / count
    residual = y - intercept - slope * x
    if not np.isfinite(residual).all():
        raise NotEstimable("NONFINITE_IBD_RESIDUAL")
    return residual


def _fit_variant(nodes, distance, fixed_predictors, y) -> dict:
    if len(nodes) == 0:
        raise NotEstimable("NO_RETAINED_EDGES")
    residual = endpoint_disjoint_residual(nodes, distance, y)
    design = np.column_stack((np.ones(len(nodes)), fixed_predictors))
    coefficients, _, rank, _ = np.linalg.lstsq(design, residual, rcond=None)
    if rank != 4:
        raise NotEstimable("RANK_DEFICIENT_PREDICTORS")
    b = float(coefficients[1])
    sd = float(residual.std(ddof=0))
    if not np.isfinite(b) or not np.isfinite(sd):
        raise NotEstimable("NONFINITE_PARTIAL_SLOPE")
    return {
        "status": "ESTIMABLE",
        "unscaled_slope": b,
        "standardized_slope": b / sd if sd > 1e-12 else None,
        "post_ibd_response_sd": sd,
        "retained_edges": int(len(nodes)),
    }


def diagnose_locality_deletions(
    edge_nodes: np.ndarray,
    chord_distance_km: np.ndarray,
    predictor_values: np.ndarray,
    p_distances: np.ndarray,
    n_localities: int,
) -> dict:
    """Run every single-locality deletion, including nonestimable variants.

    Columns are M_host, M_self, D_host_0BP. The caller must keep original
    edge-row alignment and re-create edge p-distances only in a private
    runtime from the already opened source. No original graph is mutated.
    This helper does NOT itself verify file hashes or provide an empirical
    data loader. A baseline mismatch against earlier coefficients must be
    resolved before interpreting any downstream sensitivity result.
    """
    nodes = _validate_edges(edge_nodes)
    if (isinstance(n_localities, bool) or not isinstance(n_localities, (int, np.integer))
            or n_localities < 2 or np.any(nodes >= n_localities)):
        raise ValueError("n_localities does not cover the original integer IDs")
    distance = np.asarray(chord_distance_km, dtype=float)
    y = np.asarray(p_distances, dtype=float)
    x = np.asarray(predictor_values, dtype=float)
    if x.shape != (len(nodes), 3) or not np.isfinite(x).all():
        raise ValueError("three finite aligned predictor columns are required")
    sd_x = x.std(axis=0, ddof=0)
    if np.any(sd_x <= 1e-12):
        raise ValueError("invariant full-sample predictor")
    fixed_x = (x - x.mean(axis=0)) / sd_x
    baseline = _fit_variant(nodes, distance, fixed_x, y)

    rows = []
    for locality in range(int(n_localities)):
        keep = np.all(nodes != locality, axis=1)
        row = {"removed_locality": locality, "removed_edges": int((~keep).sum())}
        try:
            row.update(_fit_variant(nodes[keep], distance[keep], fixed_x[keep], y[keep]))
        except NotEstimable as error:
            row.update({"status": "NOT_ESTIMABLE", "reason": str(error),
                        "unscaled_slope": None, "standardized_slope": None,
                        "retained_edges": int(keep.sum())})
        rows.append(row)

    valid = [r for r in rows if r["status"] == "ESTIMABLE"]
    slopes = [r["unscaled_slope"] for r in valid]
    base_sign = np.sign(baseline["unscaled_slope"])
    reversals = sum(np.sign(v) == -base_sign for v in slopes) if base_sign != 0 else None
    all_same = (all(np.sign(v) == base_sign for v in slopes)
                if len(valid) == n_localities and base_sign != 0 else None)
    return {
        "schema": "ttf_hhm_exploratory_locality_deletion_v0.1",
        "status": "EXPLORATORY_DIAGNOSTIC_NOT_CONFIRMATORY",
        "baseline": baseline,
        "n_deletions": int(n_localities),
        "n_estimable_deletions": len(valid),
        "n_nonestimable_deletions": int(n_localities) - len(valid),
        "deletions": rows,
        "n_sign_reversals": reversals,
        "min_unscaled_slope": min(slopes) if slopes else None,
        "max_unscaled_slope": max(slopes) if slopes else None,
        "all_deletions_same_nonzero_sign": all_same,
        "extrema_cover_estimable_variants_only": True,
        "predictor_units": "one ORIGINAL full-sample within-species population SD",
        "ibd_recomputed_for_every_deletion": True,
        "species_redefined": False,
        "original_graph_rewired": False,
        "primary_result_replaced": False,
        "p_value": None,
        "confidence_interval": None,
    }
