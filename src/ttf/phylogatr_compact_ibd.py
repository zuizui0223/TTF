from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit, prange

from .core import average_ranks, rank01
from .genetic_geometry import endpoint_disjoint_training_counts
from .genetic_ibd import CrossfitIBDResult


@dataclass(frozen=True)
class CompactCrossfitIBDDesign:
    """Exact endpoint-safe rank-IBD geometry without quadratic edge caches.

    The original reference design stores one training-index array and one
    geographic centered-rank array per target edge.  That is exact but O(m^2)
    memory for a species with m graph edges.  This design stores only the global
    geographic ranks plus, for each locality, the rank correction contributed by
    incident edges.  Target-specific subset ranks are reconstructed exactly at
    scoring time.
    """

    geographic_distance: np.ndarray
    edge_nodes: np.ndarray
    global_geographic_rank: np.ndarray
    vertex_rank_correction: np.ndarray
    n_training_edges: np.ndarray


def _validated_edges(edge_nodes: np.ndarray, n_edges: int) -> np.ndarray:
    nodes = np.asarray(edge_nodes, dtype=np.int64)
    if nodes.ndim != 2 or nodes.shape != (n_edges, 2):
        raise ValueError("edge_nodes must have shape (n_edges, 2)")
    if np.any(nodes < 0) or np.any(nodes[:, 0] == nodes[:, 1]):
        raise ValueError("edge nodes must be non-negative and self-edge free")
    canonical = np.sort(nodes, axis=1)
    if len(np.unique(canonical, axis=0)) != len(canonical):
        raise ValueError("duplicate undirected edges are not allowed")
    return canonical


def _vertex_rank_correction(values: np.ndarray, nodes: np.ndarray) -> np.ndarray:
    """Return less + half-tie counts contributed by each locality's incident edges."""
    x = np.asarray(values, dtype=float)
    n_vertices = int(np.max(nodes)) + 1
    correction = np.empty((n_vertices, len(x)), dtype=float)
    left_nodes = nodes[:, 0]
    right_nodes = nodes[:, 1]
    for vertex in range(n_vertices):
        incident = np.flatnonzero((left_nodes == vertex) | (right_nodes == vertex))
        if len(incident) == 0:
            correction[vertex, :] = 0.0
            continue
        ordered = np.sort(x[incident], kind="stable")
        left = np.searchsorted(ordered, x, side="left")
        right = np.searchsorted(ordered, x, side="right")
        correction[vertex, :] = 0.5 * (left + right)
    return correction


def prepare_compact_crossfit_ibd_design(
    geographic_distance: np.ndarray,
    edge_nodes: np.ndarray,
    *,
    min_training_edges: int = 5,
) -> CompactCrossfitIBDDesign:
    geographic = np.asarray(geographic_distance, dtype=float)
    if geographic.ndim != 1 or len(geographic) < 3:
        raise ValueError("geographic_distance must be a 1D array with >=3 edges")
    if not np.isfinite(geographic).all() or np.any(geographic < 0):
        raise ValueError("geographic distances must be finite and non-negative")
    if min_training_edges < 3:
        raise ValueError("min_training_edges must be >= 3")

    nodes = _validated_edges(edge_nodes, len(geographic))
    available = endpoint_disjoint_training_counts(nodes)
    if int(np.min(available)) < int(min_training_edges):
        worst = int(np.min(available))
        raise ValueError(
            "endpoint-disjoint IBD cross-fit is not estimable for every edge: "
            f"minimum available training edges={worst}, required={min_training_edges}"
        )

    return CompactCrossfitIBDDesign(
        geographic_distance=geographic.copy(),
        edge_nodes=nodes.copy(),
        global_geographic_rank=average_ranks(geographic),
        vertex_rank_correction=_vertex_rank_correction(geographic, nodes),
        n_training_edges=available.copy(),
    )


@njit(cache=True, parallel=True)
def _compact_residual_kernel(
    geographic: np.ndarray,
    genetic: np.ndarray,
    geographic_rank: np.ndarray,
    genetic_rank: np.ndarray,
    geographic_vertex_correction: np.ndarray,
    genetic_vertex_correction: np.ndarray,
    nodes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reconstruct target-specific subset ranks with fixed within-target summation order."""
    n_edges = len(geographic)
    residual = np.empty(n_edges, dtype=np.float64)
    observed = np.empty(n_edges, dtype=np.float64)
    expected = np.empty(n_edges, dtype=np.float64)
    held_geographic = np.empty(n_edges, dtype=np.float64)
    eps = np.finfo(np.float64).eps

    for target in prange(n_edges):
        left = nodes[target, 0]
        right = nodes[target, 1]
        target_x = geographic[target]
        target_y = genetic[target]
        count = 0
        sx = 0.0
        sy = 0.0
        sxx = 0.0
        sxy = 0.0
        held_x = 0.0
        held_y = 0.0

        for edge in range(n_edges):
            x_weight = 1.0 if target_x < geographic[edge] else (0.5 if target_x == geographic[edge] else 0.0)
            y_weight = 1.0 if target_y < genetic[edge] else (0.5 if target_y == genetic[edge] else 0.0)
            x_rank = (
                geographic_rank[edge]
                - geographic_vertex_correction[left, edge]
                - geographic_vertex_correction[right, edge]
                + x_weight
            )
            y_rank = (
                genetic_rank[edge]
                - genetic_vertex_correction[left, edge]
                - genetic_vertex_correction[right, edge]
                + y_weight
            )
            if edge == target:
                held_x = x_rank
                held_y = y_rank

            a = nodes[edge, 0]
            b = nodes[edge, 1]
            if a == left or b == left or a == right or b == right:
                continue
            count += 1
            sx += x_rank
            sy += y_rank
            sxx += x_rank * x_rank
            sxy += x_rank * y_rank

        n = float(count)
        x_mean = sx / n
        y_mean = sy / n
        denominator = sxx - sx * sx / n
        covariance = sxy - sx * sy / n
        beta = 0.0 if denominator <= eps else covariance / denominator
        scale = n + 1.0
        observed[target] = held_y / scale
        held_geographic[target] = held_x / scale
        expected[target] = (y_mean - beta * x_mean + beta * held_x) / scale
        residual[target] = observed[target] - expected[target]

    return residual, observed, expected, held_geographic


def crossfit_ibd_residuals_compact(
    genetic_distance: np.ndarray,
    design: CompactCrossfitIBDDesign,
) -> CrossfitIBDResult:
    genetic = np.asarray(genetic_distance, dtype=float)
    if genetic.ndim != 1 or genetic.shape != design.geographic_distance.shape:
        raise ValueError("genetic_distance shape differs from frozen compact IBD design")
    if not np.isfinite(genetic).all() or np.any(genetic < 0):
        raise ValueError("genetic distances must be finite and non-negative")

    genetic_rank = average_ranks(genetic)
    if np.array_equal(genetic_rank, design.global_geographic_rank):
        zero = np.zeros(len(genetic), dtype=float)
        held = _held_rank_fraction_for_identical_order(design)
        return CrossfitIBDResult(
            residual=zero,
            residual_turnover=np.full(len(zero), 0.5, dtype=float),
            observed_rank_fraction=held.copy(),
            expected_rank_fraction=held.copy(),
            geographic_rank_fraction=held,
            n_training_edges=design.n_training_edges.copy(),
        )

    genetic_vertex_correction = _vertex_rank_correction(genetic, design.edge_nodes)
    residual, observed, expected, held_geographic = _compact_residual_kernel(
        design.geographic_distance,
        genetic,
        design.global_geographic_rank,
        genetic_rank,
        design.vertex_rank_correction,
        genetic_vertex_correction,
        design.edge_nodes,
    )
    return CrossfitIBDResult(
        residual=residual,
        residual_turnover=rank01(residual),
        observed_rank_fraction=observed,
        expected_rank_fraction=expected,
        geographic_rank_fraction=held_geographic,
        n_training_edges=design.n_training_edges.copy(),
    )


def _held_rank_fraction_for_identical_order(design: CompactCrossfitIBDDesign) -> np.ndarray:
    x = design.geographic_distance
    nodes = design.edge_nodes
    out = np.empty(len(x), dtype=float)
    for target, (left, right) in enumerate(nodes):
        value = x[target]
        weight = (x[target] < x).astype(float) + 0.5 * (x[target] == x).astype(float)
        held_rank = (
            design.global_geographic_rank[target]
            - design.vertex_rank_correction[left, target]
            - design.vertex_rank_correction[right, target]
            + weight[target]
        )
        out[target] = held_rank / (float(design.n_training_edges[target]) + 1.0)
    return out


__all__ = [
    "CompactCrossfitIBDDesign",
    "crossfit_ibd_residuals_compact",
    "prepare_compact_crossfit_ibd_design",
]
