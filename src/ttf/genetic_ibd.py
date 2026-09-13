from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import average_ranks, rank01
from .genetic_geometry import endpoint_disjoint_training_counts


@dataclass(frozen=True)
class CrossfitIBDResult:
    """Leave-two-localities-out IBD residualization for one species."""

    residual: np.ndarray
    residual_turnover: np.ndarray
    observed_rank_fraction: np.ndarray
    expected_rank_fraction: np.ndarray
    geographic_rank_fraction: np.ndarray
    n_training_edges: np.ndarray


@dataclass(frozen=True)
class CrossfitIBDDesign:
    """Geometry-only endpoint-safe IBD design reusable across response worlds.

    Every stored object depends only on geographic edge distance and graph
    endpoints. No genetic outcome enters this design, so it can be frozen before
    any synthetic or empirical genetic response is opened.
    """

    geographic_distance: np.ndarray
    edge_nodes: np.ndarray
    training_indices: tuple[np.ndarray, ...]
    geographic_centered_ranks: tuple[np.ndarray, ...]
    geographic_rank_means: np.ndarray
    geographic_rank_denominators: np.ndarray
    held_geographic_rank_fraction: np.ndarray
    n_training_edges: np.ndarray


def _training_rank_fraction(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    return average_ranks(x) / (len(x) + 1.0)


def _rank_fraction_against(value: float, reference: np.ndarray) -> float:
    ref = np.asarray(reference, dtype=float)
    less = int(np.count_nonzero(ref < float(value)))
    equal = int(np.count_nonzero(ref == float(value)))
    return float((less + 0.5 * equal + 0.5) / (len(ref) + 1.0))


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


def prepare_crossfit_ibd_design(
    geographic_distance: np.ndarray,
    edge_nodes: np.ndarray,
    *,
    min_training_edges: int = 5,
) -> CrossfitIBDDesign:
    """Precompute every response-blind part of endpoint-safe rank IBD fitting."""
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

    indices: list[np.ndarray] = []
    centered: list[np.ndarray] = []
    means = np.empty(len(nodes), dtype=float)
    denominators = np.empty(len(nodes), dtype=float)
    held = np.empty(len(nodes), dtype=float)

    for index, (left, right) in enumerate(nodes):
        mask = (
            (nodes[:, 0] != left)
            & (nodes[:, 1] != left)
            & (nodes[:, 0] != right)
            & (nodes[:, 1] != right)
        )
        train_index = np.flatnonzero(mask).astype(np.int64, copy=False)
        x_train = geographic[train_index]
        x_rank = _training_rank_fraction(x_train)
        x_mean = float(x_rank.mean())
        x_centered = x_rank - x_mean

        indices.append(train_index)
        centered.append(x_centered)
        means[index] = x_mean
        denominators[index] = float(np.dot(x_centered, x_centered))
        held[index] = _rank_fraction_against(geographic[index], x_train)

    return CrossfitIBDDesign(
        geographic_distance=geographic.copy(),
        edge_nodes=nodes.copy(),
        training_indices=tuple(indices),
        geographic_centered_ranks=tuple(centered),
        geographic_rank_means=means,
        geographic_rank_denominators=denominators,
        held_geographic_rank_fraction=held,
        n_training_edges=available.copy(),
    )


def crossfit_ibd_residuals_prepared(
    genetic_distance: np.ndarray,
    design: CrossfitIBDDesign,
) -> CrossfitIBDResult:
    """Score one genetic-distance vector on a frozen endpoint-safe IBD design.

    If the genetic and geographic vectors have exactly the same weak rank order
    including ties, every endpoint-disjoint subset necessarily has the same weak
    rank order as well. Under the existing rank-linear IBD estimand the fitted
    slope is therefore exactly one and every held-out residual is algebraically
    zero. That identity is returned directly, avoiding floating roundoff followed
    by rank amplification. No numerical tolerance is used.
    """
    genetic = np.asarray(genetic_distance, dtype=float)
    geographic = np.asarray(design.geographic_distance, dtype=float)
    if genetic.ndim != 1 or genetic.shape != geographic.shape:
        raise ValueError("genetic_distance shape differs from frozen IBD design")
    if not np.isfinite(genetic).all() or np.any(genetic < 0):
        raise ValueError("genetic distances must be finite and non-negative")

    if np.array_equal(average_ranks(genetic), average_ranks(geographic)):
        zero = np.zeros(len(genetic), dtype=float)
        held = np.asarray(design.held_geographic_rank_fraction, dtype=float).copy()
        return CrossfitIBDResult(
            residual=zero,
            residual_turnover=np.full(len(genetic), 0.5, dtype=float),
            observed_rank_fraction=held.copy(),
            expected_rank_fraction=held.copy(),
            geographic_rank_fraction=held,
            n_training_edges=np.asarray(design.n_training_edges, dtype=np.int64).copy(),
        )

    residual = np.empty(len(genetic), dtype=float)
    observed = np.empty(len(genetic), dtype=float)
    expected = np.empty(len(genetic), dtype=float)

    for index, train_index in enumerate(design.training_indices):
        y_train = genetic[train_index]
        y_rank = _training_rank_fraction(y_train)
        y_mean = float(y_rank.mean())
        y_centered = y_rank - y_mean
        denominator = float(design.geographic_rank_denominators[index])
        beta = (
            0.0
            if denominator <= np.finfo(float).eps
            else float(
                np.dot(design.geographic_centered_ranks[index], y_centered)
                / denominator
            )
        )
        intercept = float(
            y_mean - beta * float(design.geographic_rank_means[index])
        )
        held_y = _rank_fraction_against(genetic[index], y_train)
        held_expected = float(
            intercept + beta * float(design.held_geographic_rank_fraction[index])
        )
        observed[index] = held_y
        expected[index] = held_expected
        residual[index] = float(held_y - held_expected)

    return CrossfitIBDResult(
        residual=residual,
        residual_turnover=rank01(residual),
        observed_rank_fraction=observed,
        expected_rank_fraction=expected,
        geographic_rank_fraction=np.asarray(
            design.held_geographic_rank_fraction, dtype=float
        ).copy(),
        n_training_edges=np.asarray(design.n_training_edges, dtype=np.int64).copy(),
    )


def crossfit_ibd_residuals(
    genetic_distance: np.ndarray,
    geographic_distance: np.ndarray,
    edge_nodes: np.ndarray,
    *,
    min_training_edges: int = 5,
) -> CrossfitIBDResult:
    """Remove monotone within-species IBD without endpoint leakage.

    This convenience wrapper constructs the response-blind design and then calls
    ``crossfit_ibd_residuals_prepared``. Repeated-world qualification should
    prepare the design once and reuse it.
    """
    genetic = np.asarray(genetic_distance, dtype=float)
    geographic = np.asarray(geographic_distance, dtype=float)
    if genetic.ndim != 1 or geographic.shape != genetic.shape or len(genetic) < 3:
        raise ValueError(
            "genetic_distance and geographic_distance must be equal 1D arrays with >=3 edges"
        )
    if not np.isfinite(genetic).all() or not np.isfinite(geographic).all():
        raise ValueError("genetic and geographic distances must be finite")
    if np.any(genetic < 0) or np.any(geographic < 0):
        raise ValueError("genetic and geographic distances must be non-negative")

    design = prepare_crossfit_ibd_design(
        geographic,
        edge_nodes,
        min_training_edges=int(min_training_edges),
    )
    return crossfit_ibd_residuals_prepared(genetic, design)
