from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import average_ranks, rank01
from .genetic_geometry import endpoint_disjoint_training_counts


@dataclass(frozen=True)
class CrossfitIBDResult:
    """Leave-two-localities-out IBD residualization for one species.

    ``residual_turnover`` is the within-species rank of the cross-fit residual and
    is the genetic TTF target after removing the monotone geographic-distance
    expectation. The auxiliary arrays are retained for audit and qualification.
    """

    residual: np.ndarray
    residual_turnover: np.ndarray
    observed_rank_fraction: np.ndarray
    expected_rank_fraction: np.ndarray
    geographic_rank_fraction: np.ndarray
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


def crossfit_ibd_residuals(
    genetic_distance: np.ndarray,
    geographic_distance: np.ndarray,
    edge_nodes: np.ndarray,
    *,
    min_training_edges: int = 5,
) -> CrossfitIBDResult:
    """Remove a monotone within-species IBD expectation without endpoint leakage.

    For each scored edge ``(i, j)``, the nuisance fit excludes every other edge
    incident to locality ``i`` or ``j``. The remaining endpoint-disjoint edges
    define an out-of-pair training set. Genetic and geographic distances are
    represented only by their training-set order. A rank-linear IBD expectation
    is fitted there and evaluated for the held-out edge using order fractions
    against the training set.

    This construction has three useful properties for pairwise genetic data:

    1. the expected value for an edge never uses a genetic outcome involving
       either of that edge's localities;
    2. strictly increasing transformations of genetic or geographic distance do
       not change the result;
    3. any noiseless strictly monotone IBD relation is removed exactly, while
       departures from that relation remain available to TTF.

    The returned residual is finally ranked within species so locus-specific
    genetic-distance scale does not determine cross-species field weighting.
    This is the biological IBD layer. The current core TTF v0.11 architecture is
    applied downstream as a separate geometry/inference layer: density-scaled
    species-local graphs, training-only edge-length orthogonalization, and the
    profiled-private null. Qualification of those layers does not automatically
    qualify this genetic interface; the combined pipeline must pass its own
    genetic-specific geometry gate.
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
    if min_training_edges < 3:
        raise ValueError("min_training_edges must be >= 3")

    nodes = _validated_edges(edge_nodes, len(genetic))
    available = endpoint_disjoint_training_counts(nodes)
    if int(np.min(available)) < int(min_training_edges):
        worst = int(np.min(available))
        raise ValueError(
            "endpoint-disjoint IBD cross-fit is not estimable for every edge: "
            f"minimum available training edges={worst}, required={min_training_edges}"
        )

    residual = np.empty(len(genetic), dtype=float)
    observed = np.empty(len(genetic), dtype=float)
    expected = np.empty(len(genetic), dtype=float)
    geo_fraction = np.empty(len(genetic), dtype=float)

    for index, (left, right) in enumerate(nodes):
        train = (
            (nodes[:, 0] != left)
            & (nodes[:, 1] != left)
            & (nodes[:, 0] != right)
            & (nodes[:, 1] != right)
        )
        x_train = geographic[train]
        y_train = genetic[train]

        x_rank = _training_rank_fraction(x_train)
        y_rank = _training_rank_fraction(y_train)
        x_centered = x_rank - float(x_rank.mean())
        y_centered = y_rank - float(y_rank.mean())
        denominator = float(np.dot(x_centered, x_centered))
        beta = (
            0.0
            if denominator <= np.finfo(float).eps
            else float(np.dot(x_centered, y_centered) / denominator)
        )
        intercept = float(y_rank.mean() - beta * x_rank.mean())

        held_x = _rank_fraction_against(geographic[index], x_train)
        held_y = _rank_fraction_against(genetic[index], y_train)
        held_expected = float(intercept + beta * held_x)

        geo_fraction[index] = held_x
        observed[index] = held_y
        expected[index] = held_expected
        residual[index] = float(held_y - held_expected)

    return CrossfitIBDResult(
        residual=residual,
        residual_turnover=rank01(residual),
        observed_rank_fraction=observed,
        expected_rank_fraction=expected,
        geographic_rank_fraction=geo_fraction,
        n_training_edges=available,
    )
