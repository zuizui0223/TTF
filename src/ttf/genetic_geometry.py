from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import knn_edges


@dataclass(frozen=True)
class GeneticSamplingGeometry:
    """Response-blind locality geometry for a genetic TTF analysis.

    Exact duplicate coordinate rows are collapsed to one locality.  No genetic
    distance, sequence state, break result, or other outcome enters this object.
    The retained graph is the same species-local kNN graph used by core TTF.
    """

    coordinates: np.ndarray
    record_to_locality: np.ndarray
    records_per_locality: np.ndarray
    edge_nodes: np.ndarray
    endpoint_disjoint_training_edges: np.ndarray

    @property
    def n_records(self) -> int:
        return int(len(self.record_to_locality))

    @property
    def n_localities(self) -> int:
        return int(len(self.coordinates))

    @property
    def n_edges(self) -> int:
        return int(len(self.edge_nodes))

    @property
    def min_endpoint_disjoint_training_edges(self) -> int:
        if len(self.endpoint_disjoint_training_edges) == 0:
            return 0
        return int(np.min(self.endpoint_disjoint_training_edges))


def endpoint_disjoint_training_counts(edge_nodes: np.ndarray) -> np.ndarray:
    """Count nuisance-training edges that share neither endpoint with each edge.

    These counts depend only on graph geometry.  They therefore can be frozen
    before any genetic outcome is inspected and provide a direct check that the
    leave-two-localities-out IBD nuisance fit is estimable for every scored edge.
    """
    nodes = np.asarray(edge_nodes, dtype=np.int64)
    if nodes.ndim != 2 or nodes.shape[1] != 2 or len(nodes) == 0:
        raise ValueError("edge_nodes must be a non-empty m x 2 integer array")
    if np.any(nodes < 0) or np.any(nodes[:, 0] == nodes[:, 1]):
        raise ValueError("edge nodes must be non-negative and self-edge free")
    canonical = np.sort(nodes, axis=1)
    if len(np.unique(canonical, axis=0)) != len(canonical):
        raise ValueError("duplicate undirected edges are not allowed")

    counts = np.empty(len(canonical), dtype=np.int64)
    for index, (left, right) in enumerate(canonical):
        disjoint = (
            (canonical[:, 0] != left)
            & (canonical[:, 1] != left)
            & (canonical[:, 0] != right)
            & (canonical[:, 1] != right)
        )
        counts[index] = int(np.count_nonzero(disjoint))
    return counts


def prepare_genetic_sampling_geometry(
    coordinates: np.ndarray,
    *,
    k: int = 4,
    max_distance: float | None = None,
) -> GeneticSamplingGeometry:
    """Collapse exact duplicate coordinates and freeze species-local kNN geometry.

    The collapse is deliberately exact rather than radius based.  Near-by but
    non-identical coordinates remain distinct unless a separate, prospectively
    declared locality rule is supplied upstream.  This prevents an outcome-aware
    clustering radius from entering the genetic analysis after outcomes are seen.
    """
    coords = np.asarray(coordinates, dtype=float)
    if coords.ndim != 2 or len(coords) < 2 or coords.shape[1] < 1:
        raise ValueError("coordinates must be n x d with n >= 2 and d >= 1")
    if not np.isfinite(coords).all():
        raise ValueError("coordinates must be finite")

    unique, inverse, counts = np.unique(
        coords,
        axis=0,
        return_inverse=True,
        return_counts=True,
    )
    if len(unique) < 2:
        raise ValueError("at least two unique localities are required")

    edges = knn_edges(unique, k=k, max_distance=max_distance)
    disjoint_counts = endpoint_disjoint_training_counts(edges)
    return GeneticSamplingGeometry(
        coordinates=np.asarray(unique, dtype=float),
        record_to_locality=np.asarray(inverse, dtype=np.int64),
        records_per_locality=np.asarray(counts, dtype=np.int64),
        edge_nodes=np.asarray(edges, dtype=np.int64),
        endpoint_disjoint_training_edges=disjoint_counts,
    )
