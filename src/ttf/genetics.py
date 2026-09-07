from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .core import SpeciesEdges, _validate_edge_nodes, knn_edges, rank01


@dataclass(frozen=True)
class PairwiseDistanceSample:
    """Georeferenced observations represented by a fixed pairwise distance matrix.

    This adapter is intended for traits such as genetic distance, where the
    scientifically meaningful object is a metric/dissimilarity matrix rather
    than a scalar trait attached to each record. TTF never permutes matrix cells.
    """

    species: str
    coordinates: np.ndarray
    distances: np.ndarray
    labels: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        coords = np.asarray(self.coordinates, dtype=float)
        dist = np.asarray(self.distances, dtype=float)
        if coords.ndim != 2 or coords.shape[0] < 2:
            raise ValueError("coordinates must be n x d with n >= 2")
        n = coords.shape[0]
        if dist.shape != (n, n):
            raise ValueError("distances must be an n x n matrix matching coordinates")
        if not np.isfinite(coords).all() or not np.isfinite(dist).all():
            raise ValueError("coordinates and distances must be finite")
        if np.any(dist < -1e-12):
            raise ValueError("pairwise distances must be non-negative")
        if not np.allclose(dist, dist.T, rtol=1e-8, atol=1e-10):
            raise ValueError("pairwise distance matrix must be symmetric")
        if not np.allclose(np.diag(dist), 0.0, rtol=0.0, atol=1e-10):
            raise ValueError("pairwise distance matrix diagonal must be zero")
        labels = None if self.labels is None else tuple(map(str, self.labels))
        if labels is not None:
            if len(labels) != n or len(set(labels)) != n:
                raise ValueError("labels must be unique with one label per observation")
        if not str(self.species):
            raise ValueError("species must be non-empty")
        object.__setattr__(self, "coordinates", coords)
        object.__setattr__(self, "distances", np.maximum(dist, 0.0))
        object.__setattr__(self, "labels", labels)


def build_pairwise_distance_edges(
    sample: PairwiseDistanceSample,
    *,
    k: int = 4,
    max_distance: float | None = None,
    edge_nodes: np.ndarray | None = None,
) -> SpeciesEdges:
    """Build TTF edges from a fixed pairwise genetic/dissimilarity matrix.

    Graph geometry is constructed only from geography. Trait turnover is read
    from the corresponding matrix entries and rank-standardized within species.
    No distance-matrix cell shuffling is performed or required.
    """
    nodes = (
        knn_edges(sample.coordinates, k=k, max_distance=max_distance)
        if edge_nodes is None
        else _validate_edge_nodes(edge_nodes, len(sample.coordinates))
    )
    raw = sample.distances[nodes[:, 0], nodes[:, 1]]
    turnover = rank01(raw)
    start = sample.coordinates[nodes[:, 0]]
    end = sample.coordinates[nodes[:, 1]]
    return SpeciesEdges(
        species=sample.species,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=turnover,
    )


def sequence_p_distance(
    a: str | bytes | Sequence[str],
    b: str | bytes | Sequence[str],
    *,
    min_overlap: int = 1,
) -> float:
    """Uncorrected DNA p-distance on aligned A/C/G/T sites.

    Ambiguous bases and gaps are ignored pairwise. This utility is deliberately
    conservative and alignment-agnostic; publication analyses may substitute a
    prospectively chosen genetic distance without changing the TTF estimator.
    """
    if min_overlap < 1:
        raise ValueError("min_overlap must be >= 1")

    def normalize(value: str | bytes | Sequence[str]) -> np.ndarray:
        if isinstance(value, bytes):
            text = value.decode("ascii")
        elif isinstance(value, str):
            text = value
        else:
            text = "".join(map(str, value))
        text = text.upper().replace("U", "T")
        return np.asarray(list(text), dtype="U1")

    aa = normalize(a)
    bb = normalize(b)
    if aa.shape != bb.shape or aa.ndim != 1:
        raise ValueError("sequences must be aligned and equal length")
    valid = np.isin(aa, list("ACGT")) & np.isin(bb, list("ACGT"))
    overlap = int(np.count_nonzero(valid))
    if overlap < int(min_overlap):
        raise ValueError("sequences do not have enough unambiguous overlap")
    return float(np.count_nonzero(aa[valid] != bb[valid]) / overlap)


def pairwise_sequence_distance_matrix(
    sequences: Sequence[str | bytes | Sequence[str]],
    *,
    min_overlap: int = 1,
) -> np.ndarray:
    """Compute a symmetric p-distance matrix from aligned DNA sequences."""
    seqs = tuple(sequences)
    if len(seqs) < 2:
        raise ValueError("at least two sequences are required")
    out = np.zeros((len(seqs), len(seqs)), dtype=float)
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            d = sequence_p_distance(seqs[i], seqs[j], min_overlap=min_overlap)
            out[i, j] = d
            out[j, i] = d
    return out
