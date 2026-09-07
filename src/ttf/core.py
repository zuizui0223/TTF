from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

import numpy as np

Dissimilarity = Callable[[np.ndarray, np.ndarray], float]


@dataclass(frozen=True)
class SpeciesSample:
    """Georeferenced observations for one species.

    TTF deliberately stores one species per object so graph construction cannot
    accidentally create cross-species edges.
    """

    species: str
    coordinates: np.ndarray
    trait: np.ndarray
    blocks: np.ndarray | None = None

    def __post_init__(self) -> None:
        coords = np.asarray(self.coordinates, dtype=float)
        trait = np.asarray(self.trait)
        if coords.ndim != 2 or coords.shape[0] < 2:
            raise ValueError("coordinates must be n x d with n >= 2")
        if trait.shape[0] != coords.shape[0]:
            raise ValueError("trait and coordinates must have the same first dimension")
        if not np.isfinite(coords).all():
            raise ValueError("coordinates must be finite")
        blocks = None if self.blocks is None else np.asarray(self.blocks)
        if blocks is not None and blocks.shape != (coords.shape[0],):
            raise ValueError("blocks must have one value per observation")
        object.__setattr__(self, "coordinates", coords)
        object.__setattr__(self, "trait", trait)
        object.__setattr__(self, "blocks", blocks)
        if not str(self.species):
            raise ValueError("species must be non-empty")


@dataclass(frozen=True)
class SpeciesEdges:
    """Fixed within-species graph geometry plus rank-standardized turnover."""

    species: str
    nodes: np.ndarray
    start: np.ndarray
    end: np.ndarray
    midpoint: np.ndarray
    length: np.ndarray
    turnover: np.ndarray

    @property
    def n_edges(self) -> int:
        return int(self.nodes.shape[0])


def default_dissimilarity(a: np.ndarray, b: np.ndarray) -> float:
    """Absolute difference for scalar traits, Euclidean distance for vectors."""
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    if aa.ndim == 0 and bb.ndim == 0:
        return float(abs(float(aa) - float(bb)))
    return float(np.linalg.norm(aa - bb))


def average_ranks(values: np.ndarray) -> np.ndarray:
    """1-based average ranks with deterministic tie handling."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) == 0:
        raise ValueError("values must be a non-empty vector")
    if not np.isfinite(x).all():
        raise ValueError("rank inputs must be finite")
    order = np.argsort(x, kind="stable")
    sorted_x = x[order]
    ranks = np.empty(len(x), dtype=float)
    if len(x) == 1 or np.all(sorted_x[1:] != sorted_x[:-1]):
        ranks[order] = np.arange(1, len(x) + 1, dtype=float)
        return ranks

    starts = np.r_[0, 1 + np.flatnonzero(sorted_x[1:] != sorted_x[:-1])]
    stops = np.r_[starts[1:], len(x)]
    for start, stop in zip(starts, stops):
        ranks[order[start:stop]] = 0.5 * ((start + 1) + stop)
    return ranks


def rank01(values: np.ndarray) -> np.ndarray:
    """Map within-species edge dissimilarities to (0, 1) rank scores."""
    x = np.asarray(values, dtype=float)
    return (average_ranks(x) - 0.5) / len(x)


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman correlation; constant vectors score zero rather than disappear."""
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.ndim != 1 or yy.ndim != 1 or len(xx) != len(yy):
        raise ValueError("x and y must be equal-length vectors")
    if len(xx) < 3:
        return float("nan")
    rx = average_ranks(xx)
    ry = average_ranks(yy)
    dx = rx - float(rx.mean())
    dy = ry - float(ry.mean())
    den = float(np.sqrt(np.dot(dx, dx) * np.dot(dy, dy)))
    if den <= np.finfo(float).eps:
        return 0.0
    return float(np.dot(dx, dy) / den)


def knn_edges(
    coordinates: np.ndarray,
    *,
    k: int = 4,
    max_distance: float | None = None,
) -> np.ndarray:
    """Create a unique undirected kNN graph for one species only."""
    x = np.asarray(coordinates, dtype=float)
    if x.ndim != 2 or len(x) < 2:
        raise ValueError("coordinates must be n x d with n >= 2")
    if k < 1:
        raise ValueError("k must be >= 1")
    if max_distance is not None and max_distance <= 0:
        raise ValueError("max_distance must be positive")

    delta = x[:, None, :] - x[None, :, :]
    distance = np.sqrt(np.sum(delta * delta, axis=2))
    np.fill_diagonal(distance, np.inf)
    edges: set[tuple[int, int]] = set()
    use_k = min(int(k), len(x) - 1)
    for i in range(len(x)):
        for j in np.argsort(distance[i], kind="stable")[:use_k]:
            j = int(j)
            if max_distance is not None and distance[i, j] > float(max_distance):
                continue
            a, b = (i, j) if i < j else (j, i)
            edges.add((a, b))
    if not edges:
        raise ValueError("graph has no edges under the requested distance rule")
    return np.asarray(sorted(edges), dtype=np.int64)


def _validate_edge_nodes(nodes: np.ndarray, n: int) -> np.ndarray:
    edge_nodes = np.asarray(nodes, dtype=np.int64)
    if edge_nodes.ndim != 2 or edge_nodes.shape[1] != 2 or len(edge_nodes) == 0:
        raise ValueError("edge_nodes must be a non-empty m x 2 integer array")
    if np.any(edge_nodes < 0) or np.any(edge_nodes >= n):
        raise ValueError("edge node index out of range")
    if np.any(edge_nodes[:, 0] == edge_nodes[:, 1]):
        raise ValueError("self edges are not allowed")
    canonical = np.sort(edge_nodes, axis=1)
    if len(np.unique(canonical, axis=0)) != len(canonical):
        raise ValueError("duplicate undirected edges are not allowed")
    return canonical


def edge_turnover(
    sample: SpeciesSample,
    edge_nodes: np.ndarray,
    *,
    dissimilarity: Dissimilarity | None = None,
) -> np.ndarray:
    """Rank-standardized turnover on fixed within-species edge nodes."""
    nodes = _validate_edge_nodes(edge_nodes, len(sample.coordinates))
    if dissimilarity is None:
        trait = np.asarray(sample.trait)
        try:
            left = np.asarray(trait[nodes[:, 0]], dtype=float)
            right = np.asarray(trait[nodes[:, 1]], dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "non-numeric traits require an explicit dissimilarity function"
            ) from exc
        delta = left - right
        if delta.ndim == 1:
            raw = np.abs(delta)
        else:
            axes = tuple(range(1, delta.ndim))
            raw = np.sqrt(np.sum(delta * delta, axis=axes))
    else:
        raw = np.empty(len(nodes), dtype=float)
        for i, (a, b) in enumerate(nodes):
            raw[i] = float(dissimilarity(sample.trait[a], sample.trait[b]))
    raw = np.asarray(raw, dtype=float)
    if raw.shape != (len(nodes),) or not np.isfinite(raw).all():
        raise ValueError("dissimilarity produced non-finite edge values")
    return rank01(raw)


def build_species_edges(
    sample: SpeciesSample,
    *,
    k: int = 4,
    max_distance: float | None = None,
    dissimilarity: Dissimilarity | None = None,
    edge_nodes: np.ndarray | None = None,
) -> SpeciesEdges:
    """Measure trait turnover on a fixed graph, then rank within species."""
    nodes = (
        knn_edges(sample.coordinates, k=k, max_distance=max_distance)
        if edge_nodes is None
        else _validate_edge_nodes(edge_nodes, len(sample.coordinates))
    )
    turnover = edge_turnover(
        sample,
        nodes,
        dissimilarity=dissimilarity,
    )

    start = sample.coordinates[nodes[:, 0]]
    end = sample.coordinates[nodes[:, 1]]
    midpoint = 0.5 * (start + end)
    length = np.linalg.norm(end - start, axis=1)
    return SpeciesEdges(
        species=sample.species,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=midpoint,
        length=length,
        turnover=turnover,
    )


def balanced_schedule(
    labels: Sequence[str],
    *,
    n_realizations: int,
    per_realization: int,
    seed: int,
) -> tuple[tuple[str, ...], ...]:
    """Balanced no-replacement schedule with long-run count difference <= 1."""
    items = tuple(map(str, labels))
    if len(items) == 0 or len(set(items)) != len(items):
        raise ValueError("labels must be unique and non-empty")
    if n_realizations < 1:
        raise ValueError("n_realizations must be >= 1")
    if not 1 <= per_realization <= len(items):
        raise ValueError("per_realization must be between 1 and number of labels")

    rng = np.random.default_rng(int(seed))
    total = int(n_realizations) * int(per_realization)
    base, extra = divmod(total, len(items))
    target = np.full(len(items), base, dtype=np.int64)
    if extra:
        target[rng.permutation(len(items))[:extra]] += 1

    remaining = target.copy()
    rows: list[tuple[str, ...]] = []
    for _ in range(n_realizations):
        jitter = rng.random(len(items))
        order = np.lexsort((jitter, -remaining))
        chosen = order[:per_realization]
        if np.any(remaining[chosen] <= 0):
            raise RuntimeError("balanced schedule construction failed")
        remaining[chosen] -= 1
        rows.append(tuple(items[int(i)] for i in chosen))
    if np.any(remaining != 0):
        raise RuntimeError("balanced schedule did not meet target counts")
    return tuple(rows)


def inclusion_counts(
    schedule: Sequence[Sequence[str]],
    labels: Sequence[str],
) -> dict[str, int]:
    counts = {str(label): 0 for label in labels}
    for row in schedule:
        if len(set(row)) != len(row):
            raise ValueError("schedule row contains duplicate labels")
        for label in row:
            counts[str(label)] += 1
    return counts


def split_species(
    labels: Iterable[str],
    *,
    eval_fraction: float = 0.5,
    seed: int = 0,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Deterministic species-disjoint train/evaluation split."""
    ordered = tuple(sorted(map(str, labels)))
    if len(ordered) < 2 or len(set(ordered)) != len(ordered):
        raise ValueError("at least two unique species are required")
    if not 0.0 < eval_fraction < 1.0:
        raise ValueError("eval_fraction must lie in (0, 1)")
    rng = np.random.default_rng(int(seed))
    perm = rng.permutation(len(ordered))
    n_eval = int(round(len(ordered) * float(eval_fraction)))
    n_eval = min(max(n_eval, 1), len(ordered) - 1)
    eval_ids = tuple(ordered[int(i)] for i in perm[:n_eval])
    train_ids = tuple(ordered[int(i)] for i in perm[n_eval:])
    return train_ids, eval_ids
