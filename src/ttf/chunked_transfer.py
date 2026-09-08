from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .batch import BatchTransferResult
from .core import SpeciesEdges, spearman_rho


@dataclass(frozen=True)
class ChunkedTransferGeometry:
    """Geometry-only exact dense-kernel scorer without storing the dense operator.

    The mathematical estimator is unchanged from ``prepare_transfer``: every
    training edge contributes its Gaussian kernel weight, each training species
    has total weight one, the same prior is used at every segment quadrature
    point, and held-out edge exposure is the mean over the same segment points.
    Only the execution order is changed so the full eval-edge by train-edge
    projection is never retained across held-out species.
    """

    train_species: tuple[str, ...]
    train_slices: Mapping[str, slice]
    train_positions: np.ndarray
    train_weights: np.ndarray
    eval_species: tuple[str, ...]
    eval_start: Mapping[str, np.ndarray]
    eval_end: Mapping[str, np.ndarray]
    bandwidth: float
    prior_strength: float
    prior_mean: float
    segment_points: int


def prepare_chunked_transfer(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    *,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
) -> ChunkedTransferGeometry:
    if bandwidth <= 0 or prior_strength < 0 or segment_points < 1:
        raise ValueError("invalid transfer configuration")
    if len(train_edges) == 0 or len(eval_edges) == 0:
        raise ValueError("non-empty train and evaluation edge sets are required")
    train_names = tuple(edges.species for edges in train_edges)
    eval_names = tuple(edges.species for edges in eval_edges)
    if len(set(train_names)) != len(train_names) or len(set(eval_names)) != len(eval_names):
        raise ValueError("species labels must be unique within each split")
    if set(train_names) & set(eval_names):
        raise ValueError("train and evaluation species must be disjoint")

    train_slices: dict[str, slice] = {}
    positions: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    cursor = 0
    for edges in train_edges:
        if edges.n_edges < 1:
            raise ValueError("each training species must have at least one edge")
        train_slices[edges.species] = slice(cursor, cursor + edges.n_edges)
        cursor += edges.n_edges
        positions.append(np.asarray(edges.midpoint, dtype=float))
        weights.append(np.full(edges.n_edges, 1.0 / edges.n_edges, dtype=float))

    return ChunkedTransferGeometry(
        train_species=train_names,
        train_slices=train_slices,
        train_positions=np.vstack(positions),
        train_weights=np.concatenate(weights),
        eval_species=eval_names,
        eval_start={edges.species: np.asarray(edges.start, dtype=float) for edges in eval_edges},
        eval_end={edges.species: np.asarray(edges.end, dtype=float) for edges in eval_edges},
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=float(prior_mean),
        segment_points=int(segment_points),
    )


def _validate_train_values(
    prepared: ChunkedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], int]:
    widths: set[int] = set()
    arrays: dict[str, np.ndarray] = {}
    for species in prepared.train_species:
        if species not in train_turnover:
            raise ValueError(f"missing train turnover for {species}")
        values = np.asarray(train_turnover[species], dtype=float)
        if values.ndim != 2:
            raise ValueError(f"train turnover for {species} must be edges x worlds")
        sl = prepared.train_slices[species]
        if values.shape[0] != sl.stop - sl.start:
            raise ValueError(f"train turnover shape drift for {species}")
        if not np.isfinite(values).all():
            raise ValueError(f"non-finite train turnover for {species}")
        widths.add(int(values.shape[1]))
        arrays[species] = values
    if len(widths) != 1:
        raise ValueError("train turnover batch widths disagree")
    width = widths.pop()
    if width < 1:
        raise ValueError("batch must contain at least one world")
    return arrays, width


def _train_value_chunk(
    prepared: ChunkedTransferGeometry,
    arrays: Mapping[str, np.ndarray],
    p0: int,
    p1: int,
    width: int,
) -> np.ndarray:
    out = np.empty((p1 - p0, width), dtype=float)
    filled = np.zeros(p1 - p0, dtype=bool)
    for species in prepared.train_species:
        sl = prepared.train_slices[species]
        left = max(p0, sl.start)
        right = min(p1, sl.stop)
        if left >= right:
            continue
        src0 = left - sl.start
        src1 = right - sl.start
        dst0 = left - p0
        dst1 = right - p0
        out[dst0:dst1, :] = arrays[species][src0:src1, :]
        filled[dst0:dst1] = True
    if not filled.all():
        raise RuntimeError("train chunk packing left unfilled rows")
    return out


def _kernel_block(
    points: np.ndarray,
    positions: np.ndarray,
    weights: np.ndarray,
    *,
    bandwidth2: float,
) -> np.ndarray:
    delta = points[:, None, :] - positions[None, :, :]
    distance2 = np.sum(delta * delta, axis=2)
    kernel = np.exp(-0.5 * distance2 / bandwidth2)
    kernel *= weights[None, :]
    return kernel


def score_chunked_batch(
    prepared: ChunkedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> BatchTransferResult:
    """Score response worlds with the exact dense Gaussian estimator in chunks.

    No distance cutoff, sparse truncation, kernel approximation, dtype reduction,
    or response-dependent pruning is used. The algorithm first computes the exact
    opportunity denominator at every frozen segment point. It then reconstructs
    the same normalized edge-integrated projection in train-edge chunks and
    immediately multiplies each chunk into every response world. Thus the costly
    projection is never stored for more than one edge/train chunk.

    Floating-point summation order can differ from ``PreparedTransfer``; the
    underlying operator, quadrature points, weights and priors are identical.
    """

    if edge_chunk_size < 1 or train_chunk_size < 1:
        raise ValueError("chunk sizes must be positive")
    train_arrays, width = _validate_train_values(prepared, train_turnover)
    score_map: dict[str, np.ndarray] = {}
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    t = (np.arange(prepared.segment_points, dtype=float) + 0.5) / prepared.segment_points
    h2 = prepared.bandwidth * prepared.bandwidth
    tiny = np.finfo(float).tiny
    n_train = len(prepared.train_positions)

    for species in prepared.eval_species:
        if species not in eval_turnover:
            raise ValueError(f"missing evaluation turnover for {species}")
        start = prepared.eval_start[species]
        end = prepared.eval_end[species]
        target = np.asarray(eval_turnover[species], dtype=float)
        if target.ndim != 2 or target.shape != (len(start), width):
            raise ValueError(f"evaluation turnover shape drift for {species}")
        if not np.isfinite(target).all():
            raise ValueError(f"non-finite evaluation turnover for {species}")

        points = start[:, None, :] + t[None, :, None] * (end - start)[:, None, :]
        denominator = np.empty((len(start), prepared.segment_points), dtype=float)

        # Pass 1: exact opportunity mass at each segment point.
        for e0 in range(0, len(start), edge_chunk_size):
            e1 = min(e0 + edge_chunk_size, len(start))
            flat = points[e0:e1].reshape(-1, start.shape[1])
            opportunity = np.zeros(len(flat), dtype=float)
            for p0 in range(0, n_train, train_chunk_size):
                p1 = min(p0 + train_chunk_size, n_train)
                kernel = _kernel_block(
                    flat,
                    prepared.train_positions[p0:p1],
                    prepared.train_weights[p0:p1],
                    bandwidth2=h2,
                )
                opportunity += kernel.sum(axis=1)
            denominator[e0:e1, :] = np.maximum(
                opportunity.reshape(e1 - e0, prepared.segment_points)
                + prepared.prior_strength,
                tiny,
            )

        prior_edge = (
            prepared.prior_strength * prepared.prior_mean / denominator
        ).mean(axis=1)
        predicted = np.repeat(prior_edge[:, None], width, axis=1)

        # Pass 2: reconstruct normalized integrated projection by train chunk and
        # immediately apply it to all response worlds.
        for p0 in range(0, n_train, train_chunk_size):
            p1 = min(p0 + train_chunk_size, n_train)
            values = _train_value_chunk(prepared, train_arrays, p0, p1, width)
            positions = prepared.train_positions[p0:p1]
            weights = prepared.train_weights[p0:p1]
            for e0 in range(0, len(start), edge_chunk_size):
                e1 = min(e0 + edge_chunk_size, len(start))
                flat = points[e0:e1].reshape(-1, start.shape[1])
                kernel = _kernel_block(
                    flat,
                    positions,
                    weights,
                    bandwidth2=h2,
                )
                normalized = kernel / denominator[e0:e1, :].reshape(-1, 1)
                projection = normalized.reshape(
                    e1 - e0,
                    prepared.segment_points,
                    p1 - p0,
                ).mean(axis=1)
                predicted[e0:e1, :] += projection @ values

        scores = np.asarray(
            [spearman_rho(predicted[:, i], target[:, i]) for i in range(width)],
            dtype=float,
        )
        score_map[species] = scores
        finite = np.isfinite(scores)
        sums[finite] += scores[finite]
        counts[finite] += 1

    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite held-out species scores")
    return BatchTransferResult(
        statistics=sums / counts,
        species_scores=score_map,
        n_eval_species=counts,
    )


__all__ = [
    "ChunkedTransferGeometry",
    "prepare_chunked_transfer",
    "score_chunked_batch",
]
