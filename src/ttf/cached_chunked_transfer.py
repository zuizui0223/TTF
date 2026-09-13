from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .batch import BatchTransferResult
from .chunked_transfer import (
    ChunkedTransferGeometry,
    _kernel_block,
    _train_value_chunk,
    _validate_train_values,
)
from .core import spearman_rho


@dataclass(frozen=True)
class CachedChunkedTransferGeometry:
    """Response-blind cache for the exact chunked Gaussian transfer operator.

    Only segment points and opportunity denominators are retained. The dense
    evaluation-edge by training-edge projection is still reconstructed in train
    chunks, so the scientific estimator and memory order remain those of
    ``score_chunked_batch``.
    """

    base: ChunkedTransferGeometry
    eval_points: Mapping[str, np.ndarray]
    eval_denominator: Mapping[str, np.ndarray]
    edge_chunk_size: int
    train_chunk_size: int


def prepare_cached_chunked_transfer(
    prepared: ChunkedTransferGeometry,
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> CachedChunkedTransferGeometry:
    """Cache every response-independent quantity from pass 1 exactly once."""
    if edge_chunk_size < 1 or train_chunk_size < 1:
        raise ValueError("chunk sizes must be positive")
    t = (np.arange(prepared.segment_points, dtype=float) + 0.5) / prepared.segment_points
    h2 = prepared.bandwidth * prepared.bandwidth
    tiny = np.finfo(float).tiny
    n_train = len(prepared.train_positions)
    points_map: dict[str, np.ndarray] = {}
    denominator_map: dict[str, np.ndarray] = {}

    for species in prepared.eval_species:
        start = prepared.eval_start[species]
        end = prepared.eval_end[species]
        points = start[:, None, :] + t[None, :, None] * (end - start)[:, None, :]
        denominator = np.empty((len(start), prepared.segment_points), dtype=float)
        for e0 in range(0, len(start), int(edge_chunk_size)):
            e1 = min(e0 + int(edge_chunk_size), len(start))
            flat = points[e0:e1].reshape(-1, start.shape[1])
            opportunity = np.zeros(len(flat), dtype=float)
            for p0 in range(0, n_train, int(train_chunk_size)):
                p1 = min(p0 + int(train_chunk_size), n_train)
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
        points_map[species] = points
        denominator_map[species] = denominator

    return CachedChunkedTransferGeometry(
        base=prepared,
        eval_points=points_map,
        eval_denominator=denominator_map,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )


def score_cached_chunked_batch(
    cached: CachedChunkedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
) -> BatchTransferResult:
    """Score one or many worlds with cached exact opportunity denominators."""
    prepared = cached.base
    train_arrays, width = _validate_train_values(prepared, train_turnover)
    score_map: dict[str, np.ndarray] = {}
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    h2 = prepared.bandwidth * prepared.bandwidth
    n_train = len(prepared.train_positions)

    for species in prepared.eval_species:
        if species not in eval_turnover:
            raise ValueError(f"missing evaluation turnover for {species}")
        start = prepared.eval_start[species]
        target = np.asarray(eval_turnover[species], dtype=float)
        if target.ndim != 2 or target.shape != (len(start), width):
            raise ValueError(f"evaluation turnover shape drift for {species}")
        if not np.isfinite(target).all():
            raise ValueError(f"non-finite evaluation turnover for {species}")

        points = cached.eval_points[species]
        denominator = cached.eval_denominator[species]
        prior_edge = (
            prepared.prior_strength * prepared.prior_mean / denominator
        ).mean(axis=1)
        predicted = np.repeat(prior_edge[:, None], width, axis=1)

        for p0 in range(0, n_train, cached.train_chunk_size):
            p1 = min(p0 + cached.train_chunk_size, n_train)
            values = _train_value_chunk(prepared, train_arrays, p0, p1, width)
            positions = prepared.train_positions[p0:p1]
            weights = prepared.train_weights[p0:p1]
            for e0 in range(0, len(start), cached.edge_chunk_size):
                e1 = min(e0 + cached.edge_chunk_size, len(start))
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
    "CachedChunkedTransferGeometry",
    "prepare_cached_chunked_transfer",
    "score_cached_chunked_batch",
]
