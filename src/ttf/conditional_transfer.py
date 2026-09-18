from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .batch import BatchTransferResult
from .cached_chunked_transfer import CachedChunkedTransferGeometry, score_cached_chunked_batch
from .chunked_transfer import ChunkedTransferGeometry, score_chunked_batch
from .core import spearman_rho
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test


@dataclass(frozen=True)
class TargetSourcePoolDesign:
    """Response-blind target-specific training-species pools."""

    train_species: tuple[str, ...]
    eval_species: tuple[str, ...]
    source_pool: Mapping[str, tuple[str, ...]]
    pair_coverage: Mapping[str, Mapping[str, float]]
    eligible_eval_species: tuple[str, ...]
    unsupported_eval_species: tuple[str, ...]
    support_radius: float
    minimum_target_coverage: float
    minimum_source_species: int
    require_same_group: bool


@dataclass(frozen=True)
class CachedTargetConditionedTransferGeometry:
    """Response-blind cache for target-specific source-pool denominators."""

    base: ChunkedTransferGeometry
    pools: TargetSourcePoolDesign
    active_indices: Mapping[str, np.ndarray]
    eval_points: Mapping[str, np.ndarray]
    eval_denominator: Mapping[str, np.ndarray]
    pool_weight_scale: Mapping[str, float]
    edge_chunk_size: int
    train_chunk_size: int


@dataclass(frozen=True)
class FullyCachedTargetConditionedTransferGeometry:
    """Response-blind full projection cache for target-conditioned fields.

    This is an execution-only materialization of the exact Gaussian projection
    already defined by CachedTargetConditionedTransferGeometry. It stores no
    response values and changes neither source pools nor inference.
    """

    base: ChunkedTransferGeometry
    pools: TargetSourcePoolDesign
    eval_species: tuple[str, ...]
    active_indices: Mapping[str, np.ndarray]
    projection: Mapping[str, np.ndarray]
    prior_edge: Mapping[str, np.ndarray]
    edge_chunk_size: int
    train_chunk_size: int


@dataclass(frozen=True)
class ConditionalIncrementBatch:
    """Paired target-level improvement over a reference TTF field."""

    statistics: np.ndarray
    species_increments: Mapping[str, np.ndarray]
    baseline_species_scores: Mapping[str, np.ndarray]
    conditional_species_scores: Mapping[str, np.ndarray]
    eligible_eval_species: tuple[str, ...]


@dataclass(frozen=True)
class ConditionalIncrementInference:
    statistic: float
    species_increments: Mapping[str, float]
    bootstrap: MeanBootstrapResult
    eligible_eval_species: tuple[str, ...]


def _coverage_fraction(
    target_midpoint: np.ndarray,
    source_midpoint: np.ndarray,
    *,
    radius: float,
    chunk_size: int,
) -> float:
    target = np.asarray(target_midpoint, dtype=float)
    source = np.asarray(source_midpoint, dtype=float)
    if target.ndim != 2 or source.ndim != 2 or target.shape[1] != source.shape[1]:
        raise ValueError("midpoint arrays have incompatible shapes")
    if len(target) == 0 or len(source) == 0:
        raise ValueError("source and target must each contain at least one edge")
    if radius <= 0 or chunk_size < 1:
        raise ValueError("radius and chunk_size must be positive")
    threshold2 = float(radius) ** 2
    covered = 0
    for start in range(0, len(target), int(chunk_size)):
        stop = min(start + int(chunk_size), len(target))
        delta = target[start:stop, None, :] - source[None, :, :]
        nearest2 = np.min(np.sum(delta * delta, axis=2), axis=1)
        covered += int(np.count_nonzero(nearest2 <= threshold2))
    return float(covered / len(target))


def prepare_target_source_pools(
    prepared: ChunkedTransferGeometry,
    train_midpoint: Mapping[str, np.ndarray],
    eval_midpoint: Mapping[str, np.ndarray],
    *,
    support_radius: float = 500.0,
    minimum_target_coverage: float = 0.50,
    minimum_source_species: int = 5,
    train_group: Mapping[str, str] | None = None,
    eval_group: Mapping[str, str] | None = None,
    require_same_group: bool = False,
    distance_chunk_size: int = 128,
) -> TargetSourcePoolDesign:
    """Freeze target-specific source pools from geometry and optional labels.

    A source enters a target's pool only when the predeclared fraction of target
    edge midpoints lies within the support radius of a source edge midpoint.
    Optional same-group filtering accepts only response-blind labels such as
    taxonomic order. Targets with too few sources are prospectively unsupported;
    thresholds are never relaxed by this function.
    """
    if not 0.0 <= float(minimum_target_coverage) <= 1.0:
        raise ValueError("minimum_target_coverage must lie in [0, 1]")
    if minimum_source_species < 1:
        raise ValueError("minimum_source_species must be positive")
    train_names = tuple(prepared.train_species)
    eval_names = tuple(prepared.eval_species)
    if set(train_midpoint) != set(train_names):
        raise ValueError("train midpoint species mismatch")
    if set(eval_midpoint) != set(eval_names):
        raise ValueError("evaluation midpoint species mismatch")
    if require_same_group:
        if train_group is None or eval_group is None:
            raise ValueError("same-group pools require train and evaluation labels")
        if set(train_group) != set(train_names) or set(eval_group) != set(eval_names):
            raise ValueError("group-label species mismatch")

    pools: dict[str, tuple[str, ...]] = {}
    coverages: dict[str, dict[str, float]] = {}
    eligible: list[str] = []
    unsupported: list[str] = []
    for target in eval_names:
        row: dict[str, float] = {}
        selected: list[str] = []
        target_group = "" if eval_group is None else str(eval_group[target])
        for source in train_names:
            coverage = _coverage_fraction(
                eval_midpoint[target],
                train_midpoint[source],
                radius=float(support_radius),
                chunk_size=int(distance_chunk_size),
            )
            row[source] = coverage
            group_ok = True
            if require_same_group:
                source_group = str(train_group[source])
                group_ok = bool(target_group) and source_group == target_group
            if group_ok and coverage >= float(minimum_target_coverage):
                selected.append(source)
        selected_tuple = tuple(selected)
        pools[target] = selected_tuple
        coverages[target] = row
        if len(selected_tuple) >= int(minimum_source_species):
            eligible.append(target)
        else:
            unsupported.append(target)

    return TargetSourcePoolDesign(
        train_species=train_names,
        eval_species=eval_names,
        source_pool=pools,
        pair_coverage=coverages,
        eligible_eval_species=tuple(eligible),
        unsupported_eval_species=tuple(unsupported),
        support_radius=float(support_radius),
        minimum_target_coverage=float(minimum_target_coverage),
        minimum_source_species=int(minimum_source_species),
        require_same_group=bool(require_same_group),
    )


def _validate_batch_values(
    prepared: ChunkedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], int]:
    widths: set[int] = set()
    train: dict[str, np.ndarray] = {}
    for species in prepared.train_species:
        values = np.asarray(train_turnover[species], dtype=float)
        sl = prepared.train_slices[species]
        if values.ndim != 2 or values.shape[0] != sl.stop - sl.start:
            raise ValueError(f"train turnover shape drift for {species}")
        if not np.isfinite(values).all():
            raise ValueError(f"non-finite train turnover for {species}")
        widths.add(int(values.shape[1]))
        train[species] = values
    evaluation: dict[str, np.ndarray] = {}
    for species in prepared.eval_species:
        values = np.asarray(eval_turnover[species], dtype=float)
        expected = len(prepared.eval_start[species])
        if values.ndim != 2 or values.shape[0] != expected:
            raise ValueError(f"evaluation turnover shape drift for {species}")
        if not np.isfinite(values).all():
            raise ValueError(f"non-finite evaluation turnover for {species}")
        widths.add(int(values.shape[1]))
        evaluation[species] = values
    if len(widths) != 1:
        raise ValueError("response batch widths disagree")
    width = widths.pop()
    if width < 1:
        raise ValueError("batch must contain at least one world")
    return train, evaluation, width


def _source_edge_indices(
    prepared: ChunkedTransferGeometry,
    sources: tuple[str, ...],
) -> np.ndarray:
    pieces = [
        np.arange(prepared.train_slices[name].start, prepared.train_slices[name].stop)
        for name in sources
    ]
    return np.concatenate(pieces) if pieces else np.empty(0, dtype=np.int64)



def prepare_cached_target_conditioned_transfer(
    prepared: ChunkedTransferGeometry,
    pools: TargetSourcePoolDesign,
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> CachedTargetConditionedTransferGeometry:
    """Cache response-independent geometry for target-specific source pools."""
    if tuple(pools.train_species) != tuple(prepared.train_species):
        raise ValueError("source-pool training species drift")
    if tuple(pools.eval_species) != tuple(prepared.eval_species):
        raise ValueError("source-pool evaluation species drift")
    if edge_chunk_size < 1 or train_chunk_size < 1:
        raise ValueError("chunk sizes must be positive")
    if len(pools.eligible_eval_species) < 6:
        raise ValueError("fewer than six response-blind supported evaluation species")

    t = (np.arange(prepared.segment_points, dtype=float) + 0.5) / prepared.segment_points
    h2 = prepared.bandwidth * prepared.bandwidth
    tiny = np.finfo(float).tiny
    active_map: dict[str, np.ndarray] = {}
    points_map: dict[str, np.ndarray] = {}
    denominator_map: dict[str, np.ndarray] = {}
    scale_map: dict[str, float] = {}

    for species in pools.eligible_eval_species:
        active = _source_edge_indices(prepared, pools.source_pool[species])
        if len(active) == 0:
            raise RuntimeError("eligible target has no active source edges")
        positions = prepared.train_positions[active]
        source_count = len(pools.source_pool[species])
        if source_count < 1:
            raise RuntimeError("eligible target has no selected source species")
        # Preserve the unconditional field's total training-species kernel mass.
        # Core TTF gives each species total mass one; after source restriction we
        # rescale selected species so only composition, not total opportunity/prior
        # strength, changes in the paired comparison.
        scale = float(len(prepared.train_species)) / float(source_count)
        weights = prepared.train_weights[active] * scale
        start = prepared.eval_start[species]
        end = prepared.eval_end[species]
        points = start[:, None, :] + t[None, :, None] * (end - start)[:, None, :]
        denominator = np.empty((len(start), prepared.segment_points), dtype=float)

        for e0 in range(0, len(start), int(edge_chunk_size)):
            e1 = min(e0 + int(edge_chunk_size), len(start))
            flat = points[e0:e1].reshape(-1, start.shape[1])
            opportunity = np.zeros(len(flat), dtype=float)
            for p0 in range(0, len(active), int(train_chunk_size)):
                p1 = min(p0 + int(train_chunk_size), len(active))
                delta = flat[:, None, :] - positions[p0:p1][None, :, :]
                distance2 = np.sum(delta * delta, axis=2)
                kernel = np.exp(-0.5 * distance2 / h2) * weights[p0:p1][None, :]
                opportunity += kernel.sum(axis=1)
            denominator[e0:e1, :] = np.maximum(
                opportunity.reshape(e1 - e0, prepared.segment_points)
                + prepared.prior_strength,
                tiny,
            )

        active_map[species] = active
        points_map[species] = points
        denominator_map[species] = denominator
        scale_map[species] = scale

    return CachedTargetConditionedTransferGeometry(
        base=prepared,
        pools=pools,
        active_indices=active_map,
        eval_points=points_map,
        eval_denominator=denominator_map,
        pool_weight_scale=scale_map,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )


def prepare_fully_cached_target_conditioned_transfer(
    cached: CachedTargetConditionedTransferGeometry,
    *,
    eval_species: tuple[str, ...] | None = None,
) -> FullyCachedTargetConditionedTransferGeometry:
    """Materialize exact target-by-source Gaussian projection matrices.

    Only response-blind geometry, frozen source pools, kernel weights and prior
    denominators are read. The optional target subset must come from the
    already frozen eligible evaluation species.
    """
    prepared = cached.base
    allowed = tuple(cached.pools.eligible_eval_species)
    labels = allowed if eval_species is None else tuple(map(str, eval_species))
    if not labels or len(set(labels)) != len(labels):
        raise ValueError("full projection cache requires unique evaluation species")
    if not set(labels) <= set(allowed):
        raise ValueError("full projection cache contains unsupported evaluation species")

    h2 = prepared.bandwidth * prepared.bandwidth
    projection_map: dict[str, np.ndarray] = {}
    prior_map: dict[str, np.ndarray] = {}
    active_map: dict[str, np.ndarray] = {}

    for species in labels:
        active = np.asarray(cached.active_indices[species], dtype=np.int64)
        positions = prepared.train_positions[active]
        weights = prepared.train_weights[active] * float(
            cached.pool_weight_scale[species]
        )
        points = np.asarray(cached.eval_points[species], dtype=float)
        denominator = np.asarray(cached.eval_denominator[species], dtype=float)
        n_edges = len(prepared.eval_start[species])
        matrix = np.empty((n_edges, len(active)), dtype=float)

        for e0 in range(0, n_edges, cached.edge_chunk_size):
            e1 = min(e0 + cached.edge_chunk_size, n_edges)
            flat = points[e0:e1].reshape(-1, points.shape[-1])
            denom = denominator[e0:e1, :].reshape(-1, 1)
            for p0 in range(0, len(active), cached.train_chunk_size):
                p1 = min(p0 + cached.train_chunk_size, len(active))
                delta = flat[:, None, :] - positions[p0:p1][None, :, :]
                distance2 = np.sum(delta * delta, axis=2)
                kernel = (
                    np.exp(-0.5 * distance2 / h2)
                    * weights[p0:p1][None, :]
                )
                normalized = kernel / denom
                matrix[e0:e1, p0:p1] = normalized.reshape(
                    e1 - e0, prepared.segment_points, p1 - p0
                ).mean(axis=1)

        prior = (
            prepared.prior_strength * prepared.prior_mean / denominator
        ).mean(axis=1)
        if not np.isfinite(matrix).all() or not np.isfinite(prior).all():
            raise RuntimeError("non-finite full projection cache")
        active_map[species] = active.copy()
        projection_map[species] = matrix
        prior_map[species] = prior

    return FullyCachedTargetConditionedTransferGeometry(
        base=prepared,
        pools=cached.pools,
        eval_species=labels,
        active_indices=active_map,
        projection=projection_map,
        prior_edge=prior_map,
        edge_chunk_size=int(cached.edge_chunk_size),
        train_chunk_size=int(cached.train_chunk_size),
    )


def score_fully_cached_target_conditioned_batch(
    cached: FullyCachedTargetConditionedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
) -> BatchTransferResult:
    """Score the exact target-conditioned field using a full geometry cache."""
    prepared = cached.base
    train, evaluation, width = _validate_batch_values(
        prepared, train_turnover, eval_turnover
    )
    packed = np.empty((len(prepared.train_positions), width), dtype=float)
    for name in prepared.train_species:
        packed[prepared.train_slices[name], :] = train[name]

    score_map: dict[str, np.ndarray] = {}
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    for species in cached.eval_species:
        active = cached.active_indices[species]
        predicted = np.repeat(cached.prior_edge[species][:, None], width, axis=1)
        projection = cached.projection[species]
        values = packed[active, :]
        for p0 in range(0, len(active), cached.train_chunk_size):
            p1 = min(p0 + cached.train_chunk_size, len(active))
            chunk_values = values[p0:p1, :]
            for e0 in range(0, len(predicted), cached.edge_chunk_size):
                e1 = min(e0 + cached.edge_chunk_size, len(predicted))
                predicted[e0:e1, :] += (
                    projection[e0:e1, p0:p1] @ chunk_values
                )
        target = evaluation[species]
        scores = np.asarray(
            [spearman_rho(predicted[:, i], target[:, i]) for i in range(width)],
            dtype=float,
        )
        score_map[species] = scores
        finite = np.isfinite(scores)
        sums[finite] += scores[finite]
        counts[finite] += 1
    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite conditional species scores")
    return BatchTransferResult(
        statistics=sums / counts,
        species_scores=score_map,
        n_eval_species=counts,
    )


def score_cached_target_conditioned_batch(
    cached: CachedTargetConditionedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
) -> BatchTransferResult:
    """Score source-pool fields while reusing response-blind denominators."""
    prepared = cached.base
    pools = cached.pools
    train, evaluation, width = _validate_batch_values(
        prepared, train_turnover, eval_turnover
    )
    packed = np.empty((len(prepared.train_positions), width), dtype=float)
    for name in prepared.train_species:
        packed[prepared.train_slices[name], :] = train[name]

    h2 = prepared.bandwidth * prepared.bandwidth
    score_map: dict[str, np.ndarray] = {}
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)

    for species in pools.eligible_eval_species:
        active = cached.active_indices[species]
        positions = prepared.train_positions[active]
        weights = prepared.train_weights[active] * float(cached.pool_weight_scale[species])
        values = packed[active, :]
        start = prepared.eval_start[species]
        target = evaluation[species]
        points = cached.eval_points[species]
        denominator = cached.eval_denominator[species]

        prior_edge = (
            prepared.prior_strength * prepared.prior_mean / denominator
        ).mean(axis=1)
        predicted = np.repeat(prior_edge[:, None], width, axis=1)

        for p0 in range(0, len(active), cached.train_chunk_size):
            p1 = min(p0 + cached.train_chunk_size, len(active))
            chunk_positions = positions[p0:p1]
            chunk_weights = weights[p0:p1]
            chunk_values = values[p0:p1, :]
            for e0 in range(0, len(start), cached.edge_chunk_size):
                e1 = min(e0 + cached.edge_chunk_size, len(start))
                flat = points[e0:e1].reshape(-1, start.shape[1])
                delta = flat[:, None, :] - chunk_positions[None, :, :]
                distance2 = np.sum(delta * delta, axis=2)
                kernel = np.exp(-0.5 * distance2 / h2) * chunk_weights[None, :]
                normalized = kernel / denominator[e0:e1, :].reshape(-1, 1)
                projection = normalized.reshape(
                    e1 - e0, prepared.segment_points, p1 - p0
                ).mean(axis=1)
                predicted[e0:e1, :] += projection @ chunk_values

        scores = np.asarray(
            [spearman_rho(predicted[:, i], target[:, i]) for i in range(width)],
            dtype=float,
        )
        score_map[species] = scores
        finite = np.isfinite(scores)
        sums[finite] += scores[finite]
        counts[finite] += 1

    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite conditional species scores")
    return BatchTransferResult(
        statistics=sums / counts,
        species_scores=score_map,
        n_eval_species=counts,
    )


def score_cached_conditioning_increment_batch(
    unconditional: CachedChunkedTransferGeometry,
    conditioned: CachedTargetConditionedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
) -> ConditionalIncrementBatch:
    """Paired conditional-minus-unconditional score using frozen caches."""
    baseline = score_cached_chunked_batch(
        unconditional, train_turnover, eval_turnover
    )
    conditional = score_cached_target_conditioned_batch(
        conditioned, train_turnover, eval_turnover
    )
    eligible = tuple(conditioned.pools.eligible_eval_species)
    increments = {
        name: np.asarray(conditional.species_scores[name], dtype=float)
        - np.asarray(baseline.species_scores[name], dtype=float)
        for name in eligible
    }
    matrix = np.vstack([increments[name] for name in eligible])
    return ConditionalIncrementBatch(
        statistics=np.mean(matrix, axis=0),
        species_increments=increments,
        baseline_species_scores={
            name: np.asarray(baseline.species_scores[name], dtype=float).copy()
            for name in eligible
        },
        conditional_species_scores={
            name: np.asarray(conditional.species_scores[name], dtype=float).copy()
            for name in eligible
        },
        eligible_eval_species=eligible,
    )


def score_cached_pool_difference_batch(
    reference: CachedTargetConditionedTransferGeometry,
    conditioned: CachedTargetConditionedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
) -> ConditionalIncrementBatch:
    """Paired difference between two cached response-blind pool rules."""
    if reference.base is not conditioned.base and reference.base != conditioned.base:
        raise ValueError("cached source pools do not share the same base geometry")
    eligible = tuple(conditioned.pools.eligible_eval_species)
    if not set(eligible) <= set(reference.pools.eligible_eval_species):
        raise ValueError("reference pool does not support all conditioned targets")
    reference_score = score_cached_target_conditioned_batch(
        reference, train_turnover, eval_turnover
    )
    conditioned_score = score_cached_target_conditioned_batch(
        conditioned, train_turnover, eval_turnover
    )
    increments = {
        name: np.asarray(conditioned_score.species_scores[name], dtype=float)
        - np.asarray(reference_score.species_scores[name], dtype=float)
        for name in eligible
    }
    matrix = np.vstack([increments[name] for name in eligible])
    return ConditionalIncrementBatch(
        statistics=np.mean(matrix, axis=0),
        species_increments=increments,
        baseline_species_scores={
            name: np.asarray(reference_score.species_scores[name], dtype=float).copy()
            for name in eligible
        },
        conditional_species_scores={
            name: np.asarray(conditioned_score.species_scores[name], dtype=float).copy()
            for name in eligible
        },
        eligible_eval_species=eligible,
    )

def score_target_conditioned_batch(
    prepared: ChunkedTransferGeometry,
    pools: TargetSourcePoolDesign,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> BatchTransferResult:
    """Score exact Gaussian TTF fields using a frozen source pool per target."""
    if tuple(pools.train_species) != tuple(prepared.train_species):
        raise ValueError("source-pool training species drift")
    if tuple(pools.eval_species) != tuple(prepared.eval_species):
        raise ValueError("source-pool evaluation species drift")
    if edge_chunk_size < 1 or train_chunk_size < 1:
        raise ValueError("chunk sizes must be positive")
    train, evaluation, width = _validate_batch_values(
        prepared, train_turnover, eval_turnover
    )

    score_map: dict[str, np.ndarray] = {}
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    t = (np.arange(prepared.segment_points, dtype=float) + 0.5) / prepared.segment_points
    h2 = prepared.bandwidth * prepared.bandwidth
    tiny = np.finfo(float).tiny

    packed = np.empty((len(prepared.train_positions), width), dtype=float)
    for name in prepared.train_species:
        packed[prepared.train_slices[name], :] = train[name]

    for species in pools.eligible_eval_species:
        active = _source_edge_indices(prepared, pools.source_pool[species])
        if len(active) == 0:
            raise RuntimeError("eligible target has no active source edges")
        positions = prepared.train_positions[active]
        source_count = len(pools.source_pool[species])
        if source_count < 1:
            raise RuntimeError("eligible target has no selected source species")
        scale = float(len(prepared.train_species)) / float(source_count)
        weights = prepared.train_weights[active] * scale
        values = packed[active, :]
        start = prepared.eval_start[species]
        end = prepared.eval_end[species]
        target = evaluation[species]
        points = start[:, None, :] + t[None, :, None] * (end - start)[:, None, :]
        denominator = np.empty((len(start), prepared.segment_points), dtype=float)

        for e0 in range(0, len(start), int(edge_chunk_size)):
            e1 = min(e0 + int(edge_chunk_size), len(start))
            flat = points[e0:e1].reshape(-1, start.shape[1])
            opportunity = np.zeros(len(flat), dtype=float)
            for p0 in range(0, len(active), int(train_chunk_size)):
                p1 = min(p0 + int(train_chunk_size), len(active))
                delta = flat[:, None, :] - positions[p0:p1][None, :, :]
                distance2 = np.sum(delta * delta, axis=2)
                kernel = np.exp(-0.5 * distance2 / h2) * weights[p0:p1][None, :]
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

        for p0 in range(0, len(active), int(train_chunk_size)):
            p1 = min(p0 + int(train_chunk_size), len(active))
            chunk_positions = positions[p0:p1]
            chunk_weights = weights[p0:p1]
            chunk_values = values[p0:p1, :]
            for e0 in range(0, len(start), int(edge_chunk_size)):
                e1 = min(e0 + int(edge_chunk_size), len(start))
                flat = points[e0:e1].reshape(-1, start.shape[1])
                delta = flat[:, None, :] - chunk_positions[None, :, :]
                distance2 = np.sum(delta * delta, axis=2)
                kernel = np.exp(-0.5 * distance2 / h2) * chunk_weights[None, :]
                normalized = kernel / denominator[e0:e1, :].reshape(-1, 1)
                projection = normalized.reshape(
                    e1 - e0, prepared.segment_points, p1 - p0
                ).mean(axis=1)
                predicted[e0:e1, :] += projection @ chunk_values

        scores = np.asarray(
            [spearman_rho(predicted[:, i], target[:, i]) for i in range(width)],
            dtype=float,
        )
        score_map[species] = scores
        finite = np.isfinite(scores)
        sums[finite] += scores[finite]
        counts[finite] += 1

    if len(pools.eligible_eval_species) < 6:
        raise ValueError("fewer than six response-blind supported evaluation species")
    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite conditional species scores")
    return BatchTransferResult(
        statistics=sums / counts,
        species_scores=score_map,
        n_eval_species=counts,
    )


def score_conditioning_increment_batch(
    prepared: ChunkedTransferGeometry,
    pools: TargetSourcePoolDesign,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> ConditionalIncrementBatch:
    """Return geographically conditioned minus unconditional transfer skill."""
    baseline = score_chunked_batch(
        prepared,
        train_turnover,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    conditioned = score_target_conditioned_batch(
        prepared,
        pools,
        train_turnover,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    increments = {
        species: (
            np.asarray(conditioned.species_scores[species], dtype=float)
            - np.asarray(baseline.species_scores[species], dtype=float)
        )
        for species in pools.eligible_eval_species
    }
    matrix = np.vstack([increments[name] for name in pools.eligible_eval_species])
    return ConditionalIncrementBatch(
        statistics=np.mean(matrix, axis=0),
        species_increments=increments,
        baseline_species_scores={
            name: np.asarray(baseline.species_scores[name], dtype=float).copy()
            for name in pools.eligible_eval_species
        },
        conditional_species_scores={
            name: np.asarray(conditioned.species_scores[name], dtype=float).copy()
            for name in pools.eligible_eval_species
        },
        eligible_eval_species=tuple(pools.eligible_eval_species),
    )


def score_pool_difference_batch(
    prepared: ChunkedTransferGeometry,
    reference_pools: TargetSourcePoolDesign,
    conditioned_pools: TargetSourcePoolDesign,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> ConditionalIncrementBatch:
    """Compare two response-blind source-pool rules on the same targets."""
    if tuple(reference_pools.train_species) != tuple(conditioned_pools.train_species):
        raise ValueError("source-pool training species mismatch")
    if tuple(reference_pools.eval_species) != tuple(conditioned_pools.eval_species):
        raise ValueError("source-pool evaluation species mismatch")
    eligible = tuple(conditioned_pools.eligible_eval_species)
    if not set(eligible) <= set(reference_pools.eligible_eval_species):
        raise ValueError("reference pool does not support all conditioned targets")
    reference = score_target_conditioned_batch(
        prepared,
        reference_pools,
        train_turnover,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    conditioned = score_target_conditioned_batch(
        prepared,
        conditioned_pools,
        train_turnover,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    increments = {
        name: np.asarray(conditioned.species_scores[name], dtype=float)
        - np.asarray(reference.species_scores[name], dtype=float)
        for name in eligible
    }
    matrix = np.vstack([increments[name] for name in eligible])
    return ConditionalIncrementBatch(
        statistics=np.mean(matrix, axis=0),
        species_increments=increments,
        baseline_species_scores={
            name: np.asarray(reference.species_scores[name], dtype=float).copy()
            for name in eligible
        },
        conditional_species_scores={
            name: np.asarray(conditioned.species_scores[name], dtype=float).copy()
            for name in eligible
        },
        eligible_eval_species=eligible,
    )


def infer_conditioning_increment(
    scored: ConditionalIncrementBatch,
    *,
    n_bootstrap: int = 1999,
    seed: int = 0,
) -> ConditionalIncrementInference:
    """One-sided held-out-species bootstrap for a single empirical increment."""
    if np.asarray(scored.statistics).shape != (1,):
        raise ValueError("empirical increment inference requires exactly one response world")
    species = tuple(scored.eligible_eval_species)
    increments = {
        name: float(np.asarray(scored.species_increments[name])[0])
        for name in species
    }
    bootstrap = centered_species_bootstrap_mean_test(
        [increments[name] for name in species],
        n_bootstrap=int(n_bootstrap),
        seed=int(seed),
    )
    if not np.isclose(
        bootstrap.observed_mean,
        float(scored.statistics[0]),
        atol=1e-12,
        rtol=0.0,
    ):
        raise RuntimeError("conditioning increment mean drift")
    return ConditionalIncrementInference(
        statistic=float(scored.statistics[0]),
        species_increments=increments,
        bootstrap=bootstrap,
        eligible_eval_species=species,
    )


__all__ = [
    "CachedTargetConditionedTransferGeometry",
    "FullyCachedTargetConditionedTransferGeometry",
    "ConditionalIncrementBatch",
    "ConditionalIncrementInference",
    "TargetSourcePoolDesign",
    "infer_conditioning_increment",
    "prepare_cached_target_conditioned_transfer",
    "prepare_fully_cached_target_conditioned_transfer",
    "prepare_target_source_pools",
    "score_cached_conditioning_increment_batch",
    "score_cached_pool_difference_batch",
    "score_cached_target_conditioned_batch",
    "score_fully_cached_target_conditioned_batch",
    "score_conditioning_increment_batch",
    "score_pool_difference_batch",
    "score_target_conditioned_batch",
]
