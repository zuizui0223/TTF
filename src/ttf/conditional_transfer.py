from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .batch import BatchTransferResult
from .chunked_transfer import (
    ChunkedTransferGeometry,
    prepare_chunked_transfer,
    score_chunked_batch,
)
from .core import SpeciesEdges


@dataclass(frozen=True)
class ConditionalSourceSets:
    """Response-blind source-species sets for target-specific transfer.

    Membership is determined only from frozen sampling geometry and taxonomy.
    No turnover, genetic distance, sequence identity, or held-out response enters
    this object.
    """

    geographic: Mapping[str, tuple[str, ...]]
    same_order: Mapping[str, tuple[str, ...]]
    primary_eval_species: tuple[str, ...]
    secondary_eval_species: tuple[str, ...]
    coverage: Mapping[str, Mapping[str, float]]
    radius: float
    coverage_threshold: float
    min_sources: int


@dataclass(frozen=True)
class TargetRestrictedTransfer:
    eval_species: tuple[str, ...]
    train_species_by_target: Mapping[str, tuple[str, ...]]
    prepared_by_target: Mapping[str, ChunkedTransferGeometry]


def pairwise_geographic_coverage(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    *,
    radius: float = 500.0,
) -> dict[str, dict[str, float]]:
    """Directed species co-distribution from frozen edge-midpoint geometry.

    For target species `t` and source species `s`, the value is the fraction
    of target edge midpoints lying within `radius` of at least one source edge
    midpoint. This is deliberately target-directed: the question is how much of
    the target's sampled spatial graph is represented by the source species.
    """
    if radius <= 0:
        raise ValueError("radius must be positive")
    train = tuple(train_edges)
    evaluation = tuple(eval_edges)
    if not train or not evaluation:
        raise ValueError("non-empty train and evaluation edge sets are required")
    train_names = [x.species for x in train]
    eval_names = [x.species for x in evaluation]
    if len(set(train_names)) != len(train_names) or len(set(eval_names)) != len(eval_names):
        raise ValueError("species labels must be unique within each split")
    if set(train_names) & set(eval_names):
        raise ValueError("train and evaluation species must be disjoint")
    r2 = float(radius) * float(radius)
    out: dict[str, dict[str, float]] = {}
    for target in evaluation:
        tmid = np.asarray(target.midpoint, dtype=float)
        if tmid.ndim != 2 or len(tmid) == 0 or not np.isfinite(tmid).all():
            raise ValueError(f"invalid target geometry for {target.species}")
        row: dict[str, float] = {}
        for source in train:
            smid = np.asarray(source.midpoint, dtype=float)
            if (
                smid.ndim != 2
                or len(smid) == 0
                or smid.shape[1] != tmid.shape[1]
                or not np.isfinite(smid).all()
            ):
                raise ValueError(f"invalid source geometry for {source.species}")
            nearest2 = np.min(
                np.sum((tmid[:, None, :] - smid[None, :, :]) ** 2, axis=2),
                axis=1,
            )
            row[source.species] = float(np.mean(nearest2 <= r2))
        out[target.species] = row
    return out


def build_conditional_source_sets(
    coverage: Mapping[str, Mapping[str, float]],
    taxonomy: Mapping[str, Mapping[str, str]],
    *,
    radius: float = 500.0,
    coverage_threshold: float = 0.25,
    min_sources: int = 3,
) -> ConditionalSourceSets:
    """Freeze geographic and geography+order source sets before outcomes.

    `geographic` retains source species covering at least the declared fraction
    of target edge midpoints. `same_order` is the subset sharing the target's
    non-empty order label. Target eligibility is based only on source-set size.
    """
    if radius <= 0:
        raise ValueError("radius must be positive")
    if not 0.0 < coverage_threshold <= 1.0:
        raise ValueError("coverage_threshold must lie in (0, 1]")
    if min_sources < 1:
        raise ValueError("min_sources must be positive")
    if not coverage:
        raise ValueError("coverage must be non-empty")
    geographic: dict[str, tuple[str, ...]] = {}
    same_order: dict[str, tuple[str, ...]] = {}
    primary: list[str] = []
    secondary: list[str] = []
    normalized: dict[str, dict[str, float]] = {}
    for target in sorted(coverage):
        if target not in taxonomy:
            raise ValueError(f"missing taxonomy for target {target}")
        row: dict[str, float] = {}
        for source, raw in coverage[target].items():
            if source == target:
                raise ValueError("source and target species must be disjoint")
            if source not in taxonomy:
                raise ValueError(f"missing taxonomy for source {source}")
            value = float(raw)
            if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError("coverage values must lie in [0, 1]")
            row[str(source)] = value
        if not row:
            raise ValueError(f"target {target} has no source species")
        geo = tuple(
            sorted(
                source
                for source, value in row.items()
                if value >= coverage_threshold
            )
        )
        order = str(taxonomy[target].get("order", "")).strip()
        order_sources = tuple(
            source
            for source in geo
            if order and str(taxonomy[source].get("order", "")).strip() == order
        )
        geographic[target] = geo
        same_order[target] = order_sources
        normalized[target] = dict(sorted(row.items()))
        if len(geo) >= min_sources:
            primary.append(target)
        if len(order_sources) >= min_sources:
            secondary.append(target)
    return ConditionalSourceSets(
        geographic=geographic,
        same_order=same_order,
        primary_eval_species=tuple(primary),
        secondary_eval_species=tuple(secondary),
        coverage=normalized,
        radius=float(radius),
        coverage_threshold=float(coverage_threshold),
        min_sources=int(min_sources),
    )


def prepare_target_restricted_transfer(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    source_sets: Mapping[str, Sequence[str]],
    *,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.0,
    segment_points: int = 5,
) -> TargetRestrictedTransfer:
    """Prepare the unchanged TTF kernel on a frozen source set per target."""
    train = tuple(train_edges)
    evaluation = tuple(eval_edges)
    train_map = {x.species: x for x in train}
    eval_map = {x.species: x for x in evaluation}
    if len(train_map) != len(train) or len(eval_map) != len(evaluation):
        raise ValueError("species labels must be unique")
    if set(train_map) & set(eval_map):
        raise ValueError("train and evaluation species must be disjoint")
    prepared: dict[str, ChunkedTransferGeometry] = {}
    frozen_sources: dict[str, tuple[str, ...]] = {}
    for target in sorted(source_sets):
        if target not in eval_map:
            raise ValueError(f"unknown evaluation species {target}")
        names = tuple(map(str, source_sets[target]))
        if not names or len(set(names)) != len(names):
            raise ValueError(
                f"source set for {target} must be non-empty and unique"
            )
        missing = set(names) - set(train_map)
        if missing:
            raise ValueError(
                f"unknown source species for {target}: {sorted(missing)}"
            )
        frozen_sources[target] = names
        prepared[target] = prepare_chunked_transfer(
            [train_map[name] for name in names],
            [eval_map[target]],
            bandwidth=float(bandwidth),
            prior_strength=float(prior_strength),
            prior_mean=float(prior_mean),
            segment_points=int(segment_points),
        )
    return TargetRestrictedTransfer(
        eval_species=tuple(sorted(prepared)),
        train_species_by_target=frozen_sources,
        prepared_by_target=prepared,
    )


def score_target_restricted_batch(
    prepared: TargetRestrictedTransfer,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> BatchTransferResult:
    """Score target-specific frozen source sets with the unchanged TTF kernel."""
    if not prepared.eval_species:
        raise ValueError("no evaluation species are prepared")
    score_map: dict[str, np.ndarray] = {}
    widths: set[int] = set()
    for target in prepared.eval_species:
        source_names = prepared.train_species_by_target[target]
        result = score_chunked_batch(
            prepared.prepared_by_target[target],
            {name: train_turnover[name] for name in source_names},
            {target: eval_turnover[target]},
            edge_chunk_size=int(edge_chunk_size),
            train_chunk_size=int(train_chunk_size),
        )
        values = np.asarray(result.species_scores[target], dtype=float)
        widths.add(int(len(values)))
        score_map[target] = values
    if len(widths) != 1:
        raise ValueError("target batch widths disagree")
    width = widths.pop()
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    for values in score_map.values():
        finite = np.isfinite(values)
        sums[finite] += values[finite]
        counts[finite] += 1
    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite target scores")
    return BatchTransferResult(
        statistics=sums / counts,
        species_scores=score_map,
        n_eval_species=counts,
    )


def paired_macro_increment(
    conditioned: BatchTransferResult,
    baseline: BatchTransferResult,
    species: Sequence[str],
) -> np.ndarray:
    """Species-equal paired improvement on a prospectively fixed target set."""
    names = tuple(map(str, species))
    if not names or len(set(names)) != len(names):
        raise ValueError("species must be non-empty and unique")
    differences: list[np.ndarray] = []
    width: int | None = None
    for name in names:
        if (
            name not in conditioned.species_scores
            or name not in baseline.species_scores
        ):
            raise ValueError(f"missing paired species score for {name}")
        conditioned_score = np.asarray(
            conditioned.species_scores[name], dtype=float
        )
        baseline_score = np.asarray(baseline.species_scores[name], dtype=float)
        if (
            conditioned_score.ndim != 1
            or conditioned_score.shape != baseline_score.shape
            or not np.isfinite(conditioned_score).all()
            or not np.isfinite(baseline_score).all()
        ):
            raise ValueError(f"invalid paired scores for {name}")
        if width is None:
            width = len(conditioned_score)
        elif len(conditioned_score) != width:
            raise ValueError("paired score batch widths disagree")
        differences.append(conditioned_score - baseline_score)
    return np.mean(np.vstack(differences), axis=0)


__all__ = [
    "ConditionalSourceSets",
    "TargetRestrictedTransfer",
    "pairwise_geographic_coverage",
    "build_conditional_source_sets",
    "prepare_target_restricted_transfer",
    "score_target_restricted_batch",
    "paired_macro_increment",
]
