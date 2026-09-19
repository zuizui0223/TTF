from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Mapping, Sequence

import numpy as np

from .conditional_transfer import (
    CachedTargetConditionedTransferGeometry,
    FullyCachedTargetConditionedTransferGeometry,
    MatchedTargetSourcePoolDesign,
    prepare_cached_target_conditioned_transfer,
    prepare_fully_cached_target_conditioned_transfer,
    prepare_geometry_matched_group_source_pools,
    score_fully_cached_target_conditioned_batch,
)
from .genetic_conditional_simulate import (
    GroupedGeneticSyntheticWorld,
    simulate_grouped_genetic_distance_world,
)
from .genetic_geometry import GeneticSamplingGeometry
from .geometry_control import length_orthogonalized_turnover
from .phylogatr_compact_execution import (
    PhylogatrCompactTTFDesign,
    prepare_phylogatr_compact_ttf_design,
)
from .phylogatr_compact_ibd import crossfit_ibd_residuals_compact
from .private_null_inference import envelope_upper_pvalues


DEV_REFERENCE_SEED_TAG = "conditional-order-geometry-match-v03-dev-reference"
DEV_EVALUATION_SEED_TAG = "conditional-order-geometry-match-v03-dev-evaluation"


@dataclass(frozen=True)
class ConditionalGeometryMatchedOrderDesign:
    compact: PhylogatrCompactTTFDesign
    matched_pools: MatchedTargetSourcePoolDesign
    same_order_cache: CachedTargetConditionedTransferGeometry
    different_order_cache: CachedTargetConditionedTransferGeometry
    order_by_species: Mapping[str, str]
    eligible_eval_species: tuple[str, ...]


@dataclass(frozen=True)
class ConditionalGeometryMatchedOrderExecution:
    same_order: FullyCachedTargetConditionedTransferGeometry
    different_order: FullyCachedTargetConditionedTransferGeometry
    eligible_eval_species: tuple[str, ...]


@dataclass(frozen=True)
class ConditionalGeometryMatchedOrderBatchScore:
    statistics: np.ndarray
    species_increments: Mapping[str, np.ndarray]
    same_order_species_scores: Mapping[str, np.ndarray]
    different_order_species_scores: Mapping[str, np.ndarray]
    n_eval_species: int


@dataclass(frozen=True)
class ConditionalPrivateEnvelopeScore:
    statistics: np.ndarray
    p_values: np.ndarray
    least_favourable_private_cell: tuple[str, ...]


def frozen_v03_development_seed(
    master_seed: int,
    tag: str,
    cell: str,
    replicate: int,
) -> int:
    if tag not in {DEV_REFERENCE_SEED_TAG, DEV_EVALUATION_SEED_TAG}:
        raise ValueError("v0.3 development seed tag is not authorized")
    payload = f"{int(master_seed)}|{tag}|{cell}|{int(replicate)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big", signed=False)


def prepare_conditional_geometry_matched_order(
    geometries: Mapping[str, GeneticSamplingGeometry],
    taxonomy_order: Mapping[str, str],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    bandwidth: float = 500.0,
    support_radius: float = 500.0,
    minimum_target_coverage: float = 0.50,
    minimum_source_species: int = 5,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    min_training_edges: int = 5,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> ConditionalGeometryMatchedOrderDesign:
    """Prepare the v0.3 development estimator using matched source geometry."""
    labels = set(map(str, geometries))
    if set(map(str, taxonomy_order)) != labels:
        raise ValueError("taxonomy order labels must cover every geometry species exactly")
    order_by_species = {
        name: str(taxonomy_order[name]).strip()
        for name in sorted(labels)
    }
    compact = prepare_phylogatr_compact_ttf_design(
        geometries,
        train_species=tuple(map(str, train_species)),
        eval_species=tuple(map(str, eval_species)),
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        segment_points=int(segment_points),
        min_training_edges=int(min_training_edges),
    )
    train_midpoint = {
        name: compact.template_edges[name].midpoint
        for name in compact.train_species
    }
    eval_midpoint = {
        name: compact.template_edges[name].midpoint
        for name in compact.eval_species
    }
    matched = prepare_geometry_matched_group_source_pools(
        compact.prepared,
        train_midpoint,
        eval_midpoint,
        train_group={name: order_by_species[name] for name in compact.train_species},
        eval_group={name: order_by_species[name] for name in compact.eval_species},
        train_locality_count={
            name: geometries[name].n_localities for name in compact.train_species
        },
        eval_locality_count={
            name: geometries[name].n_localities for name in compact.eval_species
        },
        support_radius=float(support_radius),
        minimum_target_coverage=float(minimum_target_coverage),
        minimum_source_species=int(minimum_source_species),
    )
    if len(matched.eligible_eval_species) < 6:
        raise ValueError("fewer than six targets survive geometry matching")

    same_cache = prepare_cached_target_conditioned_transfer(
        compact.prepared,
        matched.same_group_pools,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    different_cache = prepare_cached_target_conditioned_transfer(
        compact.prepared,
        matched.different_group_pools,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    return ConditionalGeometryMatchedOrderDesign(
        compact=compact,
        matched_pools=matched,
        same_order_cache=same_cache,
        different_order_cache=different_cache,
        order_by_species=order_by_species,
        eligible_eval_species=tuple(matched.eligible_eval_species),
    )


def prepare_conditional_geometry_matched_order_execution(
    design: ConditionalGeometryMatchedOrderDesign,
) -> ConditionalGeometryMatchedOrderExecution:
    eligible = tuple(design.eligible_eval_species)
    same = prepare_fully_cached_target_conditioned_transfer(
        design.same_order_cache,
        eval_species=eligible,
    )
    different = prepare_fully_cached_target_conditioned_transfer(
        design.different_order_cache,
        eval_species=eligible,
    )
    if same.eval_species != eligible or different.eval_species != eligible:
        raise RuntimeError("geometry-matched execution target-set drift")
    return ConditionalGeometryMatchedOrderExecution(
        same_order=same,
        different_order=different,
        eligible_eval_species=eligible,
    )


def _responses_for_world_batch(
    design: ConditionalGeometryMatchedOrderDesign,
    worlds: Sequence[GroupedGeneticSyntheticWorld],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    batch = tuple(worlds)
    if not batch:
        raise ValueError("at least one synthetic world is required")
    compact = design.compact
    used = compact.train_species + compact.eval_species
    width = len(batch)
    train_response = {
        name: np.empty((compact.template_edges[name].n_edges, width), dtype=float)
        for name in compact.train_species
    }
    eval_response = {
        name: np.empty((compact.template_edges[name].n_edges, width), dtype=float)
        for name in compact.eval_species
    }
    for column, world in enumerate(batch):
        missing = set(used) - set(world.genetic_distance)
        if missing:
            raise ValueError(
                f"synthetic world {column} is missing species: {sorted(missing)}"
            )
        for name in used:
            result = crossfit_ibd_residuals_compact(
                world.genetic_distance[name],
                compact.ibd_designs[name],
            )
            if name in compact.train_species:
                train_response[name][:, column] = length_orthogonalized_turnover(
                    result.residual_turnover,
                    compact.template_edges[name].length,
                )
            else:
                eval_response[name][:, column] = result.residual_turnover
    return train_response, eval_response


def make_conditional_geometry_matched_order_worlds(
    design: ConditionalGeometryMatchedOrderDesign,
    *,
    cell: str,
    residual_amplitude: float,
    absolute_start: int,
    count: int,
    group_mode: str,
    seed_tag: str = DEV_EVALUATION_SEED_TAG,
    ibd_strength: float = 1.0,
    noise_sd: float = 0.10,
    transition_width: float = 0.20,
    noise_dimensions: int = 2,
    master_seed: int = 20260919,
) -> tuple[GroupedGeneticSyntheticWorld, ...]:
    """Generate development-only worlds under a disjoint frozen seed namespace."""
    if count < 1 or absolute_start < 0:
        raise ValueError("count must be positive and absolute_start non-negative")
    if group_mode not in {"private", "same_order"}:
        raise ValueError("group_mode must be private or same_order")
    if seed_tag not in {DEV_REFERENCE_SEED_TAG, DEV_EVALUATION_SEED_TAG}:
        raise ValueError("unfrozen v0.3 development seed namespace")
    if group_mode == "private":
        groups = {name: f"private:{name}" for name in design.compact.geometries}
    else:
        groups = {
            name: (
                f"order:{design.order_by_species[name]}"
                if design.order_by_species[name]
                else f"private:{name}"
            )
            for name in design.compact.geometries
        }
    return tuple(
        simulate_grouped_genetic_distance_world(
            design.compact.geometries,
            groups,
            residual_amplitude=float(residual_amplitude),
            ibd_strength=float(ibd_strength),
            noise_sd=float(noise_sd),
            transition_width=float(transition_width),
            noise_dimensions=int(noise_dimensions),
            seed=frozen_v03_development_seed(
                master_seed,
                seed_tag,
                str(cell),
                int(absolute_start) + offset,
            ),
        )
        for offset in range(int(count))
    )


def score_conditional_geometry_matched_order_world_batch(
    design: ConditionalGeometryMatchedOrderDesign,
    execution: ConditionalGeometryMatchedOrderExecution,
    worlds: Sequence[GroupedGeneticSyntheticWorld],
) -> ConditionalGeometryMatchedOrderBatchScore:
    batch = tuple(worlds)
    train_response, eval_response = _responses_for_world_batch(design, batch)
    same = score_fully_cached_target_conditioned_batch(
        execution.same_order,
        train_response,
        eval_response,
    )
    different = score_fully_cached_target_conditioned_batch(
        execution.different_order,
        train_response,
        eval_response,
    )
    eligible = tuple(design.eligible_eval_species)
    if execution.eligible_eval_species != eligible:
        raise RuntimeError("geometry-matched execution target-set drift")
    increments = {
        name: (
            np.asarray(same.species_scores[name], dtype=float)
            - np.asarray(different.species_scores[name], dtype=float)
        )
        for name in eligible
    }
    matrix = np.vstack([increments[name] for name in eligible])
    statistics = np.mean(matrix, axis=0)
    return ConditionalGeometryMatchedOrderBatchScore(
        statistics=np.asarray(statistics, dtype=float).copy(),
        species_increments={
            name: np.asarray(values, dtype=float).copy()
            for name, values in increments.items()
        },
        same_order_species_scores={
            name: np.asarray(same.species_scores[name], dtype=float).copy()
            for name in eligible
        },
        different_order_species_scores={
            name: np.asarray(different.species_scores[name], dtype=float).copy()
            for name in eligible
        },
        n_eval_species=len(eligible),
    )


def calibrate_conditional_private_envelope(
    statistics: np.ndarray,
    private_references: Mapping[str, np.ndarray],
) -> ConditionalPrivateEnvelopeScore:
    """Calibrate v0.3 statistics against independent exact-geometry private nulls."""
    values = np.asarray(statistics, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("statistics must be a finite 1D array")
    p_values, labels = envelope_upper_pvalues(values, private_references)
    return ConditionalPrivateEnvelopeScore(
        statistics=values.copy(),
        p_values=np.asarray(p_values, dtype=float).copy(),
        least_favourable_private_cell=tuple(labels),
    )


__all__ = [
    "DEV_EVALUATION_SEED_TAG",
    "DEV_REFERENCE_SEED_TAG",
    "ConditionalGeometryMatchedOrderBatchScore",
    "ConditionalGeometryMatchedOrderDesign",
    "ConditionalGeometryMatchedOrderExecution",
    "ConditionalPrivateEnvelopeScore",
    "calibrate_conditional_private_envelope",
    "frozen_v03_development_seed",
    "make_conditional_geometry_matched_order_worlds",
    "prepare_conditional_geometry_matched_order",
    "prepare_conditional_geometry_matched_order_execution",
    "score_conditional_geometry_matched_order_world_batch",
]
