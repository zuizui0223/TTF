from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .batch import BatchTransferResult
from .chunked_transfer import score_chunked_batch
from .conditional_transfer import (
    ConditionalSourceSets,
    TargetRestrictedTransfer,
    build_conditional_source_sets,
    paired_macro_increment,
    pairwise_geographic_coverage,
    prepare_target_restricted_transfer,
    score_target_restricted_batch,
)
from .genetic_gate import GeneticTTFDesign, prepare_genetic_ttf_design
from .genetic_geometry import GeneticSamplingGeometry
from .genetic_ibd import crossfit_ibd_residuals_prepared
from .geometry_control import length_orthogonalized_turnover


@dataclass(frozen=True)
class GeneticConditionalDesign:
    base: GeneticTTFDesign
    source_sets: ConditionalSourceSets
    geographic_prepared: TargetRestrictedTransfer
    same_order_prepared: TargetRestrictedTransfer
    taxonomy: Mapping[str, Mapping[str, str]]


@dataclass(frozen=True)
class GeneticConditionalWorldScore:
    global_scores: BatchTransferResult
    geographic_scores: BatchTransferResult
    same_order_scores: BatchTransferResult
    delta_geographic: np.ndarray
    delta_order_given_geography: np.ndarray


def prepare_genetic_conditional_design(
    geometries: Mapping[str, GeneticSamplingGeometry],
    taxonomy: Mapping[str, Mapping[str, str]],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    bandwidth: float = 500.0,
    radius: float = 500.0,
    coverage_threshold: float = 0.25,
    min_sources: int = 3,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    min_training_edges: int = 5,
    strength_neighbours: int = 4,
) -> GeneticConditionalDesign:
    """Freeze conditional transfer geometry without reading genetic outcomes."""
    base = prepare_genetic_ttf_design(
        geometries,
        train_species=train_species,
        eval_species=eval_species,
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        segment_points=int(segment_points),
        min_training_edges=int(min_training_edges),
        strength_neighbours=int(strength_neighbours),
    )
    required = set(base.train_species) | set(base.eval_species)
    missing_taxonomy = required - set(taxonomy)
    if missing_taxonomy:
        raise ValueError(
            f"missing taxonomy: {sorted(missing_taxonomy)}"
        )
    coverage = pairwise_geographic_coverage(
        [
            base.template_edges[name]
            for name in base.train_species
        ],
        [
            base.template_edges[name]
            for name in base.eval_species
        ],
        radius=float(radius),
    )
    source_sets = build_conditional_source_sets(
        coverage,
        taxonomy,
        radius=float(radius),
        coverage_threshold=float(coverage_threshold),
        min_sources=int(min_sources),
    )
    geographic_sources = {
        target: source_sets.geographic[target]
        for target in source_sets.primary_eval_species
    }
    same_order_sources = {
        target: source_sets.same_order[target]
        for target in source_sets.secondary_eval_species
    }
    geographic_prepared = prepare_target_restricted_transfer(
        [
            base.template_edges[name]
            for name in base.train_species
        ],
        [
            base.template_edges[name]
            for name in base.eval_species
        ],
        geographic_sources,
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.0,
        segment_points=int(segment_points),
    )
    same_order_prepared = prepare_target_restricted_transfer(
        [
            base.template_edges[name]
            for name in base.train_species
        ],
        [
            base.template_edges[name]
            for name in base.eval_species
        ],
        same_order_sources,
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.0,
        segment_points=int(segment_points),
    )
    return GeneticConditionalDesign(
        base=base,
        source_sets=source_sets,
        geographic_prepared=geographic_prepared,
        same_order_prepared=same_order_prepared,
        taxonomy={
            name: dict(taxonomy[name])
            for name in sorted(required)
        },
    )


def score_genetic_conditional_world(
    design: GeneticConditionalDesign,
    genetic_distance: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> GeneticConditionalWorldScore:
    """Apply the frozen genetic response to global and conditional source sets."""
    used = design.base.train_species + design.base.eval_species
    missing = set(used) - set(genetic_distance)
    if missing:
        raise ValueError(
            f"missing genetic distances: {sorted(missing)}"
        )

    train_response: dict[str, np.ndarray] = {}
    eval_response: dict[str, np.ndarray] = {}
    for name in used:
        result = crossfit_ibd_residuals_prepared(
            genetic_distance[name],
            design.base.ibd_designs[name],
        )
        if name in design.base.train_species:
            train_response[name] = length_orthogonalized_turnover(
                result.residual_turnover,
                design.base.template_edges[name].length,
            )
        else:
            eval_response[name] = result.residual_turnover

    train_batch = {
        name: values[:, None]
        for name, values in train_response.items()
    }
    eval_batch = {
        name: values[:, None]
        for name, values in eval_response.items()
    }

    global_scores = score_chunked_batch(
        design.base.prepared,
        train_batch,
        eval_batch,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    geographic_scores = score_target_restricted_batch(
        design.geographic_prepared,
        train_batch,
        eval_batch,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    same_order_scores = score_target_restricted_batch(
        design.same_order_prepared,
        train_batch,
        eval_batch,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )

    delta_geographic = paired_macro_increment(
        geographic_scores,
        global_scores,
        design.source_sets.primary_eval_species,
    )
    delta_order_given_geography = paired_macro_increment(
        same_order_scores,
        geographic_scores,
        design.source_sets.secondary_eval_species,
    )
    return GeneticConditionalWorldScore(
        global_scores=global_scores,
        geographic_scores=geographic_scores,
        same_order_scores=same_order_scores,
        delta_geographic=delta_geographic,
        delta_order_given_geography=delta_order_given_geography,
    )


__all__ = [
    "GeneticConditionalDesign",
    "GeneticConditionalWorldScore",
    "prepare_genetic_conditional_design",
    "score_genetic_conditional_world",
]
