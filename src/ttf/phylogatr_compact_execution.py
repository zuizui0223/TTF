from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .cached_chunked_transfer import (
    CachedChunkedTransferGeometry,
    prepare_cached_chunked_transfer,
    score_cached_chunked_batch,
)
from .chunked_transfer import ChunkedTransferGeometry, prepare_chunked_transfer
from .core import SpeciesEdges, split_species
from .genetic_batch_execution import GeneticBatchScore
from .genetic_geometry import GeneticSamplingGeometry
from .genetic_simulate import GeneticSyntheticWorld
from .geometry_control import length_orthogonalized_turnover
from .phylogatr_compact_ibd import (
    CompactCrossfitIBDDesign,
    crossfit_ibd_residuals_compact,
    prepare_compact_crossfit_ibd_design,
)
from .private_strength import (
    edge_midpoint_neighbor_indices,
    training_private_strength_from_indices,
)


@dataclass(frozen=True)
class PhylogatrCompactTTFDesign:
    """Fresh-phylogatR execution design with exact compact endpoint-safe IBD."""

    geometries: Mapping[str, GeneticSamplingGeometry]
    template_edges: Mapping[str, SpeciesEdges]
    ibd_designs: Mapping[str, CompactCrossfitIBDDesign]
    train_species: tuple[str, ...]
    eval_species: tuple[str, ...]
    strength_indices: Mapping[str, np.ndarray]
    prepared: ChunkedTransferGeometry
    min_training_edges: int


def _template_species_edges(
    species: str,
    geometry: GeneticSamplingGeometry,
) -> SpeciesEdges:
    nodes = np.asarray(geometry.edge_nodes, dtype=np.int64)
    start = geometry.coordinates[nodes[:, 0]]
    end = geometry.coordinates[nodes[:, 1]]
    return SpeciesEdges(
        species=str(species),
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=np.zeros(len(nodes), dtype=float),
    )


def prepare_phylogatr_compact_ttf_design(
    geometries: Mapping[str, GeneticSamplingGeometry],
    *,
    train_species: Sequence[str] | None = None,
    eval_species: Sequence[str] | None = None,
    eval_fraction: float = 0.5,
    split_seed: int = 0,
    bandwidth: float,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    min_training_edges: int = 5,
    strength_neighbours: int = 4,
) -> PhylogatrCompactTTFDesign:
    """Mirror the frozen generic genetic design, changing only IBD storage/execution."""
    if not geometries:
        raise ValueError("at least one genetic geometry is required")
    if min_training_edges < 3:
        raise ValueError("min_training_edges must be >= 3")
    labels = tuple(sorted(map(str, geometries.keys())))
    if len(labels) != len(geometries):
        raise ValueError("species labels must be unique")
    for name in labels:
        if geometries[name].min_endpoint_disjoint_training_edges < int(min_training_edges):
            raise ValueError(
                f"{name} is not endpoint-safe for IBD cross-fitting: "
                f"minimum={geometries[name].min_endpoint_disjoint_training_edges}, "
                f"required={min_training_edges}"
            )

    if (train_species is None) != (eval_species is None):
        raise ValueError("train_species and eval_species must be supplied together")
    if train_species is None:
        train, evaluation = split_species(
            labels,
            eval_fraction=eval_fraction,
            seed=int(split_seed),
        )
    else:
        train = tuple(map(str, train_species))
        evaluation = tuple(map(str, eval_species or ()))
        if not train or not evaluation or set(train) & set(evaluation):
            raise ValueError("non-empty species-disjoint train/evaluation sets required")
        missing = (set(train) | set(evaluation)) - set(labels)
        if missing:
            raise ValueError(f"unknown species in split: {sorted(missing)}")
    if len(evaluation) < 6:
        raise ValueError("genetic qualification requires at least six held-out species")

    templates = {
        name: _template_species_edges(name, geometries[name])
        for name in train + evaluation
    }
    ibd_designs = {
        name: prepare_compact_crossfit_ibd_design(
            templates[name].length,
            templates[name].nodes,
            min_training_edges=int(min_training_edges),
        )
        for name in train + evaluation
    }
    strength_indices = {
        name: edge_midpoint_neighbor_indices(
            templates[name].midpoint,
            k=int(strength_neighbours),
        )
        for name in train
    }
    prepared = prepare_chunked_transfer(
        [templates[name] for name in train],
        [templates[name] for name in evaluation],
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.0,
        segment_points=int(segment_points),
    )
    return PhylogatrCompactTTFDesign(
        geometries={name: geometries[name] for name in labels},
        template_edges=templates,
        ibd_designs=ibd_designs,
        train_species=train,
        eval_species=evaluation,
        strength_indices=strength_indices,
        prepared=prepared,
        min_training_edges=int(min_training_edges),
    )


def prepare_phylogatr_compact_cached_transfer(
    design: PhylogatrCompactTTFDesign,
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> CachedChunkedTransferGeometry:
    return prepare_cached_chunked_transfer(
        design.prepared,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )


def score_phylogatr_compact_world_batch(
    design: PhylogatrCompactTTFDesign,
    worlds: Sequence[GeneticSyntheticWorld],
    cached_transfer: CachedChunkedTransferGeometry,
) -> GeneticBatchScore:
    """Mirror generic batch scoring with the exact compact IBD scorer only."""
    batch = tuple(worlds)
    if not batch:
        raise ValueError("at least one world is required")
    width = len(batch)
    used = design.train_species + design.eval_species
    for index, world in enumerate(batch):
        missing = set(used) - set(world.genetic_distance)
        if missing:
            raise ValueError(f"world {index} is missing species: {sorted(missing)}")

    train_response = {
        name: np.empty((design.template_edges[name].n_edges, width), dtype=float)
        for name in design.train_species
    }
    eval_response = {
        name: np.empty((design.template_edges[name].n_edges, width), dtype=float)
        for name in design.eval_species
    }
    strengths = np.empty(width, dtype=float)

    for column, world in enumerate(batch):
        train_column: dict[str, np.ndarray] = {}
        for name in used:
            result = crossfit_ibd_residuals_compact(
                world.genetic_distance[name],
                design.ibd_designs[name],
            )
            if name in design.train_species:
                response = length_orthogonalized_turnover(
                    result.residual_turnover,
                    design.template_edges[name].length,
                )
                train_response[name][:, column] = response
                train_column[name] = response
            else:
                eval_response[name][:, column] = result.residual_turnover
        strength, _ = training_private_strength_from_indices(
            train_column,
            design.strength_indices,
            design.train_species,
        )
        strengths[column] = float(strength)

    scored = score_cached_chunked_batch(
        cached_transfer,
        train_response,
        eval_response,
    )
    return GeneticBatchScore(
        statistics=np.asarray(scored.statistics, dtype=float).copy(),
        training_strengths=strengths,
        species_scores={
            name: np.asarray(scored.species_scores[name], dtype=float).copy()
            for name in design.eval_species
        },
        n_eval_species=np.asarray(scored.n_eval_species, dtype=np.int64).copy(),
    )


__all__ = [
    "PhylogatrCompactTTFDesign",
    "prepare_phylogatr_compact_cached_transfer",
    "prepare_phylogatr_compact_ttf_design",
    "score_phylogatr_compact_world_batch",
]
