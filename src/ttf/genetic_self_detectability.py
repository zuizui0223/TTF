from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import SpeciesEdges, spearman_rho
from .genetic_gate import GeneticTTFDesign
from .genetic_ibd import crossfit_ibd_residuals_prepared
from .genetic_simulate import GeneticSyntheticWorld


@dataclass(frozen=True)
class SelfFieldSpeciesDesign:
    species: str
    projection: np.ndarray
    prior_offset: np.ndarray
    eligible_counts: np.ndarray


@dataclass(frozen=True)
class GeneticSelfDetectabilityDesign:
    species: tuple[str, ...]
    by_species: Mapping[str, SelfFieldSpeciesDesign]
    bandwidth: float
    prior_strength: float
    prior_mean: float
    segment_points: int


@dataclass(frozen=True)
class GeneticSelfBatchScore:
    statistics: np.ndarray
    species_scores: Mapping[str, np.ndarray]
    n_species: np.ndarray


def _prepare_species_self_field(
    edges: SpeciesEdges,
    *,
    bandwidth: float,
    prior_strength: float,
    prior_mean: float,
    segment_points: int,
    min_training_edges: int,
) -> SelfFieldSpeciesDesign:
    if bandwidth <= 0 or prior_strength < 0 or segment_points < 1:
        raise ValueError('invalid self-field configuration')
    if min_training_edges < 1:
        raise ValueError('min_training_edges must be positive')

    nodes = np.asarray(edges.nodes, dtype=np.int64)
    n_edges = edges.n_edges
    projection = np.zeros((n_edges, n_edges), dtype=float)
    prior_offset = np.empty(n_edges, dtype=float)
    eligible_counts = np.empty(n_edges, dtype=np.int64)
    t = (np.arange(segment_points, dtype=float) + 0.5) / segment_points
    h2 = float(bandwidth) ** 2
    tiny = np.finfo(float).tiny

    for target, (left, right) in enumerate(nodes):
        allowed = (
            (nodes[:, 0] != left)
            & (nodes[:, 1] != left)
            & (nodes[:, 0] != right)
            & (nodes[:, 1] != right)
        )
        index = np.flatnonzero(allowed)
        if len(index) < int(min_training_edges):
            raise ValueError(
                f'{edges.species} target edge {target} has only {len(index)} endpoint-disjoint training edges'
            )
        eligible_counts[target] = int(len(index))
        points = edges.start[target][None, :] + t[:, None] * (
            edges.end[target] - edges.start[target]
        )[None, :]
        delta = points[:, None, :] - edges.midpoint[index][None, :, :]
        distance2 = np.sum(delta * delta, axis=2)
        kernel = np.exp(-0.5 * distance2 / h2) / float(len(index))
        denominator = np.maximum(kernel.sum(axis=1) + float(prior_strength), tiny)
        projection[target, index] = np.mean(
            kernel / denominator[:, None],
            axis=0,
        )
        prior_offset[target] = float(
            np.mean(float(prior_strength) * float(prior_mean) / denominator)
        )

    return SelfFieldSpeciesDesign(
        species=edges.species,
        projection=projection,
        prior_offset=prior_offset,
        eligible_counts=eligible_counts,
    )


def prepare_genetic_self_detectability(
    design: GeneticTTFDesign,
    *,
    species: Sequence[str] | None = None,
    bandwidth: float = 500.0,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
) -> GeneticSelfDetectabilityDesign:
    labels = tuple(design.eval_species if species is None else map(str, species))
    if not labels:
        raise ValueError('at least one species is required')
    if any(name not in design.eval_species for name in labels):
        raise ValueError('self-detectability species must belong to the frozen evaluation set')
    prepared = {
        name: _prepare_species_self_field(
            design.template_edges[name],
            bandwidth=float(bandwidth),
            prior_strength=float(prior_strength),
            prior_mean=float(prior_mean),
            segment_points=int(segment_points),
            min_training_edges=int(design.min_training_edges),
        )
        for name in labels
    }
    return GeneticSelfDetectabilityDesign(
        species=labels,
        by_species=prepared,
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=float(prior_mean),
        segment_points=int(segment_points),
    )


def score_genetic_self_world_batch(
    design: GeneticTTFDesign,
    self_design: GeneticSelfDetectabilityDesign,
    worlds: Sequence[GeneticSyntheticWorld],
) -> GeneticSelfBatchScore:
    batch = tuple(worlds)
    if not batch:
        raise ValueError('at least one world is required')
    width = len(batch)
    species_scores = {
        name: np.empty(width, dtype=float) for name in self_design.species
    }

    for name in self_design.species:
        response = np.empty((design.template_edges[name].n_edges, width), dtype=float)
        for column, world in enumerate(batch):
            if name not in world.genetic_distance:
                raise ValueError(f'world {column} is missing species {name}')
            result = crossfit_ibd_residuals_prepared(
                world.genetic_distance[name],
                design.ibd_designs[name],
            )
            response[:, column] = result.residual_turnover
        prepared = self_design.by_species[name]
        predicted = prepared.prior_offset[:, None] + prepared.projection @ response
        species_scores[name][:] = np.asarray(
            [spearman_rho(predicted[:, j], response[:, j]) for j in range(width)],
            dtype=float,
        )

    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    for values in species_scores.values():
        finite = np.isfinite(values)
        sums[finite] += values[finite]
        counts[finite] += 1
    if np.any(counts == 0):
        raise ValueError('one or more worlds had no finite within-species detectability score')
    return GeneticSelfBatchScore(
        statistics=sums / counts,
        species_scores={name: values.copy() for name, values in species_scores.items()},
        n_species=counts,
    )


__all__ = [
    'SelfFieldSpeciesDesign',
    'GeneticSelfDetectabilityDesign',
    'GeneticSelfBatchScore',
    'prepare_genetic_self_detectability',
    'score_genetic_self_world_batch',
]
