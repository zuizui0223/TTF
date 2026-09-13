from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .calibration import seed_for
from .chunked_transfer import (
    ChunkedTransferGeometry,
    prepare_chunked_transfer,
    score_chunked_batch,
)
from .core import SpeciesEdges, split_species
from .genetic_geometry import GeneticSamplingGeometry
from .genetic_ibd import CrossfitIBDResult, crossfit_ibd_residuals
from .genetic_simulate import GeneticSyntheticWorld, simulate_genetic_distance_world
from .geometry_control import length_orthogonalized_turnover
from .private_strength import (
    edge_midpoint_neighbor_indices,
    training_private_strength_from_indices,
)
from .profiled_private_null import profiled_private_pvalue


GENETIC_PRIVATE_REFERENCE_CONFIG: dict[str, tuple[float, float]] = {
    "A0": (0.0, 0.10),
    "A0p5": (0.5, 0.10),
    "A1": (1.0, 0.10),
    "A2": (2.0, 0.10),
    "A3": (3.0, 0.10),
    "A5": (5.0, 0.10),
    "A10": (10.0, 0.10),
    "infinite_snr": (1.0, 0.0),
}


@dataclass(frozen=True)
class GeneticTTFDesign:
    """Frozen outcome-blind design for genetic TTF qualification."""

    geometries: Mapping[str, GeneticSamplingGeometry]
    template_edges: Mapping[str, SpeciesEdges]
    train_species: tuple[str, ...]
    eval_species: tuple[str, ...]
    strength_indices: Mapping[str, np.ndarray]
    prepared: ChunkedTransferGeometry
    min_training_edges: int


@dataclass(frozen=True)
class GeneticWorldScore:
    statistic: float
    training_strength: float
    species_scores: Mapping[str, float]
    ibd: Mapping[str, CrossfitIBDResult]


@dataclass(frozen=True)
class GeneticGateCell:
    shared_fraction: float
    residual_amplitude: float
    n_worlds: int
    alpha: float
    rejection_rate: float
    mean_statistic: float
    median_p_value: float
    selected_pair_counts: Mapping[str, int]


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


def prepare_genetic_ttf_design(
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
) -> GeneticTTFDesign:
    """Freeze the response-blind genetic transfer operator and nuisance design."""
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
    strength_index = {
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
    return GeneticTTFDesign(
        geometries={name: geometries[name] for name in labels},
        template_edges=templates,
        train_species=train,
        eval_species=evaluation,
        strength_indices=strength_index,
        prepared=prepared,
        min_training_edges=int(min_training_edges),
    )


def score_genetic_world(
    design: GeneticTTFDesign,
    world: GeneticSyntheticWorld,
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> GeneticWorldScore:
    """Apply endpoint-safe IBD residualization then current TTF geometry control."""
    used = design.train_species + design.eval_species
    missing = set(used) - set(world.genetic_distance)
    if missing:
        raise ValueError(f"synthetic world is missing species: {sorted(missing)}")

    ibd: dict[str, CrossfitIBDResult] = {}
    train_response: dict[str, np.ndarray] = {}
    eval_response: dict[str, np.ndarray] = {}

    for name in used:
        edges = design.template_edges[name]
        result = crossfit_ibd_residuals(
            world.genetic_distance[name],
            edges.length,
            edges.nodes,
            min_training_edges=design.min_training_edges,
        )
        ibd[name] = result
        if name in design.train_species:
            train_response[name] = length_orthogonalized_turnover(
                result.residual_turnover,
                edges.length,
            )
        else:
            eval_response[name] = result.residual_turnover

    strength, _ = training_private_strength_from_indices(
        train_response,
        design.strength_indices,
        design.train_species,
    )
    scored = score_chunked_batch(
        design.prepared,
        {name: train_response[name][:, None] for name in design.train_species},
        {name: eval_response[name][:, None] for name in design.eval_species},
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    species_scores = {
        name: float(scored.species_scores[name][0])
        for name in design.eval_species
    }
    return GeneticWorldScore(
        statistic=float(scored.statistics[0]),
        training_strength=float(strength),
        species_scores=species_scores,
        ibd=ibd,
    )


def build_genetic_private_references(
    design: GeneticTTFDesign,
    *,
    configurations: Mapping[str, tuple[float, float]] = GENETIC_PRIVATE_REFERENCE_CONFIG,
    n_worlds: int = 1999,
    ibd_strength: float | Mapping[str, float] = 1.0,
    transition_width: float = 0.20,
    noise_dimensions: int = 2,
    seed: int = 20260913,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Generate independent private-null reference families for profiled inference."""
    if n_worlds < 21:
        raise ValueError("n_worlds must be at least 21")
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for label, (amplitude, noise_sd) in configurations.items():
        strengths = np.empty(int(n_worlds), dtype=float)
        statistics = np.empty(int(n_worlds), dtype=float)
        for replicate in range(int(n_worlds)):
            world = simulate_genetic_distance_world(
                design.geometries,
                shared_fraction=0.0,
                residual_amplitude=float(amplitude),
                ibd_strength=ibd_strength,
                noise_sd=float(noise_sd),
                transition_width=float(transition_width),
                noise_dimensions=int(noise_dimensions),
                seed=seed_for(seed, "genetic_private_reference", label, replicate),
            )
            score = score_genetic_world(
                design,
                world,
                edge_chunk_size=edge_chunk_size,
                train_chunk_size=train_chunk_size,
            )
            strengths[replicate] = score.training_strength
            statistics[replicate] = score.statistic
        out[str(label)] = (strengths, statistics)
    return out


def run_genetic_profiled_gate_cell(
    design: GeneticTTFDesign,
    references: Mapping[str, tuple[np.ndarray, np.ndarray]],
    *,
    shared_fraction: float,
    residual_amplitude: float,
    n_worlds: int,
    ibd_strength: float | Mapping[str, float] = 1.0,
    noise_sd: float = 0.10,
    transition_width: float = 0.20,
    noise_dimensions: int = 2,
    profile_draws: int = 999,
    selected_configs: int = 2,
    alpha: float = 0.05,
    seed: int = 20260913,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> GeneticGateCell:
    """Evaluate one genetic Gate-D cell against independent private references."""
    if n_worlds < 1:
        raise ValueError("n_worlds must be positive")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1)")

    p_values = np.empty(int(n_worlds), dtype=float)
    statistics = np.empty(int(n_worlds), dtype=float)
    selected_pairs: Counter[str] = Counter()

    for replicate in range(int(n_worlds)):
        world = simulate_genetic_distance_world(
            design.geometries,
            shared_fraction=float(shared_fraction),
            residual_amplitude=float(residual_amplitude),
            ibd_strength=ibd_strength,
            noise_sd=float(noise_sd),
            transition_width=float(transition_width),
            noise_dimensions=int(noise_dimensions),
            seed=seed_for(
                seed,
                "genetic_observed",
                shared_fraction,
                residual_amplitude,
                replicate,
            ),
        )
        score = score_genetic_world(
            design,
            world,
            edge_chunk_size=edge_chunk_size,
            train_chunk_size=train_chunk_size,
        )
        p, selected, _, _ = profiled_private_pvalue(
            score.statistic,
            score.training_strength,
            references,
            profile_draws=int(profile_draws),
            selected_configs=int(selected_configs),
        )
        p_values[replicate] = p
        statistics[replicate] = score.statistic
        selected_pairs["+".join(selected)] += 1

    return GeneticGateCell(
        shared_fraction=float(shared_fraction),
        residual_amplitude=float(residual_amplitude),
        n_worlds=int(n_worlds),
        alpha=float(alpha),
        rejection_rate=float(np.mean(p_values <= float(alpha))),
        mean_statistic=float(np.mean(statistics)),
        median_p_value=float(np.median(p_values)),
        selected_pair_counts=dict(sorted(selected_pairs.items())),
    )
