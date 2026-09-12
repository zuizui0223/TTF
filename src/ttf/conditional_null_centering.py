from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from .core import SpeciesEdges, knn_edges, rank01
from .geometry_control import length_orthogonalized_turnover
from .heterogeneous_inference import density_scaled_k
from .mismatch import PairedSpeciesSample, build_coupling_edges


@dataclass(frozen=True)
class ConditionalNullCenteringResult:
    """Frozen Q5.1 TTF-C scores over private-relation phase redraws."""

    statistics: np.ndarray
    species_scores: Mapping[str, np.ndarray]
    phases: Mapping[str, np.ndarray]
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def identical_uniform_circle_samples(
    *,
    seed: int,
    n_systems: int = 40,
    n_points: int = 60,
) -> tuple[PairedSpeciesSample, ...]:
    """Geometry-only homogeneous negative control.

    Every system receives the same evenly spaced circular coordinates after one
    world-level global rotation. Dummy state arrays are deliberately ignored by
    the conditional-null scorer.
    """

    if n_systems < 2 or n_points < 8:
        raise ValueError("homogeneous control requires >=2 systems and >=8 points")
    rng = np.random.default_rng(int(seed))
    rotation = float(rng.uniform(0.0, 2.0 * np.pi))
    theta = np.mod(
        rotation + np.linspace(0.0, 2.0 * np.pi, int(n_points), endpoint=False),
        2.0 * np.pi,
    )
    coordinates = np.column_stack((np.cos(theta), np.sin(theta)))
    zeros = np.zeros(int(n_points), dtype=float)
    return tuple(
        PairedSpeciesSample(
            species=f"sp_{i:03d}",
            coordinates=coordinates.copy(),
            state_a=zeros.copy(),
            state_b=zeros.copy(),
        )
        for i in range(int(n_systems))
    )


def _fixed_coupling_geometry(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    graph_fraction: float,
) -> tuple[
    dict[str, SpeciesEdges],
    dict[str, int],
    dict[str, int],
]:
    sample_map = {str(sample.species): sample for sample in samples}
    if len(sample_map) != len(samples):
        raise ValueError("system labels must be unique")
    train = tuple(map(str, train_species))
    evaluation = tuple(map(str, eval_species))
    if not train or not evaluation or set(train) & set(evaluation):
        raise ValueError("non-empty system-disjoint train/evaluation sets required")
    missing = (set(train) | set(evaluation)) - set(sample_map)
    if missing:
        raise ValueError(f"unknown systems in split: {sorted(missing)}")

    edges: dict[str, SpeciesEdges] = {}
    graph_k: dict[str, int] = {}
    effective_n: dict[str, int] = {}
    for name in train + evaluation:
        sample = sample_map[name]
        n = int(len(sample.coordinates))
        k = density_scaled_k(n, float(graph_fraction))
        nodes = knn_edges(sample.coordinates, k=k)
        # The returned turnover values are ignored. This call is used only to
        # construct the exact frozen Q5/Q5.1 edge geometry.
        edges[name] = build_coupling_edges(sample, edge_nodes=nodes)
        graph_k[name] = int(k)
        effective_n[name] = n
    return edges, graph_k, effective_n


def _phase_turnover_matrix(
    sample: PairedSpeciesSample,
    edges: SpeciesEdges,
    phases: np.ndarray,
) -> np.ndarray:
    """Exact raw coupling-turnover ranks for private relation rotations.

    Q5's exact relation rotation assigns [A,0] vs [0,A] according to the sign
    of sin(theta-phase). Therefore Euclidean relation dissimilarity is binary:
    zero within a side and a positive constant across sides. Rank-standardizing
    that binary crossing indicator exactly reproduces build_coupling_edges.
    """

    phi = np.asarray(phases, dtype=float)
    if phi.ndim != 1 or len(phi) < 1 or not np.isfinite(phi).all():
        raise ValueError("phases must be a finite non-empty vector")
    theta = np.mod(
        np.arctan2(sample.coordinates[:, 1], sample.coordinates[:, 0]),
        2.0 * np.pi,
    )
    side = np.sin(theta[:, None] - phi[None, :]) >= 0.0
    nodes = np.asarray(edges.nodes, dtype=np.int64)
    crossing = side[nodes[:, 0], :] != side[nodes[:, 1], :]
    out = np.empty((len(nodes), len(phi)), dtype=float)
    for j in range(len(phi)):
        out[:, j] = rank01(crossing[:, j].astype(float))
    return out


def q51_private_relation_scores_for_phases(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    phases_by_species: Mapping[str, np.ndarray],
    graph_fraction: float,
    bandwidth: float,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> ConditionalNullCenteringResult:
    """Score many independent private-phase worlds on one fixed geometry.

    Geometry, graph construction, kernel bandwidth, prior, and quadrature are
    frozen once. Only independent private relation phases vary across columns.
    The mathematical estimator for every column is exactly Q5.1 TTF-C.
    """

    train = tuple(map(str, train_species))
    evaluation = tuple(map(str, eval_species))
    sample_map = {str(sample.species): sample for sample in samples}
    edges, graph_k, effective_n = _fixed_coupling_geometry(
        samples,
        train_species=train,
        eval_species=evaluation,
        graph_fraction=float(graph_fraction),
    )

    widths = set()
    phase_map: dict[str, np.ndarray] = {}
    for name in train + evaluation:
        if name not in phases_by_species:
            raise ValueError(f"missing private phases for {name}")
        phi = np.asarray(phases_by_species[name], dtype=float)
        if phi.ndim != 1 or len(phi) < 1 or not np.isfinite(phi).all():
            raise ValueError(f"invalid private phases for {name}")
        widths.add(int(len(phi)))
        phase_map[name] = phi
    if len(widths) != 1:
        raise ValueError("all systems must have the same number of phase draws")

    raw_turnover = {
        name: _phase_turnover_matrix(sample_map[name], edges[name], phase_map[name])
        for name in train + evaluation
    }
    train_turnover = {
        name: np.column_stack(
            [
                length_orthogonalized_turnover(
                    raw_turnover[name][:, j],
                    edges[name].length,
                )
                for j in range(raw_turnover[name].shape[1])
            ]
        )
        for name in train
    }
    eval_turnover = {name: raw_turnover[name] for name in evaluation}

    prepared = prepare_chunked_transfer(
        [edges[name] for name in train],
        [edges[name] for name in evaluation],
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.0,
        segment_points=int(segment_points),
    )
    scored = score_chunked_batch(
        prepared,
        train_turnover,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    return ConditionalNullCenteringResult(
        statistics=np.asarray(scored.statistics, dtype=float),
        species_scores={
            str(name): np.asarray(values, dtype=float)
            for name, values in scored.species_scores.items()
        },
        phases=phase_map,
        graph_k=graph_k,
        effective_n=effective_n,
    )


def q51_private_relation_phase_redraw_scores(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    n_phase_draws: int,
    phase_seed: int,
    graph_fraction: float,
    bandwidth: float,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> ConditionalNullCenteringResult:
    if n_phase_draws < 2:
        raise ValueError("phase-redraw diagnostic requires at least two draws")
    names = tuple(map(str, train_species)) + tuple(map(str, eval_species))
    rng = np.random.default_rng(int(phase_seed))
    phases = {
        name: rng.uniform(0.0, 2.0 * np.pi, size=int(n_phase_draws))
        for name in names
    }
    return q51_private_relation_scores_for_phases(
        samples,
        train_species=train_species,
        eval_species=eval_species,
        phases_by_species=phases,
        graph_fraction=float(graph_fraction),
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        segment_points=int(segment_points),
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )


__all__ = [
    "ConditionalNullCenteringResult",
    "identical_uniform_circle_samples",
    "q51_private_relation_scores_for_phases",
    "q51_private_relation_phase_redraw_scores",
]
