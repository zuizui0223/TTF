from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from .conditional_null_centering import _fixed_coupling_geometry, _phase_turnover_matrix
from .core import SpeciesEdges, rank01
from .mismatch import PairedSpeciesSample
from .normalization_bias_diagnostic import length_residual_unscaled


@dataclass(frozen=True)
class PhasePropensityResult:
    """Paired unscaled-Q5.1 and analytic phase-propensity-centered scores."""

    baseline_statistics: np.ndarray
    propensity_centered_statistics: np.ndarray
    baseline_species_scores: Mapping[str, np.ndarray]
    propensity_centered_species_scores: Mapping[str, np.ndarray]
    expected_training_residual: Mapping[str, np.ndarray]
    phases: Mapping[str, np.ndarray]
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def edge_private_phase_crossing_probability(
    sample: PairedSpeciesSample,
    edges: SpeciesEdges,
) -> np.ndarray:
    """Exact P(edge crosses the Q5 private relation boundary) for uniform phase.

    Q5 defines the two relation states by the sign of ``sin(theta - phi)``.
    For an edge whose endpoint angles have shortest circular separation
    ``delta in [0, pi]``, a uniform phase places the endpoints on opposite sides
    with probability exactly ``delta / pi``.
    """

    theta = np.mod(
        np.arctan2(sample.coordinates[:, 1], sample.coordinates[:, 0]),
        2.0 * np.pi,
    )
    nodes = np.asarray(edges.nodes, dtype=np.int64)
    raw = theta[nodes[:, 0]] - theta[nodes[:, 1]]
    delta = np.abs(np.angle(np.exp(1j * raw)))
    p = delta / np.pi
    if p.shape != (edges.n_edges,) or not np.isfinite(p).all():
        raise RuntimeError("edge crossing probability shape drift")
    if np.any(p < -1e-15) or np.any(p > 1.0 + 1e-15):
        raise RuntimeError("edge crossing probability outside [0,1]")
    return np.clip(p, 0.0, 1.0)


def phase_expected_centered_rank_turnover(
    sample: PairedSpeciesSample,
    edges: SpeciesEdges,
) -> np.ndarray:
    """Exact phase expectation of centered binary relation-turnover rank.

    If ``I_e(phi)`` is the binary edge-crossing indicator and ``q(phi)`` is its
    edge mean, ``rank01(I_e) = 0.5 - 0.5*q + 0.5*I_e``.  Averaging over a
    uniform private phase therefore yields

        E[rank01(I_e) - 0.5] = 0.5 * (p_e - mean(p)).
    """

    p = edge_private_phase_crossing_probability(sample, edges)
    return 0.5 * (p - float(np.mean(p)))


def phase_expected_unscaled_q51_residual(
    sample: PairedSpeciesSample,
    edges: SpeciesEdges,
) -> np.ndarray:
    """Exact geometry-only phase mean of the unscaled Q5.1 training response.

    The Q5.1 centered rank-length projection is linear in the centered turnover
    rank.  Applying that same projection to its exact phase expectation gives
    the exact phase-expected unscaled response.  No realized response, held-out
    response, fitted coefficient across worlds, or tuning parameter enters.
    """

    expected = phase_expected_centered_rank_turnover(sample, edges)
    length_rank = rank01(np.asarray(edges.length, dtype=float))
    zc = length_rank - float(np.mean(length_rank))
    den = float(np.dot(zc, zc))
    if den <= np.finfo(float).eps:
        residual = expected.copy()
    else:
        beta = float(np.dot(zc, expected) / den)
        residual = expected - beta * zc
    if not np.isfinite(residual).all():
        raise RuntimeError("non-finite phase-expected residual")
    return residual


def phase_propensity_scores_for_phases(
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
) -> PhasePropensityResult:
    """Pair baseline and analytic phase-propensity-centered scores.

    The baseline is the unscaled rank-length residual frozen by the immediately
    preceding normalization-bias diagnosis.  The diagnostic arm subtracts only
    the exact geometry-only phase expectation of that training response.
    Evaluation turnover, graph geometry, Gaussian field, prior, quadrature, and
    held-out Spearman scoring are unchanged.
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
    baseline_train = {
        name: np.column_stack(
            [
                length_residual_unscaled(
                    raw_turnover[name][:, j],
                    edges[name].length,
                )
                for j in range(raw_turnover[name].shape[1])
            ]
        )
        for name in train
    }
    expected = {
        name: phase_expected_unscaled_q51_residual(sample_map[name], edges[name])
        for name in train
    }
    centered_train = {
        name: baseline_train[name] - expected[name][:, None]
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
    baseline = score_chunked_batch(
        prepared,
        baseline_train,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    centered = score_chunked_batch(
        prepared,
        centered_train,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )

    base_stats = np.asarray(baseline.statistics, dtype=float)
    centered_stats = np.asarray(centered.statistics, dtype=float)
    if base_stats.shape != centered_stats.shape:
        raise RuntimeError("phase-propensity paired score width drift")
    return PhasePropensityResult(
        baseline_statistics=base_stats,
        propensity_centered_statistics=centered_stats,
        baseline_species_scores={
            str(name): np.asarray(values, dtype=float)
            for name, values in baseline.species_scores.items()
        },
        propensity_centered_species_scores={
            str(name): np.asarray(values, dtype=float)
            for name, values in centered.species_scores.items()
        },
        expected_training_residual={
            str(name): np.asarray(values, dtype=float) for name, values in expected.items()
        },
        phases=phase_map,
        graph_k=graph_k,
        effective_n=effective_n,
    )


def phase_propensity_phase_redraw_scores(
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
) -> PhasePropensityResult:
    if n_phase_draws < 2:
        raise ValueError("phase-propensity diagnostic requires at least two phase draws")
    names = tuple(map(str, train_species)) + tuple(map(str, eval_species))
    rng = np.random.default_rng(int(phase_seed))
    phases = {
        name: rng.uniform(0.0, 2.0 * np.pi, size=int(n_phase_draws))
        for name in names
    }
    return phase_propensity_scores_for_phases(
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
    "PhasePropensityResult",
    "edge_private_phase_crossing_probability",
    "phase_expected_centered_rank_turnover",
    "phase_expected_unscaled_q51_residual",
    "phase_propensity_scores_for_phases",
    "phase_propensity_phase_redraw_scores",
]
