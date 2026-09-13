from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from .conditional_null_centering import _fixed_coupling_geometry, _phase_turnover_matrix
from .core import rank01
from .geometry_control import length_orthogonalized_turnover
from .mismatch import PairedSpeciesSample


@dataclass(frozen=True)
class NormalizationBiasPhaseResult:
    """Paired frozen-Q5.1 and unscaled-residual scores on identical phase draws."""

    normalized_statistics: np.ndarray
    unscaled_statistics: np.ndarray
    normalized_species_scores: Mapping[str, np.ndarray]
    unscaled_species_scores: Mapping[str, np.ndarray]
    phases: Mapping[str, np.ndarray]
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def length_residual_unscaled(
    turnover: np.ndarray,
    edge_length: np.ndarray,
) -> np.ndarray:
    """Exact Q5.1 rank-length residual before the final SD normalization.

    This copies the linear stage of ``length_orthogonalized_turnover`` exactly:
    defensive rank transform, centered rank-length projection, and the same
    near-zero residual guard.  The *only* omitted operation is multiplication by
    ``target_sd / residual_sd``.
    """

    y = np.asarray(turnover, dtype=float)
    length = np.asarray(edge_length, dtype=float)
    if y.ndim != 1 or length.shape != y.shape or len(y) < 3:
        raise ValueError("turnover and edge_length must be equal 1D arrays with >=3 edges")
    if not np.isfinite(y).all() or not np.isfinite(length).all():
        raise ValueError("turnover and edge_length must be finite")
    if np.any(length < 0):
        raise ValueError("edge lengths cannot be negative")

    y_rank = rank01(y)
    length_rank = rank01(length)
    yc = y_rank - float(y_rank.mean())
    zc = length_rank - float(length_rank.mean())
    den = float(np.dot(zc, zc))
    if den <= np.finfo(float).eps:
        residual = yc
    else:
        beta = float(np.dot(zc, yc) / den)
        residual = yc - beta * zc

    residual_sd = float(np.std(residual, ddof=0))
    if residual_sd <= np.sqrt(np.finfo(float).eps):
        return np.zeros_like(residual)
    return residual


def normalization_bias_scores_for_phases(
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
) -> NormalizationBiasPhaseResult:
    """Pair normalized Q5.1 and unscaled-residual scores on fixed geometry.

    The graph, phase draws, raw relation turnover, held-out responses, transfer
    operator, bandwidth, prior and quadrature are shared exactly.  Only the final
    per-training-system residual-SD rescaling differs.
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
    normalized_train = {
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
    unscaled_train = {
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
    eval_turnover = {name: raw_turnover[name] for name in evaluation}

    prepared = prepare_chunked_transfer(
        [edges[name] for name in train],
        [edges[name] for name in evaluation],
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.0,
        segment_points=int(segment_points),
    )
    normalized = score_chunked_batch(
        prepared,
        normalized_train,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    unscaled = score_chunked_batch(
        prepared,
        unscaled_train,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )

    norm_stats = np.asarray(normalized.statistics, dtype=float)
    raw_stats = np.asarray(unscaled.statistics, dtype=float)
    if norm_stats.shape != raw_stats.shape:
        raise RuntimeError("paired normalization diagnostic width drift")
    return NormalizationBiasPhaseResult(
        normalized_statistics=norm_stats,
        unscaled_statistics=raw_stats,
        normalized_species_scores={
            str(name): np.asarray(values, dtype=float)
            for name, values in normalized.species_scores.items()
        },
        unscaled_species_scores={
            str(name): np.asarray(values, dtype=float)
            for name, values in unscaled.species_scores.items()
        },
        phases=phase_map,
        graph_k=graph_k,
        effective_n=effective_n,
    )


def normalization_bias_phase_redraw_scores(
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
) -> NormalizationBiasPhaseResult:
    if n_phase_draws < 2:
        raise ValueError("normalization-bias diagnostic requires at least two phase draws")
    names = tuple(map(str, train_species)) + tuple(map(str, eval_species))
    rng = np.random.default_rng(int(phase_seed))
    phases = {
        name: rng.uniform(0.0, 2.0 * np.pi, size=int(n_phase_draws))
        for name in names
    }
    return normalization_bias_scores_for_phases(
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
    "NormalizationBiasPhaseResult",
    "length_residual_unscaled",
    "normalization_bias_scores_for_phases",
    "normalization_bias_phase_redraw_scores",
]
