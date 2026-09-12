from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import ChunkedTransferGeometry, _kernel_block, _validate_train_values, prepare_chunked_transfer
from .core import knn_edges, spearman_rho
from .geometry_control import length_orthogonalized_turnover
from .heterogeneous_inference import density_scaled_k
from .mismatch import PairedSpeciesSample, build_coupling_edges, build_mismatch_edges


@dataclass(frozen=True)
class TrainingLeverageEstimandDiagnostic:
    raw_statistics: np.ndarray
    loo_score_sd: np.ndarray
    loo_score_mean: np.ndarray
    loo_score_min: np.ndarray
    mean_abs_loo_score_shift: np.ndarray
    mean_prediction_loo_sd: np.ndarray
    loo_statistics_by_training_system: Mapping[str, np.ndarray]
    species_raw_scores: Mapping[str, np.ndarray]


@dataclass(frozen=True)
class PairedTrainingLeverageDiagnostic:
    mismatch: TrainingLeverageEstimandDiagnostic
    coupling: TrainingLeverageEstimandDiagnostic
    mean_effective_training_system_count: float
    mean_dominant_training_system_share: float
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def _score_with_training_leverage(
    prepared: ChunkedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
) -> tuple[TrainingLeverageEstimandDiagnostic, float, float]:
    """Reproduce the frozen score and expose exact leave-one-training-system leverage.

    The local kernel numerator and opportunity denominator are decomposed by
    training system. Removing one system is then exact subtraction, not a
    refitted approximation. Frozen bandwidth, prior, quadrature and remaining
    per-system edge weights are unchanged.
    """

    train_arrays, width = _validate_train_values(prepared, train_turnover)
    train_names = tuple(prepared.train_species)
    n_train_systems = len(train_names)
    if n_train_systems < 3:
        raise ValueError("training-leverage diagnostic requires at least three training systems")

    t = (np.arange(prepared.segment_points, dtype=float) + 0.5) / prepared.segment_points
    h2 = prepared.bandwidth * prepared.bandwidth
    tiny = np.finfo(float).tiny

    full_sum = np.zeros(width, dtype=float)
    full_count = np.zeros(width, dtype=np.int64)
    loo_sum = np.zeros((n_train_systems, width), dtype=float)
    loo_count = np.zeros((n_train_systems, width), dtype=np.int64)
    prediction_loo_sd_species: list[np.ndarray] = []
    effective_species: list[float] = []
    dominant_species: list[float] = []
    species_raw_scores: dict[str, np.ndarray] = {}

    for eval_name in prepared.eval_species:
        if eval_name not in eval_turnover:
            raise ValueError(f"missing evaluation turnover for {eval_name}")
        start = prepared.eval_start[eval_name]
        end = prepared.eval_end[eval_name]
        target = np.asarray(eval_turnover[eval_name], dtype=float)
        if target.ndim != 2 or target.shape != (len(start), width):
            raise ValueError(f"evaluation turnover shape drift for {eval_name}")
        if not np.isfinite(target).all():
            raise ValueError(f"non-finite evaluation turnover for {eval_name}")

        points = start[:, None, :] + t[None, :, None] * (end - start)[:, None, :]
        flat = points.reshape(-1, start.shape[1])
        n_flat = len(flat)

        system_mass = np.empty((n_flat, n_train_systems), dtype=float)
        system_num = np.empty((n_flat, n_train_systems, width), dtype=float)
        for j, train_name in enumerate(train_names):
            sl = prepared.train_slices[train_name]
            positions = prepared.train_positions[sl]
            weights = prepared.train_weights[sl]
            kernel = _kernel_block(flat, positions, weights, bandwidth2=h2)
            values = np.asarray(train_arrays[train_name], dtype=float)
            system_mass[:, j] = kernel.sum(axis=1)
            system_num[:, j, :] = kernel @ values

        total_mass = system_mass.sum(axis=1)
        total_num = system_num.sum(axis=1)
        denominator = np.maximum(total_mass + prepared.prior_strength, tiny)
        full_flat = (
            prepared.prior_strength * prepared.prior_mean + total_num
        ) / denominator[:, None]
        full_pred = full_flat.reshape(len(start), prepared.segment_points, width).mean(axis=1)

        shares = system_mass / np.maximum(total_mass[:, None], tiny)
        effective = 1.0 / np.maximum(np.sum(shares * shares, axis=1), tiny)
        dominant = np.max(shares, axis=1)
        effective_species.append(float(effective.mean()))
        dominant_species.append(float(dominant.mean()))

        full_scores = np.asarray(
            [spearman_rho(full_pred[:, k], target[:, k]) for k in range(width)],
            dtype=float,
        )
        species_raw_scores[eval_name] = full_scores
        finite = np.isfinite(full_scores)
        full_sum[finite] += full_scores[finite]
        full_count[finite] += 1

        loo_pred_stack = np.empty((n_train_systems, len(start), width), dtype=float)
        for j in range(n_train_systems):
            mass_j = system_mass[:, j]
            num_j = system_num[:, j, :]
            loo_den = np.maximum(total_mass - mass_j + prepared.prior_strength, tiny)
            loo_flat = (
                prepared.prior_strength * prepared.prior_mean + total_num - num_j
            ) / loo_den[:, None]
            loo_pred = loo_flat.reshape(len(start), prepared.segment_points, width).mean(axis=1)
            loo_pred_stack[j] = loo_pred
            scores = np.asarray(
                [spearman_rho(loo_pred[:, k], target[:, k]) for k in range(width)],
                dtype=float,
            )
            ok = np.isfinite(scores)
            loo_sum[j, ok] += scores[ok]
            loo_count[j, ok] += 1

        pred_sd = np.std(loo_pred_stack, axis=0, ddof=0).mean(axis=0)
        prediction_loo_sd_species.append(np.asarray(pred_sd, dtype=float))

    if np.any(full_count == 0) or np.any(loo_count == 0):
        raise ValueError("one or more worlds had no finite training-leverage scores")

    full = full_sum / full_count
    loo = loo_sum / loo_count
    loo_mean = loo.mean(axis=0)
    loo_sd = loo.std(axis=0, ddof=0)
    loo_min = loo.min(axis=0)
    mean_abs_shift = np.mean(np.abs(loo - full[None, :]), axis=0)
    mean_prediction_loo_sd = np.mean(np.vstack(prediction_loo_sd_species), axis=0)

    return (
        TrainingLeverageEstimandDiagnostic(
            raw_statistics=full,
            loo_score_sd=loo_sd,
            loo_score_mean=loo_mean,
            loo_score_min=loo_min,
            mean_abs_loo_score_shift=mean_abs_shift,
            mean_prediction_loo_sd=mean_prediction_loo_sd,
            loo_statistics_by_training_system={
                name: np.asarray(loo[j], dtype=float) for j, name in enumerate(train_names)
            },
            species_raw_scores=species_raw_scores,
        ),
        float(np.mean(effective_species)),
        float(np.mean(dominant_species)),
    )


def paired_training_leverage_diagnostic(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    graph_fraction: float,
    bandwidth: float,
    prior_strength: float = 0.25,
    segment_points: int = 5,
) -> PairedTrainingLeverageDiagnostic:
    """Diagnose exact training-system concentration/leverage for Q5 M and Q5.1 C."""

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

    mismatch_edges = {}
    coupling_edges = {}
    graph_k: dict[str, int] = {}
    effective_n: dict[str, int] = {}
    for name in train + evaluation:
        sample = sample_map[name]
        n = int(len(sample.coordinates))
        k = density_scaled_k(n, float(graph_fraction))
        nodes = knn_edges(sample.coordinates, k=k)
        mismatch_edges[name] = build_mismatch_edges(sample, edge_nodes=nodes)
        coupling_edges[name] = build_coupling_edges(sample, edge_nodes=nodes)
        graph_k[name] = int(k)
        effective_n[name] = n
        if not np.allclose(
            mismatch_edges[name].midpoint,
            coupling_edges[name].midpoint,
            atol=0.0,
            rtol=0.0,
        ):
            raise RuntimeError("mismatch/coupling geometry drift")

    m_prepared = prepare_chunked_transfer(
        [mismatch_edges[name] for name in train],
        [mismatch_edges[name] for name in evaluation],
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.5,
        segment_points=int(segment_points),
    )
    c_prepared = prepare_chunked_transfer(
        [coupling_edges[name] for name in train],
        [coupling_edges[name] for name in evaluation],
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.0,
        segment_points=int(segment_points),
    )

    m_train = {
        name: np.asarray(mismatch_edges[name].turnover, dtype=float)[:, None]
        for name in train
    }
    m_eval = {
        name: np.asarray(mismatch_edges[name].turnover, dtype=float)[:, None]
        for name in evaluation
    }
    c_train = {
        name: length_orthogonalized_turnover(
            coupling_edges[name].turnover,
            coupling_edges[name].length,
        )[:, None]
        for name in train
    }
    c_eval = {
        name: np.asarray(coupling_edges[name].turnover, dtype=float)[:, None]
        for name in evaluation
    }

    mismatch, m_eff, m_dom = _score_with_training_leverage(m_prepared, m_train, m_eval)
    coupling, c_eff, c_dom = _score_with_training_leverage(c_prepared, c_train, c_eval)
    if not np.isclose(m_eff, c_eff, atol=1e-12, rtol=0.0):
        raise RuntimeError("M/C effective training-system count drifted")
    if not np.isclose(m_dom, c_dom, atol=1e-12, rtol=0.0):
        raise RuntimeError("M/C dominant training-system share drifted")

    return PairedTrainingLeverageDiagnostic(
        mismatch=mismatch,
        coupling=coupling,
        mean_effective_training_system_count=m_eff,
        mean_dominant_training_system_share=m_dom,
        graph_k=graph_k,
        effective_n=effective_n,
    )


__all__ = [
    "TrainingLeverageEstimandDiagnostic",
    "PairedTrainingLeverageDiagnostic",
    "paired_training_leverage_diagnostic",
]
