from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import (
    ChunkedTransferGeometry,
    _kernel_block,
    _train_value_chunk,
    _validate_train_values,
    prepare_chunked_transfer,
)
from .core import knn_edges, spearman_rho
from .geometry_control import length_orthogonalized_turnover
from .heterogeneous_inference import density_scaled_k
from .mismatch import (
    PairedSpeciesSample,
    build_coupling_edges,
    build_mismatch_edges,
)
from .one_sided_geometry_control import residual_rank_correlation


@dataclass(frozen=True)
class SupportEstimandDiagnostic:
    """Held-out transfer plus response-blind support-overlap diagnostics."""

    raw_statistics: np.ndarray
    opportunity_partial_statistics: np.ndarray
    alignment_statistics: np.ndarray
    mean_prior_fraction: float
    low_support_fraction: float
    mean_opportunity: float
    species_raw_scores: Mapping[str, np.ndarray]
    species_partial_scores: Mapping[str, np.ndarray]
    species_prediction_opportunity_rho: Mapping[str, np.ndarray]
    species_target_opportunity_rho: Mapping[str, np.ndarray]
    species_alignment_product: Mapping[str, np.ndarray]
    species_mean_opportunity: Mapping[str, float]
    species_mean_prior_fraction: Mapping[str, float]
    species_low_support_fraction: Mapping[str, float]


@dataclass(frozen=True)
class PairedSupportOverlapDiagnostic:
    """Frozen Q5/Q5.1 estimands evaluated against the same support geometry."""

    mismatch: SupportEstimandDiagnostic
    coupling: SupportEstimandDiagnostic
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def _score_with_support(
    prepared: ChunkedTransferGeometry,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> SupportEstimandDiagnostic:
    """Reproduce exact chunked scoring while exposing kernel opportunity.

    This function is diagnostic only. The raw score is the frozen TTF score.
    ``opportunity_partial_statistics`` is a counterfactual diagnostic summary,
    not a qualified estimator and not a repair candidate.
    """

    if edge_chunk_size < 1 or train_chunk_size < 1:
        raise ValueError("chunk sizes must be positive")
    train_arrays, width = _validate_train_values(prepared, train_turnover)
    t = (np.arange(prepared.segment_points, dtype=float) + 0.5) / prepared.segment_points
    h2 = prepared.bandwidth * prepared.bandwidth
    tiny = np.finfo(float).tiny
    n_train = len(prepared.train_positions)

    raw_map: dict[str, np.ndarray] = {}
    partial_map: dict[str, np.ndarray] = {}
    pred_opp_map: dict[str, np.ndarray] = {}
    target_opp_map: dict[str, np.ndarray] = {}
    alignment_map: dict[str, np.ndarray] = {}
    mean_opp_map: dict[str, float] = {}
    prior_fraction_map: dict[str, float] = {}
    low_support_map: dict[str, float] = {}

    raw_sum = np.zeros(width, dtype=float)
    partial_sum = np.zeros(width, dtype=float)
    alignment_sum = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)

    for species in prepared.eval_species:
        if species not in eval_turnover:
            raise ValueError(f"missing evaluation turnover for {species}")
        start = prepared.eval_start[species]
        end = prepared.eval_end[species]
        target = np.asarray(eval_turnover[species], dtype=float)
        if target.ndim != 2 or target.shape != (len(start), width):
            raise ValueError(f"evaluation turnover shape drift for {species}")
        if not np.isfinite(target).all():
            raise ValueError(f"non-finite evaluation turnover for {species}")

        points = start[:, None, :] + t[None, :, None] * (end - start)[:, None, :]
        opportunity = np.empty((len(start), prepared.segment_points), dtype=float)

        # Exact response-blind training-kernel mass at the frozen quadrature points.
        for e0 in range(0, len(start), edge_chunk_size):
            e1 = min(e0 + edge_chunk_size, len(start))
            flat = points[e0:e1].reshape(-1, start.shape[1])
            mass = np.zeros(len(flat), dtype=float)
            for p0 in range(0, n_train, train_chunk_size):
                p1 = min(p0 + train_chunk_size, n_train)
                kernel = _kernel_block(
                    flat,
                    prepared.train_positions[p0:p1],
                    prepared.train_weights[p0:p1],
                    bandwidth2=h2,
                )
                mass += kernel.sum(axis=1)
            opportunity[e0:e1, :] = mass.reshape(
                e1 - e0, prepared.segment_points
            )

        denominator = np.maximum(opportunity + prepared.prior_strength, tiny)
        edge_opportunity = opportunity.mean(axis=1)
        edge_prior_fraction = (
            prepared.prior_strength / denominator
        ).mean(axis=1)
        mean_opp_map[species] = float(edge_opportunity.mean())
        prior_fraction_map[species] = float(edge_prior_fraction.mean())
        low_support_map[species] = float(
            np.mean(edge_opportunity <= prepared.prior_strength)
        )

        prior_edge = (
            prepared.prior_strength * prepared.prior_mean / denominator
        ).mean(axis=1)
        predicted = np.repeat(prior_edge[:, None], width, axis=1)

        for p0 in range(0, n_train, train_chunk_size):
            p1 = min(p0 + train_chunk_size, n_train)
            values = _train_value_chunk(prepared, train_arrays, p0, p1, width)
            positions = prepared.train_positions[p0:p1]
            weights = prepared.train_weights[p0:p1]
            for e0 in range(0, len(start), edge_chunk_size):
                e1 = min(e0 + edge_chunk_size, len(start))
                flat = points[e0:e1].reshape(-1, start.shape[1])
                kernel = _kernel_block(
                    flat,
                    positions,
                    weights,
                    bandwidth2=h2,
                )
                normalized = kernel / denominator[e0:e1, :].reshape(-1, 1)
                projection = normalized.reshape(
                    e1 - e0,
                    prepared.segment_points,
                    p1 - p0,
                ).mean(axis=1)
                predicted[e0:e1, :] += projection @ values

        raw = np.empty(width, dtype=float)
        partial = np.empty(width, dtype=float)
        pred_opp = np.empty(width, dtype=float)
        target_opp = np.empty(width, dtype=float)
        for column in range(width):
            raw[column] = spearman_rho(predicted[:, column], target[:, column])
            partial[column] = residual_rank_correlation(
                predicted[:, column],
                target[:, column],
                prediction_nuisances=(edge_opportunity,),
                target_nuisances=(edge_opportunity,),
            )
            pred_opp[column] = spearman_rho(predicted[:, column], edge_opportunity)
            target_opp[column] = spearman_rho(target[:, column], edge_opportunity)

        alignment = pred_opp * target_opp
        raw_map[species] = raw
        partial_map[species] = partial
        pred_opp_map[species] = pred_opp
        target_opp_map[species] = target_opp
        alignment_map[species] = alignment

        finite = np.isfinite(raw) & np.isfinite(partial) & np.isfinite(alignment)
        raw_sum[finite] += raw[finite]
        partial_sum[finite] += partial[finite]
        alignment_sum[finite] += alignment[finite]
        counts[finite] += 1

    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite held-out support diagnostics")

    return SupportEstimandDiagnostic(
        raw_statistics=raw_sum / counts,
        opportunity_partial_statistics=partial_sum / counts,
        alignment_statistics=alignment_sum / counts,
        mean_prior_fraction=float(np.mean(list(prior_fraction_map.values()))),
        low_support_fraction=float(np.mean(list(low_support_map.values()))),
        mean_opportunity=float(np.mean(list(mean_opp_map.values()))),
        species_raw_scores=raw_map,
        species_partial_scores=partial_map,
        species_prediction_opportunity_rho=pred_opp_map,
        species_target_opportunity_rho=target_opp_map,
        species_alignment_product=alignment_map,
        species_mean_opportunity=mean_opp_map,
        species_mean_prior_fraction=prior_fraction_map,
        species_low_support_fraction=low_support_map,
    )


def paired_support_overlap_diagnostic(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    graph_fraction: float,
    bandwidth: float,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> PairedSupportOverlapDiagnostic:
    """Diagnose support overlap for frozen Q5 TTF-M and Q5.1 TTF-C."""

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
            raise RuntimeError("mismatch/coupling support geometry drift")

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

    mismatch = _score_with_support(
        m_prepared,
        m_train,
        m_eval,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    coupling = _score_with_support(
        c_prepared,
        c_train,
        c_eval,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )

    # Geometry-only support summaries must be exactly shared by M and C.
    if not np.isclose(
        mismatch.mean_prior_fraction,
        coupling.mean_prior_fraction,
        atol=1e-12,
        rtol=0.0,
    ) or not np.isclose(
        mismatch.low_support_fraction,
        coupling.low_support_fraction,
        atol=1e-12,
        rtol=0.0,
    ):
        raise RuntimeError("M/C support-overlap summaries drifted on shared geometry")

    return PairedSupportOverlapDiagnostic(
        mismatch=mismatch,
        coupling=coupling,
        graph_k=graph_k,
        effective_n=effective_n,
    )


__all__ = [
    "SupportEstimandDiagnostic",
    "PairedSupportOverlapDiagnostic",
    "paired_support_overlap_diagnostic",
]
