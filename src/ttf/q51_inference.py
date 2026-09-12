from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from .core import knn_edges
from .geometry_control import length_orthogonalized_turnover
from .heterogeneous_inference import (
    density_scaled_k,
    heterogeneous_paired_heldout_species_bootstrap_test,
)
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test
from .mismatch import PairedSpeciesSample, build_coupling_edges


@dataclass(frozen=True)
class Q51PairedInferenceResult:
    """Q5.1 candidate: preserve TTF-M and geometry-control TTF-C training only.

    Q5.1 is deliberately a minimal successor to the frozen Q5 failure. The
    mismatch path is exactly the existing heterogeneous Q5 estimator. Only the
    coupling path changes: within each training system, coupling-turnover ranks
    are orthogonalized against that system's kNN edge-length ranks before the
    transferable field is learned. Held-out coupling turnover remains raw and
    the held-out score remains the ordinary within-system Spearman correlation.

    This is the same response-blind, training-only edge-length correction that
    survived the core TTF v0.3 development ablation. Q5 remains immutable.
    """

    mismatch_statistic: float
    coupling_statistic: float
    raw_coupling_statistic: float
    mismatch_species_scores: Mapping[str, float]
    coupling_species_scores: Mapping[str, float]
    raw_coupling_species_scores: Mapping[str, float]
    mismatch_bootstrap: MeanBootstrapResult
    coupling_bootstrap: MeanBootstrapResult
    raw_coupling_bootstrap: MeanBootstrapResult
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def q51_paired_heldout_species_bootstrap_test(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    graph_fraction: float,
    bandwidth: float,
    prior_strength: float = 0.25,
    mismatch_prior_mean: float = 0.5,
    segment_points: int = 5,
    n_bootstrap: int = 1999,
    seed: int = 0,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> Q51PairedInferenceResult:
    """Run the prospectively specified Q5.1 coupling-geometry candidate.

    The TTF-M branch delegates to the frozen heterogeneous Q5 implementation so
    its estimand and inference are unchanged. TTF-C is rebuilt on the exact same
    density-scaled graphs, but only training coupling turnover is projected off
    edge-length rank. The coupling field uses prior mean zero because the
    orthogonalized training response is centered. Held-out coupling turnover is
    not residualized, filtered, or used to choose a nuisance correction.
    """

    base = heterogeneous_paired_heldout_species_bootstrap_test(
        samples,
        train_species=train_species,
        eval_species=eval_species,
        graph_fraction=float(graph_fraction),
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=float(mismatch_prior_mean),
        segment_points=int(segment_points),
        n_bootstrap=int(n_bootstrap),
        seed=int(seed),
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )

    sample_map = {sample.species: sample for sample in samples}
    if len(sample_map) != len(samples):
        raise ValueError("system labels must be unique")
    train = tuple(map(str, train_species))
    evaluation = tuple(map(str, eval_species))
    if not train or not evaluation or set(train) & set(evaluation):
        raise ValueError("non-empty system-disjoint train/evaluation sets required")
    missing = (set(train) | set(evaluation)) - set(sample_map)
    if missing:
        raise ValueError(f"unknown systems in split: {sorted(missing)}")

    coupling_edges = {}
    for name in train + evaluation:
        sample = sample_map[name]
        k = density_scaled_k(len(sample.coordinates), float(graph_fraction))
        nodes = knn_edges(sample.coordinates, k=k)
        coupling_edges[name] = build_coupling_edges(sample, edge_nodes=nodes)

    prepared = prepare_chunked_transfer(
        [coupling_edges[name] for name in train],
        [coupling_edges[name] for name in evaluation],
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.0,
        segment_points=int(segment_points),
    )

    train_turnover = {
        name: length_orthogonalized_turnover(
            coupling_edges[name].turnover,
            coupling_edges[name].length,
        )[:, None]
        for name in train
    }
    eval_turnover = {
        name: np.asarray(coupling_edges[name].turnover, dtype=float)[:, None]
        for name in evaluation
    }
    scored = score_chunked_batch(
        prepared,
        train_turnover,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )

    coupling_scores = {
        name: float(scored.species_scores[name][0]) for name in evaluation
    }
    c_vec = np.asarray([coupling_scores[name] for name in evaluation], dtype=float)
    c_boot = centered_species_bootstrap_mean_test(
        c_vec,
        n_bootstrap=int(n_bootstrap),
        seed=int(seed) + 1,
    )
    corrected_statistic = float(scored.statistics[0])
    if not np.isclose(c_boot.observed_mean, corrected_statistic, atol=1e-12, rtol=0.0):
        raise RuntimeError("Q5.1 TTF-C statistic drift")

    return Q51PairedInferenceResult(
        mismatch_statistic=float(base.mismatch_statistic),
        coupling_statistic=corrected_statistic,
        raw_coupling_statistic=float(base.coupling_statistic),
        mismatch_species_scores=dict(base.mismatch_species_scores),
        coupling_species_scores=coupling_scores,
        raw_coupling_species_scores=dict(base.coupling_species_scores),
        mismatch_bootstrap=base.mismatch_bootstrap,
        coupling_bootstrap=c_boot,
        raw_coupling_bootstrap=base.coupling_bootstrap,
        graph_k=dict(base.graph_k),
        effective_n=dict(base.effective_n),
    )


__all__ = ["Q51PairedInferenceResult", "q51_paired_heldout_species_bootstrap_test"]
