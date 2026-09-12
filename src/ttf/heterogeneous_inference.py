from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from .core import knn_edges
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test
from .mismatch import PairedSpeciesSample, build_coupling_edges, build_mismatch_edges


@dataclass(frozen=True)
class HeterogeneousPairedInferenceResult:
    mismatch_statistic: float
    coupling_statistic: float
    mismatch_species_scores: Mapping[str, float]
    coupling_species_scores: Mapping[str, float]
    mismatch_bootstrap: MeanBootstrapResult
    coupling_bootstrap: MeanBootstrapResult
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def density_scaled_k(n: int, fraction: float = 0.15) -> int:
    """Predeclared density-scaled graph degree used by Q5."""
    n = int(n)
    if n < 3:
        raise ValueError("at least three observations are required")
    if not (0.0 < float(fraction) < 1.0):
        raise ValueError("fraction must lie in (0,1)")
    return min(n - 1, max(2, int(np.floor(float(fraction) * n + 0.5))))


def heterogeneous_paired_heldout_species_bootstrap_test(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    graph_fraction: float,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
    n_bootstrap: int = 1999,
    seed: int = 0,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> HeterogeneousPairedInferenceResult:
    """Run TTF-M/TTF-C with a frozen density-scaled k chosen per system.

    Q5 varies observation counts by system. The graph rule is therefore
    ``k_s=max(2, round(graph_fraction*n_s))`` after response-blind observation
    loss. Train/evaluation systems remain disjoint and every system still
    contributes equal total training weight through core TTF.
    """

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
    if len(evaluation) < 6:
        raise ValueError("at least six held-out systems are required")

    all_used = train + evaluation
    k_map: dict[str, int] = {}
    n_map: dict[str, int] = {}
    mismatch_edges = {}
    coupling_edges = {}
    for name in all_used:
        sample = sample_map[name]
        n = int(len(sample.coordinates))
        k = density_scaled_k(n, float(graph_fraction))
        n_map[name] = n
        k_map[name] = k
        nodes = knn_edges(sample.coordinates, k=k)
        mismatch_edges[name] = build_mismatch_edges(sample, edge_nodes=nodes)
        coupling_edges[name] = build_coupling_edges(sample, edge_nodes=nodes)

    prepared = prepare_chunked_transfer(
        [mismatch_edges[name] for name in train],
        [mismatch_edges[name] for name in evaluation],
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=float(prior_mean),
        segment_points=int(segment_points),
    )
    train_turnover = {
        name: np.column_stack((mismatch_edges[name].turnover, coupling_edges[name].turnover))
        for name in train
    }
    eval_turnover = {
        name: np.column_stack((mismatch_edges[name].turnover, coupling_edges[name].turnover))
        for name in evaluation
    }
    scored = score_chunked_batch(
        prepared,
        train_turnover,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )

    mismatch_scores = {name: float(scored.species_scores[name][0]) for name in evaluation}
    coupling_scores = {name: float(scored.species_scores[name][1]) for name in evaluation}
    m_vec = np.asarray([mismatch_scores[name] for name in evaluation], dtype=float)
    c_vec = np.asarray([coupling_scores[name] for name in evaluation], dtype=float)
    m_boot = centered_species_bootstrap_mean_test(
        m_vec, n_bootstrap=int(n_bootstrap), seed=int(seed)
    )
    c_boot = centered_species_bootstrap_mean_test(
        c_vec, n_bootstrap=int(n_bootstrap), seed=int(seed) + 1
    )
    if not np.isclose(m_boot.observed_mean, float(scored.statistics[0]), atol=1e-12, rtol=0.0):
        raise RuntimeError("TTF-M statistic drift under heterogeneous inference")
    if not np.isclose(c_boot.observed_mean, float(scored.statistics[1]), atol=1e-12, rtol=0.0):
        raise RuntimeError("TTF-C statistic drift under heterogeneous inference")

    return HeterogeneousPairedInferenceResult(
        mismatch_statistic=float(scored.statistics[0]),
        coupling_statistic=float(scored.statistics[1]),
        mismatch_species_scores=mismatch_scores,
        coupling_species_scores=coupling_scores,
        mismatch_bootstrap=m_boot,
        coupling_bootstrap=c_boot,
        graph_k=k_map,
        effective_n=n_map,
    )


__all__ = [
    "HeterogeneousPairedInferenceResult",
    "density_scaled_k",
    "heterogeneous_paired_heldout_species_bootstrap_test",
]
