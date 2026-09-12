from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from .core import knn_edges
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test
from .mismatch import PairedSpeciesSample, build_coupling_edges, build_mismatch_edges


@dataclass(frozen=True)
class PairedHeldoutInferenceResult:
    """Joint TTF-M / TTF-C inference on one frozen paired-state geometry.

    The two estimands share exactly the same within-system graph, train/eval
    split, quadrature points and Gaussian kernel operator.  Only the edge
    turnover response column differs.  This is an execution optimization, not
    a multivariate test: mismatch and coupling retain separate species scores,
    bootstrap nulls and p-values.
    """

    mismatch_statistic: float
    coupling_statistic: float
    mismatch_species_scores: Mapping[str, float]
    coupling_species_scores: Mapping[str, float]
    mismatch_bootstrap: MeanBootstrapResult
    coupling_bootstrap: MeanBootstrapResult


def paired_heldout_species_bootstrap_test(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    k: int,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
    n_bootstrap: int = 1999,
    seed: int = 0,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> PairedHeldoutInferenceResult:
    """Run core held-out-system inference for TTF-M and TTF-C together.

    This function intentionally does not alter the core TTF statistic.  It
    projects paired states into mismatch magnitude M and relation state R,
    constructs one fixed kNN graph per system, and feeds both turnover columns
    through the exact chunked Gaussian transfer operator.
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
    mismatch_edges = {}
    coupling_edges = {}
    for name in all_used:
        sample = sample_map[name]
        nodes = knn_edges(sample.coordinates, k=int(k))
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
        name: np.column_stack(
            (mismatch_edges[name].turnover, coupling_edges[name].turnover)
        )
        for name in train
    }
    eval_turnover = {
        name: np.column_stack(
            (mismatch_edges[name].turnover, coupling_edges[name].turnover)
        )
        for name in evaluation
    }
    scored = score_chunked_batch(
        prepared,
        train_turnover,
        eval_turnover,
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )

    mismatch_scores = {
        name: float(scored.species_scores[name][0]) for name in evaluation
    }
    coupling_scores = {
        name: float(scored.species_scores[name][1]) for name in evaluation
    }
    m_vec = np.asarray([mismatch_scores[name] for name in evaluation], dtype=float)
    c_vec = np.asarray([coupling_scores[name] for name in evaluation], dtype=float)

    mismatch_bootstrap = centered_species_bootstrap_mean_test(
        m_vec,
        n_bootstrap=int(n_bootstrap),
        seed=int(seed),
    )
    coupling_bootstrap = centered_species_bootstrap_mean_test(
        c_vec,
        n_bootstrap=int(n_bootstrap),
        seed=int(seed) + 1,
    )

    if not np.isclose(
        mismatch_bootstrap.observed_mean,
        float(scored.statistics[0]),
        atol=1e-12,
        rtol=0.0,
    ):
        raise RuntimeError("TTF-M species-score mean drifted from joint transfer statistic")
    if not np.isclose(
        coupling_bootstrap.observed_mean,
        float(scored.statistics[1]),
        atol=1e-12,
        rtol=0.0,
    ):
        raise RuntimeError("TTF-C species-score mean drifted from joint transfer statistic")

    return PairedHeldoutInferenceResult(
        mismatch_statistic=float(scored.statistics[0]),
        coupling_statistic=float(scored.statistics[1]),
        mismatch_species_scores=mismatch_scores,
        coupling_species_scores=coupling_scores,
        mismatch_bootstrap=mismatch_bootstrap,
        coupling_bootstrap=coupling_bootstrap,
    )
