from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping, Sequence

import numpy as np

from .chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from .core import knn_edges
from .geometry_control import length_orthogonalized_turnover
from .geometry_transport import (
    GeometryTransportMode,
    geometry_transport_train_weights,
    transport_effective_edge_count,
)
from .heterogeneous_inference import density_scaled_k
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test
from .mismatch import PairedSpeciesSample, build_coupling_edges


@dataclass(frozen=True)
class Q52FormalResult:
    candidate: str
    statistic: float
    species_scores: Mapping[str, float]
    bootstrap: MeanBootstrapResult
    mean_train_effective_edge_count: float
    min_train_effective_edge_count: float
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def q52_selected_transport_test(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    candidate: GeometryTransportMode,
    graph_fraction: float,
    bandwidth: float,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    n_bootstrap: int = 1999,
    seed: int = 0,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 2048,
    density_chunk_size: int = 256,
) -> Q52FormalResult:
    """Score exactly one prospectively selected Q5.2 transport candidate.

    Formal qualification must never evaluate development reserve candidates on
    confirmatory worlds.  The candidate affects only response-blind training-edge
    field mass.  Q5.1's training response, raw held-out target, raw held-out
    Spearman score, graph rule, bandwidth and priors remain unchanged.
    """
    allowed = {
        "pooled_inverse_density",
        "self_inverse_density",
        "target_density_ratio",
    }
    if candidate not in allowed:
        raise ValueError(f"candidate is not an allowed frozen Q5.2 mode: {candidate}")

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

    coupling_edges = {}
    graph_k: dict[str, int] = {}
    effective_n: dict[str, int] = {}
    for name in train + evaluation:
        sample = sample_map[name]
        n = int(len(sample.coordinates))
        k = density_scaled_k(n, float(graph_fraction))
        nodes = knn_edges(sample.coordinates, k=k)
        coupling_edges[name] = build_coupling_edges(sample, edge_nodes=nodes)
        graph_k[name] = int(k)
        effective_n[name] = n

    train_edges = [coupling_edges[name] for name in train]
    eval_edges = [coupling_edges[name] for name in evaluation]
    prepared = prepare_chunked_transfer(
        train_edges,
        eval_edges,
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=0.0,
        segment_points=int(segment_points),
    )
    weights = geometry_transport_train_weights(
        train_edges,
        eval_edges,
        mode=candidate,
        bandwidth=float(bandwidth),
        density_chunk_size=int(density_chunk_size),
    )
    if weights.shape != prepared.train_weights.shape:
        raise RuntimeError("Q5.2 formal transport-weight shape drift")
    prepared = replace(prepared, train_weights=weights)

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
    species_scores = {
        name: float(scored.species_scores[name][0]) for name in evaluation
    }
    score_vector = np.asarray([species_scores[name] for name in evaluation], dtype=float)
    bootstrap = centered_species_bootstrap_mean_test(
        score_vector,
        n_bootstrap=int(n_bootstrap),
        seed=int(seed) + 1,
    )
    statistic = float(scored.statistics[0])
    if not np.isclose(bootstrap.observed_mean, statistic, atol=1e-12, rtol=0.0):
        raise RuntimeError("Q5.2 formal statistic drift")

    effective = transport_effective_edge_count(weights, train_edges)
    eff_values = np.asarray(list(effective.values()), dtype=float)
    return Q52FormalResult(
        candidate=str(candidate),
        statistic=statistic,
        species_scores=species_scores,
        bootstrap=bootstrap,
        mean_train_effective_edge_count=float(eff_values.mean()),
        min_train_effective_edge_count=float(eff_values.min()),
        graph_k=graph_k,
        effective_n=effective_n,
    )


__all__ = ["Q52FormalResult", "q52_selected_transport_test"]
