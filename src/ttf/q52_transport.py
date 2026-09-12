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


DEVELOPMENT_MODES: tuple[str, ...] = (
    "q5_1_uniform",
    "pooled_inverse_density",
    "self_inverse_density",
    "target_density_ratio",
)


@dataclass(frozen=True)
class Q52ModeResult:
    statistic: float
    species_scores: Mapping[str, float]
    bootstrap: MeanBootstrapResult
    mean_train_effective_edge_count: float
    min_train_effective_edge_count: float


@dataclass(frozen=True)
class Q52DevelopmentResult:
    modes: Mapping[str, Q52ModeResult]
    graph_k: Mapping[str, int]
    effective_n: Mapping[str, int]


def q52_geometry_transport_development_test(
    samples: Sequence[PairedSpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    graph_fraction: float,
    bandwidth: float,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    n_bootstrap: int = 999,
    seed: int = 0,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 2048,
    density_chunk_size: int = 256,
) -> Q52DevelopmentResult:
    """Evaluate frozen Q5.2 geometry-only field-weight candidates on one world.

    Q5.2 keeps the Q5.1 response estimand fixed: training coupling turnover is
    orthogonalized against within-system edge length, held-out coupling turnover
    remains raw, and held-out scoring remains ordinary Spearman correlation.
    The only candidate dimension is the response-blind mass assigned to training
    coupling edges when the transferable field is constructed.

    ``q5_1_uniform`` is audit-only and reproduces Q5.1's uniform within-system
    training-edge mass. The three transport modes are development candidates.
    """
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

    mode_results: dict[str, Q52ModeResult] = {}
    for mode in DEVELOPMENT_MODES:
        if mode == "q5_1_uniform":
            candidate_prepared = prepared
            effective = {
                name: float(coupling_edges[name].n_edges) for name in train
            }
        else:
            weights = geometry_transport_train_weights(
                train_edges,
                eval_edges,
                mode=mode,  # type: ignore[arg-type]
                bandwidth=float(bandwidth),
                density_chunk_size=int(density_chunk_size),
            )
            if weights.shape != prepared.train_weights.shape:
                raise RuntimeError("Q5.2 packed transport-weight shape drift")
            candidate_prepared = replace(prepared, train_weights=weights)
            effective = transport_effective_edge_count(weights, train_edges)

        scored = score_chunked_batch(
            candidate_prepared,
            train_turnover,
            eval_turnover,
            edge_chunk_size=int(edge_chunk_size),
            train_chunk_size=int(train_chunk_size),
        )
        species_scores = {
            name: float(scored.species_scores[name][0]) for name in evaluation
        }
        score_vector = np.asarray(
            [species_scores[name] for name in evaluation], dtype=float
        )
        bootstrap = centered_species_bootstrap_mean_test(
            score_vector,
            n_bootstrap=int(n_bootstrap),
            seed=int(seed) + 1,
        )
        statistic = float(scored.statistics[0])
        if not np.isclose(
            bootstrap.observed_mean,
            statistic,
            atol=1e-12,
            rtol=0.0,
        ):
            raise RuntimeError(f"Q5.2 statistic drift for {mode}")
        eff_values = np.asarray(list(effective.values()), dtype=float)
        mode_results[mode] = Q52ModeResult(
            statistic=statistic,
            species_scores=species_scores,
            bootstrap=bootstrap,
            mean_train_effective_edge_count=float(eff_values.mean()),
            min_train_effective_edge_count=float(eff_values.min()),
        )

    return Q52DevelopmentResult(
        modes=mode_results,
        graph_k=graph_k,
        effective_n=effective_n,
    )


__all__ = [
    "DEVELOPMENT_MODES",
    "Q52ModeResult",
    "Q52DevelopmentResult",
    "q52_geometry_transport_development_test",
]
