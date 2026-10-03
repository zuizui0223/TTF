from __future__ import annotations

from dataclasses import dataclass
import hashlib
from collections import defaultdict
from collections.abc import Mapping, Sequence

import numpy as np

from .chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from .core import SpeciesEdges
from .genetic_geometry import GeneticSamplingGeometry
from .genetic_simulate import GeneticSyntheticWorld, simulate_genetic_distance_world
from .geometry_control import length_orthogonalized_turnover
from .phylogatr_compact_ibd import (
    crossfit_ibd_residuals_compact,
    prepare_compact_crossfit_ibd_design,
)
from .private_null_inference import envelope_upper_pvalues
from .relational_dyadic import BatchPrimaryResult, PreparedDyadicRegression, batch_primary_test


@dataclass(frozen=True)
class GeometryNullCenter:
    """Frozen dyad-specific private-world expectation used only for centering."""

    mu0: np.ndarray
    amplitudes: tuple[float, ...]
    worlds_per_amplitude: tuple[int, ...]


@dataclass(frozen=True)
class CenteredDyadicBatchResult:
    """Centered dyad response plus the existing fixed-effects primary test."""

    centered_response: np.ndarray
    regression: BatchPrimaryResult


def frozen_uint64_seed(
    master_seed: int,
    namespace: str,
    cell: str,
    replicate: int,
) -> int:
    """Deterministic disjoint seed construction for prospective synthetic worlds."""
    if not str(namespace).strip() or not str(cell).strip():
        raise ValueError("namespace and cell must be non-empty")
    if int(replicate) < 0:
        raise ValueError("replicate must be non-negative")
    payload = f"{int(master_seed)}|{namespace}|{cell}|{int(replicate)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big", signed=False)


def template_edges_from_geometries(
    geometries: Mapping[str, GeneticSamplingGeometry],
) -> dict[str, SpeciesEdges]:
    """Construct response-free edge objects from exact frozen locality graphs."""
    out: dict[str, SpeciesEdges] = {}
    for name in sorted(map(str, geometries)):
        geometry = geometries[name]
        nodes = np.asarray(geometry.edge_nodes, dtype=np.int64)
        coords = np.asarray(geometry.coordinates, dtype=float)
        start = coords[nodes[:, 0]]
        end = coords[nodes[:, 1]]
        out[name] = SpeciesEdges(
            species=name,
            nodes=nodes.copy(),
            start=start.copy(),
            end=end.copy(),
            midpoint=(0.5 * (start + end)),
            length=np.linalg.norm(end - start, axis=1),
            turnover=np.zeros(len(nodes), dtype=float),
        )
    return out


def post_ibd_response_batch(
    edge_map: Mapping[str, SpeciesEdges],
    worlds: Sequence[GeneticSyntheticWorld],
    *,
    min_training_edges: int = 5,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Convert synthetic genetic-distance worlds to fixed post-IBD responses."""
    batch = tuple(worlds)
    if not batch:
        raise ValueError("at least one synthetic world is required")
    labels = tuple(sorted(map(str, edge_map)))
    width = len(batch)
    train = {
        name: np.empty((edge_map[name].n_edges, width), dtype=float)
        for name in labels
    }
    evaluation = {
        name: np.empty((edge_map[name].n_edges, width), dtype=float)
        for name in labels
    }
    ibd_design = {
        name: prepare_compact_crossfit_ibd_design(
            edge_map[name].length,
            edge_map[name].nodes,
            min_training_edges=int(min_training_edges),
        )
        for name in labels
    }
    for column, world in enumerate(batch):
        missing = set(labels) - set(map(str, world.genetic_distance))
        if missing:
            raise ValueError(f"synthetic world is missing species: {sorted(missing)}")
        for name in labels:
            residual = crossfit_ibd_residuals_compact(
                world.genetic_distance[name],
                ibd_design[name],
            ).residual_turnover
            train[name][:, column] = length_orthogonalized_turnover(
                residual,
                edge_map[name].length,
            )
            evaluation[name][:, column] = residual
    return train, evaluation


def source_only_transfer_scores_batch(
    edge_map: Mapping[str, SpeciesEdges],
    train_response: Mapping[str, np.ndarray],
    eval_response: Mapping[str, np.ndarray],
    pairs: Sequence[tuple[str, str]],
    *,
    bandwidth_km: float = 500.0,
    prior_strength: float = 0.25,
    prior_mean: float = 0.0,
    segment_points: int = 5,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> np.ndarray:
    """Score many response worlds under the exact fixed single-source TTF operator.

    Geometry kernels are reused across all response columns. Work is grouped by
    source and streamed target-by-target by the chunked scorer, avoiding a
    persistent all-dyad dense projection cache.
    """
    pair_list = tuple((str(s), str(t)) for s, t in pairs)
    if not pair_list or len(set(pair_list)) != len(pair_list):
        raise ValueError("pairs must be unique and non-empty")
    sources = {s for s, _ in pair_list}
    targets = {t for _, t in pair_list}
    if sources & targets:
        raise ValueError("source and target roles must be species-disjoint")
    missing = (sources | targets) - set(map(str, edge_map))
    if missing:
        raise ValueError(f"missing edge geometry: {sorted(missing)}")

    widths = set()
    for name in sources:
        a = np.asarray(train_response[name], dtype=float)
        if a.ndim != 2 or a.shape[0] != edge_map[name].n_edges:
            raise ValueError(f"train response shape drift for {name}")
        widths.add(int(a.shape[1]))
    for name in targets:
        a = np.asarray(eval_response[name], dtype=float)
        if a.ndim != 2 or a.shape[0] != edge_map[name].n_edges:
            raise ValueError(f"eval response shape drift for {name}")
        widths.add(int(a.shape[1]))
    if len(widths) != 1:
        raise ValueError("response batch widths disagree")
    width = widths.pop()

    by_source: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for row, (source, target) in enumerate(pair_list):
        by_source[source].append((row, target))

    out = np.empty((len(pair_list), width), dtype=float)
    for source in sorted(by_source):
        row_targets = by_source[source]
        target_names = tuple(target for _, target in row_targets)
        if len(set(target_names)) != len(target_names):
            raise ValueError("a source-target dyad appears more than once")
        prepared = prepare_chunked_transfer(
            [edge_map[source]],
            [edge_map[target] for target in target_names],
            bandwidth=float(bandwidth_km),
            prior_strength=float(prior_strength),
            prior_mean=float(prior_mean),
            segment_points=int(segment_points),
        )
        scored = score_chunked_batch(
            prepared,
            {source: np.asarray(train_response[source], dtype=float)},
            {
                target: np.asarray(eval_response[target], dtype=float)
                for target in target_names
            },
            edge_chunk_size=int(edge_chunk_size),
            train_chunk_size=int(train_chunk_size),
        )
        for row, target in row_targets:
            out[row, :] = np.asarray(scored.species_scores[target], dtype=float)
    if not np.isfinite(out).all():
        raise RuntimeError("non-finite fixed-dyad T_st score")
    return out


def simulate_fixed_dyad_tst_batch(
    geometries: Mapping[str, GeneticSamplingGeometry],
    pairs: Sequence[tuple[str, str]],
    seeds: Sequence[int],
    *,
    shared_fraction: float,
    residual_amplitude: float,
    ibd_strength: float = 1.0,
    noise_sd: float = 0.10,
    transition_width: float = 0.20,
    noise_dimensions: int = 2,
    min_training_edges: int = 5,
    bandwidth_km: float = 500.0,
    prior_strength: float = 0.25,
    prior_mean: float = 0.0,
    segment_points: int = 5,
) -> np.ndarray:
    """Run the actual post-IBD fixed-dyad TTF operator on synthetic worlds."""
    seed_tuple = tuple(int(x) for x in seeds)
    if not seed_tuple:
        raise ValueError("at least one seed is required")
    worlds = tuple(
        simulate_genetic_distance_world(
            geometries,
            shared_fraction=float(shared_fraction),
            residual_amplitude=float(residual_amplitude),
            ibd_strength=float(ibd_strength),
            noise_sd=float(noise_sd),
            transition_width=float(transition_width),
            noise_dimensions=int(noise_dimensions),
            seed=seed,
        )
        for seed in seed_tuple
    )
    edge_map = template_edges_from_geometries(geometries)
    train, evaluation = post_ibd_response_batch(
        edge_map,
        worlds,
        min_training_edges=int(min_training_edges),
    )
    return source_only_transfer_scores_batch(
        edge_map,
        train,
        evaluation,
        pairs,
        bandwidth_km=float(bandwidth_km),
        prior_strength=float(prior_strength),
        prior_mean=float(prior_mean),
        segment_points=int(segment_points),
    )


def fit_geometry_null_center(
    private_tst_by_amplitude: Mapping[float, np.ndarray],
) -> GeometryNullCenter:
    """Freeze mu0_st with equal weight per declared private-amplitude cell."""
    if not private_tst_by_amplitude:
        raise ValueError("at least one private amplitude cell is required")
    means = []
    amplitudes = []
    counts = []
    n_dyads = None
    for amplitude in sorted(private_tst_by_amplitude, key=float):
        x = np.asarray(private_tst_by_amplitude[amplitude], dtype=float)
        if x.ndim != 2 or x.shape[1] < 1 or not np.isfinite(x).all():
            raise ValueError("private T_st cells must be finite dyad x world matrices")
        if n_dyads is None:
            n_dyads = int(x.shape[0])
        elif x.shape[0] != n_dyads:
            raise ValueError("private amplitude cells disagree on dyad rows")
        means.append(np.mean(x, axis=1))
        amplitudes.append(float(amplitude))
        counts.append(int(x.shape[1]))
    mu0 = np.mean(np.vstack(means), axis=0)
    return GeometryNullCenter(
        mu0=np.asarray(mu0, dtype=float),
        amplitudes=tuple(amplitudes),
        worlds_per_amplitude=tuple(counts),
    )


def center_dyad_scores(
    t_st: np.ndarray,
    center: GeometryNullCenter,
) -> np.ndarray:
    x = np.asarray(t_st, dtype=float)
    one = x.ndim == 1
    if one:
        x = x[:, None]
    if x.ndim != 2 or x.shape[0] != len(center.mu0):
        raise ValueError("T_st must have one row per frozen dyad")
    if not np.isfinite(x).all():
        raise ValueError("T_st must be finite")
    out = x - np.asarray(center.mu0, dtype=float)[:, None]
    return out[:, 0] if one else out


def centered_dyadic_primary_test(
    prepared: PreparedDyadicRegression,
    t_st: np.ndarray,
    center: GeometryNullCenter,
) -> CenteredDyadicBatchResult:
    centered = center_dyad_scores(t_st, center)
    tested = batch_primary_test(prepared, centered)
    return CenteredDyadicBatchResult(
        centered_response=np.asarray(centered, dtype=float),
        regression=tested,
    )


def least_favourable_beta_pvalues(
    observed_beta: np.ndarray,
    private_beta_references: Mapping[str, np.ndarray],
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Primary finite-sample inference; asymptotic cluster p-values are ignored."""
    return envelope_upper_pvalues(
        np.asarray(observed_beta, dtype=float),
        private_beta_references,
    )


__all__ = [
    "CenteredDyadicBatchResult",
    "GeometryNullCenter",
    "center_dyad_scores",
    "centered_dyadic_primary_test",
    "fit_geometry_null_center",
    "frozen_uint64_seed",
    "least_favourable_beta_pvalues",
    "post_ibd_response_batch",
    "simulate_fixed_dyad_tst_batch",
    "source_only_transfer_scores_batch",
    "template_edges_from_geometries",
]
