from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import (
    Dissimilarity,
    SpeciesEdges,
    SpeciesSample,
    build_species_edges,
    knn_edges,
)
from .transfer import PreparedTransfer, TransferResult, prepare_transfer


@dataclass(frozen=True)
class PermutationResult:
    observed: TransferResult
    null_statistics: np.ndarray
    p_value: float
    null_mean: float
    null_sd: float


def permute_trait_within_species(
    sample: SpeciesSample,
    rng: np.random.Generator,
) -> SpeciesSample:
    """Shuffle trait-to-location assignment while coordinates remain untouched.

    If ``blocks`` are supplied, exchangeability is restricted to each block.
    """
    n = len(sample.coordinates)
    source = np.arange(n, dtype=np.int64)
    if sample.blocks is None:
        source = rng.permutation(source)
    else:
        blocks = np.asarray(sample.blocks)
        for block in np.unique(blocks):
            positions = np.flatnonzero(blocks == block)
            source[positions] = rng.permutation(positions)
    return SpeciesSample(
        species=sample.species,
        coordinates=sample.coordinates,
        trait=sample.trait[source],
        blocks=sample.blocks,
    )


def fixed_graphs(
    samples: Sequence[SpeciesSample],
    *,
    k: int = 4,
    max_distance: float | None = None,
) -> dict[str, np.ndarray]:
    """Build each species graph once; trait permutations never change it."""
    if len({sample.species for sample in samples}) != len(samples):
        raise ValueError("species labels must be unique")
    return {
        sample.species: knn_edges(
            sample.coordinates,
            k=k,
            max_distance=max_distance,
        )
        for sample in samples
    }


def edges_on_fixed_graphs(
    samples: Sequence[SpeciesSample],
    graphs: Mapping[str, np.ndarray],
    *,
    dissimilarity: Dissimilarity | None = None,
) -> dict[str, SpeciesEdges]:
    return {
        sample.species: build_species_edges(
            sample,
            dissimilarity=dissimilarity,
            edge_nodes=graphs[sample.species],
        )
        for sample in samples
    }


def permutation_test(
    samples: Sequence[SpeciesSample],
    *,
    train_species: Sequence[str],
    eval_species: Sequence[str],
    k: int = 4,
    max_distance: float | None = None,
    bandwidth: float = 0.1,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
    n_permutations: int = 199,
    seed: int = 0,
    dissimilarity: Dissimilarity | None = None,
) -> PermutationResult:
    """Full TTF null: shuffle traits, refit train field, rescore held-out species."""
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")
    sample_map = {sample.species: sample for sample in samples}
    if len(sample_map) != len(samples):
        raise ValueError("species labels must be unique")
    train = tuple(map(str, train_species))
    evaluation = tuple(map(str, eval_species))
    if not train or not evaluation or set(train) & set(evaluation):
        raise ValueError("non-empty species-disjoint train/evaluation sets required")
    missing = (set(train) | set(evaluation)) - set(sample_map)
    if missing:
        raise ValueError(f"unknown species in split: {sorted(missing)}")

    used_samples = [sample_map[s] for s in train + evaluation]
    graphs = fixed_graphs(used_samples, k=k, max_distance=max_distance)
    observed_edges = edges_on_fixed_graphs(
        used_samples,
        graphs,
        dissimilarity=dissimilarity,
    )
    train_geometry = [observed_edges[s] for s in train]
    eval_geometry = [observed_edges[s] for s in evaluation]
    prepared: PreparedTransfer = prepare_transfer(
        train_geometry,
        eval_geometry,
        bandwidth=bandwidth,
        prior_strength=prior_strength,
        prior_mean=prior_mean,
        segment_points=segment_points,
    )

    observed = prepared.score(
        {s: observed_edges[s].turnover for s in train},
        {s: observed_edges[s].turnover for s in evaluation},
    )

    rng = np.random.default_rng(int(seed))
    null_statistics = np.empty(n_permutations, dtype=float)
    for b in range(n_permutations):
        permuted = [
            permute_trait_within_species(sample_map[s], rng)
            for s in train + evaluation
        ]
        perm_edges = edges_on_fixed_graphs(
            permuted,
            graphs,
            dissimilarity=dissimilarity,
        )
        result = prepared.score(
            {s: perm_edges[s].turnover for s in train},
            {s: perm_edges[s].turnover for s in evaluation},
        )
        null_statistics[b] = result.statistic

    p_value = float(
        (1 + np.count_nonzero(null_statistics >= observed.statistic))
        / (n_permutations + 1)
    )
    return PermutationResult(
        observed=observed,
        null_statistics=null_statistics,
        p_value=p_value,
        null_mean=float(null_statistics.mean()),
        null_sd=float(null_statistics.std(ddof=1)) if n_permutations > 1 else 0.0,
    )
