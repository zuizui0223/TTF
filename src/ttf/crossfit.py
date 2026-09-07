from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .core import SpeciesEdges
from .inference import SpeciesBootstrapResult, centered_species_bootstrap_mean_test
from .transfer import TransferResult, prepare_transfer


def balanced_species_folds(
    labels: Sequence[str],
    *,
    n_folds: int = 3,
    seed: int = 0,
) -> tuple[tuple[str, ...], ...]:
    """Deterministically partition species into balanced disjoint evaluation folds."""
    names = tuple(sorted(map(str, labels)))
    if len(names) < 6 or len(set(names)) != len(names):
        raise ValueError("at least six unique species are required")
    if not 2 <= int(n_folds) <= len(names):
        raise ValueError("n_folds must lie between 2 and the number of species")
    rng = np.random.default_rng(int(seed))
    permuted = [names[int(i)] for i in rng.permutation(len(names))]
    folds = tuple(tuple(permuted[i:: int(n_folds)]) for i in range(int(n_folds)))
    sizes = [len(fold) for fold in folds]
    if max(sizes) - min(sizes) > 1:
        raise RuntimeError("cross-fit folds are not balanced")
    flattened = [name for fold in folds for name in fold]
    if len(flattened) != len(names) or set(flattened) != set(names):
        raise RuntimeError("cross-fit folds do not partition species exactly once")
    return folds


def crossfit_transfer_scores(
    edge_sets: Sequence[SpeciesEdges],
    *,
    n_folds: int = 3,
    fold_seed: int = 0,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
) -> tuple[TransferResult, tuple[tuple[str, ...], ...]]:
    """Score every species once using a field learned only from other species.

    For each evaluation fold the training field is re-fit from all species not
    in that fold.  The final estimand is the macro-average of one held-out score
    per species.  No species ever contributes trait turnover to the field used
    to score itself.
    """
    if len(edge_sets) < 6:
        raise ValueError("at least six species edge sets are required")
    edge_map: Mapping[str, SpeciesEdges] = {edges.species: edges for edges in edge_sets}
    if len(edge_map) != len(edge_sets):
        raise ValueError("species labels must be unique")
    folds = balanced_species_folds(tuple(edge_map), n_folds=n_folds, seed=fold_seed)

    scores: dict[str, float] = {}
    all_names = set(edge_map)
    for evaluation in folds:
        eval_set = set(evaluation)
        train = tuple(sorted(all_names - eval_set))
        if not train:
            raise ValueError("each cross-fit fold must leave training species")
        prepared = prepare_transfer(
            [edge_map[name] for name in train],
            [edge_map[name] for name in evaluation],
            bandwidth=bandwidth,
            prior_strength=prior_strength,
            prior_mean=prior_mean,
            segment_points=segment_points,
        )
        fold_result = prepared.score(
            {name: edge_map[name].turnover for name in train},
            {name: edge_map[name].turnover for name in evaluation},
        )
        for name in evaluation:
            if name in scores:
                raise RuntimeError(f"species scored more than once: {name}")
            scores[name] = float(fold_result.species_scores[name])

    finite = np.asarray([value for value in scores.values() if np.isfinite(value)], dtype=float)
    if len(finite) < 6:
        raise ValueError("fewer than six finite cross-fitted species scores")
    observed = TransferResult(
        statistic=float(finite.mean()),
        species_scores=scores,
        n_eval_species=int(len(finite)),
    )
    return observed, folds


def crossfit_species_bootstrap_test(
    edge_sets: Sequence[SpeciesEdges],
    *,
    n_folds: int = 3,
    fold_seed: int = 0,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
    n_bootstrap: int = 1999,
    seed: int = 0,
) -> tuple[SpeciesBootstrapResult, tuple[tuple[str, ...], ...]]:
    """Centered species-bootstrap inference for cross-fitted TTF scores.

    Because training sets overlap across folds, finite-sample validity is not
    assumed from theory alone.  This estimator must pass the same prospective
    actual-geometry calibration gates before empirical use.
    """
    observed, folds = crossfit_transfer_scores(
        edge_sets,
        n_folds=n_folds,
        fold_seed=fold_seed,
        bandwidth=bandwidth,
        prior_strength=prior_strength,
        prior_mean=prior_mean,
        segment_points=segment_points,
    )
    ordered_scores = np.asarray(
        [observed.species_scores[name] for name in sorted(observed.species_scores)],
        dtype=float,
    )
    ordered_scores = ordered_scores[np.isfinite(ordered_scores)]
    bootstrap = centered_species_bootstrap_mean_test(
        ordered_scores,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    if not np.isclose(bootstrap.observed_mean, observed.statistic, atol=1e-12, rtol=0.0):
        raise RuntimeError("cross-fitted species-score mean drifted from transfer statistic")
    result = SpeciesBootstrapResult(
        observed=observed,
        species_scores=ordered_scores,
        null_statistics=bootstrap.null_means,
        null_studentized=bootstrap.null_studentized,
        p_value=bootstrap.p_value,
        null_mean=float(bootstrap.null_means.mean()),
        null_sd=float(bootstrap.null_means.std(ddof=1)),
        observed_studentized=bootstrap.observed_studentized,
    )
    return result, folds
