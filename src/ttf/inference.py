from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import Dissimilarity, SpeciesSample
from .nulls import edges_on_fixed_graphs, fixed_graphs
from .transfer import PreparedTransfer, TransferResult, prepare_transfer


@dataclass(frozen=True)
class SpeciesBootstrapResult:
    """Conditional inference using held-out species as the sampling units.

    The training boundary field is learned once.  No trait values, coordinates,
    graphs, or spatial patterns are permuted in evaluation species.  Inference
    concerns the superpopulation mean held-out transfer score conditional on the
    frozen training field.
    """

    observed: TransferResult
    species_scores: np.ndarray
    null_statistics: np.ndarray
    null_studentized: np.ndarray
    p_value: float
    null_mean: float
    null_sd: float
    observed_studentized: float


@dataclass(frozen=True)
class MeanBootstrapResult:
    observed_mean: float
    observed_studentized: float
    null_means: np.ndarray
    null_studentized: np.ndarray
    p_value: float


def centered_species_bootstrap_mean_test(
    scores: Sequence[float] | np.ndarray,
    *,
    n_bootstrap: int = 1999,
    seed: int = 0,
) -> MeanBootstrapResult:
    """One-sided centered, studentized bootstrap test of E[score] <= 0.

    Centering imposes the null mean while retaining the empirical shape,
    variance, skewness, and bounded support of held-out species scores.  The
    bootstrap resamples *species*, never observations within a species.
    """
    x = np.asarray(scores, dtype=float)
    x = x[np.isfinite(x)]
    if x.ndim != 1 or len(x) < 6:
        raise ValueError("at least six finite held-out species scores are required")
    if n_bootstrap < 99:
        raise ValueError("n_bootstrap must be >= 99 for inference")

    observed_mean = float(x.mean())
    observed_sd = float(x.std(ddof=1))
    scale = observed_sd / np.sqrt(len(x))
    if scale <= np.finfo(float).tiny:
        observed_t = float("inf") if observed_mean > 0 else 0.0
    else:
        observed_t = float(observed_mean / scale)

    centered = x - observed_mean
    rng = np.random.default_rng(int(seed))
    indices = rng.integers(0, len(centered), size=(int(n_bootstrap), len(centered)))
    draws = centered[indices]
    means = draws.mean(axis=1)
    sds = draws.std(axis=1, ddof=1)
    ses = sds / np.sqrt(len(centered))
    studentized = np.zeros(int(n_bootstrap), dtype=float)
    valid = ses > np.finfo(float).tiny
    studentized[valid] = means[valid] / ses[valid]
    studentized[~valid & (means > 0)] = np.inf
    studentized[~valid & (means < 0)] = -np.inf

    p_value = float(
        (1 + np.count_nonzero(studentized >= observed_t))
        / (int(n_bootstrap) + 1)
    )
    return MeanBootstrapResult(
        observed_mean=observed_mean,
        observed_studentized=observed_t,
        null_means=means,
        null_studentized=studentized,
        p_value=p_value,
    )


def heldout_species_bootstrap_from_prepared(
    prepared: PreparedTransfer,
    *,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    n_bootstrap: int = 1999,
    seed: int = 0,
) -> SpeciesBootstrapResult:
    """Run sharedness inference when the geometry projection is already frozen.

    Gate-I reuses one empirical sampling geometry across hundreds of synthetic
    trait worlds.  The kNN graphs and kernel projection therefore belong to the
    fixed design, not to any one simulated response.  This helper keeps those
    geometry-only quantities frozen while allowing every world's turnover ranks
    to change.
    """
    observed = prepared.score(train_turnover, eval_turnover)
    scores = np.asarray(
        [
            observed.species_scores[s]
            for s in prepared.eval_species
            if np.isfinite(observed.species_scores[s])
        ],
        dtype=float,
    )
    if len(scores) < 6:
        raise ValueError("fewer than six finite held-out species scores")

    bootstrap = centered_species_bootstrap_mean_test(
        scores,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    if not np.isclose(
        bootstrap.observed_mean,
        observed.statistic,
        atol=1e-12,
        rtol=0.0,
    ):
        raise RuntimeError("species-score mean drifted from primary transfer statistic")

    return SpeciesBootstrapResult(
        observed=observed,
        species_scores=scores,
        null_statistics=bootstrap.null_means,
        null_studentized=bootstrap.null_studentized,
        p_value=bootstrap.p_value,
        null_mean=float(bootstrap.null_means.mean()),
        null_sd=float(bootstrap.null_means.std(ddof=1)),
        observed_studentized=bootstrap.observed_studentized,
    )


def heldout_species_bootstrap_test(
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
    n_bootstrap: int = 1999,
    seed: int = 0,
    dissimilarity: Dissimilarity | None = None,
) -> SpeciesBootstrapResult:
    """Test positive transferable sharedness using untouched held-out species.

    This is the sharedness-specific inference layer.  It differs deliberately
    from within-species trait permutation: the latter tests trait-location
    exchangeability and can be anti-conservative for sharedness when every
    species has a strong but private spatial transition.
    """
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
    if len(evaluation) < 6:
        raise ValueError("at least six held-out species are required")

    used_samples = [sample_map[s] for s in train + evaluation]
    graphs = fixed_graphs(used_samples, k=k, max_distance=max_distance)
    edge_map = edges_on_fixed_graphs(
        used_samples,
        graphs,
        dissimilarity=dissimilarity,
    )
    prepared = prepare_transfer(
        [edge_map[s] for s in train],
        [edge_map[s] for s in evaluation],
        bandwidth=bandwidth,
        prior_strength=prior_strength,
        prior_mean=prior_mean,
        segment_points=segment_points,
    )
    return heldout_species_bootstrap_from_prepared(
        prepared,
        train_turnover={s: edge_map[s].turnover for s in train},
        eval_turnover={s: edge_map[s].turnover for s in evaluation},
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
