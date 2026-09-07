from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import SpeciesEdges
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test


@dataclass(frozen=True)
class CutEvidence:
    species: str
    cuts: np.ndarray
    observable: np.ndarray
    log_evidence: np.ndarray
    n_edges: int


@dataclass(frozen=True)
class MesoscopicBoundaryField:
    cuts: np.ndarray
    probability: np.ndarray
    mean_log_evidence: np.ndarray
    opportunity_species: np.ndarray
    prior_strength: float


@dataclass(frozen=True)
class MesoscopicTransferResult:
    statistic: float
    species_scores: Mapping[str, float]
    n_eval_species: int
    field: MesoscopicBoundaryField


@dataclass(frozen=True)
class MesoscopicBootstrapResult:
    observed: MesoscopicTransferResult
    bootstrap: MeanBootstrapResult

    @property
    def p_value(self) -> float:
        return float(self.bootstrap.p_value)


def _logsumexp(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) == 0:
        raise ValueError("logsumexp requires a non-empty vector")
    m = float(np.max(x))
    return m + float(np.log(np.sum(np.exp(x - m))))


def _validate_cuts(cuts: Sequence[float] | np.ndarray) -> np.ndarray:
    c = np.asarray(cuts, dtype=float)
    if c.ndim != 1 or len(c) < 2 or not np.isfinite(c).all():
        raise ValueError("cuts must be a finite vector with at least two positions")
    if np.any(np.diff(c) <= 0):
        raise ValueError("cuts must be strictly increasing")
    return c


def cut_evidence(
    edges: SpeciesEdges,
    cuts: Sequence[float] | np.ndarray,
    *,
    eps: float = 1e-8,
    max_log_evidence: float = 12.0,
) -> CutEvidence:
    """Measure species-specific evidence for each candidate boundary cut.

    Evidence is an average Gaussian log-likelihood improvement per edge of a
    binary "edge crosses cut" predictor over an intercept-only model.  It is
    based only on within-species rank turnover.  Cuts outside the species'
    geographic span are unobservable and contribute no information.

    Dividing by edge count prevents densely sampled species from dominating the
    shared field.  The cap protects zero-noise/small-SSE diagnostics from
    numerical domination without changing the ordering of ordinary evidence.
    """
    c = _validate_cuts(cuts)
    y = np.asarray(edges.turnover, dtype=float)
    if y.ndim != 1 or len(y) != edges.n_edges or not np.isfinite(y).all():
        raise ValueError("edge turnover must be finite")
    start_x = np.asarray(edges.start[:, 0], dtype=float)
    end_x = np.asarray(edges.end[:, 0], dtype=float)
    node_min = float(min(np.min(start_x), np.min(end_x)))
    node_max = float(max(np.max(start_x), np.max(end_x)))

    centered = y - float(y.mean())
    sse0 = float(np.dot(centered, centered))
    evidence = np.zeros(len(c), dtype=float)
    observable = np.zeros(len(c), dtype=bool)
    if sse0 <= float(eps):
        return CutEvidence(edges.species, c, observable, evidence, edges.n_edges)

    for j, cut in enumerate(c):
        if not (node_min < float(cut) < node_max):
            continue
        crossing = ((start_x < cut) & (end_x > cut)) | ((end_x < cut) & (start_x > cut))
        n_cross = int(np.count_nonzero(crossing))
        n_non = int(len(crossing) - n_cross)
        if n_cross == 0 or n_non == 0:
            # Geometrically within range but this graph provides no contrast at
            # the cut.  Treat as uninformative rather than negative evidence.
            continue
        observable[j] = True
        mean_cross = float(y[crossing].mean())
        mean_non = float(y[~crossing].mean())
        fitted = np.where(crossing, mean_cross, mean_non)
        residual = y - fitted
        sse1 = float(np.dot(residual, residual))
        # Gaussian profile likelihood improvement, normalized per edge.
        gain = 0.5 * np.log((sse0 + float(eps)) / (sse1 + float(eps)))
        evidence[j] = float(np.clip(gain, 0.0, float(max_log_evidence)))

    return CutEvidence(
        species=edges.species,
        cuts=c,
        observable=observable,
        log_evidence=evidence,
        n_edges=edges.n_edges,
    )


def fit_mesoscopic_boundary(
    train_edges: Sequence[SpeciesEdges],
    cuts: Sequence[float] | np.ndarray,
    *,
    prior_strength: float = 1.0,
) -> MesoscopicBoundaryField:
    """Learn a species-balanced predictive distribution over boundary cuts.

    Each cut receives the mean per-species log evidence among species that can
    actually contrast that cut, shrunk toward zero by ``prior_strength`` pseudo
    species.  A softmax converts this opportunity-corrected evidence to a
    predictive boundary-location distribution.
    """
    c = _validate_cuts(cuts)
    if len(train_edges) < 2:
        raise ValueError("at least two training species are required")
    if float(prior_strength) < 0:
        raise ValueError("prior_strength must be non-negative")

    total = np.zeros(len(c), dtype=float)
    count = np.zeros(len(c), dtype=int)
    seen = set()
    for edges in train_edges:
        if edges.species in seen:
            raise ValueError("training species must be unique")
        seen.add(edges.species)
        ev = cut_evidence(edges, c)
        total[ev.observable] += ev.log_evidence[ev.observable]
        count[ev.observable] += 1

    mean = np.zeros(len(c), dtype=float)
    denom = count.astype(float) + float(prior_strength)
    active = denom > 0
    mean[active] = total[active] / denom[active]

    # A cut observed by no training species must remain at the geometry-only
    # prior rather than becoming an artificial high-evidence refuge.
    logits = mean.copy()
    logits -= float(np.max(logits))
    probability = np.exp(logits)
    probability /= float(probability.sum())
    return MesoscopicBoundaryField(
        cuts=c,
        probability=probability,
        mean_log_evidence=mean,
        opportunity_species=count,
        prior_strength=float(prior_strength),
    )


def mesoscopic_species_score(
    field: MesoscopicBoundaryField,
    edges: SpeciesEdges,
) -> float:
    """Held-out predictive log-score gain over a geometry-only private cut.

    The shared model uses the training boundary distribution restricted to cuts
    the held-out species can evaluate.  The null/private model is uniform over
    exactly those same observable cuts.  Therefore finite candidate-boundary
    chance commonness is part of the baseline by construction.
    """
    ev = cut_evidence(edges, field.cuts)
    mask = ev.observable
    n = int(np.count_nonzero(mask))
    if n < 2:
        return float("nan")
    log_like = ev.log_evidence[mask]
    p = np.asarray(field.probability[mask], dtype=float)
    p /= float(p.sum())
    shared_log = _logsumexp(np.log(p) + log_like)
    private_log = _logsumexp(log_like) - float(np.log(n))
    return float(shared_log - private_log)


def mesoscopic_transfer(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    cuts: Sequence[float] | np.ndarray,
    *,
    prior_strength: float = 1.0,
) -> MesoscopicTransferResult:
    field = fit_mesoscopic_boundary(train_edges, cuts, prior_strength=prior_strength)
    scores: dict[str, float] = {}
    for edges in eval_edges:
        if edges.species in scores:
            raise ValueError("evaluation species must be unique")
        scores[edges.species] = mesoscopic_species_score(field, edges)
    finite = np.asarray([v for v in scores.values() if np.isfinite(v)], dtype=float)
    if len(finite) < 6:
        raise ValueError("at least six finite held-out mesoscopic scores are required")
    return MesoscopicTransferResult(
        statistic=float(finite.mean()),
        species_scores=scores,
        n_eval_species=len(finite),
        field=field,
    )


def mesoscopic_bootstrap_test(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    cuts: Sequence[float] | np.ndarray,
    *,
    prior_strength: float = 1.0,
    n_bootstrap: int = 1999,
    seed: int = 0,
) -> MesoscopicBootstrapResult:
    observed = mesoscopic_transfer(
        train_edges,
        eval_edges,
        cuts,
        prior_strength=prior_strength,
    )
    finite = np.asarray(
        [v for v in observed.species_scores.values() if np.isfinite(v)],
        dtype=float,
    )
    bootstrap = centered_species_bootstrap_mean_test(
        finite,
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    if not np.isclose(bootstrap.observed_mean, observed.statistic, atol=1e-12, rtol=0.0):
        raise RuntimeError("mesoscopic species-score mean drift")
    return MesoscopicBootstrapResult(observed=observed, bootstrap=bootstrap)
