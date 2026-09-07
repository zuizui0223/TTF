from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import SpeciesEdges
from .inference import MeanBootstrapResult, centered_species_bootstrap_mean_test
from .mesoscopic import cut_evidence


@dataclass(frozen=True)
class CutPosterior:
    """Species-normalized posterior-like evidence over candidate mesh cuts."""

    species: str
    cuts: np.ndarray
    observable: np.ndarray
    probability: np.ndarray
    private_probability: np.ndarray
    total_log_evidence: np.ndarray
    n_edges: int


@dataclass(frozen=True)
class ConsensusBoundaryField:
    """Species-equal boundary consensus on a common discrete spatial mesh."""

    cuts: np.ndarray
    probability: np.ndarray
    logits: np.ndarray
    prior_strength: float
    n_training_species: int
    n_informative_species: int
    iterations: int


@dataclass(frozen=True)
class ConsensusTransferResult:
    statistic: float
    species_scores: Mapping[str, float]
    n_eval_species: int
    field: ConsensusBoundaryField


@dataclass(frozen=True)
class ConsensusBootstrapResult:
    observed: ConsensusTransferResult
    bootstrap: MeanBootstrapResult

    @property
    def p_value(self) -> float:
        return float(self.bootstrap.p_value)


def _softmax(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) == 0 or not np.isfinite(x).all():
        raise ValueError("softmax requires a non-empty finite vector")
    z = x - float(np.max(x))
    out = np.exp(z)
    out /= float(out.sum())
    return out


def species_cut_posterior(
    edges: SpeciesEdges,
    cuts: Sequence[float] | np.ndarray,
    *,
    evidence_temperature: float = 1.0,
) -> CutPosterior:
    """Convert cut evidence into one normalized distribution per species.

    ``cut_evidence`` stores the profile-likelihood improvement per edge.  Here
    it is multiplied by the species' edge count to recover the corresponding
    total profile-likelihood improvement before normalization.  This lets a
    well-resolved species express a sharper location distribution without ever
    giving that species more than one normalized distribution downstream.

    The private reference distribution is uniform over exactly the cuts that
    this species can contrast with its fixed graph.
    """

    if not np.isfinite(evidence_temperature) or float(evidence_temperature) <= 0:
        raise ValueError("evidence_temperature must be positive and finite")
    ev = cut_evidence(edges, cuts)
    observable = np.asarray(ev.observable, dtype=bool)
    probability = np.zeros(len(ev.cuts), dtype=float)
    private = np.zeros(len(ev.cuts), dtype=float)
    total = np.zeros(len(ev.cuts), dtype=float)
    n_observable = int(np.count_nonzero(observable))
    if n_observable:
        total[observable] = (
            float(evidence_temperature)
            * float(edges.n_edges)
            * np.asarray(ev.log_evidence[observable], dtype=float)
        )
        private[observable] = 1.0 / n_observable
        logits = total[observable] + np.log(private[observable])
        probability[observable] = _softmax(logits)
    return CutPosterior(
        species=edges.species,
        cuts=np.asarray(ev.cuts, dtype=float),
        observable=observable,
        probability=probability,
        private_probability=private,
        total_log_evidence=total,
        n_edges=int(edges.n_edges),
    )


def fit_consensus_boundary_from_posteriors(
    posteriors: Sequence[CutPosterior],
    *,
    prior_strength: float = 1.0,
    tolerance: float = 1e-11,
    max_iterations: int = 20000,
) -> ConsensusBoundaryField:
    """Fit the species-equal observability-conditioned consensus distribution.

    For species s, q_s is compared only with the global field restricted to the
    cuts O_s that species can observe.  The objective is

        sum_s sum_{c in O_s} q_sc log r_sc(p)
        + lambda sum_c (1/C) log p_c,

    where r_s is p renormalized on O_s.  The prior term is one fractional
    uniform pseudo-species.  The objective is concave in global logits, and a
    deterministic gradient ascent is sufficient because the mesh is small.

    Crucially, if q_s is uniform on O_s for every species, p=uniform is an exact
    stationary solution regardless of the missing-cut pattern.
    """

    if len(posteriors) < 2:
        raise ValueError("at least two training posteriors are required")
    if not np.isfinite(prior_strength) or float(prior_strength) <= 0:
        raise ValueError("prior_strength must be positive and finite")
    if not np.isfinite(tolerance) or float(tolerance) <= 0:
        raise ValueError("tolerance must be positive and finite")
    if int(max_iterations) < 1:
        raise ValueError("max_iterations must be positive")

    cuts = np.asarray(posteriors[0].cuts, dtype=float)
    n_cuts = len(cuts)
    if n_cuts < 2:
        raise ValueError("at least two cuts are required")
    seen: set[str] = set()
    informative: list[CutPosterior] = []
    for posterior in posteriors:
        if posterior.species in seen:
            raise ValueError("training species must be unique")
        seen.add(posterior.species)
        if posterior.cuts.shape != cuts.shape or not np.allclose(
            posterior.cuts, cuts, atol=0.0, rtol=0.0
        ):
            raise ValueError("all species must use the same candidate cuts")
        mask = np.asarray(posterior.observable, dtype=bool)
        if mask.shape != (n_cuts,):
            raise ValueError("posterior observability shape mismatch")
        if int(np.count_nonzero(mask)) < 2:
            continue
        q = np.asarray(posterior.probability[mask], dtype=float)
        if not np.isfinite(q).all() or np.any(q < 0) or not np.isclose(
            q.sum(), 1.0, atol=1e-10, rtol=0.0
        ):
            raise ValueError("invalid species cut posterior")
        informative.append(posterior)
    if len(informative) < 2:
        raise ValueError("at least two informative training species are required")

    theta = np.zeros(n_cuts, dtype=float)
    uniform = np.full(n_cuts, 1.0 / n_cuts, dtype=float)
    denominator = float(len(informative)) + float(prior_strength)
    iterations = 0
    for iteration in range(1, int(max_iterations) + 1):
        global_probability = _softmax(theta)
        gradient = float(prior_strength) * (uniform - global_probability)
        for posterior in informative:
            mask = posterior.observable
            restricted = _softmax(theta[mask])
            gradient[mask] += posterior.probability[mask] - restricted
        # Logits are identifiable only up to an additive constant.  Every
        # objective component has zero-sum gradient analytically; recentering
        # removes numerical drift without changing the fitted probabilities.
        gradient -= float(gradient.mean())
        updated = theta + gradient / denominator
        updated -= float(updated.mean())
        iterations = iteration
        if float(np.max(np.abs(updated - theta))) <= float(tolerance):
            theta = updated
            break
        theta = updated
    else:
        raise RuntimeError("mesoscopic consensus optimizer did not converge")

    return ConsensusBoundaryField(
        cuts=cuts,
        probability=_softmax(theta),
        logits=theta,
        prior_strength=float(prior_strength),
        n_training_species=len(posteriors),
        n_informative_species=len(informative),
        iterations=int(iterations),
    )


def fit_consensus_boundary(
    train_edges: Sequence[SpeciesEdges],
    cuts: Sequence[float] | np.ndarray,
    *,
    prior_strength: float = 1.0,
    evidence_temperature: float = 1.0,
) -> ConsensusBoundaryField:
    posteriors = [
        species_cut_posterior(
            edges,
            cuts,
            evidence_temperature=evidence_temperature,
        )
        for edges in train_edges
    ]
    return fit_consensus_boundary_from_posteriors(
        posteriors,
        prior_strength=prior_strength,
    )


def consensus_species_score_from_posterior(
    field: ConsensusBoundaryField,
    posterior: CutPosterior,
) -> float:
    """Expected held-out log-score gain over its geometry-private prior."""

    if posterior.cuts.shape != field.cuts.shape or not np.allclose(
        posterior.cuts, field.cuts, atol=0.0, rtol=0.0
    ):
        raise ValueError("posterior and field cuts differ")
    mask = np.asarray(posterior.observable, dtype=bool)
    n_observable = int(np.count_nonzero(mask))
    if n_observable < 2:
        return float("nan")
    q = np.asarray(posterior.probability[mask], dtype=float)
    u = np.asarray(posterior.private_probability[mask], dtype=float)
    predicted = np.asarray(field.probability[mask], dtype=float)
    predicted /= float(predicted.sum())
    if np.any(predicted <= 0) or np.any(u <= 0):
        raise RuntimeError("consensus score requires positive predictive probabilities")
    return float(np.sum(q * (np.log(predicted) - np.log(u))))


def consensus_species_score(
    field: ConsensusBoundaryField,
    edges: SpeciesEdges,
    *,
    evidence_temperature: float = 1.0,
) -> float:
    posterior = species_cut_posterior(
        edges,
        field.cuts,
        evidence_temperature=evidence_temperature,
    )
    return consensus_species_score_from_posterior(field, posterior)


def consensus_transfer(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    cuts: Sequence[float] | np.ndarray,
    *,
    prior_strength: float = 1.0,
    evidence_temperature: float = 1.0,
) -> ConsensusTransferResult:
    field = fit_consensus_boundary(
        train_edges,
        cuts,
        prior_strength=prior_strength,
        evidence_temperature=evidence_temperature,
    )
    scores: dict[str, float] = {}
    for edges in eval_edges:
        if edges.species in scores:
            raise ValueError("evaluation species must be unique")
        scores[edges.species] = consensus_species_score(
            field,
            edges,
            evidence_temperature=evidence_temperature,
        )
    finite = np.asarray([v for v in scores.values() if np.isfinite(v)], dtype=float)
    if len(finite) < 6:
        raise ValueError("at least six finite held-out consensus scores are required")
    return ConsensusTransferResult(
        statistic=float(finite.mean()),
        species_scores=scores,
        n_eval_species=int(len(finite)),
        field=field,
    )


def consensus_bootstrap_test(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    cuts: Sequence[float] | np.ndarray,
    *,
    prior_strength: float = 1.0,
    evidence_temperature: float = 1.0,
    n_bootstrap: int = 1999,
    seed: int = 0,
) -> ConsensusBootstrapResult:
    observed = consensus_transfer(
        train_edges,
        eval_edges,
        cuts,
        prior_strength=prior_strength,
        evidence_temperature=evidence_temperature,
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
        raise RuntimeError("consensus species-score mean drift")
    return ConsensusBootstrapResult(observed=observed, bootstrap=bootstrap)
