from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .core import SpeciesEdges
from .mesoscopic import (
    CutEvidence,
    MesoscopicBoundaryField,
    _logsumexp,
    _validate_cuts,
    cut_evidence,
)


@dataclass(frozen=True)
class MesoscopicAlignmentRandomizationResult:
    """Conditional test of cross-species boundary-location alignment.

    Held-out species evidence is frozen.  Within every training species, the
    multiset of boundary-evidence strengths is retained but its labels are
    independently permuted among exactly the cuts that species can observe.
    The training field is re-fit for every randomization.
    """

    observed_statistic: float
    observed_species_scores: np.ndarray
    null_statistics: np.ndarray
    p_value: float
    null_mean: float
    null_sd: float
    observed_field: MesoscopicBoundaryField
    n_eval_species: int


def fit_mesoscopic_boundary_from_evidence(
    evidences: Sequence[CutEvidence],
    *,
    prior_strength: float = 1.0,
) -> MesoscopicBoundaryField:
    if len(evidences) < 2:
        raise ValueError("at least two training species are required")
    if not np.isfinite(prior_strength) or float(prior_strength) < 0:
        raise ValueError("prior_strength must be finite and non-negative")

    cuts = _validate_cuts(evidences[0].cuts)
    total = np.zeros(len(cuts), dtype=float)
    count = np.zeros(len(cuts), dtype=int)
    seen: set[str] = set()
    for evidence in evidences:
        if evidence.species in seen:
            raise ValueError("training species must be unique")
        seen.add(evidence.species)
        if evidence.cuts.shape != cuts.shape or not np.allclose(
            evidence.cuts, cuts, atol=0.0, rtol=0.0
        ):
            raise ValueError("all species must use the same candidate cuts")
        mask = np.asarray(evidence.observable, dtype=bool)
        values = np.asarray(evidence.log_evidence, dtype=float)
        if mask.shape != cuts.shape or values.shape != cuts.shape:
            raise ValueError("cut evidence shape mismatch")
        if not np.isfinite(values).all():
            raise ValueError("cut evidence must be finite")
        total[mask] += values[mask]
        count[mask] += 1

    mean = np.zeros(len(cuts), dtype=float)
    denom = count.astype(float) + float(prior_strength)
    active = denom > 0
    mean[active] = total[active] / denom[active]
    logits = mean - float(np.max(mean))
    probability = np.exp(logits)
    probability /= float(probability.sum())
    return MesoscopicBoundaryField(
        cuts=cuts,
        probability=probability,
        mean_log_evidence=mean,
        opportunity_species=count,
        prior_strength=float(prior_strength),
    )


def permute_cut_evidence_labels(
    evidence: CutEvidence,
    rng: np.random.Generator,
) -> CutEvidence:
    """Relabel evidence only within a species' observable boundary elements."""

    mask = np.asarray(evidence.observable, dtype=bool)
    ids = np.flatnonzero(mask)
    values = np.zeros_like(np.asarray(evidence.log_evidence, dtype=float))
    if len(ids):
        values[ids] = np.asarray(evidence.log_evidence, dtype=float)[rng.permutation(ids)]
    return CutEvidence(
        species=evidence.species,
        cuts=np.asarray(evidence.cuts, dtype=float),
        observable=mask.copy(),
        log_evidence=values,
        n_edges=int(evidence.n_edges),
    )


def mesoscopic_species_score_from_evidence(
    field: MesoscopicBoundaryField,
    evidence: CutEvidence,
) -> float:
    if evidence.cuts.shape != field.cuts.shape or not np.allclose(
        evidence.cuts, field.cuts, atol=0.0, rtol=0.0
    ):
        raise ValueError("evidence and field cuts differ")
    mask = np.asarray(evidence.observable, dtype=bool)
    n = int(np.count_nonzero(mask))
    if n < 2:
        return float("nan")
    log_like = np.asarray(evidence.log_evidence[mask], dtype=float)
    p = np.asarray(field.probability[mask], dtype=float)
    p /= float(p.sum())
    shared_log = _logsumexp(np.log(p) + log_like)
    private_log = _logsumexp(log_like) - float(np.log(n))
    return float(shared_log - private_log)


def _score_eval_evidence(
    field: MesoscopicBoundaryField,
    eval_evidence: Sequence[CutEvidence],
) -> tuple[float, np.ndarray]:
    scores = np.asarray(
        [mesoscopic_species_score_from_evidence(field, evidence) for evidence in eval_evidence],
        dtype=float,
    )
    scores = scores[np.isfinite(scores)]
    if len(scores) < 6:
        raise ValueError("at least six finite held-out mesoscopic scores are required")
    return float(scores.mean()), scores


def mesoscopic_alignment_randomization_test(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    cuts: Sequence[float] | np.ndarray,
    *,
    prior_strength: float = 1.0,
    n_randomizations: int = 999,
    seed: int = 0,
) -> MesoscopicAlignmentRandomizationResult:
    """Test transfer beyond finite-mesh chance commonness.

    Null assumption: conditional on each species' fixed graph, observability set,
    and multiset of boundary-evidence strengths, boundary-element labels are
    exchangeable within each training species under no cross-species shared
    location.  Evaluation evidence is never permuted.
    """

    cuts_array = _validate_cuts(cuts)
    if int(n_randomizations) < 99:
        raise ValueError("n_randomizations must be >= 99")
    train_evidence = [cut_evidence(edges, cuts_array) for edges in train_edges]
    eval_evidence = [cut_evidence(edges, cuts_array) for edges in eval_edges]

    observed_field = fit_mesoscopic_boundary_from_evidence(
        train_evidence,
        prior_strength=prior_strength,
    )
    observed_statistic, observed_scores = _score_eval_evidence(
        observed_field,
        eval_evidence,
    )

    rng = np.random.default_rng(int(seed))
    null = np.empty(int(n_randomizations), dtype=float)
    for b in range(int(n_randomizations)):
        randomized = [
            permute_cut_evidence_labels(evidence, rng)
            for evidence in train_evidence
        ]
        field = fit_mesoscopic_boundary_from_evidence(
            randomized,
            prior_strength=prior_strength,
        )
        null[b], _ = _score_eval_evidence(field, eval_evidence)

    p_value = float(
        (1 + np.count_nonzero(null >= observed_statistic))
        / (int(n_randomizations) + 1)
    )
    return MesoscopicAlignmentRandomizationResult(
        observed_statistic=float(observed_statistic),
        observed_species_scores=observed_scores,
        null_statistics=null,
        p_value=p_value,
        null_mean=float(null.mean()),
        null_sd=float(null.std(ddof=1)),
        observed_field=observed_field,
        n_eval_species=int(len(observed_scores)),
    )
