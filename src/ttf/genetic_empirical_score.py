from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .chunked_transfer import score_chunked_batch
from .core import rank01, spearman_rho
from .genetic_gate import GeneticTTFDesign, GeneticWorldScore
from .genetic_ibd import CrossfitIBDResult, crossfit_ibd_residuals_prepared
from .genetic_self_detectability import GeneticSelfDetectabilityDesign
from .geometry_control import length_orthogonalized_turnover
from .private_strength import training_private_strength_from_indices


@dataclass(frozen=True)
class EmpiricalSelfScore:
    statistic: float
    species_scores: Mapping[str, float]


@dataclass(frozen=True)
class DescriptiveTotalTransfer:
    """Unqualified raw-distance summary; never an input to primary inference."""

    statistic: float | None
    species_scores: Mapping[str, float | None]

    def as_dict(self) -> dict:
        undefined = sorted(name for name, score in self.species_scores.items() if score is None)
        return {
            "name": "total_genetic_distance_transfer",
            "status": "DESCRIPTIVE_ONLY" if not undefined else "NOT_EVALUABLE_DESCRIPTIVE",
            "statistic": self.statistic,
            "heldout_species_scores": dict(sorted(self.species_scores.items())),
            "frozen_evaluation_species_count": len(self.species_scores),
            "finite_evaluation_species_count": len(self.species_scores) - len(undefined),
            "undefined_species": undefined,
            "undefined_reason": "nonfinite_score" if undefined else None,
            "constant_vector_score": 0.0,
            "biological_ibd_adjustment": False,
            "training_geometry_control": "within_species_edge_length_rank_orthogonalization",
            "same_primary_geometry_split_kernel": True,
            "inferentially_qualified": False,
            "used_for_primary_decision": False,
            "claim_limit": "Descriptive only; no significance test or rescue of the primary post-IBD decision.",
        }


def _validated_distance_vector(
    design: GeneticTTFDesign,
    genetic_distance: Mapping[str, np.ndarray],
    name: str,
) -> np.ndarray:
    values = np.asarray(genetic_distance[name], dtype=float)
    expected = len(design.template_edges[name].nodes)
    if values.ndim != 1 or len(values) != expected:
        raise ValueError(
            f"genetic distance vector length drift for {name}: {values.shape!r}, expected {expected}"
        )
    if not np.isfinite(values).all() or np.any(values < 0.0):
        raise ValueError(f"genetic distances must be finite and non-negative for {name}")
    return values


def score_genetic_distance_mapping(
    design: GeneticTTFDesign,
    genetic_distance: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> GeneticWorldScore:
    """Score one empirical genetic-distance mapping using the frozen genetic TTF design.

    This is deliberately a thin adapter around the already-qualified numerical
    operations used by ``score_genetic_world``. It accepts only one vector of
    genetic distances per frozen species graph and introduces no empirical
    tuning, resplitting, edge filtering, or alternative nuisance model.
    """
    used = design.train_species + design.eval_species
    missing = set(used) - set(map(str, genetic_distance.keys()))
    extra = set(map(str, genetic_distance.keys())) - set(used)
    if missing:
        raise ValueError(f"genetic distance mapping is missing species: {sorted(missing)}")
    if extra:
        raise ValueError(f"genetic distance mapping has unexpected species: {sorted(extra)}")

    ibd: dict[str, CrossfitIBDResult] = {}
    train_response: dict[str, np.ndarray] = {}
    eval_response: dict[str, np.ndarray] = {}

    for name in used:
        values = _validated_distance_vector(design, genetic_distance, name)
        result = crossfit_ibd_residuals_prepared(values, design.ibd_designs[name])
        ibd[name] = result
        edges = design.template_edges[name]
        if name in design.train_species:
            train_response[name] = length_orthogonalized_turnover(
                result.residual_turnover,
                edges.length,
            )
        else:
            eval_response[name] = result.residual_turnover

    strength, _ = training_private_strength_from_indices(
        train_response,
        design.strength_indices,
        design.train_species,
    )
    scored = score_chunked_batch(
        design.prepared,
        {name: train_response[name][:, None] for name in design.train_species},
        {name: eval_response[name][:, None] for name in design.eval_species},
        edge_chunk_size=int(edge_chunk_size),
        train_chunk_size=int(train_chunk_size),
    )
    species_scores = {
        name: float(scored.species_scores[name][0]) for name in design.eval_species
    }
    if not np.isfinite(scored.statistics[0]) or not np.isfinite(strength):
        raise RuntimeError("non-finite empirical genetic TTF score")
    if any(not np.isfinite(value) for value in species_scores.values()):
        raise RuntimeError("non-finite empirical held-out species score")
    return GeneticWorldScore(
        statistic=float(scored.statistics[0]),
        training_strength=float(strength),
        species_scores=species_scores,
        ibd=ibd,
    )


def score_self_detectability_mapping(
    design: GeneticTTFDesign,
    self_design: GeneticSelfDetectabilityDesign,
    genetic_distance: Mapping[str, np.ndarray],
) -> EmpiricalSelfScore:
    """Score empirical within-species self-predictability on frozen eval species.

    Endpoint-safe post-IBD turnover is projected by the same precomputed self-field
    used in synthetic qualification. No empirical threshold, species filtering,
    bandwidth selection, or response-dependent geometry enters this adapter.
    """
    if tuple(self_design.species) != tuple(design.eval_species):
        raise ValueError("empirical self design must cover the exact frozen evaluation species")
    species_scores: dict[str, float] = {}
    for name in self_design.species:
        if name not in genetic_distance:
            raise ValueError(f"genetic distance mapping is missing self species: {name}")
        values = _validated_distance_vector(design, genetic_distance, name)
        result = crossfit_ibd_residuals_prepared(values, design.ibd_designs[name])
        response = np.asarray(result.residual_turnover, dtype=float)
        prepared = self_design.by_species[name]
        predicted = prepared.prior_offset + prepared.projection @ response
        score = float(spearman_rho(predicted, response))
        if not np.isfinite(score):
            raise RuntimeError(f"non-finite empirical self score for {name}")
        species_scores[name] = score
    statistic = float(np.mean(np.asarray(list(species_scores.values()), dtype=float)))
    if not np.isfinite(statistic):
        raise RuntimeError("non-finite empirical self-detectability statistic")
    return EmpiricalSelfScore(statistic=statistic, species_scores=species_scores)


def score_total_genetic_distance_mapping(
    design: GeneticTTFDesign,
    genetic_distance: Mapping[str, np.ndarray],
    *,
    edge_chunk_size: int = 32,
    train_chunk_size: int = 4096,
) -> DescriptiveTotalTransfer:
    """Describe transfer before biological IBD adjustment on the frozen design.

    Retain the primary field, split, species weights and training-only geometric
    length control. Replace post-IBD turnover with raw distance rank. No primary
    reference distribution is valid for this different response, so this adapter
    supplies no p-value. The inherited scorer assigns constant vectors zero.
    Non-finite scores retain their species as null and
    make the aggregate undefined, rather than changing the evaluation population.
    Invalid inputs still raise; only a non-finite score is reported as unavailable.
    """
    used = design.train_species + design.eval_species
    if set(genetic_distance) != set(used):
        raise ValueError("total distance mapping must cover the exact frozen species")
    ranks = {
        name: rank01(_validated_distance_vector(design, genetic_distance, name))
        for name in used
    }
    train = {
        name: length_orthogonalized_turnover(ranks[name], design.template_edges[name].length)[:, None]
        for name in design.train_species
    }
    try:
        scored = score_chunked_batch(
            design.prepared,
            train,
            {name: ranks[name][:, None] for name in design.eval_species},
            edge_chunk_size=edge_chunk_size,
            train_chunk_size=train_chunk_size,
        )
    except ValueError as exc:
        # The inherited scorer rejects an all-undefined batch. For this one-world
        # descriptive call that means every frozen evaluation score is undefined.
        if str(exc) != "one or more worlds had no finite held-out species scores":
            raise
        return DescriptiveTotalTransfer(None, {name: None for name in design.eval_species})
    scores = {
        name: float(scored.species_scores[name][0])
        if np.isfinite(scored.species_scores[name][0]) else None
        for name in design.eval_species
    }
    statistic = (
        float(scored.statistics[0])
        if all(value is not None for value in scores.values()) else None
    )
    return DescriptiveTotalTransfer(statistic, scores)


__all__ = [
    "DescriptiveTotalTransfer",
    "EmpiricalSelfScore",
    "score_genetic_distance_mapping",
    "score_self_detectability_mapping",
    "score_total_genetic_distance_mapping",
]
