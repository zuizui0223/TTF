from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .core import SpeciesEdges, spearman_rho


@dataclass(frozen=True)
class FieldEstimate:
    points: np.ndarray
    value: np.ndarray
    opportunity: np.ndarray


@dataclass(frozen=True)
class TransferResult:
    """Primary TTF estimand: macro-average held-out species transfer skill."""

    statistic: float
    species_scores: Mapping[str, float]
    n_eval_species: int


@dataclass(frozen=True)
class KernelBoundaryModel:
    """Opportunity-corrected continuous boundary field learned from train species."""

    positions: np.ndarray
    values: np.ndarray
    weights: np.ndarray
    bandwidth: float
    prior_strength: float = 0.25
    prior_mean: float = 0.5

    def __post_init__(self) -> None:
        x = np.asarray(self.positions, dtype=float)
        y = np.asarray(self.values, dtype=float)
        w = np.asarray(self.weights, dtype=float)
        if x.ndim != 2 or y.shape != (len(x),) or w.shape != (len(x),):
            raise ValueError("positions, values, and weights have incompatible shapes")
        if len(x) == 0 or not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError("model inputs must be finite and non-empty")
        if np.any(w <= 0) or not np.isfinite(w).all():
            raise ValueError("weights must be finite and positive")
        if self.bandwidth <= 0 or self.prior_strength < 0:
            raise ValueError("bandwidth must be positive and prior_strength non-negative")
        object.__setattr__(self, "positions", x)
        object.__setattr__(self, "values", y)
        object.__setattr__(self, "weights", w)

    def _kernel(self, points: np.ndarray) -> np.ndarray:
        q = np.asarray(points, dtype=float)
        if q.ndim != 2 or q.shape[1] != self.positions.shape[1]:
            raise ValueError("query points have incompatible dimensionality")
        delta = q[:, None, :] - self.positions[None, :, :]
        distance2 = np.sum(delta * delta, axis=2)
        return np.exp(-0.5 * distance2 / (self.bandwidth * self.bandwidth)) * self.weights[None, :]

    def predict_points(self, points: np.ndarray) -> FieldEstimate:
        q = np.asarray(points, dtype=float)
        kernel = self._kernel(q)
        opportunity = kernel.sum(axis=1)
        numerator = kernel @ self.values
        denominator = opportunity + float(self.prior_strength)
        value = (
            numerator + float(self.prior_strength) * float(self.prior_mean)
        ) / np.maximum(denominator, np.finfo(float).tiny)
        return FieldEstimate(points=q, value=value, opportunity=opportunity)

    def predict_segments(
        self,
        start: np.ndarray,
        end: np.ndarray,
        *,
        n_points: int = 5,
    ) -> np.ndarray:
        """Approximate mean boundary exposure along each held-out edge."""
        a = np.asarray(start, dtype=float)
        b = np.asarray(end, dtype=float)
        if a.shape != b.shape or a.ndim != 2:
            raise ValueError("start and end must be equal-shaped m x d arrays")
        if n_points < 1:
            raise ValueError("n_points must be >= 1")
        t = (np.arange(n_points, dtype=float) + 0.5) / n_points
        points = a[:, None, :] + t[None, :, None] * (b - a)[:, None, :]
        predicted = self.predict_points(points.reshape(-1, a.shape[1])).value
        return predicted.reshape(len(a), n_points).mean(axis=1)


def fit_boundary_model(
    edge_sets: Sequence[SpeciesEdges],
    *,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
) -> KernelBoundaryModel:
    """Fit with equal total weight per training species."""
    if len(edge_sets) == 0:
        raise ValueError("at least one training species is required")
    if len({edges.species for edges in edge_sets}) != len(edge_sets):
        raise ValueError("training species must be unique")

    positions: list[np.ndarray] = []
    values: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    for edges in edge_sets:
        if edges.n_edges < 1:
            raise ValueError("each training species must have at least one edge")
        positions.append(edges.midpoint)
        values.append(edges.turnover)
        weights.append(np.full(edges.n_edges, 1.0 / edges.n_edges, dtype=float))
    return KernelBoundaryModel(
        positions=np.vstack(positions),
        values=np.concatenate(values),
        weights=np.concatenate(weights),
        bandwidth=float(bandwidth),
        prior_strength=float(prior_strength),
        prior_mean=float(prior_mean),
    )


def transfer_statistic(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    *,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
) -> TransferResult:
    """Learn on train species and score edge-level predictions in unseen species."""
    if len(eval_edges) == 0:
        raise ValueError("at least one evaluation species is required")
    train_species = {edges.species for edges in train_edges}
    eval_species = {edges.species for edges in eval_edges}
    if train_species & eval_species:
        raise ValueError("train and evaluation species must be disjoint")

    model = fit_boundary_model(
        train_edges,
        bandwidth=bandwidth,
        prior_strength=prior_strength,
        prior_mean=prior_mean,
    )
    scores: dict[str, float] = {}
    for edges in eval_edges:
        predicted = model.predict_segments(
            edges.start,
            edges.end,
            n_points=segment_points,
        )
        scores[edges.species] = spearman_rho(predicted, edges.turnover)
    finite = np.asarray([score for score in scores.values() if np.isfinite(score)], dtype=float)
    if len(finite) == 0:
        raise ValueError("no evaluation species had enough edges for scoring")
    return TransferResult(
        statistic=float(finite.mean()),
        species_scores=scores,
        n_eval_species=int(len(finite)),
    )


@dataclass(frozen=True)
class PreparedTransfer:
    """Geometry-only projection reused across permutation replicates.

    This makes the null honest and cheap: geometry and kernel weights are frozen,
    but turnover ranks, train fit values, and held-out evaluation targets change
    on every permutation.
    """

    train_species: tuple[str, ...]
    train_slices: Mapping[str, slice]
    eval_species: tuple[str, ...]
    eval_projection: Mapping[str, np.ndarray]
    eval_opportunity: Mapping[str, np.ndarray]
    eval_prior_offset: Mapping[str, np.ndarray]
    prior_strength: float
    prior_mean: float

    def score(
        self,
        train_turnover: Mapping[str, np.ndarray],
        eval_turnover: Mapping[str, np.ndarray],
    ) -> TransferResult:
        train_values = np.empty(
            max(s.stop for s in self.train_slices.values()), dtype=float
        )
        for species in self.train_species:
            values = np.asarray(train_turnover[species], dtype=float)
            sl = self.train_slices[species]
            if values.shape != (sl.stop - sl.start,):
                raise ValueError(f"train turnover shape drift for {species}")
            train_values[sl] = values

        scores: dict[str, float] = {}
        for species in self.eval_species:
            projection = self.eval_projection[species]
            predicted = (
                projection @ train_values
                + self.eval_prior_offset[species]
            )
            target = np.asarray(eval_turnover[species], dtype=float)
            if target.shape != predicted.shape:
                raise ValueError(f"evaluation turnover shape drift for {species}")
            scores[species] = spearman_rho(predicted, target)

        finite = np.asarray([v for v in scores.values() if np.isfinite(v)], dtype=float)
        if len(finite) == 0:
            raise ValueError("no finite held-out species scores")
        return TransferResult(
            statistic=float(finite.mean()),
            species_scores=scores,
            n_eval_species=int(len(finite)),
        )


def prepare_transfer(
    train_edges: Sequence[SpeciesEdges],
    eval_edges: Sequence[SpeciesEdges],
    *,
    bandwidth: float,
    prior_strength: float = 0.25,
    prior_mean: float = 0.5,
    segment_points: int = 5,
) -> PreparedTransfer:
    """Precompute geometry-only kernels for fast full-pipeline permutation tests."""
    if bandwidth <= 0 or segment_points < 1:
        raise ValueError("bandwidth and segment_points must be positive")
    if len(train_edges) == 0 or len(eval_edges) == 0:
        raise ValueError("non-empty train and evaluation edge sets are required")
    train_names = tuple(edges.species for edges in train_edges)
    eval_names = tuple(edges.species for edges in eval_edges)
    if len(set(train_names)) != len(train_names) or len(set(eval_names)) != len(eval_names):
        raise ValueError("species labels must be unique within each split")
    if set(train_names) & set(eval_names):
        raise ValueError("train and evaluation species must be disjoint")

    train_slices: dict[str, slice] = {}
    positions: list[np.ndarray] = []
    base_weights: list[np.ndarray] = []
    cursor = 0
    for edges in train_edges:
        train_slices[edges.species] = slice(cursor, cursor + edges.n_edges)
        cursor += edges.n_edges
        positions.append(edges.midpoint)
        base_weights.append(np.full(edges.n_edges, 1.0 / edges.n_edges, dtype=float))
    train_positions = np.vstack(positions)
    train_weights = np.concatenate(base_weights)

    eval_projection: dict[str, np.ndarray] = {}
    eval_opportunity: dict[str, np.ndarray] = {}
    eval_prior_offset: dict[str, np.ndarray] = {}
    t = (np.arange(segment_points, dtype=float) + 0.5) / segment_points
    for edges in eval_edges:
        points = (
            edges.start[:, None, :]
            + t[None, :, None] * (edges.end - edges.start)[:, None, :]
        )
        flat = points.reshape(-1, train_positions.shape[1])
        delta = flat[:, None, :] - train_positions[None, :, :]
        distance2 = np.sum(delta * delta, axis=2)
        kernel = (
            np.exp(-0.5 * distance2 / (float(bandwidth) * float(bandwidth)))
            * train_weights[None, :]
        )
        point_opportunity = kernel.sum(axis=1)
        denominator = np.maximum(
            point_opportunity + float(prior_strength),
            np.finfo(float).tiny,
        )
        normalized_kernel = kernel / denominator[:, None]
        normalized_kernel = normalized_kernel.reshape(
            edges.n_edges, segment_points, -1
        ).mean(axis=1)
        prior_offset = (
            float(prior_strength)
            * float(prior_mean)
            / denominator
        ).reshape(edges.n_edges, segment_points).mean(axis=1)
        eval_projection[edges.species] = normalized_kernel
        eval_opportunity[edges.species] = point_opportunity.reshape(
            edges.n_edges, segment_points
        ).mean(axis=1)
        eval_prior_offset[edges.species] = prior_offset

    return PreparedTransfer(
        train_species=train_names,
        train_slices=train_slices,
        eval_species=eval_names,
        eval_projection=eval_projection,
        eval_opportunity=eval_opportunity,
        eval_prior_offset=eval_prior_offset,
        prior_strength=float(prior_strength),
        prior_mean=float(prior_mean),
    )
