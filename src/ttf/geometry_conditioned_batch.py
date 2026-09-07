from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .geometry_control import (
    length_orthogonalized_turnover,
    partial_spearman_rho,
    spearman_length_association,
)
from .transfer import PreparedTransfer


@dataclass(frozen=True)
class GeometryConditionedBatchResult:
    """Exact transfer scores after prospectively defined edge-length control."""

    statistics: np.ndarray
    species_scores: Mapping[str, np.ndarray]
    n_eval_species: np.ndarray
    eval_turnover_length_rho: Mapping[str, np.ndarray]
    eval_prediction_length_rho: Mapping[str, np.ndarray]


def _width(values: Mapping[str, np.ndarray], names: tuple[str, ...]) -> int:
    widths: set[int] = set()
    for name in names:
        array = np.asarray(values[name], dtype=float)
        if array.ndim != 2:
            raise ValueError(f"turnover for {name} must be edges x worlds")
        widths.add(int(array.shape[1]))
    if len(widths) != 1:
        raise ValueError("turnover batch widths disagree")
    width = widths.pop()
    if width < 1:
        raise ValueError("batch must contain at least one world")
    return width


def score_prepared_length_conditioned_batch(
    prepared: PreparedTransfer,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    *,
    train_edge_length: Mapping[str, np.ndarray],
    eval_edge_length: Mapping[str, np.ndarray],
) -> GeometryConditionedBatchResult:
    """Score worlds after removing within-species monotone edge-length nuisance.

    The training field is fit to length-orthogonalized turnover. Held-out species
    are scored by partial Spearman correlation between field exposure and raw
    turnover, controlling the same species' edge-length rank. Geometry and the
    prepared kernel projection remain fixed exactly; only the response-side
    estimand changes, so this is a new v0.3 development statistic rather than a
    reinterpretation of v0.2 Gate-I failures.

    ``prepared`` should be built with ``prior_mean=0`` because orthogonalized
    training turnover is centered at zero.
    """
    width = _width(train_turnover, prepared.train_species)
    if _width(eval_turnover, prepared.eval_species) != width:
        raise ValueError("train and evaluation batch widths disagree")

    n_train_edges = max(s.stop for s in prepared.train_slices.values())
    train_values = np.empty((n_train_edges, width), dtype=float)
    for species in prepared.train_species:
        values = np.asarray(train_turnover[species], dtype=float)
        length = np.asarray(train_edge_length[species], dtype=float)
        sl = prepared.train_slices[species]
        expected = (sl.stop - sl.start, width)
        if values.shape != expected or length.shape != (expected[0],):
            raise ValueError(f"training geometry/turnover shape drift for {species}")
        for column in range(width):
            train_values[sl, column] = length_orthogonalized_turnover(
                values[:, column],
                length,
            )

    score_map: dict[str, np.ndarray] = {}
    target_length_map: dict[str, np.ndarray] = {}
    prediction_length_map: dict[str, np.ndarray] = {}
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)

    for species in prepared.eval_species:
        projection = prepared.eval_projection[species]
        predicted = projection @ train_values + prepared.eval_prior_offset[species][:, None]
        target = np.asarray(eval_turnover[species], dtype=float)
        length = np.asarray(eval_edge_length[species], dtype=float)
        if target.shape != predicted.shape or length.shape != (predicted.shape[0],):
            raise ValueError(f"evaluation geometry/turnover shape drift for {species}")

        scores = np.empty(width, dtype=float)
        target_length = np.empty(width, dtype=float)
        prediction_length = np.empty(width, dtype=float)
        for column in range(width):
            scores[column] = partial_spearman_rho(
                predicted[:, column],
                target[:, column],
                length,
            )
            target_length[column] = spearman_length_association(
                target[:, column],
                length,
            )
            prediction_length[column] = spearman_length_association(
                predicted[:, column],
                length,
            )

        score_map[species] = scores
        target_length_map[species] = target_length
        prediction_length_map[species] = prediction_length
        finite = np.isfinite(scores)
        sums[finite] += scores[finite]
        counts[finite] += 1

    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite held-out species scores")

    return GeometryConditionedBatchResult(
        statistics=sums / counts,
        species_scores=score_map,
        n_eval_species=counts,
        eval_turnover_length_rho=target_length_map,
        eval_prediction_length_rho=prediction_length_map,
    )
