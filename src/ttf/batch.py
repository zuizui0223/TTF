from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .core import spearman_rho
from .transfer import PreparedTransfer


@dataclass(frozen=True)
class BatchTransferResult:
    """Exact prepared-transfer scores for several response worlds at once."""

    statistics: np.ndarray
    species_scores: Mapping[str, np.ndarray]
    n_eval_species: np.ndarray


def _batch_width(
    values: Mapping[str, np.ndarray],
    species: tuple[str, ...],
    *,
    label: str,
) -> int:
    widths: set[int] = set()
    for name in species:
        if name not in values:
            raise ValueError(f"missing {label} turnover for {name}")
        array = np.asarray(values[name], dtype=float)
        if array.ndim != 2:
            raise ValueError(f"{label} turnover for {name} must be edges x worlds")
        widths.add(int(array.shape[1]))
    if len(widths) != 1:
        raise ValueError(f"{label} turnover batch widths disagree")
    width = widths.pop()
    if width < 1:
        raise ValueError("batch must contain at least one world")
    return width


def score_prepared_batch(
    prepared: PreparedTransfer,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
) -> BatchTransferResult:
    """Score many trait worlds while reusing one exact dense geometry projection.

    This changes only the linear-algebra execution shape.  Every column is the
    same calculation as ``PreparedTransfer.score`` for that response world;
    opportunity correction and the edge-integrated kernel projection are not
    approximated or truncated.
    """

    train_width = _batch_width(
        train_turnover,
        prepared.train_species,
        label="train",
    )
    eval_width = _batch_width(
        eval_turnover,
        prepared.eval_species,
        label="evaluation",
    )
    if train_width != eval_width:
        raise ValueError("train and evaluation batch widths disagree")
    width = train_width

    n_train_edges = max(s.stop for s in prepared.train_slices.values())
    train_values = np.empty((n_train_edges, width), dtype=float)
    for species in prepared.train_species:
        values = np.asarray(train_turnover[species], dtype=float)
        sl = prepared.train_slices[species]
        expected = (sl.stop - sl.start, width)
        if values.shape != expected:
            raise ValueError(
                f"train turnover shape drift for {species}: {values.shape} != {expected}"
            )
        train_values[sl, :] = values

    score_map: dict[str, np.ndarray] = {}
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    for species in prepared.eval_species:
        projection = prepared.eval_projection[species]
        predicted = (
            projection @ train_values
            + prepared.eval_prior_offset[species][:, None]
        )
        target = np.asarray(eval_turnover[species], dtype=float)
        expected = (predicted.shape[0], width)
        if target.shape != expected:
            raise ValueError(
                f"evaluation turnover shape drift for {species}: {target.shape} != {expected}"
            )
        scores = np.asarray(
            [spearman_rho(predicted[:, i], target[:, i]) for i in range(width)],
            dtype=float,
        )
        score_map[species] = scores
        finite = np.isfinite(scores)
        sums[finite] += scores[finite]
        counts[finite] += 1

    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite held-out species scores")
    return BatchTransferResult(
        statistics=sums / counts,
        species_scores=score_map,
        n_eval_species=counts,
    )
