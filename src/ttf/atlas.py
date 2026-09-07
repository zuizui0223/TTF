from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RecurrenceEstimate:
    probability: np.ndarray
    valid_realizations: np.ndarray
    top_fraction: float


def recurrence_probability(
    fields: np.ndarray,
    *,
    top_fraction: float = 0.10,
) -> RecurrenceEstimate:
    """Descriptive recurrence probability across balanced realizations.

    This is intentionally not the TTF sharedness test. It answers how often a
    location belongs to the high-turnover band under repeated resampling.
    """
    x = np.asarray(fields, dtype=float)
    if x.ndim != 2 or x.shape[0] < 1 or x.shape[1] < 1:
        raise ValueError("fields must be realization x cell")
    if not 0.0 < top_fraction < 1.0:
        raise ValueError("top_fraction must lie in (0, 1)")

    hits = np.zeros(x.shape[1], dtype=float)
    valid = np.zeros(x.shape[1], dtype=np.int64)
    for row in x:
        keep = np.isfinite(row)
        if not np.any(keep):
            continue
        threshold = float(np.quantile(row[keep], 1.0 - top_fraction))
        hits[keep] += (row[keep] >= threshold).astype(float)
        valid[keep] += 1
    probability = np.full(x.shape[1], np.nan, dtype=float)
    mask = valid > 0
    probability[mask] = hits[mask] / valid[mask]
    return RecurrenceEstimate(
        probability=probability,
        valid_realizations=valid,
        top_fraction=float(top_fraction),
    )
