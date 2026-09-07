from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable
import math

import numpy as np

from .calibration import CalibrationCell


@dataclass(frozen=True)
class BinomialInterval:
    successes: int
    trials: int
    estimate: float
    low: float
    high: float

    def to_dict(self) -> dict[str, float | int]:
        return asdict(self)


@dataclass(frozen=True)
class PrecisionQualificationReport:
    type1_pass: bool
    power_pass: bool
    type1_upper_ceiling: float
    power_lower_floor: float
    max_zero_shared_estimate: float
    max_zero_shared_upper95: float
    full_shared_moderate_estimate: float
    full_shared_moderate_lower95: float
    moderate_amplitude: float
    zero_shared_intervals: tuple[dict[str, float | int], ...]
    full_shared_interval: dict[str, float | int]

    @property
    def passed(self) -> bool:
        return bool(self.type1_pass and self.power_pass)

    def to_dict(self) -> dict:
        out = asdict(self)
        out["passed"] = self.passed
        return out


def wilson_interval(
    successes: int,
    trials: int,
    *,
    z: float = 1.959963984540054,
) -> BinomialInterval:
    """Two-sided Wilson score interval for a binomial proportion."""
    successes = int(successes)
    trials = int(trials)
    if trials < 1 or successes < 0 or successes > trials:
        raise ValueError("invalid binomial counts")
    p = successes / trials
    z2 = float(z) ** 2
    denom = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / denom
    half = float(z) * math.sqrt(
        (p * (1.0 - p) + z2 / (4.0 * trials)) / trials
    ) / denom
    return BinomialInterval(
        successes=successes,
        trials=trials,
        estimate=float(p),
        low=float(max(0.0, center - half)),
        high=float(min(1.0, center + half)),
    )


def _cell_interval(cell: CalibrationCell) -> BinomialInterval:
    raw = float(cell.rejection_rate) * int(cell.n_replicates)
    count = int(round(raw))
    if not np.isclose(raw, count, atol=1e-8, rtol=0.0):
        raise ValueError("calibration rejection rate is not compatible with replicate count")
    return wilson_interval(count, int(cell.n_replicates))


def qualify_calibration_precision(
    cells: Iterable[CalibrationCell],
    *,
    moderate_amplitude: float,
    type1_upper_ceiling: float = 0.10,
    power_lower_floor: float = 0.80,
) -> PrecisionQualificationReport:
    """Qualification that requires confidence bounds, not point estimates alone.

    Every positive-amplitude zero-shared arm must have a Wilson 95% upper bound
    no greater than ``type1_upper_ceiling``.  The fully shared moderate-signal
    arm must have a Wilson 95% lower bound at least ``power_lower_floor``.
    """
    data = tuple(cells)
    zero_shared = sorted(
        (
            cell
            for cell in data
            if np.isclose(cell.shared_fraction, 0.0) and cell.amplitude > 0
        ),
        key=lambda cell: cell.amplitude,
    )
    positive = [
        cell
        for cell in data
        if np.isclose(cell.shared_fraction, 1.0)
        and np.isclose(cell.amplitude, float(moderate_amplitude))
    ]
    if not zero_shared:
        raise ValueError("precision calibration lacks zero-shared positive-amplitude arms")
    if len(positive) != 1:
        raise ValueError("precision calibration requires exactly one full-shared moderate arm")

    zero_rows: list[dict[str, float | int]] = []
    zero_intervals: list[BinomialInterval] = []
    for cell in zero_shared:
        interval = _cell_interval(cell)
        zero_intervals.append(interval)
        zero_rows.append(
            {
                "shared_fraction": float(cell.shared_fraction),
                "amplitude": float(cell.amplitude),
                **interval.to_dict(),
            }
        )

    positive_interval = _cell_interval(positive[0])
    positive_row = {
        "shared_fraction": float(positive[0].shared_fraction),
        "amplitude": float(positive[0].amplitude),
        **positive_interval.to_dict(),
    }
    max_upper = max(interval.high for interval in zero_intervals)
    max_estimate = max(interval.estimate for interval in zero_intervals)

    return PrecisionQualificationReport(
        type1_pass=bool(max_upper <= float(type1_upper_ceiling)),
        power_pass=bool(positive_interval.low >= float(power_lower_floor)),
        type1_upper_ceiling=float(type1_upper_ceiling),
        power_lower_floor=float(power_lower_floor),
        max_zero_shared_estimate=float(max_estimate),
        max_zero_shared_upper95=float(max_upper),
        full_shared_moderate_estimate=float(positive_interval.estimate),
        full_shared_moderate_lower95=float(positive_interval.low),
        moderate_amplitude=float(moderate_amplitude),
        zero_shared_intervals=tuple(zero_rows),
        full_shared_interval=positive_row,
    )
