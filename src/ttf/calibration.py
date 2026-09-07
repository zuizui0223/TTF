from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from typing import Iterable, Sequence

import numpy as np

from .core import split_species
from .inference import heldout_species_bootstrap_test
from .nulls import permutation_test
from .simulate import simulate_circular_boundary_world


def seed_for(master_seed: int, *parts: object) -> int:
    payload = "|".join(map(str, (int(master_seed),) + parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


@dataclass(frozen=True)
class CalibrationCell:
    shared_fraction: float
    amplitude: float
    n_replicates: int
    alpha: float
    rejection_rate: float
    mean_statistic: float
    mean_null_statistic: float
    median_p_value: float

    def to_dict(self) -> dict[str, float | int]:
        return asdict(self)


@dataclass(frozen=True)
class QualificationReport:
    type1_pass: bool
    power_pass: bool
    max_zero_shared_rejection: float
    full_shared_moderate_power: float
    type1_ceiling: float
    power_floor: float
    moderate_amplitude: float

    @property
    def passed(self) -> bool:
        return bool(self.type1_pass and self.power_pass)

    def to_dict(self) -> dict[str, float | bool]:
        out = asdict(self)
        out["passed"] = self.passed
        return out


def run_calibration(
    *,
    shared_fractions: Sequence[float],
    amplitudes: Sequence[float],
    n_replicates: int,
    n_permutations: int,
    n_species: int = 40,
    records_per_species: int = 60,
    eval_fraction: float = 0.5,
    k: int = 4,
    bandwidth: float = 0.2,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    noise_sd: float = 0.8,
    transition_width: float = 0.2,
    alpha: float = 0.05,
    seed: int = 20260907,
    inference: str = "trait_permutation",
) -> tuple[CalibrationCell, ...]:
    """Estimate type-I and power over the sharedness x amplitude plane.

    ``trait_permutation`` is retained as a diagnostic test of trait-location
    exchangeability. ``heldout_species_bootstrap`` is the sharedness-specific
    conditional test: it leaves every species' spatial structure untouched and
    resamples held-out species scores under a centered null mean.
    """
    if n_replicates < 1 or n_permutations < 1:
        raise ValueError("n_replicates and n_permutations must be >= 1")
    if inference not in {"trait_permutation", "heldout_species_bootstrap"}:
        raise ValueError("unknown inference method")
    if inference == "heldout_species_bootstrap" and n_permutations < 99:
        raise ValueError("heldout species bootstrap requires at least 99 resamples")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1)")

    cells: list[CalibrationCell] = []
    for shared in map(float, shared_fractions):
        for amplitude in map(float, amplitudes):
            p_values: list[float] = []
            observed: list[float] = []
            null_means: list[float] = []
            for replicate in range(n_replicates):
                world_seed = seed_for(seed, "world", shared, amplitude, replicate)
                split_seed = seed_for(seed, "split", shared, amplitude, replicate)
                null_seed = seed_for(seed, inference, shared, amplitude, replicate)
                world = simulate_circular_boundary_world(
                    n_species=n_species,
                    records_per_species=records_per_species,
                    shared_fraction=shared,
                    amplitude=amplitude,
                    noise_sd=noise_sd,
                    transition_width=transition_width,
                    seed=world_seed,
                )
                train, evaluation = split_species(
                    [sample.species for sample in world.samples],
                    eval_fraction=eval_fraction,
                    seed=split_seed,
                )
                common = dict(
                    samples=world.samples,
                    train_species=train,
                    eval_species=evaluation,
                    k=k,
                    bandwidth=bandwidth,
                    prior_strength=prior_strength,
                    segment_points=segment_points,
                    seed=null_seed,
                )
                if inference == "trait_permutation":
                    result = permutation_test(
                        **common,
                        n_permutations=n_permutations,
                    )
                else:
                    result = heldout_species_bootstrap_test(
                        **common,
                        n_bootstrap=n_permutations,
                    )
                p_values.append(result.p_value)
                observed.append(result.observed.statistic)
                null_means.append(result.null_mean)

            p = np.asarray(p_values, dtype=float)
            cells.append(
                CalibrationCell(
                    shared_fraction=shared,
                    amplitude=amplitude,
                    n_replicates=n_replicates,
                    alpha=float(alpha),
                    rejection_rate=float(np.mean(p <= alpha)),
                    mean_statistic=float(np.mean(observed)),
                    mean_null_statistic=float(np.mean(null_means)),
                    median_p_value=float(np.median(p)),
                )
            )
    return tuple(cells)


def qualify_calibration(
    cells: Iterable[CalibrationCell],
    *,
    moderate_amplitude: float,
    type1_ceiling: float = 0.10,
    power_floor: float = 0.80,
) -> QualificationReport:
    """Prospective gates required before an empirical TTF claim is interpretable."""
    data = tuple(cells)
    zero_shared = [
        cell
        for cell in data
        if np.isclose(cell.shared_fraction, 0.0) and cell.amplitude > 0
    ]
    full_shared = [
        cell
        for cell in data
        if np.isclose(cell.shared_fraction, 1.0)
        and np.isclose(cell.amplitude, float(moderate_amplitude))
    ]
    if not zero_shared:
        raise ValueError("calibration lacks mandatory shared=0, amplitude>0 cells")
    if len(full_shared) != 1:
        raise ValueError("calibration must contain one full-shared moderate-amplitude cell")

    max_type1 = max(cell.rejection_rate for cell in zero_shared)
    power = full_shared[0].rejection_rate
    return QualificationReport(
        type1_pass=bool(max_type1 <= float(type1_ceiling)),
        power_pass=bool(power >= float(power_floor)),
        max_zero_shared_rejection=float(max_type1),
        full_shared_moderate_power=float(power),
        type1_ceiling=float(type1_ceiling),
        power_floor=float(power_floor),
        moderate_amplitude=float(moderate_amplitude),
    )
