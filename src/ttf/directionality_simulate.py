from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping

import numpy as np

from .directionality import OrientedPairedSpeciesSample
from .mismatch import PairedSpeciesSample

DirectionalWorldMode = Literal[
    "shared_breakdown",
    "shared_recoupling",
    "shared_neutral_rotation",
    "private_directional_change",
    "shared_sign_flip",
]


@dataclass(frozen=True)
class DirectionalSyntheticWorld:
    samples: tuple[OrientedPairedSpeciesSample, ...]
    mode: str
    front_by_species: Mapping[str, float]
    transition_width: float
    low_mismatch: float
    high_mismatch: float


def _transition(x: np.ndarray, center: float, width: float) -> np.ndarray:
    return 0.5 * (1.0 + np.tanh((x - float(center)) / float(width)))


def _relation_rotation(h: np.ndarray, amplitude: float) -> np.ndarray:
    side = h >= 0.5
    state = np.zeros((len(h), 2), dtype=float)
    state[~side, 0] = float(amplitude)
    state[side, 1] = float(amplitude)
    return state


def simulate_directional_world(
    *,
    mode: DirectionalWorldMode,
    n_species: int = 40,
    records_per_species: int = 60,
    low_mismatch: float = 1.0,
    high_mismatch: float = 3.0,
    noise_sd: float = 0.25,
    transition_width: float = 0.12,
    seed: int = 0,
) -> DirectionalSyntheticWorld:
    """Generate Q4 directionality worlds on a prospectively oriented axis.

    The signed orientation coordinate is x itself and is generated before state
    outcomes. Shared worlds place the relation front at x=0. Private worlds use
    independent front centers. The sign-flip world alternates deterioration and
    recoupling across the same front, preserving a transferable direction-free
    relation front while cancelling shared directional mismatch change.
    """

    allowed = {
        "shared_breakdown",
        "shared_recoupling",
        "shared_neutral_rotation",
        "private_directional_change",
        "shared_sign_flip",
    }
    if mode not in allowed:
        raise ValueError(f"unknown directional world mode: {mode}")
    if n_species < 8 or records_per_species < 20:
        raise ValueError("directional worlds require >=8 systems and >=20 records/system")
    if not (0 < low_mismatch < high_mismatch):
        raise ValueError("require 0 < low_mismatch < high_mismatch")
    if transition_width <= 0 or noise_sd < 0:
        raise ValueError("transition_width positive and noise_sd non-negative required")

    rng = np.random.default_rng(int(seed))
    samples: list[OrientedPairedSpeciesSample] = []
    fronts: dict[str, float] = {}
    private_fronts = rng.uniform(-0.75, 0.75, size=n_species)

    for i in range(n_species):
        name = f"sp_{i:03d}"
        x = np.sort(rng.uniform(-1.0, 1.0, size=records_per_species))
        coordinates = x[:, None]
        center = float(private_fronts[i]) if mode == "private_directional_change" else 0.0
        fronts[name] = center
        h = _transition(x, center, transition_width)

        if mode == "shared_neutral_rotation":
            state_a = _relation_rotation(h, amplitude=float(high_mismatch))
            state_b = np.zeros_like(state_a)
        else:
            if mode == "shared_breakdown":
                left, right = float(low_mismatch), float(high_mismatch)
            elif mode == "shared_recoupling":
                left, right = float(high_mismatch), float(low_mismatch)
            elif mode == "private_directional_change":
                left, right = float(low_mismatch), float(high_mismatch)
            elif mode == "shared_sign_flip":
                if i % 2 == 0:
                    left, right = float(low_mismatch), float(high_mismatch)
                else:
                    left, right = float(high_mismatch), float(low_mismatch)
            else:  # pragma: no cover
                raise AssertionError(mode)
            latent = left + (right - left) * h
            state_a = latent + rng.normal(0.0, float(noise_sd), size=records_per_species)
            state_b = np.zeros(records_per_species, dtype=float)

        paired = PairedSpeciesSample(
            species=name,
            coordinates=coordinates,
            state_a=state_a,
            state_b=state_b,
        )
        samples.append(OrientedPairedSpeciesSample(sample=paired, orientation=x))

    return DirectionalSyntheticWorld(
        samples=tuple(samples),
        mode=str(mode),
        front_by_species=fronts,
        transition_width=float(transition_width),
        low_mismatch=float(low_mismatch),
        high_mismatch=float(high_mismatch),
    )


__all__ = ["DirectionalSyntheticWorld", "simulate_directional_world"]
