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
    "strict_private_periodic_directional_change",
    "shared_front_zone_directional_change",
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
    strict_private_arc_halfwidth: float = 1.0,
    seed: int = 0,
) -> DirectionalSyntheticWorld:
    """Generate Q4 directionality worlds on prospectively oriented supports.

    Existing Q4 v0.1 modes are preserved exactly. ``private_directional_change``
    remains the frozen bounded-center construction used by the immutable failed
    Q4 v0.1 gate.

    Q4.1 adds two explicitly different worlds:

    ``strict_private_periodic_directional_change``
        Each system receives a phase anchor uniformly on the full unit circle.
        Absolute TTF coordinates are a local circular arc centered on that
        independently sampled phase, while the separate outcome-independent
        local orientation coordinate spans [-1, 1] and places the low-to-high
        mismatch transition at orientation zero. The common directional sign is
        retained but the absolute population turnover-intensity field is
        rotationally invariant.

    ``shared_front_zone_directional_change``
        Retains the old finite-line bounded-center construction but labels it as
        a shared front-density zone alternative rather than a strict-private
        type-I null.
    """

    allowed = {
        "shared_breakdown",
        "shared_recoupling",
        "shared_neutral_rotation",
        "private_directional_change",
        "shared_sign_flip",
        "strict_private_periodic_directional_change",
        "shared_front_zone_directional_change",
    }
    if mode not in allowed:
        raise ValueError(f"unknown directional world mode: {mode}")
    if n_species < 8 or records_per_species < 20:
        raise ValueError("directional worlds require >=8 systems and >=20 records/system")
    if not (0 < low_mismatch < high_mismatch):
        raise ValueError("require 0 < low_mismatch < high_mismatch")
    if transition_width <= 0 or noise_sd < 0:
        raise ValueError("transition_width positive and noise_sd non-negative required")
    if not (0.0 < float(strict_private_arc_halfwidth) < np.pi):
        raise ValueError("strict_private_arc_halfwidth must lie in (0, pi)")

    rng = np.random.default_rng(int(seed))
    samples: list[OrientedPairedSpeciesSample] = []
    fronts: dict[str, float] = {}

    # Preserve the RNG draw used by every Q4 v0.1 mode. This is deliberately
    # unconditional so adding Q4.1 modes cannot alter frozen Q4 v0.1 worlds.
    private_fronts = rng.uniform(-0.75, 0.75, size=n_species)
    strict_phases = (
        rng.uniform(0.0, 2.0 * np.pi, size=n_species)
        if mode == "strict_private_periodic_directional_change"
        else None
    )

    for i in range(n_species):
        name = f"sp_{i:03d}"

        if mode == "strict_private_periodic_directional_change":
            # The absolute spatial front is private and uniformly rotated. The
            # orientation coordinate is a separate predeclared local axis.
            g = np.sort(rng.uniform(-1.0, 1.0, size=records_per_species))
            phase = float(strict_phases[i])  # type: ignore[index]
            angle = phase + float(strict_private_arc_halfwidth) * g
            coordinates = np.column_stack((np.cos(angle), np.sin(angle)))
            fronts[name] = phase
            h = _transition(g, 0.0, transition_width)
            latent = float(low_mismatch) + (
                float(high_mismatch) - float(low_mismatch)
            ) * h
            state_a = latent + rng.normal(
                0.0, float(noise_sd), size=records_per_species
            )
            state_b = np.zeros(records_per_species, dtype=float)
            orientation = g

        else:
            # Keep the original Q4 v0.1 construction unchanged for all old
            # modes. shared_front_zone is intentionally the same bounded-center
            # geometry as the old PRIVATE cell, but with corrected semantics.
            x = np.sort(rng.uniform(-1.0, 1.0, size=records_per_species))
            coordinates = x[:, None]
            center = (
                float(private_fronts[i])
                if mode in {
                    "private_directional_change",
                    "shared_front_zone_directional_change",
                }
                else 0.0
            )
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
                elif mode in {
                    "private_directional_change",
                    "shared_front_zone_directional_change",
                }:
                    left, right = float(low_mismatch), float(high_mismatch)
                elif mode == "shared_sign_flip":
                    if i % 2 == 0:
                        left, right = float(low_mismatch), float(high_mismatch)
                    else:
                        left, right = float(high_mismatch), float(low_mismatch)
                else:  # pragma: no cover
                    raise AssertionError(mode)
                latent = left + (right - left) * h
                state_a = latent + rng.normal(
                    0.0, float(noise_sd), size=records_per_species
                )
                state_b = np.zeros(records_per_species, dtype=float)
            orientation = x

        paired = PairedSpeciesSample(
            species=name,
            coordinates=coordinates,
            state_a=state_a,
            state_b=state_b,
        )
        samples.append(
            OrientedPairedSpeciesSample(sample=paired, orientation=orientation)
        )

    return DirectionalSyntheticWorld(
        samples=tuple(samples),
        mode=str(mode),
        front_by_species=fronts,
        transition_width=float(transition_width),
        low_mismatch=float(low_mismatch),
        high_mismatch=float(high_mismatch),
    )


__all__ = ["DirectionalSyntheticWorld", "simulate_directional_world"]
