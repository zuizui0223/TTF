from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping

import numpy as np

from .mismatch import PairedSpeciesSample

PairedWorldMode = Literal[
    "coupled_shared_transition",
    "private_mismatch_transition",
    "shared_mismatch_transition",
    "shared_relational_rotation",
]


@dataclass(frozen=True)
class PairedSyntheticWorld:
    samples: tuple[PairedSpeciesSample, ...]
    mode: str
    boundary_phase: Mapping[str, float]
    transition_width: float
    amplitude: float


def _sigmoid_transition(theta: np.ndarray, phase: float, width: float) -> np.ndarray:
    return 0.5 * (1.0 + np.tanh(np.sin(theta - float(phase)) / float(width)))


def simulate_paired_transition_world(
    *,
    mode: PairedWorldMode,
    n_species: int = 40,
    records_per_species: int = 60,
    amplitude: float = 2.0,
    noise_sd: float = 0.0,
    transition_width: float = 0.2,
    seed: int = 0,
) -> PairedSyntheticWorld:
    """Generate structural qualification worlds for TTF-M and TTF-C.

    Modes
    -----
    coupled_shared_transition
        A and B undergo exactly the same shared transition. Component turnover
        is strong and transferable, but mismatch and relative coupling states
        are constant. This is the mandatory false-positive control.

    private_mismatch_transition
        A is a baseline and B departs across a strong system-private boundary.
        Every system has structured mismatch, but boundary phase is independent.

    shared_mismatch_transition
        A is a baseline and B departs across one common boundary. Mismatch
        magnitude and relative coupling both carry a shared transition.

    shared_relational_rotation
        The relative vector A-B rotates from e1 to e2 across one common
        boundary while keeping unit norm. Mismatch magnitude is exactly
        constant when noise_sd=0, but the coupling relation changes.
    """

    allowed = {
        "coupled_shared_transition",
        "private_mismatch_transition",
        "shared_mismatch_transition",
        "shared_relational_rotation",
    }
    if mode not in allowed:
        raise ValueError(f"unknown paired world mode: {mode}")
    if n_species < 4:
        raise ValueError("n_species must be >= 4")
    if records_per_species < 8:
        raise ValueError("records_per_species must be >= 8")
    if amplitude < 0 or noise_sd < 0 or transition_width <= 0:
        raise ValueError("amplitude/noise must be non-negative and width positive")

    rng = np.random.default_rng(int(seed))
    samples: list[PairedSpeciesSample] = []
    phases: dict[str, float] = {}

    private_phases = rng.uniform(0.0, 2.0 * np.pi, n_species)
    for i in range(n_species):
        name = f"sp_{i:03d}"
        theta = rng.uniform(0.0, 2.0 * np.pi, records_per_species)
        coordinates = np.column_stack((np.cos(theta), np.sin(theta)))

        if mode == "private_mismatch_transition":
            phase = float(private_phases[i])
        else:
            phase = 0.0
        phases[name] = phase
        h = _sigmoid_transition(theta, phase, transition_width)

        if mode == "coupled_shared_transition":
            shared_noise = rng.normal(0.0, float(noise_sd), records_per_species)
            state = float(amplitude) * h + shared_noise
            state_a = state
            state_b = state.copy()

        elif mode in {"private_mismatch_transition", "shared_mismatch_transition"}:
            state_a = np.zeros(records_per_species, dtype=float)
            state_b = (
                float(amplitude) * h
                + rng.normal(0.0, float(noise_sd), records_per_species)
            )

        else:  # shared_relational_rotation
            # phi moves from 0 to pi/2 across the shared boundary. Angular
            # noise preserves ||A-B||=amplitude exactly.
            phi = 0.5 * np.pi * h
            if noise_sd:
                phi = phi + rng.normal(0.0, float(noise_sd), records_per_species)
            state_a = float(amplitude) * np.column_stack((np.cos(phi), np.sin(phi)))
            state_b = np.zeros_like(state_a)

        samples.append(
            PairedSpeciesSample(
                species=name,
                coordinates=coordinates,
                state_a=state_a,
                state_b=state_b,
            )
        )

    return PairedSyntheticWorld(
        samples=tuple(samples),
        mode=str(mode),
        boundary_phase=phases,
        transition_width=float(transition_width),
        amplitude=float(amplitude),
    )
