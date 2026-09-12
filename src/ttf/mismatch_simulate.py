from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping

import numpy as np

from .mismatch import PairedSpeciesSample

PairedWorldMode = Literal[
    "coupled_shared_transition",
    "coupled_shared_transition_independent_noise",
    "private_mismatch_transition",
    "shared_mismatch_transition",
    "private_relational_rotation",
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


def _exact_relation_rotation(
    h: np.ndarray,
    *,
    amplitude: float,
) -> np.ndarray:
    """Two orthogonal relation states with exactly equal Euclidean norm."""
    side = h >= 0.5
    state = np.zeros((len(h), 2), dtype=float)
    state[~side, 0] = float(amplitude)
    state[side, 1] = float(amplitude)
    return state


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
    """Generate structural and qualification worlds for TTF-M and TTF-C.

    Modes
    -----
    coupled_shared_transition
        A and B undergo exactly the same shared transition and share the same
        noise realization. Component turnover is strong and transferable, but
        mismatch and relative coupling states are exactly constant. This is the
        original structural false-positive control.

    coupled_shared_transition_independent_noise
        A and B share the same strong latent spatial transition but receive
        independent observation noise. There is still no mismatch/coupling
        boundary in expectation. This is the harder formal false-positive
        control because both projected states contain spatially unstructured
        measurement variation while the component fields themselves share a
        strong boundary.

    private_mismatch_transition
        A is a baseline and B departs across a strong system-private boundary.
        Every system has structured mismatch, but boundary phase is independent.

    shared_mismatch_transition
        A is a baseline and B departs across one common boundary. Mismatch
        magnitude and relative coupling both carry a shared transition.

    private_relational_rotation
        The relative vector A-B switches between two orthogonal directions at
        an independently phased boundary in every system. Euclidean mismatch
        magnitude is exactly constant but each system has strong private
        coupling change. This is a TTF-C-specific private-structure null.

    shared_relational_rotation
        The relative vector A-B switches between the same two orthogonal
        directions across one common boundary while keeping identical norm.
        With ``noise_sd=0`` this construction is numerically exact, so mismatch
        magnitude is constant and only the coupling relation changes.
    """

    allowed = {
        "coupled_shared_transition",
        "coupled_shared_transition_independent_noise",
        "private_mismatch_transition",
        "shared_mismatch_transition",
        "private_relational_rotation",
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

        if mode in {"private_mismatch_transition", "private_relational_rotation"}:
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

        elif mode == "coupled_shared_transition_independent_noise":
            latent = float(amplitude) * h
            state_a = latent + rng.normal(0.0, float(noise_sd), records_per_species)
            state_b = latent + rng.normal(0.0, float(noise_sd), records_per_species)

        elif mode in {"private_mismatch_transition", "shared_mismatch_transition"}:
            state_a = np.zeros(records_per_species, dtype=float)
            state_b = (
                float(amplitude) * h
                + rng.normal(0.0, float(noise_sd), records_per_species)
            )

        elif mode in {"private_relational_rotation", "shared_relational_rotation"}:
            state_a = _exact_relation_rotation(h, amplitude=float(amplitude))
            if noise_sd:
                # Rotation noise preserves the relation norm analytically but
                # can introduce floating-point norm variation. Formal
                # constant-magnitude controls therefore freeze noise_sd=0.
                angle = rng.normal(0.0, float(noise_sd), records_per_species)
                c = np.cos(angle)
                s = np.sin(angle)
                x = state_a[:, 0].copy()
                y = state_a[:, 1].copy()
                state_a[:, 0] = c * x - s * y
                state_a[:, 1] = s * x + c * y
            state_b = np.zeros_like(state_a)

        else:  # pragma: no cover - guarded by allowed set
            raise AssertionError(mode)

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
