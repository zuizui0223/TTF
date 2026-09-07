from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .core import SpeciesSample


@dataclass(frozen=True)
class SyntheticWorld:
    samples: tuple[SpeciesSample, ...]
    shared_species: tuple[str, ...]
    boundary_phase: Mapping[str, float]
    shared_fraction: float
    amplitude: float


def simulate_circular_boundary_world(
    *,
    n_species: int = 40,
    records_per_species: int = 60,
    shared_fraction: float,
    amplitude: float,
    noise_sd: float = 0.8,
    transition_width: float = 0.2,
    seed: int = 0,
) -> SyntheticWorld:
    """Generate an orthogonal sharedness x amplitude calibration world.

    Each species is sampled on the same circular geographic domain. Every
    species may have a strong spatial transition. ``shared_fraction`` controls
    only whether the transition *location* is common across species; amplitude
    controls within-species strength.

    At shared_fraction=0 and amplitude>0, species have strong but randomly
    phased boundaries. Because phase is uniform on the circle, the population
    has no privileged geographic boundary location. This is the mandatory
    adversarial type-I arm.
    """
    if n_species < 4:
        raise ValueError("n_species must be >= 4")
    if records_per_species < 8:
        raise ValueError("records_per_species must be >= 8")
    if not 0.0 <= shared_fraction <= 1.0:
        raise ValueError("shared_fraction must lie in [0, 1]")
    if amplitude < 0 or noise_sd < 0 or transition_width <= 0:
        raise ValueError("amplitude/noise must be non-negative and width positive")

    rng = np.random.default_rng(int(seed))
    names = np.asarray([f"sp_{i:03d}" for i in range(n_species)], dtype=object)
    n_shared = int(round(n_species * float(shared_fraction)))
    shared_index = set(map(int, rng.permutation(n_species)[:n_shared]))
    private_phase = rng.uniform(0.0, np.pi, n_species)

    samples: list[SpeciesSample] = []
    phases: dict[str, float] = {}
    shared_names: list[str] = []
    for i, name_obj in enumerate(names):
        name = str(name_obj)
        phase = 0.0 if i in shared_index else float(private_phase[i])
        theta = rng.uniform(0.0, 2.0 * np.pi, records_per_species)
        coordinates = np.column_stack((np.cos(theta), np.sin(theta)))
        latent = float(amplitude) * np.tanh(
            np.sin(theta - phase) / float(transition_width)
        )
        trait = latent + rng.normal(0.0, float(noise_sd), records_per_species)
        samples.append(
            SpeciesSample(
                species=name,
                coordinates=coordinates,
                trait=trait,
            )
        )
        phases[name] = phase
        if i in shared_index:
            shared_names.append(name)

    return SyntheticWorld(
        samples=tuple(samples),
        shared_species=tuple(sorted(shared_names)),
        boundary_phase=phases,
        shared_fraction=float(shared_fraction),
        amplitude=float(amplitude),
    )
