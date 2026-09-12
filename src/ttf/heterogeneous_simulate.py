from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping

import numpy as np

from .mismatch import PairedSpeciesSample

Q5ResponseMode = Literal[
    "null_component",
    "null_private_mismatch",
    "null_private_relation",
    "power_shared_mismatch",
    "power_shared_relation",
]
GeometryProfile = Literal["matched", "shifted"]


@dataclass(frozen=True)
class SystemGeometryRecord:
    species: str
    split: str
    nominal_n: int
    effective_n: int
    clustering: str
    truncation: str
    missingness: float


@dataclass(frozen=True)
class HeterogeneousPairedWorld:
    samples: tuple[PairedSpeciesSample, ...]
    geometry: tuple[SystemGeometryRecord, ...]
    response_mode: str
    geometry_profile: str
    shared_phase: float
    private_phase: Mapping[str, float]


def _transition(theta: np.ndarray, phase: float, width: float) -> np.ndarray:
    return 0.5 * (1.0 + np.tanh(np.sin(theta - float(phase)) / float(width)))


def _exact_relation_rotation(h: np.ndarray, amplitude: float) -> np.ndarray:
    side = h >= 0.5
    state = np.zeros((len(h), 2), dtype=float)
    state[~side, 0] = float(amplitude)
    state[side, 1] = float(amplitude)
    return state


def _expanded_labels(values: list, counts: list[int], rng: np.random.Generator) -> list:
    out = []
    for value, count in zip(values, counts):
        out.extend([value] * int(count))
    if len(out) != 20:
        raise RuntimeError("geometry allocation must produce exactly 20 systems per split")
    rng.shuffle(out)
    return out


def _geometry_assignments(profile: GeometryProfile, rng: np.random.Generator) -> dict[str, list]:
    n_levels = [30, 45, 60, 90, 120]
    if profile == "matched":
        n_train = _expanded_labels(n_levels, [4, 4, 4, 4, 4], rng)
        n_eval = _expanded_labels(n_levels, [4, 4, 4, 4, 4], rng)
        cl_train = _expanded_labels(["uniform", "moderate", "strong"], [8, 7, 5], rng)
        cl_eval = _expanded_labels(["uniform", "moderate", "strong"], [8, 7, 5], rng)
    elif profile == "shifted":
        n_train = _expanded_labels(n_levels, [2, 3, 4, 5, 6], rng)
        n_eval = _expanded_labels(n_levels, [6, 5, 4, 3, 2], rng)
        cl_train = _expanded_labels(["uniform", "moderate", "strong"], [12, 6, 2], rng)
        cl_eval = _expanded_labels(["uniform", "moderate", "strong"], [3, 7, 10], rng)
    else:  # pragma: no cover
        raise AssertionError(profile)

    trunc_train = _expanded_labels(["full", "one_sided", "interior_gap"], [10, 5, 5], rng)
    trunc_eval = _expanded_labels(["full", "one_sided", "interior_gap"], [10, 5, 5], rng)
    miss_train = _expanded_labels([0.0, 0.10, 0.25], [8, 7, 5], rng)
    miss_eval = _expanded_labels([0.0, 0.10, 0.25], [8, 7, 5], rng)
    return {
        "n": n_train + n_eval,
        "clustering": cl_train + cl_eval,
        "truncation": trunc_train + trunc_eval,
        "missingness": miss_train + miss_eval,
    }


def _base_angles(n: int, clustering: str, rng: np.random.Generator) -> np.ndarray:
    if clustering == "uniform":
        return rng.uniform(0.0, 2.0 * np.pi, size=n)
    center = rng.uniform(0.0, 2.0 * np.pi)
    if clustering == "moderate":
        use_cluster = rng.random(n) < 0.60
        out = rng.uniform(0.0, 2.0 * np.pi, size=n)
        out[use_cluster] = rng.vonmises(center, 2.0, size=int(np.count_nonzero(use_cluster)))
        return np.mod(out, 2.0 * np.pi)
    if clustering == "strong":
        use_cluster = rng.random(n) < 0.85
        out = rng.uniform(0.0, 2.0 * np.pi, size=n)
        out[use_cluster] = rng.vonmises(center, 8.0, size=int(np.count_nonzero(use_cluster)))
        return np.mod(out, 2.0 * np.pi)
    raise ValueError(f"unknown clustering mode: {clustering}")


def _allowed(theta: np.ndarray, truncation: str, anchor: float) -> np.ndarray:
    rel = np.angle(np.exp(1j * (theta - float(anchor))))
    if truncation == "full":
        return np.ones(len(theta), dtype=bool)
    if truncation == "one_sided":
        return rel <= (0.75 * np.pi)
    if truncation == "interior_gap":
        return np.abs(rel) >= 0.35
    raise ValueError(f"unknown truncation mode: {truncation}")


def _sample_geometry(
    n: int,
    clustering: str,
    truncation: str,
    rng: np.random.Generator,
) -> np.ndarray:
    anchor = rng.uniform(0.0, 2.0 * np.pi)
    accepted: list[float] = []
    attempts = 0
    while len(accepted) < n:
        attempts += 1
        if attempts > 1000:
            raise RuntimeError("geometry rejection sampler failed")
        draw = _base_angles(max(32, 2 * (n - len(accepted))), clustering, rng)
        ok = _allowed(draw, truncation, anchor)
        accepted.extend(map(float, draw[ok]))
    return np.asarray(accepted[:n], dtype=float)


def simulate_q5_world(
    *,
    response_mode: Q5ResponseMode,
    geometry_profile: GeometryProfile,
    amplitude: float = 2.0,
    noise_sd: float = 0.8,
    transition_width: float = 0.2,
    seed: int = 0,
) -> HeterogeneousPairedWorld:
    """Generate Q5 worlds with geometry drawn before response states.

    Forty fixed system labels are used so the first 20 are training systems and
    the remaining 20 are held out. Geometry assignments and coordinates are
    sampled before shared/private response phases and state values are created.
    Missingness is a response-blind row mask applied to geometry before states
    are evaluated on the retained rows.
    """

    allowed = {
        "null_component",
        "null_private_mismatch",
        "null_private_relation",
        "power_shared_mismatch",
        "power_shared_relation",
    }
    if response_mode not in allowed:
        raise ValueError(f"unknown Q5 response mode: {response_mode}")
    if geometry_profile not in {"matched", "shifted"}:
        raise ValueError(f"unknown geometry profile: {geometry_profile}")
    if amplitude < 0 or noise_sd < 0 or transition_width <= 0:
        raise ValueError("invalid response parameters")

    rng = np.random.default_rng(int(seed))
    assignments = _geometry_assignments(geometry_profile, rng)
    retained_theta: list[np.ndarray] = []
    geometry_rows: list[SystemGeometryRecord] = []
    labels = [f"sp_{i:03d}" for i in range(40)]

    # Geometry, truncation and missingness are fully realized before response
    # phases are drawn. This preserves the Q5 response-blind geometry contract.
    for i, name in enumerate(labels):
        nominal_n = int(assignments["n"][i])
        clustering = str(assignments["clustering"][i])
        truncation = str(assignments["truncation"][i])
        missingness = float(assignments["missingness"][i])
        theta = _sample_geometry(nominal_n, clustering, truncation, rng)
        keep = rng.random(nominal_n) >= missingness
        # Deterministically protect the method's minimal graph support without
        # selecting on response values; this is still geometry-only.
        if int(np.count_nonzero(keep)) < 8:
            keep[:] = False
            keep[:8] = True
        theta = theta[keep]
        retained_theta.append(theta)
        geometry_rows.append(
            SystemGeometryRecord(
                species=name,
                split="train" if i < 20 else "eval",
                nominal_n=nominal_n,
                effective_n=int(len(theta)),
                clustering=clustering,
                truncation=truncation,
                missingness=missingness,
            )
        )

    shared_phase = float(rng.uniform(0.0, 2.0 * np.pi))
    private_phases = rng.uniform(0.0, 2.0 * np.pi, size=40)
    samples: list[PairedSpeciesSample] = []
    phase_map: dict[str, float] = {}

    for i, (name, theta) in enumerate(zip(labels, retained_theta)):
        phase = (
            float(private_phases[i])
            if response_mode in {"null_private_mismatch", "null_private_relation"}
            else shared_phase
        )
        phase_map[name] = phase
        h = _transition(theta, phase, transition_width)
        coordinates = np.column_stack((np.cos(theta), np.sin(theta)))

        if response_mode == "null_component":
            latent = float(amplitude) * h
            state_a = latent + rng.normal(0.0, float(noise_sd), size=len(theta))
            state_b = latent + rng.normal(0.0, float(noise_sd), size=len(theta))
        elif response_mode in {"null_private_mismatch", "power_shared_mismatch"}:
            state_a = np.zeros(len(theta), dtype=float)
            state_b = float(amplitude) * h + rng.normal(0.0, float(noise_sd), size=len(theta))
        elif response_mode in {"null_private_relation", "power_shared_relation"}:
            state_a = _exact_relation_rotation(h, amplitude=float(amplitude))
            state_b = np.zeros_like(state_a)
        else:  # pragma: no cover
            raise AssertionError(response_mode)

        samples.append(
            PairedSpeciesSample(
                species=name,
                coordinates=coordinates,
                state_a=state_a,
                state_b=state_b,
            )
        )

    return HeterogeneousPairedWorld(
        samples=tuple(samples),
        geometry=tuple(geometry_rows),
        response_mode=str(response_mode),
        geometry_profile=str(geometry_profile),
        shared_phase=shared_phase,
        private_phase=phase_map,
    )


__all__ = [
    "SystemGeometryRecord",
    "HeterogeneousPairedWorld",
    "simulate_q5_world",
]
