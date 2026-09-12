from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .heterogeneous_simulate import (
    HeterogeneousPairedWorld,
    SystemGeometryRecord,
    _exact_relation_rotation,
    _sample_geometry,
    _transition,
)
from .mismatch import PairedSpeciesSample

FactorProfile = Literal[
    "balanced",
    "n_only",
    "clustering_only",
    "truncation_only",
    "missingness_only",
    "all_matched",
    "n_shift_only",
    "clustering_shift_only",
    "n_clustering_shift",
    "all_shifted",
]
ResponseMode = Literal[
    "null_private_mismatch",
    "null_private_relation",
    "power_shared_mismatch",
]


def _expanded(values: list, counts: list[int], rng: np.random.Generator) -> list:
    out: list = []
    for value, count in zip(values, counts):
        out.extend([value] * int(count))
    if len(out) != 20:
        raise RuntimeError("factor profile must allocate 20 systems per split")
    rng.shuffle(out)
    return out


def _profile_assignments(profile: FactorProfile, rng: np.random.Generator) -> dict[str, list]:
    n_levels = [30, 45, 60, 90, 120]
    cl_levels = ["uniform", "moderate", "strong"]
    tr_levels = ["full", "one_sided", "interior_gap"]
    miss_levels = [0.0, 0.10, 0.25]

    n_train = [60] * 20
    n_eval = [60] * 20
    cl_train = ["uniform"] * 20
    cl_eval = ["uniform"] * 20
    tr_train = ["full"] * 20
    tr_eval = ["full"] * 20
    miss_train = [0.0] * 20
    miss_eval = [0.0] * 20

    if profile in {"n_only", "all_matched"}:
        n_train = _expanded(n_levels, [4, 4, 4, 4, 4], rng)
        n_eval = _expanded(n_levels, [4, 4, 4, 4, 4], rng)
    if profile in {"clustering_only", "all_matched"}:
        cl_train = _expanded(cl_levels, [8, 7, 5], rng)
        cl_eval = _expanded(cl_levels, [8, 7, 5], rng)
    if profile in {"truncation_only", "all_matched", "all_shifted"}:
        tr_train = _expanded(tr_levels, [10, 5, 5], rng)
        tr_eval = _expanded(tr_levels, [10, 5, 5], rng)
    if profile in {"missingness_only", "all_matched", "all_shifted"}:
        miss_train = _expanded(miss_levels, [8, 7, 5], rng)
        miss_eval = _expanded(miss_levels, [8, 7, 5], rng)

    if profile in {"n_shift_only", "n_clustering_shift", "all_shifted"}:
        n_train = _expanded(n_levels, [2, 3, 4, 5, 6], rng)
        n_eval = _expanded(n_levels, [6, 5, 4, 3, 2], rng)
    if profile in {"clustering_shift_only", "n_clustering_shift", "all_shifted"}:
        cl_train = _expanded(cl_levels, [12, 6, 2], rng)
        cl_eval = _expanded(cl_levels, [3, 7, 10], rng)

    if profile not in {
        "balanced", "n_only", "clustering_only", "truncation_only",
        "missingness_only", "all_matched", "n_shift_only",
        "clustering_shift_only", "n_clustering_shift", "all_shifted",
    }:
        raise ValueError(f"unknown factor profile: {profile}")

    return {
        "n": n_train + n_eval,
        "clustering": cl_train + cl_eval,
        "truncation": tr_train + tr_eval,
        "missingness": miss_train + miss_eval,
    }


def simulate_factor_world(
    *,
    profile: FactorProfile,
    response_mode: ResponseMode,
    amplitude: float,
    noise_sd: float,
    transition_width: float,
    seed: int,
) -> HeterogeneousPairedWorld:
    """Fresh Q5-development world isolating one or more geometry factors.

    Geometry is realized before response phases. The response families are the
    three Q5 cases most informative about the frozen TTF-C failure. This helper
    is diagnostic only and does not modify the frozen Q5 simulator.
    """
    if response_mode not in {
        "null_private_mismatch", "null_private_relation", "power_shared_mismatch"
    }:
        raise ValueError(f"unsupported response mode: {response_mode}")
    rng = np.random.default_rng(int(seed))
    assignments = _profile_assignments(profile, rng)
    labels = [f"sp_{i:03d}" for i in range(40)]
    retained_theta: list[np.ndarray] = []
    geometry_rows: list[SystemGeometryRecord] = []

    for i, name in enumerate(labels):
        nominal_n = int(assignments["n"][i])
        clustering = str(assignments["clustering"][i])
        truncation = str(assignments["truncation"][i])
        missingness = float(assignments["missingness"][i])
        theta = _sample_geometry(nominal_n, clustering, truncation, rng)
        keep = rng.random(nominal_n) >= missingness
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
        phase = float(private_phases[i]) if response_mode.startswith("null_private") else shared_phase
        phase_map[name] = phase
        h = _transition(theta, phase, float(transition_width))
        coordinates = np.column_stack((np.cos(theta), np.sin(theta)))
        if response_mode in {"null_private_mismatch", "power_shared_mismatch"}:
            state_a = np.zeros(len(theta), dtype=float)
            state_b = float(amplitude) * h + rng.normal(0.0, float(noise_sd), size=len(theta))
        else:
            state_a = _exact_relation_rotation(h, float(amplitude))
            state_b = np.zeros_like(state_a)
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
        geometry_profile=str(profile),
        shared_phase=shared_phase,
        private_phase=phase_map,
    )


__all__ = ["FactorProfile", "ResponseMode", "simulate_factor_world"]
