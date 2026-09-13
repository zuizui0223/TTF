from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry


@dataclass(frozen=True)
class GeneticSyntheticWorld:
    """Synthetic pairwise genetic distances on frozen locality graphs."""

    genetic_distance: Mapping[str, np.ndarray]
    residual_edge_signal: Mapping[str, np.ndarray]
    shared_species: tuple[str, ...]
    boundary_normal: Mapping[str, np.ndarray]
    boundary_offset: Mapping[str, float]
    shared_normal: np.ndarray
    shared_offset: float
    residual_amplitude: float
    noise_sd: float


def _unit_normal(rng: np.random.Generator, active: np.ndarray) -> np.ndarray:
    vector = np.zeros(len(active), dtype=float)
    draw = rng.normal(size=int(np.count_nonzero(active)))
    norm = float(np.linalg.norm(draw))
    while norm <= np.finfo(float).tiny:
        draw = rng.normal(size=int(np.count_nonzero(active)))
        norm = float(np.linalg.norm(draw))
    vector[active] = draw / norm
    return vector


def _ibd_strength_map(
    labels: tuple[str, ...],
    value: float | Mapping[str, float],
) -> dict[str, float]:
    if isinstance(value, Mapping):
        missing = set(labels) - set(map(str, value.keys()))
        if missing:
            raise ValueError(f"missing IBD strength for species: {sorted(missing)}")
        out = {name: float(value[name]) for name in labels}
    else:
        out = {name: float(value) for name in labels}
    if any((not np.isfinite(v)) or v < 0.0 for v in out.values()):
        raise ValueError("IBD strengths must be finite and non-negative")
    return out


def simulate_genetic_distance_world(
    geometries: Mapping[str, GeneticSamplingGeometry],
    *,
    shared_fraction: float,
    residual_amplitude: float,
    ibd_strength: float | Mapping[str, float] = 1.0,
    noise_sd: float = 0.10,
    transition_width: float = 0.20,
    noise_dimensions: int = 2,
    seed: int = 0,
) -> GeneticSyntheticWorld:
    """Generate outcome-blind calibration worlds for the genetic TTF interface.

    Nontrivial worlds combine an IBD component with orthogonal endpoint
    transition and locality-noise components. The exactly noise-free,
    zero-residual-amplitude arm is special only in numerical representation:
    because the nuisance estimand is rank based, every strictly positive scalar
    multiple of geographic distance is mathematically equivalent. We therefore
    use canonical positive scale one and copy geographic edge distance exactly
    whenever ``ibd_strength > 0``. This preserves weak ordering bit-for-bit and
    avoids floating rescaling of near-tied edges. With zero IBD strength the pure
    arm remains an exactly zero distance vector.

    No empirical genetic value enters the simulator.
    """
    if not 0.0 <= float(shared_fraction) <= 1.0:
        raise ValueError("shared_fraction must lie in [0, 1]")
    if residual_amplitude < 0 or noise_sd < 0 or transition_width <= 0:
        raise ValueError("amplitude/noise must be non-negative and width positive")
    if noise_dimensions < 1:
        raise ValueError("noise_dimensions must be >= 1")
    if not geometries:
        raise ValueError("at least one species geometry is required")

    labels = tuple(sorted(map(str, geometries.keys())))
    if len(labels) != len(geometries):
        raise ValueError("species labels must be unique")
    dimensions = {geometries[name].coordinates.shape[1] for name in labels}
    if len(dimensions) != 1:
        raise ValueError("all genetic geometries must use the same coordinate dimension")

    strengths = _ibd_strength_map(labels, ibd_strength)
    pooled = np.vstack([geometries[name].coordinates for name in labels])
    center = np.median(pooled, axis=0)
    radial_scale = float(np.sqrt(np.sum(np.var(pooled, axis=0))))
    if radial_scale <= np.sqrt(np.finfo(float).eps):
        raise ValueError("pooled sampling geometry has no spatial extent")

    axis_scale = np.std(pooled, axis=0)
    active = axis_scale > np.sqrt(np.finfo(float).eps)
    if not np.any(active):
        raise ValueError("pooled sampling geometry has no varying coordinate dimension")
    safe_axis_scale = axis_scale.copy()
    safe_axis_scale[~active] = 1.0
    z = {
        name: (geometries[name].coordinates - center) / safe_axis_scale
        for name in labels
    }

    rng = np.random.default_rng(int(seed))
    n_shared = int(round(len(labels) * float(shared_fraction)))
    shared_index = set(map(int, rng.permutation(len(labels))[:n_shared]))
    shared_normal = _unit_normal(rng, active)
    pooled_z = np.vstack([z[name] for name in labels])
    shared_offset = float(np.median(pooled_z @ shared_normal))

    genetic: dict[str, np.ndarray] = {}
    residual_truth: dict[str, np.ndarray] = {}
    normals: dict[str, np.ndarray] = {}
    offsets: dict[str, float] = {}
    shared_names: list[str] = []

    exact_ibd_only = float(residual_amplitude) == 0.0 and float(noise_sd) == 0.0

    for index, name in enumerate(labels):
        geometry = geometries[name]
        coords = geometry.coordinates
        zz = z[name]
        if index in shared_index:
            normal = shared_normal.copy()
            offset = shared_offset
            shared_names.append(name)
        else:
            normal = _unit_normal(rng, active)
            offset = float(np.median(zz @ normal))

        boundary_state = np.tanh(((zz @ normal) - offset) / float(transition_width))
        nodes = geometry.edge_nodes
        geographic = np.linalg.norm(
            coords[nodes[:, 0]] - coords[nodes[:, 1]],
            axis=1,
        )
        boundary_delta = float(residual_amplitude) * (
            boundary_state[nodes[:, 0]] - boundary_state[nodes[:, 1]]
        )

        if exact_ibd_only:
            if strengths[name] > 0.0:
                genetic[name] = geographic.copy()
            else:
                genetic[name] = np.zeros_like(geographic)
        else:
            ibd_component = strengths[name] * geographic / radial_scale
            noise = rng.normal(
                0.0,
                float(noise_sd),
                size=(len(coords), int(noise_dimensions)),
            )
            noise_delta = noise[nodes[:, 0]] - noise[nodes[:, 1]]
            genetic[name] = np.sqrt(
                ibd_component * ibd_component
                + boundary_delta * boundary_delta
                + np.sum(noise_delta * noise_delta, axis=1)
            )
        residual_truth[name] = np.abs(boundary_delta)
        normals[name] = normal.copy()
        offsets[name] = float(offset)

    return GeneticSyntheticWorld(
        genetic_distance=genetic,
        residual_edge_signal=residual_truth,
        shared_species=tuple(sorted(shared_names)),
        boundary_normal=normals,
        boundary_offset=offsets,
        shared_normal=shared_normal.copy(),
        shared_offset=float(shared_offset),
        residual_amplitude=float(residual_amplitude),
        noise_sd=float(noise_sd),
    )
