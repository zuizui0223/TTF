from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry


@dataclass(frozen=True)
class GeneticSyntheticWorld:
    """Synthetic pairwise genetic distances on frozen locality graphs.

    Genetic distances are generated from a latent Euclidean state so they remain
    symmetric, non-negative, and endpoint-dependent rather than independent edge
    draws.  The base latent coordinates reproduce strict monotone IBD when the
    residual amplitude and noise are zero.  A separate latent boundary dimension
    adds private or shared place-specific differentiation.
    """

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

    For locality ``i`` in species ``s`` the latent genetic state contains three
    independent pieces:

    1. a scaled copy of the observed geographic coordinate (ordinary IBD);
    2. one private/shared transition coordinate;
    3. optional independent locality noise coordinates.

    Genetic edge distance is Euclidean distance between the latent states at the
    two edge endpoints.  Consequently, when ``residual_amplitude=noise_sd=0``,
    genetic distance is exactly a positive scalar multiple of geographic edge
    length and the endpoint-safe rank IBD residualizer should remove it exactly.

    Shared residual species use one common hyperplane in the pooled standardized
    geographic frame.  Private species use independent hyperplanes whose offsets
    pass through the median projection of that species.  No empirical genetic
    value enters the simulator.
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

    # A single scalar scale preserves Euclidean geographic-distance ordering
    # exactly in the IBD-only arm.  Per-axis standardization is used separately
    # only to define the orientation of synthetic transition fields.
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
        base = strengths[name] * (coords - center) / radial_scale
        boundary = (float(residual_amplitude) * boundary_state)[:, None]
        noise = rng.normal(
            0.0,
            float(noise_sd),
            size=(len(coords), int(noise_dimensions)),
        )
        latent = np.hstack([base, boundary, noise])

        nodes = geometry.edge_nodes
        delta = latent[nodes[:, 0]] - latent[nodes[:, 1]]
        genetic[name] = np.sqrt(np.sum(delta * delta, axis=1))
        residual_truth[name] = float(residual_amplitude) * np.abs(
            boundary_state[nodes[:, 0]] - boundary_state[nodes[:, 1]]
        )
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
