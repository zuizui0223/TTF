from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry


@dataclass(frozen=True)
class GroupedGeneticSyntheticWorld:
    genetic_distance: Mapping[str, np.ndarray]
    residual_edge_signal: Mapping[str, np.ndarray]
    group_by_species: Mapping[str, str]
    group_normal: Mapping[str, np.ndarray]
    group_offset: Mapping[str, float]
    residual_amplitude: float
    noise_sd: float


def _unit_normal(
    rng: np.random.Generator,
    active: np.ndarray,
) -> np.ndarray:
    vector = np.zeros(len(active), dtype=float)
    draw = rng.normal(size=int(np.count_nonzero(active)))
    norm = float(np.linalg.norm(draw))
    while norm <= np.finfo(float).tiny:
        draw = rng.normal(size=int(np.count_nonzero(active)))
        norm = float(np.linalg.norm(draw))
    vector[active] = draw / norm
    return vector


def simulate_grouped_genetic_distance_world(
    geometries: Mapping[str, GeneticSamplingGeometry],
    groups: Mapping[str, str],
    *,
    residual_amplitude: float,
    ibd_strength: float = 1.0,
    noise_sd: float = 0.10,
    transition_width: float = 0.20,
    noise_dimensions: int = 2,
    seed: int = 0,
) -> GroupedGeneticSyntheticWorld:
    """Simulate spatial genetic structure shared only within frozen groups.

    Group labels must be fixed without genetic outcomes, for example taxonomic
    order. Species in the same group share one boundary field; different groups
    receive independent fields. A private-null world is obtained by assigning
    every species a unique group label.
    """
    if not geometries:
        raise ValueError("at least one geometry is required")
    if (
        residual_amplitude < 0
        or ibd_strength < 0
        or noise_sd < 0
        or transition_width <= 0
    ):
        raise ValueError("invalid simulation parameter")
    if noise_dimensions < 1:
        raise ValueError("noise_dimensions must be positive")
    labels = tuple(sorted(map(str, geometries)))
    if set(labels) != set(map(str, groups)):
        raise ValueError("groups must label every species exactly once")
    group_by_species = {
        name: str(groups[name]).strip()
        for name in labels
    }
    if any(not value for value in group_by_species.values()):
        raise ValueError("group labels must be non-empty")
    dimensions = {
        np.asarray(geometries[name].coordinates).shape[1]
        for name in labels
    }
    if len(dimensions) != 1:
        raise ValueError("all geometries must have the same dimensionality")

    pooled = np.vstack([
        geometries[name].coordinates
        for name in labels
    ])
    center = np.median(pooled, axis=0)
    radial_scale = float(np.sqrt(np.sum(np.var(pooled, axis=0))))
    if radial_scale <= np.sqrt(np.finfo(float).eps):
        raise ValueError("pooled geometry has no spatial extent")
    axis_scale = np.std(pooled, axis=0)
    active = axis_scale > np.sqrt(np.finfo(float).eps)
    if not np.any(active):
        raise ValueError("pooled geometry has no varying coordinate dimension")
    safe_axis_scale = axis_scale.copy()
    safe_axis_scale[~active] = 1.0
    standardized = {
        name: (geometries[name].coordinates - center) / safe_axis_scale
        for name in labels
    }

    rng = np.random.default_rng(int(seed))
    unique_groups = tuple(sorted(set(group_by_species.values())))
    normals = {
        group: _unit_normal(rng, active)
        for group in unique_groups
    }
    offsets: dict[str, float] = {}
    for group in unique_groups:
        members = [
            name
            for name in labels
            if group_by_species[name] == group
        ]
        pooled_group = np.vstack([
            standardized[name]
            for name in members
        ])
        offsets[group] = float(
            np.median(pooled_group @ normals[group])
        )

    genetic: dict[str, np.ndarray] = {}
    residual: dict[str, np.ndarray] = {}
    for name in labels:
        geometry = geometries[name]
        nodes = np.asarray(
            geometry.edge_nodes,
            dtype=np.int64,
        )
        group = group_by_species[name]
        state = np.tanh(
            (
                (standardized[name] @ normals[group])
                - offsets[group]
            )
            / float(transition_width)
        )
        geographic = np.linalg.norm(
            geometry.coordinates[nodes[:, 0]]
            - geometry.coordinates[nodes[:, 1]],
            axis=1,
        )
        boundary_delta = float(residual_amplitude) * (
            state[nodes[:, 0]]
            - state[nodes[:, 1]]
        )
        noise = rng.normal(
            0.0,
            float(noise_sd),
            size=(
                geometry.n_localities,
                int(noise_dimensions),
            ),
        )
        noise_delta = noise[nodes[:, 0]] - noise[nodes[:, 1]]
        ibd_component = (
            float(ibd_strength)
            * geographic
            / radial_scale
        )
        genetic[name] = np.sqrt(
            ibd_component * ibd_component
            + boundary_delta * boundary_delta
            + np.sum(noise_delta * noise_delta, axis=1)
        )
        residual[name] = np.abs(boundary_delta)

    return GroupedGeneticSyntheticWorld(
        genetic_distance=genetic,
        residual_edge_signal=residual,
        group_by_species=group_by_species,
        group_normal={
            name: normal.copy()
            for name, normal in normals.items()
        },
        group_offset=dict(offsets),
        residual_amplitude=float(residual_amplitude),
        noise_sd=float(noise_sd),
    )


__all__ = [
    "GroupedGeneticSyntheticWorld",
    "simulate_grouped_genetic_distance_world",
]
