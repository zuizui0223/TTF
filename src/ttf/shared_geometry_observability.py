from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .geometry import SpeciesGeometry
from .private_geometry_control import pooled_affine_frame


@dataclass(frozen=True)
class SharedGeometryObservability:
    directions: int
    equal_species_mean_contrast_by_direction: np.ndarray
    species_mean_contrast: Mapping[str, float]
    species_median_contrast: Mapping[str, float]
    species_nontrivial_fraction: Mapping[str, float]

    def summary(self) -> dict[str, object]:
        x = self.equal_species_mean_contrast_by_direction
        species = np.asarray(list(self.species_mean_contrast.values()), dtype=float)
        return {
            "directions": self.directions,
            "equal_species_mean_contrast": {
                "mean": float(np.mean(x)),
                "median": float(np.median(x)),
                "q05": float(np.quantile(x, 0.05)),
                "q25": float(np.quantile(x, 0.25)),
                "q75": float(np.quantile(x, 0.75)),
                "q95": float(np.quantile(x, 0.95)),
                "min": float(np.min(x)),
                "max": float(np.max(x)),
            },
            "species_mean_contrast": {
                "mean": float(np.mean(species)),
                "median": float(np.median(species)),
                "q10": float(np.quantile(species, 0.10)),
                "q90": float(np.quantile(species, 0.90)),
                "min": float(np.min(species)),
                "max": float(np.max(species)),
            },
        }


def _random_normals(
    *,
    dimension: int,
    active: np.ndarray,
    n_directions: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(int(seed))
    use = np.asarray(active, dtype=bool)
    draws = rng.normal(size=(int(n_directions), int(np.count_nonzero(use))))
    norms = np.linalg.norm(draws, axis=1)
    bad = norms <= np.finfo(float).tiny
    while np.any(bad):
        draws[bad] = rng.normal(size=(int(np.count_nonzero(bad)), draws.shape[1]))
        norms = np.linalg.norm(draws, axis=1)
        bad = norms <= np.finfo(float).tiny
    draws /= norms[:, None]
    normals = np.zeros((int(n_directions), int(dimension)), dtype=float)
    normals[:, use] = draws
    return normals


def shared_transition_observability(
    geometries: Sequence[SpeciesGeometry],
    graphs: Mapping[str, np.ndarray],
    *,
    transition_width: float = 0.2,
    n_directions: int = 2048,
    seed: int = 20260912,
    nontrivial_edge_contrast: float = 0.25,
) -> SharedGeometryObservability:
    """Measure outcome-free observability of the simulator's shared transition.

    For each random shared normal, this reproduces the Gate-I shared hyperplane:
    coordinates are pooled-affine standardized and the hyperplane passes through
    the pooled median projection.  It then measures the noise-free absolute tanh
    contrast on each species' already-fixed graph.  No trait outcome, estimator,
    bandwidth, training/evaluation split, or candidate score enters this audit.

    The quantity is diagnostic only.  It is not a pass/fail criterion and must
    not be used to discard difficult geometries post hoc.
    """
    data = tuple(geometries)
    if not data:
        raise ValueError("at least one geometry is required")
    if transition_width <= 0 or n_directions < 32:
        raise ValueError("transition_width must be positive and directions >= 32")
    if not 0 <= nontrivial_edge_contrast <= 2:
        raise ValueError("nontrivial edge contrast must lie in [0,2]")
    names = [item.species for item in data]
    if set(names) != set(graphs):
        raise ValueError("graphs must match geometry species exactly")

    center, scale, active = pooled_affine_frame(data)
    normals = _random_normals(
        dimension=data[0].coordinates.shape[1],
        active=active,
        n_directions=n_directions,
        seed=seed,
    )
    standardized = {
        item.species: (item.coordinates - center) / scale for item in data
    }
    pooled_z = np.vstack([standardized[item.species] for item in data])
    pooled_projection = pooled_z @ normals.T
    offsets = np.median(pooled_projection, axis=0)

    species_by_direction = []
    species_mean: dict[str, float] = {}
    species_median: dict[str, float] = {}
    species_nontrivial: dict[str, float] = {}
    for item in data:
        z = standardized[item.species]
        nodes = np.asarray(graphs[item.species], dtype=np.int64)
        if nodes.ndim != 2 or nodes.shape[1] != 2 or len(nodes) < 1:
            raise ValueError(f"invalid graph for {item.species}")
        projection = z @ normals.T
        latent = np.tanh((projection - offsets[None, :]) / float(transition_width))
        contrast = np.abs(latent[nodes[:, 0], :] - latent[nodes[:, 1], :])
        mean_by_direction = contrast.mean(axis=0)
        species_by_direction.append(mean_by_direction)
        species_mean[item.species] = float(np.mean(mean_by_direction))
        species_median[item.species] = float(np.median(mean_by_direction))
        species_nontrivial[item.species] = float(
            np.mean(np.max(contrast, axis=0) >= float(nontrivial_edge_contrast))
        )

    equal_species = np.mean(np.vstack(species_by_direction), axis=0)
    return SharedGeometryObservability(
        directions=int(n_directions),
        equal_species_mean_contrast_by_direction=equal_species,
        species_mean_contrast=species_mean,
        species_median_contrast=species_median,
        species_nontrivial_fraction=species_nontrivial,
    )
