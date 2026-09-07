from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .core import SpeciesSample, knn_edges


@dataclass(frozen=True)
class GeometryNormalization:
    center: np.ndarray
    scale: float
    coordinates: Mapping[str, np.ndarray]


@dataclass(frozen=True)
class SemiSyntheticWorld:
    """Synthetic trait transitions placed on frozen empirical species geometry."""

    samples: tuple[SpeciesSample, ...]
    shared_species: tuple[str, ...]
    shared_normal: np.ndarray | None
    shared_offset: float | None
    crossing_species: tuple[str, ...]
    bandwidth: float


def normalize_geometry(
    geometry: Mapping[str, np.ndarray],
) -> GeometryNormalization:
    """Center and scale a multi-species geometry using only coordinates.

    The pooled RMS radial distance is one after normalization.  This makes a
    bandwidth rule portable across empirical data sets without looking at trait
    values or benchmark outcomes.
    """
    if len(geometry) < 2:
        raise ValueError("at least two species are required")
    clean: dict[str, np.ndarray] = {}
    for species, coords in geometry.items():
        x = np.asarray(coords, dtype=float)
        if x.ndim != 2 or x.shape[1] != 2 or len(x) < 2:
            raise ValueError("each geometry must be n x 2 with n >= 2")
        if not np.isfinite(x).all():
            raise ValueError("coordinates must be finite")
        clean[str(species)] = x
    if len(set(clean)) != len(clean):
        raise ValueError("species labels must be unique")

    pooled = np.vstack(tuple(clean.values()))
    center = pooled.mean(axis=0)
    scale = float(np.sqrt(np.mean(np.sum((pooled - center) ** 2, axis=1))))
    if scale <= np.finfo(float).tiny:
        raise ValueError("geometry has zero spatial scale")
    normalized = {species: (coords - center) / scale for species, coords in clean.items()}
    return GeometryNormalization(center=center, scale=scale, coordinates=normalized)


def geometry_bandwidth(
    geometry: Mapping[str, np.ndarray],
    *,
    k: int = 4,
) -> float:
    """Outcome-free kernel bandwidth: pooled median species-local kNN edge length."""
    normalized = normalize_geometry(geometry).coordinates
    lengths: list[np.ndarray] = []
    for coords in normalized.values():
        nodes = knn_edges(coords, k=k)
        lengths.append(np.linalg.norm(coords[nodes[:, 0]] - coords[nodes[:, 1]], axis=1))
    pooled = np.concatenate(lengths)
    value = float(np.median(pooled))
    if not np.isfinite(value) or value <= 0:
        raise RuntimeError("could not derive a positive geometry bandwidth")
    return value


def _random_unit_vector(rng: np.random.Generator) -> np.ndarray:
    angle = float(rng.uniform(0.0, 2.0 * np.pi))
    return np.asarray([np.cos(angle), np.sin(angle)], dtype=float)


def _draw_shared_boundary(
    geometry: Mapping[str, np.ndarray],
    *,
    rng: np.random.Generator,
    min_crossing_species: int,
    max_attempts: int = 2000,
) -> tuple[np.ndarray, float, tuple[str, ...]]:
    names = tuple(sorted(geometry))
    min_crossing = min(max(int(min_crossing_species), 1), len(names))
    pooled = np.vstack([geometry[name] for name in names])
    for _ in range(int(max_attempts)):
        normal = _random_unit_vector(rng)
        pooled_projection = pooled @ normal
        # Keep the shared boundary away from the extreme geographic tails while
        # allowing its exact position to vary among worlds.
        q = float(rng.uniform(0.35, 0.65))
        offset = float(np.quantile(pooled_projection, q))
        crossing = tuple(
            name
            for name in names
            if float(np.min(geometry[name] @ normal)) < offset
            < float(np.max(geometry[name] @ normal))
        )
        if len(crossing) >= min_crossing:
            return normal, offset, crossing
    raise RuntimeError(
        f"could not draw a shared boundary crossed by {min_crossing} species"
    )


def simulate_on_geometry(
    geometry: Mapping[str, np.ndarray],
    *,
    shared_fraction: float,
    amplitude: float,
    noise_sd: float = 0.8,
    transition_width: float = 0.20,
    min_shared_crossing_species: int = 6,
    k: int = 4,
    seed: int = 0,
) -> SemiSyntheticWorld:
    """Generate private/shared transition worlds on frozen empirical geometry.

    ``shared_fraction=0`` is deliberately adversarial: every species receives a
    strong spatial transition, but its orientation and location are private.
    ``shared_fraction=1`` gives all species one common boundary.  The common
    boundary is drawn using coordinates only and must cross at least the
    requested number of species, defining an explicit empirical identifiability
    condition rather than silently choosing an easy benchmark after seeing
    outcomes.
    """
    if not 0.0 <= float(shared_fraction) <= 1.0:
        raise ValueError("shared_fraction must lie in [0, 1]")
    if float(amplitude) < 0 or float(noise_sd) < 0 or float(transition_width) <= 0:
        raise ValueError("amplitude/noise must be non-negative and width positive")

    normalized = normalize_geometry(geometry).coordinates
    names = tuple(sorted(normalized))
    rng = np.random.default_rng(int(seed))
    n_shared = int(round(float(shared_fraction) * len(names)))

    shared_normal: np.ndarray | None = None
    shared_offset: float | None = None
    crossing: tuple[str, ...] = ()
    shared_names: tuple[str, ...] = ()
    if n_shared > 0:
        shared_normal, shared_offset, crossing = _draw_shared_boundary(
            normalized,
            rng=rng,
            min_crossing_species=min(min_shared_crossing_species, n_shared),
        )
        crossing_order = list(crossing)
        rng.shuffle(crossing_order)
        selected = crossing_order[: min(n_shared, len(crossing_order))]
        if len(selected) < n_shared:
            remaining = [name for name in names if name not in selected]
            rng.shuffle(remaining)
            selected.extend(remaining[: n_shared - len(selected)])
        shared_names = tuple(sorted(selected))

    samples: list[SpeciesSample] = []
    shared_set = set(shared_names)
    for name in names:
        coords = normalized[name]
        if name in shared_set:
            assert shared_normal is not None and shared_offset is not None
            normal = shared_normal
            offset = shared_offset
        else:
            normal = _random_unit_vector(rng)
            projection = coords @ normal
            offset = float(np.median(projection))
        latent = float(amplitude) * np.tanh(
            ((coords @ normal) - float(offset)) / float(transition_width)
        )
        trait = latent + rng.normal(0.0, float(noise_sd), size=len(coords))
        samples.append(SpeciesSample(species=name, coordinates=coords, trait=trait))

    return SemiSyntheticWorld(
        samples=tuple(samples),
        shared_species=shared_names,
        shared_normal=shared_normal,
        shared_offset=shared_offset,
        crossing_species=crossing,
        bandwidth=geometry_bandwidth(normalized, k=k),
    )
