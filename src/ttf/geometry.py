from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Mapping, Sequence

import numpy as np

from .core import SpeciesSample


@dataclass(frozen=True)
class SpeciesGeometry:
    """Empirical sampling geometry with no empirical trait values attached.

    Gate-I semi-synthetic qualification deliberately separates the observed
    sampling frame (species identities, coordinates, record counts and optional
    dependence blocks) from synthetic trait generation.
    """

    species: str
    coordinates: np.ndarray
    blocks: np.ndarray | None = None

    def __post_init__(self) -> None:
        coords = np.asarray(self.coordinates, dtype=float)
        if coords.ndim != 2 or coords.shape[0] < 2 or coords.shape[1] < 1:
            raise ValueError("coordinates must be n x d with n >= 2 and d >= 1")
        if not np.isfinite(coords).all():
            raise ValueError("coordinates must be finite")
        blocks = None if self.blocks is None else np.asarray(self.blocks)
        if blocks is not None and blocks.shape != (coords.shape[0],):
            raise ValueError("blocks must have one value per observation")
        if not str(self.species):
            raise ValueError("species must be non-empty")
        object.__setattr__(self, "coordinates", coords)
        object.__setattr__(self, "blocks", blocks)

    @classmethod
    def from_sample(cls, sample: SpeciesSample) -> "SpeciesGeometry":
        """Drop empirical trait values while retaining the sampling frame."""
        return cls(
            species=sample.species,
            coordinates=sample.coordinates.copy(),
            blocks=None if sample.blocks is None else sample.blocks.copy(),
        )


@dataclass(frozen=True)
class GeometrySyntheticWorld:
    samples: tuple[SpeciesSample, ...]
    shared_species: tuple[str, ...]
    boundary_normal: Mapping[str, np.ndarray]
    boundary_offset: Mapping[str, float]
    shared_normal: np.ndarray
    shared_offset: float
    shared_fraction: float
    amplitude: float
    geometry_fingerprint: str


def _validated_geometry(
    geometries: Sequence[SpeciesGeometry],
) -> tuple[SpeciesGeometry, ...]:
    data = tuple(sorted(geometries, key=lambda item: str(item.species)))
    if len(data) < 4:
        raise ValueError("at least four species geometries are required")
    labels = [item.species for item in data]
    if len(set(labels)) != len(labels):
        raise ValueError("species labels must be unique")
    dimensions = {item.coordinates.shape[1] for item in data}
    if len(dimensions) != 1:
        raise ValueError("all species must use the same coordinate dimension")
    return data


def geometry_fingerprint(geometries: Sequence[SpeciesGeometry]) -> str:
    """Stable SHA-256 fingerprint of the fixed Gate-I sampling frame."""
    data = _validated_geometry(geometries)
    digest = hashlib.sha256()
    for item in data:
        label = item.species.encode("utf-8")
        digest.update(len(label).to_bytes(8, "little"))
        digest.update(label)
        coords = np.asarray(item.coordinates, dtype="<f8", order="C")
        digest.update(np.asarray(coords.shape, dtype="<i8").tobytes())
        digest.update(coords.tobytes(order="C"))
        if item.blocks is None:
            digest.update(b"blocks:none")
        else:
            digest.update(b"blocks:present")
            for value in item.blocks:
                encoded = str(value).encode("utf-8")
                digest.update(len(encoded).to_bytes(8, "little"))
                digest.update(encoded)
    return digest.hexdigest()


def _unit_normal(rng: np.random.Generator, active: np.ndarray) -> np.ndarray:
    vector = np.zeros(len(active), dtype=float)
    draw = rng.normal(size=int(np.count_nonzero(active)))
    norm = float(np.linalg.norm(draw))
    while norm <= np.finfo(float).tiny:
        draw = rng.normal(size=int(np.count_nonzero(active)))
        norm = float(np.linalg.norm(draw))
    vector[active] = draw / norm
    return vector


def simulate_fixed_geometry_boundary_world(
    geometries: Sequence[SpeciesGeometry],
    *,
    shared_fraction: float,
    amplitude: float,
    noise_sd: float = 0.8,
    transition_width: float = 0.2,
    seed: int = 0,
) -> GeometrySyntheticWorld:
    """Generate synthetic traits on an unchanged empirical sampling geometry.

    Coordinates, species counts, record counts and optional blocks are retained
    exactly.  Only trait values are simulated.  A globally shared transition is
    a common hyperplane in a pooled affine-standardized coordinate system;
    private transitions receive species-specific hyperplanes whose offsets pass
    through the median projection of that species.  The standardization affects
    synthetic trait generation only: TTF graphs and boundary estimation still
    use the original coordinates supplied by the caller.

    ``transition_width`` is measured in pooled standardized-coordinate units.
    This keeps the synthetic signal definition invariant to coordinate units
    while preserving the raw empirical geometry for the estimator itself.
    """
    data = _validated_geometry(geometries)
    if not 0.0 <= shared_fraction <= 1.0:
        raise ValueError("shared_fraction must lie in [0, 1]")
    if amplitude < 0 or noise_sd < 0 or transition_width <= 0:
        raise ValueError("amplitude/noise must be non-negative and width positive")

    pooled = np.vstack([item.coordinates for item in data])
    center = np.median(pooled, axis=0)
    raw_scale = np.std(pooled, axis=0)
    active = raw_scale > np.sqrt(np.finfo(float).eps)
    if not np.any(active):
        raise ValueError("sampling geometry has no varying coordinate dimension")
    scale = raw_scale.copy()
    scale[~active] = 1.0
    standardized = {
        item.species: (item.coordinates - center) / scale for item in data
    }

    rng = np.random.default_rng(int(seed))
    n_species = len(data)
    n_shared = int(round(n_species * float(shared_fraction)))
    shared_index = set(map(int, rng.permutation(n_species)[:n_shared]))

    shared_normal = _unit_normal(rng, active)
    pooled_z = np.vstack([standardized[item.species] for item in data])
    shared_offset = float(np.median(pooled_z @ shared_normal))

    samples: list[SpeciesSample] = []
    normals: dict[str, np.ndarray] = {}
    offsets: dict[str, float] = {}
    shared_names: list[str] = []
    for index, item in enumerate(data):
        z = standardized[item.species]
        if index in shared_index:
            normal = shared_normal.copy()
            offset = shared_offset
            shared_names.append(item.species)
        else:
            normal = _unit_normal(rng, active)
            offset = float(np.median(z @ normal))
        latent = float(amplitude) * np.tanh(
            ((z @ normal) - offset) / float(transition_width)
        )
        trait = latent + rng.normal(0.0, float(noise_sd), len(z))
        samples.append(
            SpeciesSample(
                species=item.species,
                coordinates=item.coordinates.copy(),
                trait=trait,
                blocks=None if item.blocks is None else item.blocks.copy(),
            )
        )
        normals[item.species] = normal.copy()
        offsets[item.species] = float(offset)

    return GeometrySyntheticWorld(
        samples=tuple(samples),
        shared_species=tuple(sorted(shared_names)),
        boundary_normal=normals,
        boundary_offset=offsets,
        shared_normal=shared_normal.copy(),
        shared_offset=float(shared_offset),
        shared_fraction=float(shared_fraction),
        amplitude=float(amplitude),
        geometry_fingerprint=geometry_fingerprint(data),
    )
