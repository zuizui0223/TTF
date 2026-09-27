from __future__ import annotations

import math

import numpy as np

from .relational_environment import nearest_coverage


EARTH_RADIUS_KM = 6371.0088


def latlon_to_ecef_km(
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> np.ndarray:
    lat = np.radians(np.asarray(latitude, dtype=float))
    lon = np.radians(np.asarray(longitude, dtype=float))
    if lat.shape != lon.shape or lat.ndim != 1 or len(lat) == 0:
        raise ValueError("latitude and longitude must be aligned non-empty vectors")
    if not (np.isfinite(lat).all() and np.isfinite(lon).all()):
        raise ValueError("coordinates must be finite")
    cos_lat = np.cos(lat)
    return EARTH_RADIUS_KM * np.column_stack((
        cos_lat * np.cos(lon),
        cos_lat * np.sin(lon),
        np.sin(lat),
    ))


def chord_radius_km(arc_radius_km: float) -> float:
    """Convert a great-circle arc radius to the equivalent 3-D chord radius."""
    radius = float(arc_radius_km)
    if not 0 < radius < math.pi * EARTH_RADIUS_KM:
        raise ValueError("arc radius must lie in (0, pi * Earth radius)")
    return 2.0 * EARTH_RADIUS_KM * math.sin(
        radius / (2.0 * EARTH_RADIUS_KM)
    )


def directed_occurrence_coverage(
    target_latlon: np.ndarray,
    source_latlon: np.ndarray,
    *,
    radius_km: float,
) -> float:
    """Fraction of target occurrence points within radius of any source point.

    This is a response-blind directed range co-opportunity diagnostic. It uses
    great-circle distance on the spherical Earth through an equivalent ECEF
    chord threshold.
    """
    target = np.asarray(target_latlon, dtype=float)
    source = np.asarray(source_latlon, dtype=float)
    if (
        target.ndim != 2
        or source.ndim != 2
        or target.shape[1] != 2
        or source.shape[1] != 2
        or len(target) == 0
        or len(source) == 0
    ):
        raise ValueError("occurrence coordinates must be non-empty n x 2 arrays")
    target_xyz = latlon_to_ecef_km(target[:, 0], target[:, 1])
    source_xyz = latlon_to_ecef_km(source[:, 0], source[:, 1])
    return nearest_coverage(
        target_xyz,
        source_xyz,
        radius=chord_radius_km(radius_km),
    )
