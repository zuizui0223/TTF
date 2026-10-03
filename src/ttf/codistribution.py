from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class GeographicCoopportunity:
    """Response-blind directed source-target geographic support.

    The primary relation is symmetric_min_coverage.  Direction-specific
    coverages are retained because they diagnose nested or asymmetric range
    support without using any genetic response.
    """

    source: str
    target: str
    source_to_target_coverage: float
    target_to_source_coverage: float
    symmetric_min_coverage: float
    centroid_distance: float
    source_edges: int
    target_edges: int


def _validated_points(points: np.ndarray, *, label: str) -> np.ndarray:
    values = np.asarray(points, dtype=float)
    if values.ndim != 2 or values.shape[0] < 1 or values.shape[1] < 1:
        raise ValueError(f"{label} must be a non-empty n x d array")
    if not np.isfinite(values).all():
        raise ValueError(f"{label} must be finite")
    return values


def coverage_fraction(
    query_points: np.ndarray,
    reference_points: np.ndarray,
    *,
    radius: float = 500.0,
    chunk_size: int = 128,
) -> float:
    """Fraction of query points within radius of at least one reference point."""
    query = _validated_points(query_points, label="query_points")
    reference = _validated_points(reference_points, label="reference_points")
    if query.shape[1] != reference.shape[1]:
        raise ValueError("query and reference point dimensionality must match")
    if not np.isfinite(radius) or float(radius) <= 0.0:
        raise ValueError("radius must be finite and positive")
    if int(chunk_size) < 1:
        raise ValueError("chunk_size must be positive")

    threshold2 = float(radius) ** 2
    covered = 0
    for start in range(0, len(query), int(chunk_size)):
        stop = min(start + int(chunk_size), len(query))
        delta = query[start:stop, None, :] - reference[None, :, :]
        nearest2 = np.min(np.sum(delta * delta, axis=2), axis=1)
        covered += int(np.count_nonzero(nearest2 <= threshold2))
    return float(covered / len(query))


def geographic_coopportunity(
    source: str,
    target: str,
    source_midpoints: np.ndarray,
    target_midpoints: np.ndarray,
    *,
    radius: float = 500.0,
    chunk_size: int = 128,
) -> GeographicCoopportunity:
    """Compute the frozen symmetric geographic co-opportunity relation.

    Coverage is evaluated in both directions because a small range can be
    completely nested within a much larger range.  The primary relation takes
    the minimum of the two directional coverages, requiring reciprocal support.
    """
    source_name = str(source)
    target_name = str(target)
    if not source_name or not target_name:
        raise ValueError("source and target labels must be non-empty")

    source_points = _validated_points(source_midpoints, label="source_midpoints")
    target_points = _validated_points(target_midpoints, label="target_midpoints")
    if source_points.shape[1] != target_points.shape[1]:
        raise ValueError("source and target point dimensionality must match")

    source_to_target = coverage_fraction(
        source_points,
        target_points,
        radius=float(radius),
        chunk_size=int(chunk_size),
    )
    target_to_source = coverage_fraction(
        target_points,
        source_points,
        radius=float(radius),
        chunk_size=int(chunk_size),
    )
    centroid_distance = float(
        np.linalg.norm(np.mean(source_points, axis=0) - np.mean(target_points, axis=0))
    )
    return GeographicCoopportunity(
        source=source_name,
        target=target_name,
        source_to_target_coverage=source_to_target,
        target_to_source_coverage=target_to_source,
        symmetric_min_coverage=float(min(source_to_target, target_to_source)),
        centroid_distance=centroid_distance,
        source_edges=int(len(source_points)),
        target_edges=int(len(target_points)),
    )


def geographic_coopportunity_table(
    source_midpoints: Mapping[str, np.ndarray],
    target_midpoints: Mapping[str, np.ndarray],
    *,
    source_species: Sequence[str] | None = None,
    target_species: Sequence[str] | None = None,
    radius: float = 500.0,
    chunk_size: int = 128,
    exclude_self: bool = True,
) -> tuple[GeographicCoopportunity, ...]:
    """Build the response-blind directed dyad geometry table.

    The table is geometry only.  It does not select source pools and it accepts
    no genetic response.  All eligible source-target dyads are retained.
    """
    source_labels = (
        tuple(map(str, source_species))
        if source_species is not None
        else tuple(sorted(map(str, source_midpoints)))
    )
    target_labels = (
        tuple(map(str, target_species))
        if target_species is not None
        else tuple(sorted(map(str, target_midpoints)))
    )
    if not source_labels or not target_labels:
        raise ValueError("at least one source and target species are required")
    if len(set(source_labels)) != len(source_labels):
        raise ValueError("source species must be unique")
    if len(set(target_labels)) != len(target_labels):
        raise ValueError("target species must be unique")
    if not set(source_labels) <= set(map(str, source_midpoints)):
        raise ValueError("source_species contains labels without source midpoints")
    if not set(target_labels) <= set(map(str, target_midpoints)):
        raise ValueError("target_species contains labels without target midpoints")

    out: list[GeographicCoopportunity] = []
    for target in target_labels:
        for source in source_labels:
            if exclude_self and source == target:
                continue
            out.append(
                geographic_coopportunity(
                    source,
                    target,
                    np.asarray(source_midpoints[source], dtype=float),
                    np.asarray(target_midpoints[target], dtype=float),
                    radius=float(radius),
                    chunk_size=int(chunk_size),
                )
            )
    return tuple(out)


__all__ = [
    "GeographicCoopportunity",
    "coverage_fraction",
    "geographic_coopportunity",
    "geographic_coopportunity_table",
]
