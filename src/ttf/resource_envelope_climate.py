from __future__ import annotations

import hashlib
import math
from typing import Iterable, Sequence

import numpy as np


CLIMATE_SPLIT_TAG = "butterfly-resource-climate-crossfit-v0.1"


def observed_unit_split(
    species: str,
    observed_units: Iterable[str],
    *,
    train_fraction: float = 0.5,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    units = sorted({str(unit).strip() for unit in observed_units if str(unit).strip()})
    if not 0.0 < train_fraction < 1.0:
        raise ValueError("train_fraction must lie strictly between zero and one")
    if len(units) < 2:
        return tuple(units), ()
    ranked = sorted(
        units,
        key=lambda unit: (
            hashlib.sha256(
                f"{CLIMATE_SPLIT_TAG}|{species}|{unit}".encode("utf-8")
            ).hexdigest(),
            unit,
        ),
    )
    cut = max(1, min(len(ranked) - 1, int(round(len(ranked) * train_fraction))))
    return tuple(sorted(ranked[:cut])), tuple(sorted(ranked[cut:]))


def climate_mismatch(
    unit_vectors: np.ndarray,
    training_vectors: np.ndarray,
) -> np.ndarray:
    units = np.asarray(unit_vectors, dtype=float)
    train = np.asarray(training_vectors, dtype=float)
    if units.ndim != 2 or train.ndim != 2 or units.shape[1] != train.shape[1]:
        raise ValueError("unit and training climate arrays must be n x p with shared p")
    if len(train) < 2:
        raise ValueError("at least two training climate rows are required")
    mean = np.mean(train, axis=0)
    sd = np.std(train, axis=0)
    usable = sd > np.sqrt(np.finfo(float).eps)
    if not np.any(usable):
        raise ValueError("training climate has no varying axis")
    z = (units[:, usable] - mean[usable]) / sd[usable]
    return np.sqrt(np.mean(z * z, axis=1))


def probability_greater(
    left: Sequence[float],
    right: Sequence[float],
) -> float | None:
    """Return P(left > right) + 0.5 P(tie) across all cross-pairs."""
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return None
    greater = 0
    ties = 0
    for x in a:
        greater += int(np.count_nonzero(x > b))
        ties += int(np.count_nonzero(x == b))
    return float((greater + 0.5 * ties) / (len(a) * len(b)))


def radical_inverse(index: int, base: int) -> float:
    if index < 1 or base < 2:
        raise ValueError("index must be >=1 and base >=2")
    value = 0.0
    factor = 1.0 / base
    n = int(index)
    while n:
        value += factor * (n % base)
        n //= base
        factor /= base
    return value


def deterministic_interior_points(geometry, *, maximum_points: int = 16):
    """Generate bounded deterministic lon/lat samples covering polygon components."""
    from shapely.geometry import Point

    if maximum_points < 1:
        raise ValueError("maximum_points must be positive")
    if geometry.is_empty:
        return ()

    if geometry.geom_type == "Polygon":
        components = [geometry]
    elif geometry.geom_type == "MultiPolygon":
        components = [g for g in geometry.geoms if not g.is_empty]
    else:
        components = [geometry]

    components = sorted(
        components,
        key=lambda g: (-float(g.area), tuple(map(float, g.bounds))),
    )
    total_area = sum(max(0.0, float(g.area)) for g in components)
    points: list[tuple[float, float]] = []

    def add_point(point):
        xy = (float(point.x), float(point.y))
        if xy not in points:
            points.append(xy)

    # Guarantee representation of the largest components first.
    for component in components[:maximum_points]:
        add_point(component.representative_point())
        if len(points) >= maximum_points:
            return tuple(points)

    for component_index, component in enumerate(components):
        if len(points) >= maximum_points:
            break
        area = max(0.0, float(component.area))
        if total_area > 0:
            target = max(
                1,
                int(round(maximum_points * area / total_area)),
            )
        else:
            target = 1
        existing_before = len(points)
        minx, miny, maxx, maxy = map(float, component.bounds)
        if not (
            math.isfinite(minx)
            and math.isfinite(miny)
            and math.isfinite(maxx)
            and math.isfinite(maxy)
            and maxx >= minx
            and maxy >= miny
        ):
            continue
        for k in range(1, 4001):
            if len(points) - existing_before >= target or len(points) >= maximum_points:
                break
            u = radical_inverse(k + 17 * component_index, 2)
            v = radical_inverse(k + 31 * component_index, 3)
            candidate = Point(minx + u * (maxx - minx), miny + v * (maxy - miny))
            if component.covers(candidate):
                add_point(candidate)

    return tuple(points[:maximum_points])


__all__ = [
    "CLIMATE_SPLIT_TAG",
    "climate_mismatch",
    "deterministic_interior_points",
    "observed_unit_split",
    "probability_greater",
    "radical_inverse",
]
