from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np


def occurrence_priority(species: str, source_key: int) -> str:
    return hashlib.sha256(
        f"relational-env-occ-v0.1|{species}|{int(source_key)}".encode("utf-8")
    ).hexdigest()


def deterministic_page_offsets(total_count: int, *, page_size: int = 300, maximum_pages: int = 20) -> tuple[int, ...]:
    """Return the frozen GBIF page offsets without using any ecological response."""
    total = max(0, min(int(total_count), 100_000))
    if total <= 0:
        return ()
    if page_size < 1 or maximum_pages < 1:
        raise ValueError("page_size and maximum_pages must be positive")
    pages = int(math.ceil(total / page_size))
    if pages <= maximum_pages:
        return tuple(i * page_size for i in range(pages))
    max_offset = max(0, total - page_size)
    raw = np.linspace(0.0, float(max_offset), num=maximum_pages)
    snapped = [
        min(max_offset, int(round(float(value) / page_size)) * page_size)
        for value in raw
    ]
    # For counts > maximum_pages*page_size, 20 distinct page bins should exist.
    offsets: list[int] = []
    for value in snapped:
        if value not in offsets:
            offsets.append(value)
    if len(offsets) < maximum_pages:
        candidates = list(range(0, max_offset + 1, page_size))
        if candidates[-1] != max_offset:
            candidates.append(max_offset)
        for value in candidates:
            if value not in offsets:
                offsets.append(value)
            if len(offsets) == maximum_pages:
                break
        offsets.sort()
    return tuple(offsets[:maximum_pages])


def haversine_km(a: tuple[float, float], b: tuple[float, float]) -> float:
    radius = 6371.0088
    lat1, lon1 = map(math.radians, a)
    lat2, lon2 = map(math.radians, b)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    h = (
        math.sin(dlat / 2.0) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    )
    return 2.0 * radius * math.asin(min(1.0, math.sqrt(h)))


def filter_and_thin_occurrences(
    species: str,
    records: Iterable[Mapping[str, object]],
    *,
    minimum_distance_km: float = 10.0,
    maximum_retained: int = 200,
) -> list[dict[str, object]]:
    """Exact-coordinate deduplication, frozen hash priority, then greedy thinning."""
    if minimum_distance_km <= 0 or maximum_retained < 1:
        raise ValueError("invalid thinning parameters")
    by_coordinate: dict[tuple[float, float], dict[str, object]] = {}
    for record in records:
        try:
            key = int(record["source_key"])
            lat = float(record["latitude"])
            lon = float(record["longitude"])
        except (KeyError, TypeError, ValueError):
            continue
        if key <= 0 or not np.isfinite(lat) or not np.isfinite(lon):
            continue
        if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
            continue
        # Literal 0,0 is frozen as invalid for this environmental occurrence layer.
        if lat == 0.0 and lon == 0.0:
            continue
        coord = (lat, lon)
        current = by_coordinate.get(coord)
        if current is None or key < int(current["source_key"]):
            by_coordinate[coord] = {
                "species": species,
                "source_key": key,
                "latitude": lat,
                "longitude": lon,
            }

    ordered = sorted(
        by_coordinate.values(),
        key=lambda row: (
            occurrence_priority(species, int(row["source_key"])),
            int(row["source_key"]),
        ),
    )
    retained: list[dict[str, object]] = []
    retained_coords: list[tuple[float, float]] = []
    for row in ordered:
        coord = (float(row["latitude"]), float(row["longitude"]))
        if all(haversine_km(coord, other) >= minimum_distance_km for other in retained_coords):
            retained.append(dict(row))
            retained_coords.append(coord)
            if len(retained) >= maximum_retained:
                break
    for index, row in enumerate(retained):
        row["priority_rank"] = int(index)
        row["priority_sha256"] = occurrence_priority(species, int(row["source_key"]))
    return retained


def canonical_pca(environment: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Standardize pooled environment and return canonicalized eigenvectors/scores."""
    x = np.asarray(environment, dtype=float)
    if x.ndim != 2 or x.shape[0] < 3 or x.shape[1] < 2 or not np.isfinite(x).all():
        raise ValueError("environment must be a finite n x p matrix")
    mean = x.mean(axis=0)
    sd = x.std(axis=0, ddof=0)
    if np.any(sd <= np.finfo(float).eps):
        raise ValueError("environment contains a constant variable")
    z = (x - mean) / sd
    covariance = (z.T @ z) / len(z)
    eigval, eigvec = np.linalg.eigh(covariance)
    order = np.argsort(eigval)[::-1]
    eigval = eigval[order]
    eigvec = eigvec[:, order]
    for axis in range(eigvec.shape[1]):
        column = eigvec[:, axis]
        pivot = int(np.argmax(np.abs(column)))
        if column[pivot] < 0:
            eigvec[:, axis] *= -1.0
    score = z @ eigvec
    return mean, sd, eigval, eigvec


def project_whiten(
    environment: np.ndarray,
    mean: np.ndarray,
    sd: np.ndarray,
    eigval: np.ndarray,
    eigvec: np.ndarray,
    *,
    axes: int = 2,
) -> np.ndarray:
    x = np.asarray(environment, dtype=float)
    z = (x - np.asarray(mean, dtype=float)) / np.asarray(sd, dtype=float)
    score = z @ np.asarray(eigvec, dtype=float)[:, :axes]
    scale = np.sqrt(np.asarray(eigval, dtype=float)[:axes])
    if np.any(scale <= np.finfo(float).eps):
        raise ValueError("cannot whiten a zero-variance PC axis")
    return score / scale


def frozen_grid_bounds(pca_sample_whitened: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(pca_sample_whitened, dtype=float)
    if x.ndim != 2 or x.shape[1] != 2 or not np.isfinite(x).all():
        raise ValueError("whitened PCA sample must be finite n x 2")
    low = np.percentile(x, 0.5, axis=0)
    high = np.percentile(x, 99.5, axis=0)
    span = high - low
    if np.any(span <= np.finfo(float).eps):
        raise ValueError("degenerate PCA grid extent")
    return low - 0.10 * span, high + 0.10 * span


def grid_centres(low: np.ndarray, high: np.ndarray, *, cells: int = 50) -> np.ndarray:
    if cells < 2:
        raise ValueError("cells must be >=2")
    lo = np.asarray(low, dtype=float)
    hi = np.asarray(high, dtype=float)
    edges0 = np.linspace(lo[0], hi[0], cells + 1)
    edges1 = np.linspace(lo[1], hi[1], cells + 1)
    c0 = 0.5 * (edges0[:-1] + edges0[1:])
    c1 = 0.5 * (edges1[:-1] + edges1[1:])
    a, b = np.meshgrid(c0, c1, indexing="ij")
    return np.column_stack((a.ravel(), b.ravel()))


def normalized_kde_grid(points: np.ndarray, centres: np.ndarray) -> np.ndarray:
    x = np.asarray(points, dtype=float)
    q = np.asarray(centres, dtype=float)
    if x.ndim != 2 or x.shape[1] != 2 or len(x) < 2:
        raise ValueError("KDE needs at least two 2-D points")
    h = float(len(x) ** (-1.0 / 6.0))
    delta = q[:, None, :] - x[None, :, :]
    density = np.exp(-0.5 * np.sum(delta * delta, axis=2) / (h * h)).mean(axis=1)
    total = float(density.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("KDE grid has zero/non-finite mass")
    return density / total


def schoener_d(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError("density grids must be equal-length vectors")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("density grids must be finite")
    if not np.isclose(a.sum(), 1.0, atol=1e-9) or not np.isclose(b.sum(), 1.0, atol=1e-9):
        raise ValueError("density grids must sum to one")
    value = 1.0 - 0.5 * float(np.abs(a - b).sum())
    return float(min(1.0, max(0.0, value)))


def relation_pairs(
    density: np.ndarray,
    source_indices: Sequence[int],
    target_indices: Sequence[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.asarray(density, dtype=float)
    src = np.asarray(source_indices, dtype=np.int64)
    tgt = np.asarray(target_indices, dtype=np.int64)
    source_out = np.repeat(src, len(tgt))
    target_out = np.tile(tgt, len(src))
    relation = np.empty(len(source_out), dtype=float)
    cursor = 0
    for s in src:
        block = 1.0 - 0.5 * np.abs(d[int(s)][None, :] - d[tgt]).sum(axis=1)
        relation[cursor : cursor + len(tgt)] = np.clip(block, 0.0, 1.0)
        cursor += len(tgt)
    return source_out, target_out, relation
