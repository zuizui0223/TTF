from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import numpy as np

from .codistributed_geometry_null import GeometryNullCenter, frozen_uint64_seed


def frozen_seed_range(
    master_seed: int,
    namespace: str,
    cell: str,
    start: int,
    stop: int,
) -> tuple[int, ...]:
    """Return frozen seeds for original replicate indices [start, stop)."""
    start = int(start)
    stop = int(stop)
    if start < 0 or stop <= start:
        raise ValueError("seed range must satisfy 0 <= start < stop")
    return tuple(
        frozen_uint64_seed(master_seed, namespace, cell, replicate)
        for replicate in range(start, stop)
    )


def validate_exact_ranges(
    observed: Iterable[tuple[int, int]],
    expected: Iterable[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    """Require exactly the predeclared shard ranges, with no gaps or duplicates."""
    got = tuple(sorted((int(a), int(b)) for a, b in observed))
    want = tuple(sorted((int(a), int(b)) for a, b in expected))
    if got != want:
        raise RuntimeError(f"shard range drift: {got} != {want}")
    if not want:
        raise RuntimeError("empty shard contract")
    cursor = want[0][0]
    if cursor != 0:
        raise RuntimeError("shard contract must start at replicate zero")
    for start, stop in want:
        if start != cursor or stop <= start:
            raise RuntimeError("shard contract has a gap, overlap, or empty range")
        cursor = stop
    return want


def aggregate_equal_amplitude_center(
    cell_means: Mapping[float, np.ndarray],
    amplitudes: Sequence[float],
    *,
    worlds_per_amplitude: int,
) -> GeometryNullCenter:
    """Aggregate already-computed per-amplitude means in frozen amplitude order.

    Each input vector must itself be the exact single-batch mean for that
    amplitude. This helper only performs the final equal-amplitude vstack/mean
    used by the monolithic runner.
    """
    ordered=tuple(float(a) for a in amplitudes)
    if not ordered or len(set(ordered))!=len(ordered):
        raise ValueError("amplitudes must be unique and non-empty")
    if set(map(float,cell_means))!=set(ordered):
        raise ValueError("center amplitude cells do not match frozen amplitudes")
    n=int(worlds_per_amplitude)
    if n<1:
        raise ValueError("worlds_per_amplitude must be positive")
    vectors=[]
    width=None
    for amplitude in ordered:
        x=np.asarray(cell_means[amplitude],dtype=float)
        if x.ndim!=1 or len(x)==0 or not np.isfinite(x).all():
            raise ValueError(f"invalid center mean for amplitude {amplitude:g}")
        if width is None:
            width=len(x)
        elif len(x)!=width:
            raise ValueError("center amplitude dyad-count drift")
        vectors.append(x)
    return GeometryNullCenter(
        mu0=np.mean(np.vstack(vectors),axis=0),
        amplitudes=ordered,
        worlds_per_amplitude=tuple([n]*len(ordered)),
    )


__all__ = [
    "aggregate_equal_amplitude_center",
    "frozen_seed_range",
    "validate_exact_ranges",
]
