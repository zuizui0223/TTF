from __future__ import annotations

from collections.abc import Iterable

from .codistributed_geometry_null import frozen_uint64_seed


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


__all__ = ["frozen_seed_range", "validate_exact_ranges"]
