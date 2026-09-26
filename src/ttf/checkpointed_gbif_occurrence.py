from __future__ import annotations

import hashlib
import math
import time
from typing import Iterable


PILOT_OCCURRENCE_TAG = "butterfly-resource-envelope-occurrence-pilot-v0.1"


def deterministic_page_offsets(
    total_count: int,
    *,
    page_size: int = 300,
    maximum_pages: int = 12,
) -> tuple[int, ...]:
    """Spread a bounded number of pages across the full GBIF result index."""
    total = max(0, min(int(total_count), 100_000))
    if total <= 0:
        return ()
    if page_size < 1 or maximum_pages < 1:
        raise ValueError("page_size and maximum_pages must be positive")
    pages = int(math.ceil(total / page_size))
    if pages <= maximum_pages:
        return tuple(i * page_size for i in range(pages))
    max_offset = max(0, total - page_size)
    if maximum_pages == 1:
        return (0,)
    raw = [
        i * max_offset / (maximum_pages - 1)
        for i in range(maximum_pages)
    ]
    snapped = [
        min(max_offset, int(round(value / page_size)) * page_size)
        for value in raw
    ]
    offsets: list[int] = []
    for value in snapped:
        if value not in offsets:
            offsets.append(value)
    if len(offsets) < maximum_pages:
        for value in range(0, max_offset + 1, page_size):
            if value not in offsets:
                offsets.append(value)
            if len(offsets) == maximum_pages:
                break
    return tuple(sorted(offsets[:maximum_pages]))


def species_state_key(species: str) -> str:
    name = str(species).strip()
    if not name:
        raise ValueError("species must be non-empty")
    digest = hashlib.sha256(
        f"{PILOT_OCCURRENCE_TAG}|{name}".encode("utf-8")
    ).hexdigest()
    return f"{digest[:16]}-{name.replace(' ', '_')}"


def missing_page_offsets(
    planned_offsets: Iterable[int],
    completed_offsets: Iterable[int],
) -> tuple[int, ...]:
    completed = {int(x) for x in completed_offsets}
    return tuple(int(x) for x in planned_offsets if int(x) not in completed)


def request_timeout_seconds(
    deadline_monotonic: float,
    *,
    per_request_cap: float = 30.0,
    minimum: float = 1.0,
    now_monotonic: float | None = None,
) -> float:
    now = time.monotonic() if now_monotonic is None else float(now_monotonic)
    remaining = float(deadline_monotonic) - now
    if remaining < minimum:
        raise TimeoutError("species total deadline exhausted")
    return float(min(per_request_cap, remaining))


__all__ = [
    "PILOT_OCCURRENCE_TAG",
    "deterministic_page_offsets",
    "missing_page_offsets",
    "request_timeout_seconds",
    "species_state_key",
]
