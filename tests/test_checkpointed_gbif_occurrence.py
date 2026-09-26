from __future__ import annotations

import pytest

from ttf.checkpointed_gbif_occurrence import (
    deterministic_page_offsets,
    missing_page_offsets,
    occurrence_window_chunks,
    request_timeout_seconds,
    species_state_key,
)


def test_deterministic_page_offsets_cover_bounded_index():
    offsets = deterministic_page_offsets(100_000, page_size=300, maximum_pages=12)
    assert len(offsets) == 12
    assert offsets == tuple(sorted(set(offsets)))
    assert offsets[0] == 0
    assert offsets[-1] <= 99_700


def test_small_result_uses_every_page_once():
    assert deterministic_page_offsets(601, page_size=300, maximum_pages=12) == (
        0,
        300,
        600,
    )


def test_missing_page_offsets_only_returns_unfinished_pages():
    assert missing_page_offsets((0, 300, 600, 900), (300, 900)) == (0, 600)


def test_request_timeout_respects_species_deadline_and_request_cap():
    assert request_timeout_seconds(
        150.0,
        now_monotonic=100.0,
        per_request_cap=30.0,
    ) == 30.0
    assert request_timeout_seconds(
        110.0,
        now_monotonic=100.0,
        per_request_cap=30.0,
    ) == 10.0
    with pytest.raises(TimeoutError):
        request_timeout_seconds(
            100.5,
            now_monotonic=100.0,
            minimum=1.0,
        )


def test_species_state_key_is_stable_and_name_specific():
    assert species_state_key("Papilio machaon") == species_state_key("Papilio machaon")
    assert species_state_key("Papilio machaon") != species_state_key("Pieris rapae")
    assert species_state_key("Papilio machaon").endswith("-Papilio_machaon")


def test_occurrence_window_chunks_preserve_exact_ordinal_window():
    chunks = occurrence_window_chunks(79_800, 300, chunk_size=50)
    assert chunks == (
        (79_800, 50),
        (79_850, 50),
        (79_900, 50),
        (79_950, 50),
        (80_000, 50),
        (80_050, 50),
    )
    covered = [
        index
        for offset, limit in chunks
        for index in range(offset, offset + limit)
    ]
    assert covered == list(range(79_800, 80_100))


def test_occurrence_window_chunks_handle_short_final_window():
    assert occurrence_window_chunks(600, 1, chunk_size=50) == ((600, 1),)
