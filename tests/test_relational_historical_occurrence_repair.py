from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path("scripts/acquire_relational_historical_occurrences_v02.py")


def load():
    spec = importlib.util.spec_from_file_location("hist_occ_repair", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_repair_partition_covers_exact_1000_once():
    m = load()
    names = [f"Species {i:04d}" for i in range(1000)]
    seen = []
    for batch in range(m.FROZEN_BATCH_COUNT):
        block = m.batch_species(names, batch, m.FROZEN_BATCH_SIZE)
        assert len(block) == 4
        seen.extend(block)
    assert seen == names
    assert len(set(seen)) == 1000


def test_repair_partition_rejects_candidate_or_batch_drift():
    m = load()
    names = [f"Species {i:04d}" for i in range(1000)]
    with pytest.raises(RuntimeError, match="1000-species"):
        m.batch_species(names[:-1], 0, 4)
    with pytest.raises(RuntimeError, match="batch size drift"):
        m.batch_species(names, 0, 5)
    with pytest.raises(ValueError, match="batch index"):
        m.batch_species(names, 250, 4)


def test_timeout_is_terminal_request_error_not_ecological_failure():
    m = load()
    ledger = m.timeout_ledger("Alpha beta", m.FROZEN_SPECIES_TIMEOUT_SECONDS)
    assert ledger["species"] == "Alpha beta"
    assert ledger["status"] == "REQUEST_ERROR"
    assert "3300" in ledger["error"]


def test_frozen_repair_constants_leave_job_margin():
    m = load()
    assert m.FROZEN_BATCH_SIZE == 4
    assert m.FROZEN_BATCH_COUNT == 250
    assert m.FROZEN_SPECIES_TIMEOUT_SECONDS == 3300
    assert m.FROZEN_BATCH_SIZE * m.FROZEN_SPECIES_TIMEOUT_SECONDS == 13200
    assert 13200 < 240 * 60
