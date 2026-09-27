from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "analyze_butterfly_host_plant_drivers.py"
)
SPEC = importlib.util.spec_from_file_location(
    "butterfly_host_plant_drivers",
    SCRIPT,
)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_concentration_summary_tracks_top_hosts():
    rows = [
        {"total_fractional_butterfly_unit_credit": 50.0},
        {"total_fractional_butterfly_unit_credit": 30.0},
        {"total_fractional_butterfly_unit_credit": 10.0},
        {"total_fractional_butterfly_unit_credit": 10.0},
    ]
    result = mod.concentration_summary(rows)
    assert result["total_fractional_butterfly_unit_credit"] == 100.0
    assert result["host_plants_with_positive_credit"] == 4
    assert result["share_top_1_host"] == 0.5
    assert result["share_top_5_hosts"] == 1.0
    assert result["hosts_to_reach_25_percent"] == 1
    assert result["hosts_to_reach_50_percent"] == 1
    assert result["hosts_to_reach_75_percent"] == 2


def test_concentration_summary_rejects_empty_credit():
    import pytest

    with pytest.raises(RuntimeError):
        mod.concentration_summary([])
