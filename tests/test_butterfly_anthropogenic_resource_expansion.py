from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "analyze_butterfly_anthropogenic_resource_expansion.py"
)
SPEC = importlib.util.spec_from_file_location(
    "butterfly_anthropogenic_resource_expansion",
    SCRIPT,
)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_host_breadth_bands():
    assert mod.host_band(1) == "1_family"
    assert mod.host_band(2) == "2_families"
    assert mod.host_band(3) == "3_to_5_families"
    assert mod.host_band(5) == "3_to_5_families"
    assert mod.host_band(6) == "6plus_families"


def test_summary_detects_expansion_and_preserves_species_count():
    rows = [
        {
            "species": "A a",
            "host_family_count": 1.0,
            "host_breadth_stratum": "1_family",
            "resolved_host_species": 2,
            "native_resource_units": 10,
            "contemporary_resource_units": 20,
            "introduced_added_units": 10,
            "contemporary_over_native_ratio": 2.0,
            "log_resource_expansion": 0.6466271649250525,
            "introduced_share_of_contemporary": 0.5,
        },
        {
            "species": "B b",
            "host_family_count": 6.0,
            "host_breadth_stratum": "6plus_families",
            "resolved_host_species": 8,
            "native_resource_units": 100,
            "contemporary_resource_units": 100,
            "introduced_added_units": 0,
            "contemporary_over_native_ratio": 1.0,
            "log_resource_expansion": 0.0,
            "introduced_share_of_contemporary": 0.0,
        },
        {
            "species": "C c",
            "host_family_count": 2.0,
            "host_breadth_stratum": "2_families",
            "resolved_host_species": 3,
            "native_resource_units": 30,
            "contemporary_resource_units": 60,
            "introduced_added_units": 30,
            "contemporary_over_native_ratio": 2.0,
            "log_resource_expansion": 0.6768866596881653,
            "introduced_share_of_contemporary": 0.5,
        },
    ]
    summary = mod.summarize_rows(rows)
    assert summary["species"] == 3
    assert summary["species_expanded"] == 2
    assert summary["total_native_species_units"] == 140
    assert summary["total_contemporary_species_units"] == 180
    assert summary["total_introduced_added_species_units"] == 40


def test_spearman_returns_none_for_constant_vector():
    assert mod.spearman([1, 1, 1], [1, 2, 3]) is None
