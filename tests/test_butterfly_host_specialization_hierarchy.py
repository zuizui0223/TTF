from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "analyze_butterfly_host_specialization_hierarchy.py"
)
SPEC = importlib.util.spec_from_file_location(
    "butterfly_host_specialization_hierarchy",
    SCRIPT,
)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_summarize_stratum_detects_species_level_portfolio_gradient():
    rows = [
        {
            "species": "a",
            "resolved_host_species": 1,
            "effective_contributor_number": 1.0,
            "maximum_single_host_fractional_share": 1.0,
            "introduced_added_units": 20,
        },
        {
            "species": "b",
            "resolved_host_species": 2,
            "effective_contributor_number": 1.5,
            "maximum_single_host_fractional_share": 0.8,
            "introduced_added_units": 30,
        },
        {
            "species": "c",
            "resolved_host_species": 5,
            "effective_contributor_number": 3.0,
            "maximum_single_host_fractional_share": 0.4,
            "introduced_added_units": 40,
        },
        {
            "species": "d",
            "resolved_host_species": 10,
            "effective_contributor_number": 5.0,
            "maximum_single_host_fractional_share": 0.2,
            "introduced_added_units": 60,
        },
        {
            "species": "e",
            "resolved_host_species": 20,
            "effective_contributor_number": 8.0,
            "maximum_single_host_fractional_share": 0.1,
            "introduced_added_units": 80,
        },
        {
            "species": "f",
            "resolved_host_species": 30,
            "effective_contributor_number": 10.0,
            "maximum_single_host_fractional_share": 0.08,
            "introduced_added_units": 100,
        },
    ]
    result = mod.summarize_stratum(rows)
    assert result["species"] == 6
    assert result["spearman"][
        "resolved_host_species_vs_effective_contributor_number"
    ] > 0.99
    assert result["spearman"][
        "resolved_host_species_vs_maximum_single_host_fractional_share"
    ] < -0.99
    assert set(result["resolved_host_species_tertiles"]) == {"low", "mid", "high"}
