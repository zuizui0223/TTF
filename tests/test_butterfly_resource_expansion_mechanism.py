from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "analyze_butterfly_resource_expansion_mechanism.py"
)
SPEC = importlib.util.spec_from_file_location(
    "butterfly_resource_expansion_mechanism",
    SCRIPT,
)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_fractional_credit_sums_exactly_over_added_union():
    result = mod.fractional_host_contributions(
        {
            "host1": {"A", "B"},
            "host2": {"B", "C"},
        },
        {"A", "B", "C"},
    )
    assert result["introduced_added_units"] == 3
    assert result["contributing_host_species"] == 2
    assert abs(
        result["maximum_single_host_fractional_share"] - 0.5
    ) < 1e-12
    assert abs(result["top_two_host_fractional_share"] - 1.0) < 1e-12
    assert abs(result["effective_contributor_number"] - 2.0) < 1e-12
    assert abs(
        result["fraction_added_units_with_multiple_contributing_hosts"]
        - 1 / 3
    ) < 1e-12
    assert abs(sum(result["fractional_credits"].values()) - 3.0) < 1e-12


def test_dominant_host_is_detected_with_overlap():
    result = mod.fractional_host_contributions(
        {
            "dominant": {"A", "B", "C"},
            "minor": {"C"},
        },
        {"A", "B", "C"},
    )
    assert abs(
        result["maximum_single_host_fractional_share"] - (2.5 / 3)
    ) < 1e-12
    assert abs(result["top_two_host_fractional_share"] - 1.0) < 1e-12
    assert result["effective_contributor_number"] < 2.0


def test_no_added_units_returns_null_contribution_metrics():
    result = mod.fractional_host_contributions(
        {"host1": set()},
        set(),
    )
    assert result["introduced_added_units"] == 0
    assert result["contributing_host_species"] == 0
    assert result["maximum_single_host_fractional_share"] is None
    assert result["effective_contributor_number"] is None
