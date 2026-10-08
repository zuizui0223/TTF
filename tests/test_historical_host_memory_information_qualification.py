import copy

import pytest

from scripts.qualify_historical_host_memory_predictor import validate_and_qualify

RULE = {
    "external_admissibility": {"minimum_species_with_complete_frozen_predictor": 500},
    "predictor_information": {
        "minimum_species_fraction_with_unique_fraction_at_least_0.05": 0.70,
        "minimum_panel_median_unique_fraction": 0.10,
        "maximum_standardized_predictor_condition_number": 30.0,
        "minimum_species_fraction_passing_condition_number": 0.90,
    },
}


def _census(unique="0.3", condition="2.0"):
    names = [f"Species t{i:03}" for i in range(642)]
    rows = [
        {
            "species": sp, "status": "complete", "edges": "20",
            "unique_fraction_M_host": unique,
            "predictor_condition_number": condition,
        } for sp in names
    ]
    return rows, names


def test_qualification_uses_both_deterministic_panels_and_frozen_thresholds():
    rows, names = _census()
    result, assignments = validate_and_qualify(rows, names, RULE)
    assert result["decision"] == "PASS_TO_SYNTHETIC_QUALIFICATION"
    assert result["complete_species"] == 642
    assert result["split_counts"] == {"development": 321, "confirmatory": 321}
    assert len({p["species"] for p in assignments}) == 642
    assert assignments == validate_and_qualify(list(reversed(rows)), list(reversed(names)), RULE)[1]


def test_qualification_fails_without_retuning_on_low_unique_information():
    rows, names = _census(unique="0.01", condition="2.0")
    result, _ = validate_and_qualify(rows, names, RULE)
    assert result["decision"] == "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_PREDICTOR_INFORMATION"
    assert "INSUFFICIENT_NONREDUNDANT_HOST_VARIANCE" in result["qualification_failures"]
    assert "PANEL_UNIQUE_HOST_VARIANCE_TOO_SMALL" in result["qualification_failures"]


def test_external_species_attrition_is_separate_from_biological_null():
    rows, names = _census()
    for row in rows[:200]:
        row["status"] = "nonfinite_edge_climate"
        row["unique_fraction_M_host"] = ""
        row["predictor_condition_number"] = ""
    result, assignments = validate_and_qualify(rows, names, RULE)
    assert result["complete_species"] == 442
    assert len(assignments) == 442
    assert result["decision"] == "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_PREDICTOR_INFORMATION"
    assert "INSUFFICIENT_COMPLETE_EXTERNAL_PREDICTOR" in result["qualification_failures"]


def test_source_identity_and_finite_diagnostics_fail_closed():
    rows, names = _census()
    with pytest.raises(RuntimeError, match="candidate species identity"):
        validate_and_qualify(rows, names[:-1] + ["Another species"], RULE)
    rows[0]["unique_fraction_M_host"] = "nan"
    with pytest.raises(RuntimeError, match="non-finite"):
        validate_and_qualify(rows, names, RULE)


def test_frozen_threshold_cannot_change():
    rows, names = _census()
    alt = copy.deepcopy(RULE)
    alt["predictor_information"]["minimum_panel_median_unique_fraction"] = 0.09
    with pytest.raises(RuntimeError, match="thresholds drift"):
        validate_and_qualify(rows, names, alt)
