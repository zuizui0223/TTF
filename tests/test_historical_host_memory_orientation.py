import pytest
from scripts.describe_historical_host_memory_orientation import describe


def test_orientation_counts_positive_and_negative_species():
    roles={"A":"development","B":"confirmatory"}
    rows=[
        {"species":"A","M_host":"1","M_self":"1"},
        {"species":"A","M_host":"2","M_self":"2"},
        {"species":"A","M_host":"3","M_self":"3"},
        {"species":"B","M_host":"-3","M_self":"3"},
        {"species":"B","M_host":"-2","M_self":"2"},
        {"species":"B","M_host":"-1","M_self":"1"},
    ]
    results=describe(rows,roles)
    assert results["all"]["species"]==2
    assert results["development"]["species_mean_host_memory_positive"]==1
    assert results["confirmatory"]["species_mean_host_memory_negative"]==1
    assert results["all"]["median_within_species_correlation_M_host_M_self"]==pytest.approx(0.0,abs=1e-12)


def test_orientation_does_not_treat_missing_species_as_zero():
    with pytest.raises(RuntimeError,match="do not equal"):
        describe([{"species":"A","M_host":"1","M_self":"1"}],{"A":"development","B":"confirmatory"})


def test_nonfinite_climate_predictor_abstains():
    with pytest.raises(RuntimeError,match="nonfinite predictor"):
        describe([{"species":"A","M_host":"nan","M_self":"1"}],{"A":"development"})
