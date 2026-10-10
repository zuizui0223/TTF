"""Combinatorial masks-unknown robustness assertions; no empirical sequence file."""
import pytest
from scripts.audit_historical_host_memory_information_attrition_bound import worst_case_bound


def fixture(*, low=16, medium=0, ill=0):
    roles=(
        [{"species":f"C{i:03d}","panel":"confirmatory"} for i in range(321)]
        + [{"species":f"D{i:03d}","panel":"development"} for i in range(320)]
    )
    diag=[]
    for i,r in enumerate(roles):
        is_confirm=i<321
        n=i if is_confirm else -1
        val=(0.04 if is_confirm and n<low else
             0.06 if is_confirm and n<low+medium else 0.24)
        cond=35 if is_confirm and n<ill else 2
        diag.append({"species":r["species"],"status":"complete",
                     "unique_fraction_M_host":str(val),
                     "predictor_condition_number":str(cond)})
    diag.append({"species":"EXCLUDED","status":"nonfinite_edge_climate",
                 "unique_fraction_M_host":"","predictor_condition_number":""})
    rule={
        "post_mask_requalification":{"minimum_confirmatory_survivor_species":200},
        "predictor_information":{
            "minimum_species_fraction_with_unique_fraction_at_least_0.05":0.70,
            "minimum_panel_median_unique_fraction":0.10,
            "maximum_standardized_predictor_condition_number":30,
            "minimum_species_fraction_passing_condition_number":0.90,
        },
    }
    return diag,roles,rule


def test_arbitrary_mask_survivor_set_information_guaranteed_with_sixteen_low():
    v=worst_case_bound(*fixture())
    assert v["confirmatory_species"]==321
    assert v["unique_fraction_below_0_05_in_full_confirmatory"]==16
    assert v["worst_case_at_200"]["minimum_possible_fraction_unique_ge_0_05"]==pytest.approx(0.92)
    assert v["worst_case_at_200"]["minimum_possible_median_unique_fraction"]==pytest.approx(0.24)
    assert v["all_possible_confirmatory_survivor_subsets_ge_200_pass_information"] is True


def test_too_many_low_unique_information_species_breaks_guarantee():
    v=worst_case_bound(*fixture(low=81))
    assert v["worst_case_at_200"]["minimum_possible_fraction_unique_ge_0_05"]==pytest.approx(119/200)
    assert v["threshold_checks"]["fraction_unique_ge_0_70"] is False
    assert v["all_possible_confirmatory_survivor_subsets_ge_200_pass_information"] is False


def test_median_can_fail_even_when_fraction_above_point_zero_five_passes():
    v=worst_case_bound(*fixture(low=0,medium=150))
    assert v["worst_case_at_200"]["minimum_possible_median_unique_fraction"]==pytest.approx(0.06)
    assert v["threshold_checks"]["fraction_unique_ge_0_70"] is True
    assert v["threshold_checks"]["median_unique_ge_0_10"] is False


def test_condition_number_can_fail_independently():
    v=worst_case_bound(*fixture(ill=21))
    assert v["worst_case_at_200"]["minimum_possible_fraction_condition_le_30"]==pytest.approx(179/200)
    assert v["threshold_checks"]["fraction_condition_ge_0_90"] is False


def test_source_role_and_threshold_drift_fail_closed():
    d,roles,rule=fixture()
    roles[0]["panel"]="development"
    with pytest.raises(RuntimeError,match="assignment drift"):
        worst_case_bound(d,roles,rule)
    d,roles,rule=fixture()
    rule["post_mask_requalification"]["minimum_confirmatory_survivor_species"]=199
    with pytest.raises(RuntimeError,match="thresholds drift"):
        worst_case_bound(d,roles,rule)
