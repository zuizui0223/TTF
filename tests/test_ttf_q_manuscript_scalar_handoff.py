import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_ttf_q_manuscript_scalar_handoff_matches_frozen_receipts():
    handoff=json.loads(
        (ROOT/"manuscript/generated/ttf_q_scalar_results_v0.1.json").read_text()
    )
    c1=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_c1_result_receipt_v0.1.json").read_text()
    )
    c2=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_c2_result_receipt_v0.1.json").read_text()
    )
    ladder=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_b1_c1_information_ladder_v0.1.json").read_text()
    )
    bc=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_bc_calibrated_detectability_result_v0.2.json").read_text()
    )
    known=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_known_truth_decomposition_v0.1.json").read_text()
    )
    kd=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_known_truth_detectability_v0.2.json").read_text()
    )
    final=json.loads(
        (ROOT/"benchmarks/frozen/relational_program_final_state_v0.3.json").read_text()
    )

    assert handoff["information_decomposition"]["known_truth_scenarios"] == known["decomposition_grid"]["scenarios"]
    assert handoff["information_decomposition"]["known_truth_max_absolute_error"] == known["exact_recovery"]["maximum_absolute_error_across_all_nested_information_quantities"]

    assert handoff["climate_design"]["admissible_species"] == c1["design"]["admissible_species"]
    assert handoff["climate_design"]["development_dyads"] == c1["design"]["development_dyads"]
    assert handoff["climate_design"]["confirmatory_dyads"] == c1["design"]["confirmatory_dyads"]

    assert handoff["current_climate"]["development_endpoint_survival"] == ladder["development"]["B1_current_climate"]["source_target_fe_retained_variance_fraction"]
    assert handoff["current_climate"]["development_unique_before_geography"] == ladder["development"]["B1_current_climate"]["total_unique_variance_fraction"]
    assert handoff["current_climate"]["development_unique_after_geography"] == bc["current_climate"]["total_unique_variance_fraction_after_geography"]

    hist=handoff["historical_climate"]
    assert hist["development_endpoint_survival"] == c1["development"]["source_target_fe_retained_variance_fraction"]
    assert hist["development_unique_before_geography"] == c1["development"]["total_unique_variance_fraction"]
    assert hist["confirmatory_unique_before_geography"] == c1["confirmatory"]["total_unique_variance_fraction"]
    assert hist["development_unique_after_geography"] == c2["development"]["C2_total_unique_variance_fraction"]
    assert hist["confirmatory_unique_after_geography"] == c2["confirmatory"]["C2_total_unique_variance_fraction"]
    assert hist["development_geography_retained_fraction"] == c2["development"]["geography_retained_fraction"]
    assert hist["confirmatory_geography_retained_fraction"] == c2["confirmatory"]["geography_retained_fraction"]

    current_env=bc["current_climate"]["calibrated_evaluable_envelope"]
    historical_env=bc["historical_given_current"]["calibrated_evaluable_envelope"]
    assert handoff["current_climate"]["development_evaluable_mde_A0_A1_A2_A3"] == [
        current_env[k]["evaluable_grid_mde"] for k in ("A0","A1","A2","A3")
    ]
    assert hist["development_evaluable_mde_A0_A1_A2_A3"] == [
        historical_env[k]["evaluable_grid_mde"] for k in ("A0","A1","A2","A3")
    ]

    for name in ("reference","endpoint_loss","control_redundancy","source_concentration"):
        source=kd["scenarios"][name]
        target=handoff["known_truth_detectability"][name]
        assert target["unique_fraction"] == source["observed_total_unique_variance_fraction"]
        if name == "reference":
            assert target["signal_effective_sources"] == source["signal_effective_sources"]
        expected=[
            source["evaluable_grid_mde_by_A"][k]
            for k in ("A0","A1","A2","A3")
        ]
        assert target["evaluable_mde_A0_A1_A2_A3"] == expected

    assert handoff["boundary"]["relational_v03_final_state"] == final["final_state"]
    assert handoff["boundary"]["genetic_response_used"] is False
