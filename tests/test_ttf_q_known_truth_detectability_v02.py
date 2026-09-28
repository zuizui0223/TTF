import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_calibrated_known_truth_result_separates_power_from_calibration():
    contract=json.loads(
        (ROOT/"docs/supporting/ttf_q_known_truth_detectability_v0.2.json").read_text()
    )
    result=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_known_truth_detectability_v0.2.json").read_text()
    )

    assert contract["status"]=="FROZEN_AFTER_V01_DIAGNOSTIC_BEFORE_V02_RESULTS"
    assert result["status"]=="FROZEN_COMPLETE_CALIBRATED_KNOWN_TRUTH_DETECTABILITY_BENCHMARK"
    assert result["single_cell_reference"]["all_four_binary_qualified"] is False
    assert result["single_cell_reference"]["same_binary_classification_pairs"]==6
    assert result["single_cell_reference"]["same_binary_but_different_evaluable_mde_pairs"]==5

    reference=result["scenarios"]["reference"]
    endpoint=result["scenarios"]["endpoint_loss"]
    redundant=result["scenarios"]["control_redundancy"]
    concentrated=result["scenarios"]["source_concentration"]

    assert all(cell["pass"] for cell in reference["type1_calibration_by_A"].values())
    assert all(cell["pass"] for cell in endpoint["type1_calibration_by_A"].values())
    assert all(cell["pass"] for cell in redundant["type1_calibration_by_A"].values())
    assert not any(cell["pass"] for cell in concentrated["type1_calibration_by_A"].values())

    assert reference["evaluable_grid_mde_by_A"]=={
        "A0":0.05,"A1":0.08,"A2":0.08,"A3":0.10
    }
    assert endpoint["evaluable_grid_mde_by_A"]=={
        "A0":0.08,"A1":0.10,"A2":0.15,"A3":0.20
    }
    assert redundant["evaluable_grid_mde_by_A"]==endpoint["evaluable_grid_mde_by_A"]
    assert all(value is None for value in concentrated["evaluable_grid_mde_by_A"].values())

    assert reference["observed_total_unique_variance_fraction"] > endpoint["observed_total_unique_variance_fraction"]
    assert reference["observed_total_unique_variance_fraction"] > redundant["observed_total_unique_variance_fraction"]
    assert concentrated["observed_total_unique_variance_fraction"] == reference["observed_total_unique_variance_fraction"]
    assert concentrated["signal_effective_sources"] < 0.25 * reference["signal_effective_sources"]

    conclusion=result["methodological_conclusion"]
    assert conclusion["single_binary_label_conflates_distinct_design_states"] is True
    assert conclusion["effective_endpoint_counts_are_not_effective_sample_sizes"] is True
    assert result["boundary"]==contract["boundary"]
