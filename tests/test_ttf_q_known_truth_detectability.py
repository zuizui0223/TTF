import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_known_truth_detectability_shows_binary_gate_information_loss():
    contract=json.loads(
        (ROOT/"docs/supporting/ttf_q_known_truth_detectability_v0.1.json").read_text()
    )
    result=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_known_truth_detectability_v0.1.json").read_text()
    )

    assert contract["status"]=="FROZEN_BEFORE_DETECTABILITY_RESULTS"
    assert result["single_cell_reference"]["all_four_scenarios_binary_qualified"] is False
    assert result["single_cell_reference"]["same_binary_classification_pair_count"]==6
    assert result["single_cell_reference"]["same_binary_but_different_envelope_pair_count"]==6

    reference=result["scenarios"]["reference"]
    endpoint=result["scenarios"]["endpoint_loss"]
    redundant=result["scenarios"]["control_redundancy"]
    concentrated=result["scenarios"]["source_concentration"]

    assert endpoint["total_unique_variance_fraction"] < reference["total_unique_variance_fraction"]
    assert redundant["total_unique_variance_fraction"] < reference["total_unique_variance_fraction"]
    assert endpoint["minimum_detectable_effect_by_private_amplitude"]["A3"] > reference["minimum_detectable_effect_by_private_amplitude"]["A3"]
    assert redundant["minimum_detectable_effect_by_private_amplitude"]["A0"] > reference["minimum_detectable_effect_by_private_amplitude"]["A0"]

    assert concentrated["total_unique_variance_fraction"] == reference["total_unique_variance_fraction"]
    assert concentrated["signal_effective_sources"] < 0.25 * reference["signal_effective_sources"]
    assert result["methodological_findings"]["endpoint_concentration_is_not_a_monotone_proxy_for_power"] is True
    assert result["methodological_findings"]["effective_endpoint_counts_are_independent_signal_concentration_diagnostics_not_effective_sample_sizes"] is True
    assert result["boundary"]==contract["boundary"]
