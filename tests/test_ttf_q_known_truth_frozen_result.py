import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_known_truth_frozen_receipt_matches_pre_result_contract():
    contract = json.loads(
        (ROOT / "docs/supporting/ttf_q_known_truth_benchmark_v0.1.json").read_text()
    )
    result = json.loads(
        (ROOT / "benchmarks/frozen/ttf_q_known_truth_decomposition_v0.1.json").read_text()
    )

    assert contract["status"] == "FROZEN_BEFORE_BENCHMARK_RESULTS"
    assert result["status"] == "FROZEN_COMPLETE_KNOWN_TRUTH_INFORMATION_DECOMPOSITION"
    assert result["decomposition_grid"]["scenarios"] == 243
    assert (
        result["exact_recovery"]["maximum_absolute_error_across_all_nested_information_quantities"]
        < result["exact_recovery"]["tolerance_declared_before_result"]
    )
    assert result["exact_recovery"]["all_scenarios_within_tolerance"] is True
    assert result["signal_concentration"]["source_concentration_detection_fraction"] == 1.0
    assert result["signal_concentration"]["target_concentration_detection_fraction"] == 1.0
    assert result["signal_concentration"]["median_source_effective_count_ratio_concentrated_vs_broad"] < 0.30
    assert result["signal_concentration"]["median_target_effective_count_ratio_concentrated_vs_broad"] < 0.30
    assert result["boundary"] == contract["boundary"]
