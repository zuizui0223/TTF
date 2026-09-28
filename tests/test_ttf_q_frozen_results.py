import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path: str):
    return json.loads((ROOT / path).read_text())


def test_ttf_q_frozen_results_are_response_blind_and_v03_stays_closed():
    final_state = load("benchmarks/frozen/relational_program_final_state_v0.3.json")
    c1 = load("benchmarks/frozen/ttf_q_c1_result_receipt_v0.1.json")
    c2 = load("benchmarks/frozen/ttf_q_c2_result_receipt_v0.1.json")
    ladder = load("benchmarks/frozen/ttf_q_b1_c1_information_ladder_v0.1.json")
    detect = load("benchmarks/frozen/ttf_q_bc_detectability_result_v0.1.json")

    assert final_state["final_state"] == "CLOSED_NO_EVALUABLE_TEST"

    for payload in (c1, c2, ladder, detect):
        boundary = payload.get("boundary", {})
        if "relational_v03_final_state_unchanged" in boundary:
            assert boundary["relational_v03_final_state_unchanged"] == final_state["final_state"]
        assert boundary.get("relational_v03_reopened", False) is False
        assert boundary.get("genetic_response_used", False) is False

    assert c1["design"]["admissible_species"] == 653
    assert c1["development"]["total_unique_variance_fraction"] == 0.3647140731871479
    assert c1["confirmatory"]["total_unique_variance_fraction"] == 0.3812840593549338

    assert c2["development"]["geography_retained_fraction"] == 0.819650642967475
    assert c2["confirmatory"]["geography_retained_fraction"] == 0.8295950394914621
    assert 0 < c2["development"]["geography_attributable_fraction_of_C1_unique_information"] < 0.25
    assert 0 < c2["confirmatory"]["geography_attributable_fraction_of_C1_unique_information"] < 0.25

    assert ladder["development"]["B1_current_climate"]["total_unique_variance_fraction"] > 0.40
    assert ladder["confirmatory"]["B1_current_climate"]["total_unique_variance_fraction"] > 0.44
    assert ladder["development"]["C1_history_given_current"]["total_unique_variance_fraction"] > 0.36
    assert ladder["confirmatory"]["C1_history_given_current"]["total_unique_variance_fraction"] > 0.38

    assert detect["B_current_climate"]["minimum_detectable_effect_for_80pct_power_by_private_amplitude"] == {
        "A0": 0.02,
        "A1": 0.02,
        "A2": 0.05,
        "A3": 0.05,
    }
    assert detect["C_historical_given_current"]["minimum_detectable_effect_for_80pct_power_by_private_amplitude"] == {
        "A0": 0.02,
        "A1": 0.03,
        "A2": 0.05,
        "A3": 0.08,
    }


def test_ttf_q_information_is_broadly_distributed_not_single_endpoint_dominated():
    c2 = load("benchmarks/frozen/ttf_q_c2_result_receipt_v0.1.json")

    for panel in ("development", "confirmatory"):
        result = c2[panel]
        assert result["signal_effective_sources_after_geography"] > 50
        assert result["signal_effective_targets_after_geography"] > 50
        assert result["max_source_signal_share_after_geography"] < 0.05
        assert result["max_target_signal_share_after_geography"] < 0.05
