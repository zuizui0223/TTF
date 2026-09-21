import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load(path):
    return json.loads((ROOT / path).read_text())


def test_v03_family_is_exactly_two_fixed_slots():
    program = load("docs/supporting/relational_program_v0.3.json")
    family = program["prospective_family"]
    assert family["familywise_alpha"] == 0.05
    assert family["maximum_confirmatory_openings"] == 2
    assert family["slots_fixed_before_first_v03_qualification"] == ["B", "C"]
    assert family["allocation"] == {
        "B_environmental_niche_similarity": 0.025,
        "C_historical_climate_exposure_similarity": 0.025,
    }
    assert family["no_alpha_recycling"] is True
    assert family["no_replacement_slots"] is True
    assert program["descriptive_only"]["conditional_order_contrast_v0_3"]["inferential_alpha"] == 0


def test_study_b_chain_is_bound_to_alpha_0p025_and_v03_freshness():
    program = load("docs/supporting/relational_program_v0.3.json")
    relation = load("docs/supporting/relational_environment_relation_rule_v0.3.json")
    opportunity = load("docs/supporting/relational_environment_opportunity_rule_v0.2.json")
    qualification = load("docs/supporting/relational_environment_qualification_rule_v0.2.json")
    freshness = load("docs/supporting/relational_future_family_freshness_amendment_v0.1.json")

    assert program["slot_B"]["execution_chain"] == [
        "docs/supporting/relational_environment_relation_rule_v0.3.json",
        "docs/supporting/relational_environment_opportunity_rule_v0.2.json",
        "docs/supporting/relational_environment_qualification_rule_v0.2.json",
    ]
    assert relation["parent_program"] == "docs/supporting/relational_program_v0.3.json"
    assert opportunity["parent_relation_rule"] == "docs/supporting/relational_environment_relation_rule_v0.3.json"
    assert qualification["parent_opportunity_rule"] == "docs/supporting/relational_environment_opportunity_rule_v0.2.json"
    assert qualification["inference"]["alpha"] == 0.025
    assert qualification["gate_numeric"] == {
        "p_value_cutoff": 0.025,
        "private_type1_wilson95_upper_max": 0.05,
        "relational_power_wilson95_lower_min": 0.8,
    }
    assert freshness["exact_overlap_with_study_B_candidates"]["union_overlap_species"] == 70
    assert freshness["execution_rule"]["backfill"] is False
    assert freshness["execution_rule"]["replacement_candidates"] is False


def test_study_c_is_frozen_before_b_and_disjoint_from_entire_b_candidate_universe():
    program = load("docs/supporting/relational_program_v0.3.json")
    history = load("docs/supporting/relational_historical_climate_exposure_rule_v0.1.json")
    freshness = load("docs/supporting/relational_future_family_freshness_amendment_v0.1.json")

    assert program["slot_C"]["alpha_one_sided"] == 0.025
    assert program["slot_C"]["status_at_freeze"] == (
        "FULL_RELATION_PANEL_AND_QUALIFICATION_CONTRACT_FROZEN_BEFORE_B_QUALIFICATION"
    )
    assert history["alpha_one_sided"] == 0.025
    assert history["synthetic_qualification"]["alpha"] == 0.025
    assets = history["historical_climate"]["asset_contract"]
    assert assets["version"] == "1.0"
    assert assets["time_index"] == {"LGM_21ka_BP": -190, "present_0_BP": 20}
    assert len(assets["required_basenames"]) == 8
    assert set(name.split("_")[2] for name in assets["required_basenames"]) == {
        "bio01", "bio07", "bio12", "bio15"
    }
    assert "entire 1000-species Study-B candidate universe" in freshness["study_C_reservation"]["rule"]
    assert history["independent_species_domain"]["exclude_entire_study_B_candidate_universe"] == (
        "benchmarks/frozen/relational_fresh_candidate_species_v0.1.csv"
    )


def test_terminal_states_distinguish_two_nulls_from_not_evaluable():
    program = load("docs/supporting/relational_program_v0.3.json")
    states = program["terminal_states"]
    assert "CLOSED_TWO_NULLS" in states
    assert "CLOSED_PARTIAL_NOT_EVALUABLE" in states
    assert "CLOSED_NO_EVALUABLE_TEST" in states
    assert "both evaluable qualified nulls" in states["CLOSED_TWO_NULLS"]
    assert "not a two-null biological conclusion" in states["CLOSED_PARTIAL_NOT_EVALUABLE"]
