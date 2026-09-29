import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_palearctic_lgm_route_is_frozen_before_subpanel_response():
    rule=json.loads(
        (ROOT/"docs/supporting/genetic_palearctic_lgm_subpanel_rule_v0.1.json").read_text()
    )
    assert rule["status"].startswith("FROZEN_")
    chronology=rule["chronology_and_claim_boundary"]
    assert chronology["global_211_species_phase4_result_already_opened"] is True
    assert chronology["independent_confirmatory_claim_for_this_same_archive"] is False
    assert chronology["subgroup_specific_pairwise_transfer_response_inspected_for_this_route"] is False
    assert chronology["species_level_genetic_scores_inspected_for_subgroup_selection"] is False

    filt=rule["authoritative_subpanel_filter"]
    assert filt["minimum_species"]==40
    assert filt["minimum_source_clusters"]==20
    assert filt["minimum_target_clusters"]==20
    assert filt["palearctic"]["minimum_fraction_of_selected_localities_in_palearctic"]==0.80

    lgm=rule["lgm_climatic_suitability_model"]
    assert lgm["primary_relation"].startswith("R_LGM")
    assert lgm["nuisance_current_relation"].startswith("R_present")
    assert lgm["thresholded_refugial_components_primary"] is False

    opening=rule["ttf_q_pre_response_characterization"]["opening_rule"]
    assert opening["minimum_total_unique_variance_fraction"]==0.15
    assert opening["maximum_evaluable_grid_mde_at_A2"]==0.10
