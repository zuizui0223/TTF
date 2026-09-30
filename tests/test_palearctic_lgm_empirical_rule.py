import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_palearctic_lgm_v03_empirical_rule_is_frozen_before_outcome():
    p=json.loads(
        (ROOT/"docs/supporting/genetic_palearctic_lgm_empirical_opening_rule_v0.3.json").read_text()
    )
    assert p["status"]=="FROZEN_BEFORE_V03_R_LGM_TTF_Q_RESULT_OR_SUBPANEL_GENETIC_RESPONSE"
    assert p["opening_condition"]["required_ttf_q_status"]=="PASS_TO_ONE_SHOT_V03_SUBPANEL_GENETIC_RESPONSE"
    assert p["opening_condition"]["genetic_identity_may_be_opened_before_pass"] is False
    assert p["opening_condition"]["one_shot_only"] is True

    g=p["empirical_species_geometry"]
    assert g["neighbor_fraction"]==0.15
    assert g["minimum_palearctic_localities"]==12
    assert g["minimum_endpoint_disjoint_training_edges"]==5
    assert g["no_backfill"] is True
    assert g["no_role_resplit"] is True
    assert "rerun the frozen v0.3 TTF-Q" in g["survivor_rule"]

    r=p["genetic_response"]
    assert r["minimum_comparable_fraction"]==0.5
    assert r["minimum_ibd_training_edges"]==5
    assert r["kernel"]=={
        "bandwidth_km":500,
        "prior_strength":0.25,
        "prior_mean":0,
        "segment_points":5,
        "training_species_per_field":1,
    }

    m=p["relational_model"]
    assert m["predictors_in_order"]==[
        "z_R_LGM",
        "z_R_present",
        "z_geographic_coverage",
        "centered_same_order",
        "centered_same_family",
        "z_abs_log_retained_occurrence_count_ratio",
    ]
    assert m["source_fixed_effect"] is True
    assert m["target_fixed_effect"] is True
    assert m["primary_index"]==0
    assert m["alpha_one_sided"]==0.025
    assert "two-way source/target cluster-robust" in m["covariance"]

    assert p["serialization"]["sequence_identity"] is False
    assert p["serialization"]["edge_level_genetic_distance_vectors"] is False
    assert all(v is False for v in p["response_firewall"].values())
