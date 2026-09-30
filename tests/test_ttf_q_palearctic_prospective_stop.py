import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_ttf_q_palearctic_prospective_stop_is_response_blind():
    case=json.loads(
        (ROOT/"benchmarks/frozen/ttf_q_palearctic_prospective_stop_case_v0.1.json").read_text()
    )
    authority=json.loads(
        (ROOT/"benchmarks/frozen/genetic_palearctic_lgm_v03_qualification_authority_correction_v0.1.json").read_text()
    )

    assert case["status"]=="RESPONSE_BLIND_STOP_BEFORE_SUBPANEL_GENETIC_RESPONSE"
    assert case["realized_external_design"]["admissible_species"]==21
    assert case["realized_external_design"]["source_clusters"]==10
    assert case["realized_external_design"]["target_clusters"]==11
    assert case["realized_external_design"]["directed_dyads"]==110

    det=case["authoritative_deterministic_ttf_q"]
    assert det["unique_information"]["pass"] is False
    assert det["source_signal_breadth"]["pass"] is False
    assert det["target_signal_breadth"]["pass"] is False
    assert det["overall_pass"] is False
    assert det["decision"]=="NOT_EVALUABLE_NONLEPIDOPTERA_PALEARCTIC_LGM_RELATION"

    assert authority["authoritative_decision"]["status"]==det["decision"]
    assert case["diagnostic_full_v03_synthetic_run"]["authority_for_stop"] is False

    firewall=case["biological_response_firewall"]
    assert firewall["subpanel_nucleotide_identity_opened"] is False
    assert firewall["pairwise_T_st_computed"] is False
    assert firewall["beta_LGM_computed"] is False
    assert firewall["genetic_response_authorized"] is False

    assert case["interpretation"]["biological_null_claim_allowed"] is False
