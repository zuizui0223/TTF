import hashlib
import json
import runpy
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/screen_palearctic_holometabola_panel.py"
NS=runpy.run_path(str(SCRIPT))


def test_palearctic_lgm_contract_is_frozen_before_pairwise_response():
    contract=json.loads(
        (ROOT/"docs/supporting/genetic_palearctic_holometabola_lgm_program_v0.1.json").read_text()
    )
    assert contract["status"]=="FROZEN_BEFORE_PAIRWISE_GENETIC_RESPONSE_CONSTRUCTION"
    assert contract["provenance_boundary"]["species_pair_genetic_concordance_opened_for_this_program"] is False
    assert contract["panel_rule"]["taxonomy"]["allowed_orders"] == [
        "Coleoptera","Hymenoptera","Lepidoptera"
    ]
    assert contract["panel_rule"]["palearctic_membership"]["minimum_fraction_inside_palearctic"]==0.8
    assert contract["panel_rule"]["minimum_species_to_continue"]==30
    assert contract["ttf_q_qualification"]["must_precede_pairwise_genetic_response_construction"] is True


def test_deterministic_split_is_species_disjoint_and_balanced():
    species=[f"Species {i:02d}" for i in range(48)]
    a,b=NS["deterministic_split"](species,"palearctic-holometabola-lgm-v0.1")
    assert len(a)==24 and len(b)==24
    assert not set(a)&set(b)
    assert set(a)|set(b)==set(species)


def test_feasibility_receipt_is_response_blind_and_non_authoritative():
    receipt=json.loads(
        (ROOT/"benchmarks/frozen/genetic_palearctic_holometabola_feasibility_prescreen_v0.1.json").read_text()
    )
    assert receipt["status"]=="DESCRIPTIVE_RESPONSE_BLIND_FEASIBILITY_ONLY"
    assert receipt["firewall"]["species_level_phase4_scores_read"] is False
    assert receipt["firewall"]["pairwise_genetic_concordance_constructed"] is False
    assert receipt["provisional_core_palearctic_screen"]["authority"]=="NON_AUTHORITATIVE_FEASIBILITY_ONLY"
    assert receipt["provisional_core_palearctic_screen"]["eligible_species"]==48
