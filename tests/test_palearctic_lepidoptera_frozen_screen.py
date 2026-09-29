import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def test_palearctic_lepidoptera_screen_is_response_blind_and_large_enough():
    r=json.loads((ROOT/"benchmarks/frozen/palearctic_lepidoptera_response_blind_screen_v0.1.json").read_text())
    assert r["status"]=="PASS_TO_LGM_PREDICTOR_CONSTRUCTION"
    assert r["census"]["survivors_total"]==211
    assert r["census"]["insecta"]==183
    assert r["census"]["lepidoptera"]==90
    assert r["census"]["primary_strict_species"]==26
    assert r["census"]["primary_ordered_nonself_dyads"]==650
    assert r["census"]["primary_families"]==13
    assert len(r["primary_species"])==26
    assert all(x["core_pal_frac_30"]==1.0 for x in r["primary_species"])
    assert all(v is False for v in r["response_firewall"].values())

def test_palearctic_screen_contract_precedes_any_subpanel_response():
    rule=json.loads((ROOT/"docs/supporting/palearctic_lepidoptera_screen_rule_v0.1.json").read_text())
    assert rule["status"]=="FROZEN_BEFORE_SUBPANEL_GENETIC_RESPONSE_READ"
    assert rule["primary_screen"]["minimum_species_to_continue"]==20
    assert rule["response_firewall"]["subpanel_transfer_response_used"] is False
