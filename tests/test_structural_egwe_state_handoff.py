from __future__ import annotations

import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
C=ROOT/"docs/supporting/structural_egwe_state_handoff_v0.1.json"

def load():
    return json.loads(C.read_text())

def test_predictor_increment_is_R_plus_C_over_R():
    x=load()
    p=x["predictor_spaces"]["primary_increment"]
    assert p["name"]=="C|R"
    assert p["larger_space"]=="RC"
    assert p["baseline_space"]=="R"
    assert p["estimand"]=="T_RC - T_R"

def test_species_cannot_be_selected_on_favorable_structural_result():
    x=load()
    assert "select species because within-species C-minus-R was favorable" in x["eligibility"]["forbidden"]

def test_burned_pilot_never_enters_transfer_denominator():
    x=load()
    assert x["burned_pilot_rule"]["enters_transfer_denominator"] is False

def test_genetic_result_is_not_structural_replication():
    x=load()
    assert x["current_boundaries"]["counts_as_Structural_replication"] is False
    assert x["current_boundaries"]["PNW_RMNP_enter_transfer_denominator"] is False
