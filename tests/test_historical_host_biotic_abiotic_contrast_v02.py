"""Contract boundary: biotic memory requires an insect-climate counterfactual."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"docs/supporting/genetic_historical_host_biotic_abiotic_contrast_v0.2.json"

def test_biotic_abiotic_contract_cannot_claim_positive_from_host_only():
    p=json.loads(CONTRACT.read_text())
    assert p["status"]=="SCIENTIFIC_SUCCESSOR_FROZEN_BEFORE_ANY_LGM_EXPOSURE_OR_GENETIC_RESPONSE"
    assert p["phase_status"].startswith("DESIGN_ONLY")
    controls=p["exposure"]["mandatory_controls"]
    assert "LGM insect abiotic edge resistance" in controls
    assert "present insect abiotic edge resistance" in controls
    assert "present host edge resistance" in controls
    assert p["pre_response_requirements"]["must_rerun_synthetic_qualification_for_this_new_multivariable_estimand"] is True
    assert p["pre_response_requirements"]["v0_1_synthetic_qualification_not_transferable"] is True
    assert p["pre_response_requirements"]["no_interaction_specific_positive_claim_without_all_controls"] is True
    assert p["panel"]["candidate_species"]==135
    assert all(value is False for value in p["response_firewall"].values())
