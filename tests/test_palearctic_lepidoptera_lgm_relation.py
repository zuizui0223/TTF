from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def test_lgm_relation_builder_has_no_genetic_response_dependency():
    source=(ROOT/"scripts/build_palearctic_lepidoptera_lgm_relation.py").read_text()
    assert "species_scores" not in source
    assert "genetic_empirical" not in source
    assert "phase4_empirical" not in source
    assert "R_LGM" in source
    assert "R_current" in source
    assert "geographic_coverage" in source

def test_lgm_rule_fixes_model_before_response():
    import json
    rule=json.loads((ROOT/"docs/supporting/palearctic_lepidoptera_lgm_sdm_rule_v0.1.json").read_text())
    assert rule["status"]=="FROZEN_BEFORE_SUBPANEL_GENETIC_RESPONSE_READ"
    assert rule["sdm"]["ridge_logistic_C"]==1.0
    assert rule["sdm"]["maximum_background_cells"]==2000
    assert rule["sdm"]["minimum_background_cells"]==100
    assert rule["pairwise_primary_relation"]["grid_resolution_degrees"]==0.5
