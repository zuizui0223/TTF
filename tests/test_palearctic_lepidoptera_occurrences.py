from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def test_palearctic_occurrence_wrapper_reuses_frozen_fetch_species():
    source=(ROOT/"scripts/acquire_palearctic_lepidoptera_occurrences.py").read_text()
    assert 'scripts/acquire_relational_environment_occurrences.py' in source
    assert 'fetch_species' in source
    assert 'species_scores' not in source
    assert 'genetic_empirical' not in source

def test_palearctic_occurrence_aggregator_requires_20_species_before_lgm():
    source=(ROOT/"docs/supporting/palearctic_lepidoptera_lgm_sdm_rule_v0.1.json").read_text()
    assert '"minimum_predictor_admissible_species": 20' in source
