import ast
import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_palearctic_occurrence_rule_is_frozen_and_response_blind():
    rule=json.loads(
        (ROOT/"docs/supporting/genetic_palearctic_holometabola_occurrence_acquisition_rule_v0.1.json").read_text()
    )
    assert rule["status"]=="FROZEN_BEFORE_FULL_OCCURRENCE_ASSET"
    assert rule["gbif"]["year"]=="1990,2026"
    assert rule["gbif"]["maximum_pages"]==20
    assert rule["retention"]["minimum_distance_km"]==10.0
    assert rule["retention"]["maximum_retained_per_species"]==1000
    assert rule["retention"]["minimum_required_per_species"]==50
    assert rule["firewall"]["pairwise_genetic_concordance_opened"] is False


def test_palearctic_occurrence_acquirer_uses_camel_case_gbif_search_params():
    source=(ROOT/"scripts/acquire_palearctic_holometabola_occurrences.py").read_text()
    ast.parse(source)
    for token in ("taxonKey","hasCoordinate","hasGeospatialIssue","occurrenceStatus"):
        assert token in source
    assert "species_scores" not in source
