import json
from pathlib import Path


RULE = Path("benchmarks/frozen/relational_environment_rate_limit_recovery_v0.1.json")
CORE = Path("scripts/acquire_relational_environment_occurrences.py")
RECOVERY = Path("scripts/recover_relational_environment_rate_limit.py")


def test_rate_limit_recovery_is_http_only_and_pre_result():
    p = json.loads(RULE.read_text())
    assert p["schema"] == "ttf_relational_environment_rate_limit_recovery_v0.1"
    assert p["status"] == "FROZEN_HTTP_ONLY_RECOVERY_BEFORE_RELATION_RESULT"
    assert p["source_request_error_species"] == [
        "Pyrrhosoma nymphula",
        "Lithobates clamitans",
        "Diarsia rubi",
        "Mythimna impura",
    ]
    assert p["source_failure_class"] == "GBIF HTTP 429 REQUEST_ERROR only"
    assert p["relation_result_seen"] is False
    assert p["genetic_response_used"] is False
    assert all(v is False for v in p["response_firewall"].values())


def test_rate_limit_recovery_preserves_frozen_query_and_offsets():
    p = json.loads(RULE.read_text())
    contract = p["unchanged_scientific_contract"]
    assert contract == {
        "exact_taxon_match": True,
        "year": "2010,2026",
        "hasCoordinate": True,
        "hasGeospatialIssue": False,
        "occurrenceStatus": "PRESENT",
        "deterministic_page_offsets": True,
        "page_size": 300,
        "maximum_pages": 20,
        "exact_coordinate_dedup": True,
        "thinning_km": 10.0,
        "maximum_retained": 200,
        "minimum_occurrences": 30,
    }
    core = CORE.read_text()
    assert '"year": "2010,2026"' in core
    assert '"hasCoordinate": "true"' in core
    assert '"hasGeospatialIssue": "false"' in core
    assert '"occurrenceStatus": "PRESENT"' in core
    assert "deterministic_page_offsets(total, page_size=300, maximum_pages=20)" in core
    assert "minimum_distance_km=10.0" in core
    assert "maximum_retained=200" in core


def test_recovery_changes_transport_function_only():
    text = RECOVERY.read_text()
    assert "acquisition.get_json = slow_get_json" in text
    assert "acquisition.fetch_species(species)" in text
    assert "time.sleep(2.0)" in text
    assert "time.sleep(float(rule[\"transport_only_changes\"][\"initial_cooldown_seconds\"]))" in text
