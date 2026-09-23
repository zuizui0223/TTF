import importlib.util
from pathlib import Path


SCRIPT = Path("scripts/acquire_relational_environment_occurrences.py")
COMPLETE_SCRIPT = Path("scripts/complete_relational_environment_transport.py")


def load_complete():
    spec = importlib.util.spec_from_file_location("envcomplete", COMPLETE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load():
    spec = importlib.util.spec_from_file_location("envocc", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_environment_acquisition_uses_current_gbif_search_parameter_names(monkeypatch):
    m = load()
    calls = []

    def fake_get_json(path, params, retries=8):
        calls.append((path, dict(params)))
        if path == "species/match":
            return {
                "usageKey": 1,
                "matchType": "EXACT",
                "rank": "SPECIES",
                "canonicalName": "Test species",
            }
        if path == "occurrence/search" and params.get("limit") == 1:
            return {"count": 1}
        if path == "occurrence/search":
            return {
                "results": [{
                    "key": 123,
                    "decimalLatitude": 35.0,
                    "decimalLongitude": 139.0,
                }]
            }
        raise AssertionError(path)

    monkeypatch.setattr(m, "get_json", fake_get_json)
    m.fetch_species("Test species")

    occurrence_calls = [params for path, params in calls if path == "occurrence/search"]
    assert occurrence_calls
    for params in occurrence_calls:
        assert params["taxonKey"] == 1
        assert params["hasCoordinate"] == "true"
        assert params["hasGeospatialIssue"] == "false"
        assert params["occurrenceStatus"] == "PRESENT"
        assert params["year"] == "2010,2026"
        assert "taxon_key" not in params
        assert "has_coordinate" not in params
        assert "has_geospatial_issue" not in params
        assert "occurrence_status" not in params


def test_frozen_transport_retry_uses_cutover_parallelism():
    m = load_complete()
    assert m.RETRY_WORKERS == 20
    source = COMPLETE_SCRIPT.read_text()
    assert "ThreadPoolExecutor(max_workers=RETRY_WORKERS)" in source
    assert "for name in pending:" in source
