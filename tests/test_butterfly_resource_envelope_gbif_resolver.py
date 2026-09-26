from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "acquire_butterfly_resource_envelope_occurrences.py"
SPEC = importlib.util.spec_from_file_location("butterfly_resource_occurrence_acquire", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
acquire = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(acquire)


def test_resolver_keeps_exact_strict_species(monkeypatch):
    calls = []

    def fake_get_json(path, params, *, deadline, request_seconds=30.0, attempts=3):
        calls.append((path, params))
        assert path == "species/match"
        return {
            "usageKey": 123,
            "canonicalName": "Papilio machaon",
            "rank": "SPECIES",
            "matchType": "EXACT",
        }

    monkeypatch.setattr(acquire, "get_json", fake_get_json)
    result = acquire.resolve_species_metadata("Papilio machaon", deadline=999.0, request_seconds=30.0)
    assert result["status"] == "MATCHED"
    assert result["usage_key"] == 123
    assert result["resolver_route"] == "species_match_strict"
    assert len(calls) == 1


def test_resolver_recovers_unique_exact_accepted_species_from_search(monkeypatch):
    def fake_get_json(path, params, *, deadline, request_seconds=30.0, attempts=3):
        if path == "species/match":
            return {
                "usageKey": 0,
                "canonicalName": None,
                "rank": None,
                "matchType": "NONE",
            }
        assert path == "species/search"
        return {
            "results": [
                {
                    "key": 456,
                    "canonicalName": "Graphium sarpedon",
                    "rank": "SPECIES",
                    "taxonomicStatus": "ACCEPTED",
                },
                {
                    "key": 999,
                    "canonicalName": "Graphium sarpedonides",
                    "rank": "SPECIES",
                    "taxonomicStatus": "ACCEPTED",
                },
            ]
        }

    monkeypatch.setattr(acquire, "get_json", fake_get_json)
    result = acquire.resolve_species_metadata("Graphium sarpedon", deadline=999.0, request_seconds=30.0)
    assert result["status"] == "MATCHED"
    assert result["usage_key"] == 456
    assert result["match_type"] == "EXACT_SEARCH_FALLBACK"
    assert result["resolver_route"] == "species_search_exact_accepted_unique"


def test_resolver_rejects_ambiguous_exact_accepted_search(monkeypatch):
    def fake_get_json(path, params, *, deadline, request_seconds=30.0, attempts=3):
        if path == "species/match":
            return {
                "usageKey": 0,
                "canonicalName": None,
                "rank": None,
                "matchType": "NONE",
            }
        return {
            "results": [
                {
                    "key": 1,
                    "canonicalName": "Alpha beta",
                    "rank": "SPECIES",
                    "taxonomicStatus": "ACCEPTED",
                },
                {
                    "key": 2,
                    "canonicalName": "Alpha beta",
                    "rank": "SPECIES",
                    "taxonomicStatus": "ACCEPTED",
                },
            ]
        }

    monkeypatch.setattr(acquire, "get_json", fake_get_json)
    result = acquire.resolve_species_metadata("Alpha beta", deadline=999.0, request_seconds=30.0)
    assert result["status"] == "REJECTED_GBIF_TAXON_MATCH"
    assert result["exact_accepted_search_candidates"] == 2
