import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_palearctic_insect_lgm_contract_and_parent_are_frozen_response_blind():
    contract = json.loads(
        (ROOT / "docs/supporting/palearctic_insect_lgm_refugia_v0.1.json").read_text()
    )
    candidates = json.loads(
        (ROOT / "benchmarks/frozen/palearctic_insect_lgm_taxonomic_candidates_v0.1.json").read_text()
    )
    prescreen = json.loads(
        (ROOT / "benchmarks/frozen/palearctic_insect_lgm_response_blind_prescreen_v0.1.json").read_text()
    )

    assert contract["status"].startswith("FROZEN_RESPONSE_BLIND")
    assert candidates["status"].startswith("FROZEN_RESPONSE_BLIND")
    assert len(candidates["species"]) == 143
    assert len(set(candidates["species"])) == 143
    assert candidates["source"]["parent_survivors"] == 211
    assert candidates["response_firewall"] == {
        "species_level_Phase4_scores_used": False,
        "pairwise_genetic_distances_used": False,
        "subpanel_T_st_computed": False,
        "subpanel_beta_LGM_computed": False,
    }
    assert prescreen["response_blind_counts"]["Insecta_among_211"] == 183
    assert prescreen["response_blind_counts"]["coarse_proxy_0p80_terrestrial_core_orders"] == 47
    assert prescreen["genetic_response_used"] is False
    assert contract["response_blind_panel_eligibility"]["palearctic_realm"][
        "minimum_fraction_retained_occurrences_in_realm"
    ] == 0.80


def test_palearctic_role_assignment_is_deterministic_and_response_blind():
    module = load_script(
        ROOT / "scripts/classify_palearctic_insect_lgm_panel.py",
        "palearctic_panel",
    )
    species = ["Gamma species", "Alpha species", "Beta species", "Delta species"]
    first = module.role_order(species)
    second = module.role_order(list(reversed(species)))
    assert first == second
    assert set(first) == set(species)


def test_palearctic_gbif_acquisition_reuses_frozen_search_contract(monkeypatch):
    module = load_script(
        ROOT / "scripts/acquire_palearctic_insect_lgm_occurrences.py",
        "palearctic_acquire",
    )
    calls = []

    def fake_get_json(path, params, retries=8):
        calls.append((path, dict(params)))
        if path == "species/match":
            return {
                "usageKey": 1,
                "matchType": "EXACT",
                "rank": "SPECIES",
                "canonicalName": "Test species",
                "kingdom": "Animalia",
                "phylum": "Arthropoda",
                "class": "Insecta",
                "order": "Lepidoptera",
                "family": "Testidae",
                "genus": "Test",
            }
        if path == "occurrence/search" and params.get("limit") == 1:
            return {"count": 1}
        if path == "occurrence/search":
            return {
                "results": [{
                    "key": 123,
                    "decimalLatitude": 50.0,
                    "decimalLongitude": 10.0,
                }]
            }
        raise AssertionError(path)

    monkeypatch.setattr(module, "get_json", fake_get_json)
    retained, ledger = module.fetch_species("Test species")
    assert ledger["taxonomy"]["class"] == "Insecta"
    assert ledger["taxonomy"]["order"] == "Lepidoptera"
    assert len(retained) == 1

    occurrence_calls = [params for path, params in calls if path == "occurrence/search"]
    assert occurrence_calls
    for params in occurrence_calls:
        assert params["taxonKey"] == 1
        assert params["hasCoordinate"] == "true"
        assert params["hasGeospatialIssue"] == "false"
        assert params["occurrenceStatus"] == "PRESENT"
        assert params["year"] == "2010,2026"
