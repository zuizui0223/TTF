import csv
import importlib.util
import json
from pathlib import Path


SCRIPT = Path("scripts/freeze_relational_historical_candidate_census.py")


def load():
    spec = importlib.util.spec_from_file_location("historycensus", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_history_species_hash_is_source_bound():
    m = load()
    a = m.species_hash("a" * 64, "Alpha beta")
    assert a == m.species_hash("a" * 64, "Alpha beta")
    assert a != m.species_hash("b" * 64, "Alpha beta")


def test_verified_prior_exclusions_union_three_frozen_universes(tmp_path):
    m = load()

    s1_names = ["Alpha beta", "Beta gamma"]
    s2_names = ["Gamma delta", "Delta epsilon"]
    b_names = ["Delta epsilon", "Epsilon zeta"]

    b = tmp_path / "b.csv"
    with b.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["species"])
        writer.writeheader()
        writer.writerows([{"species": x} for x in b_names])

    manifest = tmp_path / "prior.json"
    manifest.write_text(json.dumps({
        "schema": "ttf_relational_historical_prior_exclusion_manifest_v0.1",
        "S1_butterfly_trait": {
            "species": s1_names,
            "source_full_design_sha256_verified": "1" * 64,
            "full_design_sha_match": True,
        },
        "S2_host_resource_geography": {
            "species": s2_names,
            "source_full_census_sha256_verified": "2" * 64,
            "selected_species_sha256_verified": m.species_digest(s2_names),
            "selected_species_digest_match": True,
        },
        "B_environment_candidate_universe": {
            "species": b_names,
        },
        "membership_checks": {
            "union_species": 5,
        },
        "response_firewall": {
            "prior_genetic_values_read": False,
            "Study_B_genetic_response_used": False,
            "Study_C_sequence_identity_opened": False,
            "Study_C_T_st_opened": False,
        },
    }))

    history = {
        "independent_species_domain": {
            "prior_universe_reproduction": {
                "S1_butterfly_trait": {
                    "expected_design_sha256": "1" * 64,
                    "expected_species": 2,
                },
                "S2_host_resource_geography": {
                    "expected_full_census_sha256": "2" * 64,
                    "expected_selected_species_sha256": m.species_digest(s2_names),
                    "expected_species": 2,
                },
                "B_environment_candidate_universe": {
                    "expected_sha256": m.sha256_path(b),
                    "expected_species": 2,
                },
            }
        }
    }

    union, receipt = m.verified_prior_exclusions(
        history,
        manifest_path=manifest,
        b_candidates_path=b,
    )
    assert union == {
        "Alpha beta", "Beta gamma", "Gamma delta",
        "Delta epsilon", "Epsilon zeta",
    }
    assert receipt["S1_species"] == 2
    assert receipt["S2_species"] == 2
    assert receipt["B_species"] == 2
    assert receipt["union_species"] == 5
