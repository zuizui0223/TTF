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
    s1 = tmp_path / "s1.json"
    s1.write_text(json.dumps({
        "schema": "ttf_relational_prior_S1_species_exclusion_v0.1",
        "source_design_sha256": "1" * 64,
        "species": s1_names,
        "species_count": 2,
        "species_list_sha256": m.species_digest(s1_names),
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }))

    s2_names = ["Gamma delta", "Delta epsilon"]
    s2 = tmp_path / "s2.json"
    s2.write_text(json.dumps({
        "schema": "ttf_relational_prior_S2_species_exclusion_v0.1",
        "source_full_census_sha256": "2" * 64,
        "species": s2_names,
        "species_count": 2,
        "species_list_sha256": m.species_digest(s2_names),
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }))

    b = tmp_path / "b.csv"
    with b.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["species"])
        writer.writeheader()
        writer.writerows([
            {"species": "Delta epsilon"},
            {"species": "Epsilon zeta"},
        ])

    history = {
        "independent_species_domain": {
            "prior_universe_reproduction": {
                "S1_butterfly_trait": {
                    "expected_design_sha256": "1" * 64,
                    "expected_species": 2,
                    "expected_species_list_sha256": m.species_digest(s1_names),
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
        s1_manifest_path=s1,
        s2_manifest_path=s2,
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
