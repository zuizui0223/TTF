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

    s1 = tmp_path / "s1.json"
    s1.write_text(json.dumps({"species_order": ["Alpha beta", "Beta gamma"]}))

    s2_names = ["Gamma delta", "Delta epsilon"]
    s2 = tmp_path / "s2.json"
    s2.write_text(json.dumps({"species": [{"species": x} for x in s2_names]}))

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
                    "expected_design_sha256": m.sha256_path(s1),
                    "expected_species": 2,
                },
                "S2_host_resource_geography": {
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
        s1_design_path=s1,
        s2_census_path=s2,
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
