from __future__ import annotations

import hashlib
import json
from pathlib import Path


def _load(path: str) -> dict:
    return json.loads(Path(path).read_text())


def test_decker_exclusion_manifest_is_exactly_221_species() -> None:
    payload = _load("benchmarks/frozen/genetic_decker_species_exclusion_v0.1.json")
    species = list(payload["species"])
    assert payload["schema"] == "ttf_genetic_decker_species_exclusion_v0.1"
    assert len(species) == payload["species_count"] == 221
    assert species == sorted(species)
    digest = hashlib.sha256(("\n".join(species) + "\n").encode("utf-8")).hexdigest()
    assert digest == payload["species_list_sha256"] == "1f9e1e17f936e2cf21b85b8c5157c01459b69f785c500666ae19dbd352530319"
    assert all(value is False for value in payload["outcome_firewall"].values())


def test_fresh_phase1_numeric_contract_is_closed_before_download() -> None:
    execution = _load("docs/supporting/genetic_phylogatr_phase1_execution_rule_v0.1.json")
    protocol = _load("docs/supporting/genetic_phylogatr_confirmatory_protocol_v0.1.json")
    parser = _load("docs/supporting/genetic_phylogatr_phase1_parser_rule_v0.1.json")
    digest = _load("docs/supporting/genetic_phylogatr_phase1_digest_rule_v0.1.json")
    assert execution["schema"] == "ttf_genetic_phylogatr_phase1_execution_rule_v0.1"
    assert execution["phase1"] == {
        "minimum_unique_localities": 12,
        "neighbor_fraction": 0.15,
        "minimum_endpoint_disjoint_ibd_training_edges": 5,
        "maximum_species": 250,
        "minimum_panel_size": 150,
        "train_fraction_semantics": "lower floor(n/2) by frozen split hash is training; remainder evaluation",
        "support_census_bandwidth_km": 500.0,
        "support_census_is_descriptive_only": True,
    }
    assert all(value is False for value in execution["outcome_firewall"].values())
    assert all(value is False for value in protocol["outcome_firewall_at_freeze"].values())
    assert all(value is False for value in parser["outcome_firewall"].values())
    assert all(value is False for value in digest["firewall"].values())
