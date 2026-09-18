from __future__ import annotations

import json
from pathlib import Path


PROTOCOL = Path("docs/supporting/genetic_conditional_transfer_protocol_v0.1.json")
CENSUS = Path("benchmarks/frozen/genetic_conditional_response_blind_census_v0.1.json")


def test_conditional_genetic_protocol_is_response_blind_and_species_disjoint() -> None:
    protocol = json.loads(PROTOCOL.read_text())
    census = json.loads(CENSUS.read_text())

    assert protocol["schema"] == "ttf_genetic_conditional_transfer_protocol_v0.1"
    assert protocol["status"].startswith("FROZEN_RESPONSE_BLIND")
    assert protocol["parent_empirical_result"]["rerun_or_retuning_allowed"] is False

    panel = protocol["response_blind_panel"]
    assert panel["exclude_parent_phase1_species_count"] == 250
    assert panel["minimum_unique_localities"] == 12
    assert panel["neighbor_fraction"] == 0.15
    assert panel["minimum_endpoint_disjoint_ibd_training_edges"] == 5

    geo = protocol["co_distribution"]
    assert geo["radius_km"] == 500
    assert geo["minimum_coverage_fraction"] == 0.25
    assert geo["minimum_source_species_per_target"] == 3

    primary = protocol["hypotheses"]["primary"]
    assert primary["name"] == "order_given_geography"
    assert "Delta_order_given_geo" in primary["estimand"]
    assert "not a measurement of dispersal ability" in primary["interpretation_if_qualified_positive"]

    assert census["status"] == "RESPONSE_BLIND_CONDITIONAL_PANEL_FEASIBLE"
    assert census["parent_phase1_species_excluded"] == 250
    assert census["parent_new_species_overlap"] == 0
    assert census["metadata_eligible_species_after_exclusions"] == 14612
    assert census["selection"]["selected_species"] == 250
    assert census["selection"]["reserve_ranks_checked"] == 250
    assert census["selection"]["geometry_failures_before_250"] == 0
    assert census["split"] == {
        "train": 125,
        "eval": 125,
        "train_species_list_sha256": "787b396805eed23fb2e5ede53d02c7433f2ec3ad67ac97cdc100d04c96ed2a10",
        "eval_species_list_sha256": "d95b77b24ffb612576035b383737bbeade9c16eb7e198d54b2ae27bd7e2f4216",
    }

    support = census["response_blind_support"]
    assert support["geo_coverage_at_least_0_25"]["targets_with_at_least_3_sources"] == 113
    assert support["geo_coverage_at_least_0_25_same_order"]["targets_with_at_least_3_sources"] == 80
    assert support["geo_coverage_at_least_0_25_same_family"]["targets_with_at_least_3_sources"] == 11
    assert support["geo_coverage_at_least_0_25_same_family"]["decision"] == "DO_NOT_USE_AS_CONFIRMATORY_CONDITION"

    assert census["sequence_identity_opened"] is False
    assert census["pairwise_genetic_distances_opened"] is False
    assert census["empirical_ttf_opened"] is False

    firewall = protocol["outcome_firewall"]
    assert firewall and all(value is False for value in firewall.values())
    assert protocol["empirical_opening"]["status"] == "NOT_AUTHORIZED"
