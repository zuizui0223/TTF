from __future__ import annotations

import json
from pathlib import Path


CENSUS = Path("benchmarks/frozen/genetic_conditional_two_panel_response_blind_census_v0.1.json")
PROTOCOL = Path("docs/supporting/genetic_conditional_transfer_protocol_v0.1.json")


def test_conditional_successor_is_two_panel_and_response_blind() -> None:
    census = json.loads(CENSUS.read_text())
    protocol = json.loads(PROTOCOL.read_text())

    assert census["status"] == "RESPONSE_BLIND_TWO_PANEL_FEASIBILITY_COMPLETE"
    assert census["source_archive"]["sha256"] == "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
    assert census["census"]["metadata_eligible_species_after_exclusions"] == 14612
    assert census["census"]["selected_total"] == 500
    assert census["census"]["geometry_failures_before_500"] == 0
    assert census["census"]["development_species"] == 250
    assert census["census"]["confirmatory_species"] == 250
    assert census["census"]["development_confirmatory_overlap"] == 0
    assert census["exclusions"]["old_new_overlap"] == 0
    assert all(value is False for value in census["outcome_firewall"].values())

    assert census["development_panel"]["support_500km"]["coverage_ge_0_25"]["target_species_with_ge_5_sources"] == 113
    assert census["confirmatory_panel"]["support_500km"]["coverage_ge_0_25"]["target_species_with_ge_5_sources"] == 108
    assert census["development_panel"]["support_500km"]["coverage_ge_0_25_same_order"]["target_species_with_ge_5_sources"] == 79
    assert census["confirmatory_panel"]["support_500km"]["coverage_ge_0_25_same_order"]["target_species_with_ge_5_sources"] == 86

    assert protocol["status"] == "FROZEN_RESPONSE_BLIND_ESTIMAND_BEFORE_ANY_NEW_PANEL_NUCLEOTIDE_IDENTITY"
    assert protocol["geographic_conditioning"]["support_radius_km"] == 500.0
    assert protocol["geographic_conditioning"]["source_in_pool"] == "target coverage >= 0.25"
    assert protocol["geographic_conditioning"]["minimum_source_species_per_target"] == 5
    assert protocol["primary_estimand"]["name"] == "same_order_increment_beyond_geography"\n    assert protocol["primary_estimand"]["group_label"] == "exact non-empty phylogatR order"\n    assert protocol["secondary_estimand"]["name"] == "geographic_conditioning_increment_support_diagnostic"\n    assert "already spatially local" in protocol["secondary_estimand"]["role"]
    assert protocol["excluded_estimands"]["same_family"].startswith("not pursued")
    assert protocol["parent_empirical_result"]["result_sha256"] == "f5d19fa50c7c18cbf2a110c0afb9c70dcd543c44b77521c015938013843e72fb"
    assert all(value is False for value in protocol["outcome_firewall"].values())
