from __future__ import annotations

import json
from pathlib import Path


CENSUS = Path("benchmarks/frozen/genetic_conditional_two_panel_response_blind_census_v0.1.json")
PROTOCOL = Path("docs/supporting/genetic_conditional_transfer_protocol_v0.1.json")
QUALIFICATION = Path("docs/supporting/genetic_conditional_order_qualification_rule_v0.1.json")


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

    geometry = census["response_blind_geometry_artifacts"]
    assert geometry["development"]["geometry_csv_sha256"] == "d4affc28a2f9da89d0566b0f52123f4711009c99e985baff9f44767da996a9cd"
    assert geometry["development"]["manifest_sha256"] == "f78bc6fa735010270db287764be0f818ef09ce416b00e0cb0c1aa1f00f7513ac"
    assert geometry["confirmatory"]["geometry_csv_sha256"] == "fc93f1aaac2abc369f3e94c99b73e1d139f635638e197d1f71f977a49a44dbba"
    assert geometry["confirmatory"]["manifest_sha256"] == "debe44fabd10c53b2cf5fad13ad0c182e57429edd193619fa7a5c34e4dbdb880"
    assert geometry["development"]["exact_geometry_reconstructed_without_nucleotide_identity"] is True
    assert geometry["confirmatory"]["exact_geometry_reconstructed_without_nucleotide_identity"] is True

    assert census["development_panel"]["support_500km"]["coverage_ge_0_25"]["target_species_with_ge_5_sources"] == 113
    assert census["confirmatory_panel"]["support_500km"]["coverage_ge_0_25"]["target_species_with_ge_5_sources"] == 108
    assert census["development_panel"]["support_500km"]["coverage_ge_0_25_same_order"]["target_species_with_ge_5_sources"] == 79
    assert census["confirmatory_panel"]["support_500km"]["coverage_ge_0_25_same_order"]["target_species_with_ge_5_sources"] == 86

    assert protocol["status"] == "FROZEN_RESPONSE_BLIND_ESTIMAND_BEFORE_ANY_NEW_PANEL_NUCLEOTIDE_IDENTITY"
    assert protocol["geographic_conditioning"]["support_radius_km"] == 500.0
    assert protocol["geographic_conditioning"]["source_in_pool"] == "target coverage >= 0.25"
    assert protocol["geographic_conditioning"]["minimum_source_species_per_target"] == 5
    assert protocol["primary_estimand"]["name"] == "same_order_increment_beyond_geography"
    assert protocol["primary_estimand"]["group_label"] == "exact non-empty phylogatR order"
    assert protocol["secondary_estimand"]["name"] == "geographic_conditioning_increment_support_diagnostic"
    assert "already spatially local" in protocol["secondary_estimand"]["role"]
    assert protocol["excluded_estimands"]["same_family"].startswith("not pursued")
    assert protocol["parent_empirical_result"]["result_sha256"] == "f5d19fa50c7c18cbf2a110c0afb9c70dcd543c44b77521c015938013843e72fb"
    assert all(value is False for value in protocol["outcome_firewall"].values())


def test_conditional_order_qualification_is_frozen_before_outcomes() -> None:
    rule = json.loads(QUALIFICATION.read_text())

    assert rule["schema"] == "ttf_genetic_conditional_order_qualification_rule_v0.1"
    assert rule["status"] == "FROZEN_BEFORE_ANY_CONDITIONAL_SYNTHETIC_QUALIFICATION_RESULT_OR_NEW_PANEL_NUCLEOTIDE_IDENTITY"

    panel = rule["development_panel"]
    assert panel["species"] == 250
    assert panel["train_species"] == 125
    assert panel["eval_species"] == 125
    assert panel["species_digest_sha256"] == "43ddab768fca83f995b8380069646d373cb112c7cb9b3042f19ffd311a4a6583"
    assert panel["empirical_sequence_identity_opened"] is False
    assert panel["empirical_genetic_outcome_opened"] is False

    estimator = rule["frozen_estimator"]
    assert estimator["support_radius_km"] == 500.0
    assert estimator["minimum_target_coverage"] == 0.25
    assert estimator["minimum_source_species"] == 5
    assert estimator["eligible_same_order_eval_species"] == 79
    assert estimator["bootstrap_resamples"] == 1999

    private = rule["synthetic_worlds"]["private_null_cells"]
    assert [cell["residual_amplitude"] for cell in private] == [0.5, 1.0, 2.0, 3.0]
    assert all(cell["worlds"] == 500 for cell in private)
    positive = rule["synthetic_worlds"]["positive_control"]
    assert positive["cell"] == "same_order_A2"
    assert positive["residual_amplitude"] == 2.0
    assert positive["worlds"] == 500

    assert rule["qualification"]["type1_wilson95_upper_ceiling"] == 0.10
    assert rule["qualification"]["power_wilson95_lower_floor"] == 0.80
    assert all(value is False for value in rule["outcome_firewall"].values())
