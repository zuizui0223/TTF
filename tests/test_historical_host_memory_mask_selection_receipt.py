"""Public frozen aggregate-only mask-selection receipt review, no private source."""
import json
import math
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
RECEIPT=REPO/"benchmarks/frozen/historical_host_memory_mask_attrition_representativeness_v0.1.json"
MASK=REPO/"benchmarks/frozen/historical_host_memory_private_confirmatory_mask_259_survivor_replay_v0.1.json"


def test_aggregate_sha_identity_and_exact_mask_counts():
    r=json.loads(RECEIPT.read_text())
    src=json.loads(MASK.read_text())
    assert r["status"]=="EXPLORATORY_POST_MASK_REPRESENTATIVENESS_AUDIT"
    assert r["source_sha256"]["mask_ledger"]==src["confirmatory_mask"]["private_ledger_sha256"]
    assert r["source_sha256"]["candidate_table"]==src["source"]["candidate_sha256"]
    assert r["source_sha256"]["species_diagnostics"]==src["source"]["species_diagnostics_sha256"]
    assert r["source_sha256"]["roles"]==src["source"]["panel_roles_sha256"]
    assert r["confirmatory_species"]==321
    assert r["mask_pass_species"]==src["confirmatory_mask"]["valid_mask_species"]==259
    assert r["mask_fail_species"]==src["confirmatory_mask"]["failed_mask_species"]==62
    assert r["original_confirmatory_edges"]==src["confirmatory_mask"]["original_edges_in_confirmatory"]==41592
    assert r["survivor_confirmatory_edges"]==src["survivor_information"]["original_edges_surviving_exactly"]==20362
    assert r["invalid_original_edges_with_finite_masks"]==src["confirmatory_mask"]["failed_edges_lacking_jointly_canonical_sites"]==521


def test_graph_exposure_bins_and_family_rates_are_accounted_for():
    r=json.loads(RECEIPT.read_text())
    bins=r["original_graph_edge_count_bins"]
    assert len(bins)==5
    assert sum(x["original_species"] for x in bins)==321
    assert sum(x["mask_survivor_species"] for x in bins)==259
    assert all(0<=x["mask_survival_rate"]<=1 for x in bins)
    assert bins[0]["mask_survivor_species"]==154
    assert bins[0]["original_species"]==176
    assert bins[-1]["mask_survivor_species"]==4
    assert bins[-1]["original_species"]==12
    families=r["families_with_at_least_10_original_species"]
    assert len(families)==9
    assert all(x["original_species"]==x["mask_survivor_species"]+x["mask_fail_species"] for x in families)
    assert all(abs(x["mask_survival_rate"]-x["mask_survivor_species"]/x["original_species"])<1e-14 for x in families)
    by={x["family"]:x for x in families}
    assert by["Gracillariidae"]["mask_survivor_species"]==19
    assert by["Gracillariidae"]["original_species"]==39
    assert by["Noctuidae"]["mask_survivor_species"]==51
    assert by["Noctuidae"]["original_species"]==53


def test_exploratory_selection_does_not_change_preregistered_claim():
    r=json.loads(RECEIPT.read_text())
    assert math.isclose(r["surviving_original_edge_fraction"],20362/41592)
    assert math.isclose(r["mask_pass_species_fraction"],259/321)
    assert r["group_medians"]["pass"]["frozen_edge_count"]==43
    assert r["group_medians"]["fail"]["frozen_edge_count"]==78.5
    assert r["interpretation"]["no_extra_exclusions_or_backfill"] is True
    assert r["interpretation"]["families_are_not_an_authorized_subgroup_analysis"] is True
    assert r["interpretation"]["observed_genetic_distance_or_effect_unopened"] is True
    assert not any(key in r for key in ("species_rows","per_species_mask","nucleotide_identity","observed_genetic_effect"))
