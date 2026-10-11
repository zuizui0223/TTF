"""Public aggregate receipt QA. Raw sequences/masks stay in the owner's private Library."""
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RECEIPT = REPO / "benchmarks/frozen/historical_host_memory_private_confirmatory_mask_259_survivor_replay_v0.1.json"


def wilson(k: int, n: int) -> tuple[float, float]:
    z = 1.959963984540054
    p = k/n
    d = 1.0+z*z/n
    center = (p+z*z/(2*n))/d
    width = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return center-width, center+width


def test_frozen_receipt_original_sources_are_canonical():
    report=json.loads(RECEIPT.read_text())
    pre=json.loads((REPO/"benchmarks/frozen/historical_host_memory_synthetic_qualification_result_v0.1.json").read_text())
    qualified=json.loads((REPO/"benchmarks/frozen/historical_host_memory_information_qualification_result_v0.2.json").read_text())
    recovered=json.loads((REPO/"benchmarks/frozen/historical_host_memory_original_archive_provenance_correction_v0.1.json").read_text())
    assert report["source"]["archive_sha256"] == recovered["original_source"]["canonical_sha256"]
    assert len(report["source"]["archive_sha256"]) == 64
    assert report["source"]["archive_size_bytes"] == 274988692
    assert report["source"]["candidate_sha256"] == qualified["sources"]["frozen_642_candidates"]
    assert report["source"]["panel_roles_sha256"] == pre["design"]["assignments_sha256"]
    assert report["source"]["full_external_predictor_sha256"] == pre["design"]["predictor_sha256"]
    assert report["source"]["original_synthetic_result_sha256"] == pre["source_run"]["exact_inner_receipt_sha256"]


def test_exact_259_survivors_original_edge_subset_and_information_gates():
    v=json.loads(RECEIPT.read_text())
    m=v["confirmatory_mask"]
    i=v["survivor_information"]
    assert m["panel_species"] == 321
    assert m["valid_mask_species"] == i["survivor_species"] == 259
    assert m["failed_mask_species"] == 62
    assert sum(m["fail_species_reason_counts"].values()) == 62
    assert m["failed_edges_lacking_jointly_canonical_sites"] == 521
    assert i["original_edges_surviving_exactly"] == 20362
    assert i["removed_species_aggregate_edges"] + i["original_edges_surviving_exactly"] == m["original_edges_in_confirmatory"] == 41592
    assert i["species_unique_fraction_at_least_0_05"] == 244
    assert i["species_unique_fraction_at_least_0_05_fraction"] == 244/259
    assert i["median_unique_fraction_M_host"] >= 0.10
    assert i["species_unique_fraction_at_least_0_05_fraction"] >= 0.70
    assert i["species_condition_number_le_30_fraction"] >= 0.90
    assert i["survivor_edges_equal_byte_for_byte_to_frozen_71520_edge_predictor_filtered_by_survivor_species"] is True
    assert i["survivor_roles_equal_byte_for_byte_to_frozen_641_species_roles_filtered_by_survivor_species"] is True
    for sha in ("survivor_edges_sha256", "survivor_roles_sha256"):
        assert len(i[sha]) == 64


def test_independent_post_mask_simulation_wilson_and_never_empirical_opening():
    v=json.loads(RECEIPT.read_text())
    s=v["independent_survivor_synthetic"]
    assert (s["null_reference_worlds"],s["independent_null_evaluation_worlds"],s["positive_worlds"]) == (1999,500,500)
    assert s["seed_namespace"] == "historical-host-memory-survivor-v0.1"
    assert s["null_rejections"]==30 and s["positive_rejections"]==485
    nlo,nhi=wilson(s["null_rejections"],500)
    plo,phi=wilson(s["positive_rejections"],500)
    assert abs(nlo-s["null_wilson_lower"])<1e-14
    assert abs(nhi-s["null_wilson_upper"])<1e-14
    assert abs(plo-s["positive_wilson_lower"])<1e-14
    assert abs(phi-s["positive_wilson_upper"])<1e-14
    assert s["null_wilson_upper"] <= s["null_wilson_upper_gate"] == 0.10
    assert s["positive_wilson_lower"] >= s["positive_wilson_lower_gate"] == 0.80
    assert s["qualification_status"] == "PASS_TO_SEPARATE_CONFIRMATORY_IDENTITY_AUTHORIZATION_REVIEW"
    f=v["reproduction_and_firewall"]
    assert all(f[k] is False for k in (
        "raw_zip_uploaded_to_public_GitHub","nucleotide_identity_saved",
        "genetic_distance_computed","empirical_post_ibd_turnover_computed",
        "empirical_host_history_coefficient_computed"))
    assert f["separate_qualified_one_shot_identity_authorization_required"] is True
