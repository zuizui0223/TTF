"""Validate frozen one-shot P1 receipt without reopening raw VCF/genotypes."""
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1] / "research" / "prospective"
P1=json.loads((ROOT/"senecio_p1_first_opening_aggregate_receipt_20261115.json").read_text())
CROSS=json.loads((ROOT/"senecio_p0_p1_crosslayer_descriptive_20261115.json").read_text())


def test_source_execution_pin_and_correct_unit():
    assert P1["source_git_blob"]=="fe1a452b183e0690e013cef57e59bedb13c7e072"
    assert P1["source_execution"]["run_id"]==43693277926
    assert P1["source"]["vcf_header_sample_count"]==1009
    assert P1["source"]["matched_sample_count"]==986
    assert P1["source"]["source_variant_records"]==9269
    assert P1["source"]["distinct_snv_contigs"]==3431
    assert P1["heldout_pair_count"]==8
    assert P1["train_only_marker_count_per_fold"]==100
    assert P1["status"]=="EXPLORATORY_P1_GENOMIC_LOCALITY_TRANSFER_REPLAY"


def test_eight_pair_aggregates_use_pair_not_individual_weight():
    folds=P1["folds"]
    assert len(folds)==8
    assert len({f["pair"] for f in folds})==8
    assert sum(f["n_dune_original"]+f["n_headland_original"] for f in folds)==986
    assert all(f["status"]=="EVALUABLE_LOCALITY_AUC" for f in folds)
    assert all(0<=f["auc"]<=1 for f in folds)
    assert all(f["distinct_training_selected_contigs"]==100 for f in folds)
    assert all(f["n_dune_original"]==f["n_dune_test"] and f["n_headland_original"]==f["n_headland_test"] for f in folds)
    mean=sum(f["auc"] for f in folds)/8
    assert math.isclose(mean,P1["summary"]["pair_equal_mean_auc"],abs_tol=1e-12)
    assert sum(f["auc"]>.5 for f in folds)==P1["summary"]["pairs_auc_above_0_5"]==3
    assert all(f["dune_marker_call_rate"]>.9 and f["headland_marker_call_rate"]>.9 for f in folds)


def test_post_outcome_crosslayer_comparison_does_not_change_estimands():
    assert CROSS["status"]=="POSTOUTCOME_DESCRIPTIVE_CROSSLAYER_SENECIO_P0_P1_ONLY"
    assert CROSS["predeclared_primary"] is False
    assert CROSS["source_p1_run_id"]==43693277926
    assert len(CROSS["rows"])==8
    keyed={f["pair"]:f for f in P1["folds"]}
    for row in CROSS["rows"]:
        assert 0<=row["morphology_signs_agree"]<=9
        assert row["genotype_marker_count"]==100
        assert math.isclose(row["genotype_selected_snp_auc"],keyed[row["pair"]]["auc"])
    assert CROSS["summary"]["pairs_with_both"]==5
    assert sum(x["morphology_signs_agree"] for x in CROSS["rows"])==65
    assert "causal_genetic_loci" in P1["not_claimed"]
