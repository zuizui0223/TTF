"""Offline source-only TTF-Q four-origin geometry tests. No phenotype values."""
import copy
import importlib.util
import json
from pathlib import Path
import pytest

SOURCE=Path(__file__).resolve().parents[1]/"research"/"prospective"
P=SOURCE/"penstemon_four_origin_source_genome_blocks_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_four_source_blocks",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def original():
    return json.loads((SOURCE/"penstemon_habroanthus_52_ncbi_biosample_taxon_receipt_20261115.json").read_text())


def test_four_real_source_blocks_are_strictly_accession_disjoint():
    out=m.specimen_geometry(original())
    assert out["status"]=="FOUR_ORIGIN_GENOMIC_SOURCE_COVERAGE_ONLY_PHENOTYPE_JOIN_HOLD"
    assert out["genome_samples_in_four_origin_blocks"]==33
    assert out["genome_samples_outside_origin_blocks"]==19
    assert out["phenotyped_samples_in_four_origin_blocks_possible_bounds"]==[15,33]
    assert len(out["strict_origin_heldout_folds"])==4
    assert [f["heldout_genome_samples"] for f in out["strict_origin_heldout_folds"]]==[10,4,9,10]
    assert [f["training_genome_samples"] for f in out["strict_origin_heldout_folds"]]==[23,29,24,23]
    assert all(f["training_origin_block_count"]==3 for f in out["strict_origin_heldout_folds"])
    assert all(f["heldout_derived_genome_samples"]>0 and f["heldout_bee_comparator_genome_samples"]>0 for f in out["strict_origin_heldout_folds"])
    assert all(v["possible_phenotyped_count_bounds_from_34_of_52"][0]==0 for v in out["block_support"].values())
    assert out["no_flower_trait_snp_or_fitness_outcome_read"] is True
    assert out["no_TTF_cross_origin_result_claim"] is True


def test_missing_original_biosample_fails_closed():
    payload=original()
    payload["source_matched_taxonomy"]["source_biosample_taxon_crosswalk"].pop()
    with pytest.raises(ValueError,match="source specimen ID"):
        m.specimen_geometry(payload)


def test_duplicated_original_accession_fails_closed():
    payload=original()
    a=payload["source_matched_taxonomy"]["source_biosample_taxon_crosswalk"]
    a[-1]["biosample_accession"]=a[0]["biosample_accession"]
    with pytest.raises(ValueError,match="source specimen ID"):
        m.specimen_geometry(payload)


def test_unavailable_source_comparator_must_not_be_guessed():
    payload=original()
    a=payload["source_matched_taxonomy"]["source_biosample_taxon_crosswalk"]
    for item in a:
        if item["scientific_name"]=="Penstemon glaber":
            item["scientific_name"]="Penstemon hallii"
    counter=payload["source_matched_taxonomy"]["source_taxon_sample_counts"]
    counter["Penstemon hallii"]+=counter.pop("Penstemon glaber")
    payload["source_matched_taxonomy"]["source_distinct_taxon_labels"]=20
    with pytest.raises(ValueError,match="source taxon labels"):
        m.specimen_geometry(payload)


def test_github_source_blob_pin():
    data=(SOURCE/"penstemon_habroanthus_52_ncbi_biosample_taxon_receipt_20261115.json").read_bytes()
    assert m.source_blob_sha(data)==m.BLOB_SHA
