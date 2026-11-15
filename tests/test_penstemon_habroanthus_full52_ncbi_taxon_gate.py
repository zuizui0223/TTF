"""Synthetic strict all-52 Penstemon BioSample species source-coverage tests."""
import importlib.util
from pathlib import Path
import pytest

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_habroanthus_full52_ncbi_taxon_gate_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_52_spp_source",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def make_source():
    acc=[f"SAMN6099{i:04d}" for i in range(52)]
    species=[
       "Penstemon barbatus","Penstemon virgatus","Penstemon cardinalis",
       "Penstemon glaber","Penstemon eatonii","Penstemon laevis",
       "Penstemon labrosus","Penstemon putus","Penstemon deaveri",
    ]
    rows=[{"source_biosample_accession":id,
           "scientific_name":species[i%len(species)]} for i,id in enumerate(acc)]
    summary={"status":"NCBI_BIOSAMPLE_SOURCE_SUMMARY_ONLY","records":rows}
    return acc,summary


def test_all_four_sister_groups_verified_without_biological_responses():
    acc,summary=make_source()
    result=m.source_coverage(summary,acc)
    assert result["source_genome_sample_count"]==52
    assert result["source_distinct_taxon_labels"]==9
    assert result["number_of_four_primary_sister_source_pairs_covered"]==4
    assert result["sample_attribute_values_never_exported"] is True
    assert result["not_a_molecular_phenotype_or_species_transfer_result"] is True


def test_unmatched_sample_fails_taxonomic_admission():
    acc,summary=make_source()
    summary["records"][-1]["source_biosample_accession"]="SAMN00000"
    result=m.source_coverage(summary,acc)
    assert result["status"]=="NCBI_SOURCE_SAMPLE_ACCESSIONS_NOT_FULLY_JOINED"
    assert result["not_a_biological_null"] is True


def test_source_taxon_missing_prevents_full_source_label_claim():
    acc,summary=make_source()
    summary["records"][1]["scientific_name"]=""
    result=m.source_coverage(summary,acc)
    assert result["status"]=="BIOPROJECT_SOURCE_COMPLETE_TAXON_UNRESOLVED"
    assert result["blank_scientific_names"]==1


def test_fewer_than_52_source_samples_fail_closed():
    acc,summary=make_source()
    with pytest.raises(ValueError,match="52"):
        m.source_coverage(summary,acc[:50])


def test_comparison_names_are_frozen_prior_to_source_values():
    assert m.N_SOURCE==52
    assert sorted(m.SOURCE_COMPARISONS)==[
        "barbatus_to_virgatus","cardinalis_to_glaber",
        "eatonii_to_laevis","labrosus_to_deaveri","labrosus_to_putus"]
