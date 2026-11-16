"""Offline original Table S1 source recovery guards; no phenotype magnitudes."""
import importlib.util
from pathlib import Path

F=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_authentic_s1_s3_one_zip_recovery_20261116.py"
spec=importlib.util.spec_from_file_location("penstemon_source_recovery",F)
s=importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


def test_explicit_partial_source_never_becomes_52_34_admitted():
    text="Table S1\nPopulation Phenotyped\nCAW-018 P. barbatus yes\nCAW-043 P. virgatus no\n"
    t=s.original_s1_partial_source(text,["Penstemon barbatus","Penstemon virgatus"])
    assert t["complete_rows_from_pdf"]==2
    assert t["explicit_phenotyped_yes"]==1
    assert t["source_taxon_unmatched_ids"]==[]
    assert t["status"]=="INCOMPLETE_SOURCE_DO_NOT_REDEFINE_52_34_GATE"
    assert t["four_origin_partial_lower_bound_coverage"]["B1_barbatus_virgatus"]["origin_both_derived_and_any_bee_phenotyped"] is False


def test_unmatched_species_is_not_guessed():
    text="Table S1\nCAW-190 P. strictiformis yes\n"
    t=s.original_s1_partial_source(text,["Penstemon comarrhenus","Penstemon cardinalis"])
    assert t["source_taxon_unmatched_ids"]==["CAW-190"]
    assert t["explicit_phenotyped_yes"]==1


def test_original_pdf_sha_is_fixed():
    assert s.PDF_SHA256=="3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943"
    assert s.MAX_ATTEMPTS==2
