"""Synthetic original public Penstemon Table S1 52/34 presence-only guards."""
import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_original_full52_table_s1_phenotyped_gate_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_s1_full_source",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fixture():
    return """Table S1. Collection information
Population identifier Species determination Phenotyped
CAW-018* P. gl aber Midland, CO CAW KANU no
CAW-034 P. car di nalis (car-034) Bonito Lake, NM CAW KANU yes
LW-038 P. bar ba tus (bar-038) Victor, CO LW KANU yes
"""


def test_source_voucher_ids_and_yes_no_no_trait_numeric():
    d=m.parse_s1_page(fixture(),
       ["Penstemon glaber","Penstemon cardinalis","Penstemon barbatus"],
       expected_n=3,expected_yes=2)
    assert d["source_sample_count"]==3
    assert d["source_phenotyped_count"]==2
    assert d["source_taxon_matched_count"]==3
    keys={r["original_voucher_population_id"]:r for r in d["original_voucher_rows"]}
    assert keys["CAW-018"]["original_5taxon_introgression_marker_asterisk"] is True
    assert keys["CAW-018"]["phenotyped"] is False
    assert keys["CAW-034"]["source_species_taxon_match"]=="Penstemon cardinalis"
    assert keys["LW-038"]["phenotyped"] is True
    assert d["sample_to_SAMN_individual_ID_join"]=="NOT_VERIFIED"


def test_wrong_table_total_or_wrong_yes_no_count_fails_closed():
    with pytest.raises(ValueError,match="source phenotyped count"):
        m.parse_s1_page(fixture(),
          ["Penstemon glaber","Penstemon cardinalis","Penstemon barbatus"],
          expected_n=3,expected_yes=3)
    with pytest.raises(ValueError,match="row structure"):
        m.parse_s1_page(fixture(),
          ["Penstemon glaber","Penstemon cardinalis","Penstemon barbatus"],
          expected_n=4,expected_yes=2)


def test_taxonomic_synonyms_never_guessed():
    d=m.parse_s1_page(fixture(),
       ["Penstemon glaber","Penstemon cardinalis","Penstemon virgatus"],
       expected_n=3,expected_yes=2)
    assert "LW-038" in d["unmatched_original_voucher_ids_no_taxon_synonym_guess"]


def test_source_pdf_SHA_pin_and_original_NCBI_receipt_blob_pin():
    assert m.PDF_SHA256=="3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943"
    assert m.git_blob_sha(m.GENOME_SOURCE.read_bytes())==m.GENOME_SOURCE_BLOB


def test_visually_authenticated_unhyphenated_TTC026503_restores_one_source_row_only():
    source=("Table S1. Collection information\nPopulation identifier Species determination Phenotyped\n"
       "CAW-034* P. cardinalis (car-034) Bonito Lake, NM CAW KANU yes\n"
       "TTC026503 P. cardinalis Guadalupe Mtns, TX MJ TTC no\n")
    d=m.parse_s1_page(source,["Penstemon cardinalis"],expected_n=2,expected_yes=1)
    keys={r["original_voucher_population_id"]:r for r in d["original_voucher_rows"]}
    assert set(keys)=={"CAW-034","TTC026503"}
    assert keys["TTC026503"]["phenotyped"] is False
    assert keys["TTC026503"]["source_species_taxon_match"]=="Penstemon cardinalis"


def test_original_34_phenotyped_gate_remains_strict_after_52nd_record_fix():
    source=("Table S1\nPopulation Phenotyped\n"
       "CAW-034 P. cardinalis yes\n"
       "TTC026503 P. cardinalis no\n")
    with pytest.raises(ValueError,match="source phenotyped count 1 differs"):
        m.parse_s1_page(source,["Penstemon cardinalis"],expected_n=2,expected_yes=2)
