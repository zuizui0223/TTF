"""Synthetic original Table S1 text-layout/identity-only tests, no biology."""
from pathlib import Path
import importlib.util

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_original_S1_two_parser_52_51_source_diagnostic_20261116.py"
spec=importlib.util.spec_from_file_location("penstemon_s1_2engine",FILE)
s=importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)

NAMES=["Penstemon barbatus","Penstemon virgatus"]

def test_preserves_complete_asterisk_id_and_phenotyped_flag_only():
    raw=("Table S1\nPopulation Phenotyped\n"
         "CAW-001* P. barbatus yes\n"
         "CAW-002 P. virgatus no\n"
         "BWS-003 P. barbatus\n"
         "yes\n")
    p=s.text_profile(raw,"synthetic",NAMES)
    assert p["exact_rows"]==2
    assert p["exact_phenotyped_yes"]==1 and p["exact_phenotyped_no"]==1
    assert p["id_mentions_not_complete_rows"]==["BWS-003"]
    assert "Sample Phenotyped" not in repr(p)
    assert p["source_taxa_from_complete_rows"] == {"Penstemon barbatus":1,"Penstemon virgatus":1}

def test_disallows_table_s2_data_from_source_scope():
    page11="Table S1\nPhenotyped\nCAW-001 P. barbatus yes"
    page12="BWS-002 P. virgatus no\nTable S2\nCAW-003 P. barbatus yes"
    data=s.select_s1_source_text(page11,page12)
    assert "BWS-002" in data and "CAW-003" not in data

def test_any_position_detects_linebroken_source_ids_without_repair():
    text="Table S1\nPhenotyped\nCAW-001 P. barbatus yes\ncontinued BWS - 002 species no"
    p=s.text_profile(text,"synthetic",NAMES)
    assert "BWS-002" in p["id_mentions_not_complete_rows"]
    assert p["exact_rows"]==1


def test_original_unhyphenated_TTC026503_is_not_invented_synonym():
    original=("Table S1\nPopulation Phenotyped\n"
       "CAW-034 P. cardinalis yes\n"
       "TTC026503 P. cardinalis Guadalupe Mtns, TX MJ TTC no\n")
    profile=s.text_profile(original,"visual-verified",["Penstemon cardinalis"])
    assert profile["exact_rows"]==2
    assert profile["exact_phenotyped_yes"]==1
    assert profile["exact_phenotyped_no"]==1
    assert "TTC026503" in profile["exact_row_ids"]
