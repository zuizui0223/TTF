"""Bounded original Table S1 source-layout preview parser tests."""
import importlib.util
from pathlib import Path
import pytest

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_original_table_s1_specimen_header_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_s1_header",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_source_header_exposes_only_bounded_original_keys():
    text="Table S1: voucher inventory\nSample Species\nCAW-023 Penstemon barbatus\nBWS-104 Penstemon virgatus\nLW-009 P. eatonii"
    r=m.inspect_s1_layout(text)
    assert r["status"]=="ORIGINAL_TABLE_S1_HEADER_AND_LAYOUT_PREVIEW_ONLY"
    assert r["sample_author_key_prefix_occurrences"]=={"CAW":1,"BWS":1,"LW":1,"SAMN":0}
    assert r["source_table_proves_34_phenotyped"]=="NOT_ASSESSED"
    assert r["no_floral_trait_or_genotype_values_read"] is True


def test_source_without_table_title_is_not_accepted():
    with pytest.raises(ValueError,match="heading missing"):
        m.inspect_s1_layout("Table S2\nCAW-104")
