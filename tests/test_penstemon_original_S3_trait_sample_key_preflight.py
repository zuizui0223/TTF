"""Synthetic-only Table S3 original voucher-key schema tests."""
import importlib.util
from pathlib import Path
import pytest

F=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_original_S3_trait_sample_key_preflight_20261116.py"
sp=importlib.util.spec_from_file_location("s3",F)
m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)

def test_source_ids_are_identifiers_not_trait_values():
    text="Table S3\nPopulation Nectar Volume\nCAW-094 P. labrosus 22.345 8.9\nBWS-134 P. barbatus 30.1\nTable S4\nCAW-999 0.4\n"
    d=m.page_metadata(text,"mock")
    assert d["distinct_voucher_tokens"]==2
    assert d["voucher_tokens"]==["BWS-134","CAW-094"]
    assert "22.345" not in str(d) and "CAW-999" not in str(d)

def test_no_table_s3_heading_fails_closed():
    with pytest.raises(ValueError,match="Table S3"):
        m.page_metadata("Table S2\nP. labrosus", "mock")
