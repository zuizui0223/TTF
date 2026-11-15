"""Costus source-only explicit author-key and unnamed R-rowname test."""
import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_osf_source_taxon_crosswalk_20261115.py"
spec=importlib.util.spec_from_file_location("costus_osf_taxon",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_unnamed_source_rowname_is_not_a_scent_value():
    raw=b",beta_ocimene,other_scent\nCostus allenii,12.0,99\nCostus lasius,0,21\n"
    r=m.inspect(raw,[""])
    assert r["keys"][""]==["costus allenii","costus lasius"]
    assert "12.0" not in repr(r) and "99" not in repr(r)
    assert r["no_numeric_value_read"] is True


def test_exact_author_crosswalk_keys_no_species_merging():
    raw=b"Name,Name_in_info,Pollination\nCostus_allenii,allenii,Bird\nCostus_scaber_18044,scaber 18044,Bee\nCostus_scaber_18085,scaber 18085,Bee\n"
    r=m.inspect(raw,["Name","Name_in_info"])
    assert r["distinct_taxon_keys_per_source_column"]["Name_in_info"]==3
    assert r["keys"]["Name_in_info"]==["allenii","scaber 18044","scaber 18085"]


def test_missing_expected_original_key_column_fail_closed():
    with pytest.raises(ValueError,match="key fields differ"):
        m.inspect(b"Species,amount\nA,1\n",["Name"])
