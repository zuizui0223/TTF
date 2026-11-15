"""Synthetic Costus published tree-tip and categorical-source readiness tests."""
import importlib.util
from pathlib import Path
import pytest
pytest.importorskip("dendropy")  # Optional source-only tree CI provides this dependency.

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_original_tree_label_geometry_20261115.py"
spec=importlib.util.spec_from_file_location("costus_tree_source",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_tree_tip_taxa_are_not_inferred_from_scent_data():
    n=b"(Costus_allenii_98050:1,(Costus_lasius_18080:2,Costus_wilsonii_18050:3):1);"
    keys=m.original_tree_tips(n)
    assert keys==["costus_allenii_98050","costus_lasius_18080","costus_wilsonii_18050"]


def test_duplicate_tree_tip_label_rejected():
    with pytest.raises(Exception,match="Multiple occurrences|duplicate|repeated"):
        m.original_tree_tips(b"(Costus_A:1,Costus_A:2);")


def test_scent_mean_keys_are_from_first_column_only():
    payload=b",ebetaocimene,unknown\nTaxon_a,33,999\n"
    with pytest.raises(ValueError,match="count"):
        m.means_row_labels(payload)


def test_input_tree_metadata_is_locked_independently_of_scent_values():
    assert m.FILES[m.TREE_NAME]["sha256"] == "24f1cbfaa82ce4d33aad1ba80787cc611a6730aafb38f8244954c4d0d50ae7a1"
    assert m.FILES[m.INFO_NAME]["sha256"] != m.FILES[m.MEANS_NAME]["sha256"]
