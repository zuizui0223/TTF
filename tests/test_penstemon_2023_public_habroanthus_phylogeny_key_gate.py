"""Synthetic Penstemon 2023 source-only tip/key coverage test."""
import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_2023_public_habroanthus_phylogeny_key_gate_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_2023_tips",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_public_tree_tips_are_unique_and_exact():
    a=b"((CAW123_1_S454:0.2,cardBL_1_S443:0.2):0.1,CAW98_1_S441:0.5);"
    assert m.tip_labels(a)==["CAW123_1_S454","CAW98_1_S441","cardBL_1_S443"]
    with pytest.raises(Exception,match="uplicate|epeated"):
        m.tip_labels(b"(a:1,a:2);")


def test_groupmap_only_exact_sample_labels():
    rows=b"cardBL_1_S443\tcardinalis\nCAW123_1_S454\teatonii\n"
    assert m.group_map(rows)=={"cardBL_1_S443":"cardinalis","CAW123_1_S454":"eatonii"}
    with pytest.raises(ValueError,match="duplicate"):
        m.group_map(rows+rows.splitlines(keepends=True)[0])
