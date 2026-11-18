"""Synthetic Penstemon 2023 source-only tip/key coverage test."""
import importlib.util
from pathlib import Path
import pytest

# Only the source-specific workflow installs DendroPy. The frozen TTF-Q
# base test environment intentionally need not acquire a phylogeny parser.
pytest.importorskip("dendropy")

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


def test_source_iqtree_quoted_internal_annotations_keep_tip_ids():
    raw=b"((CAW123_1_S454:0.2,cardBL_1_S443:0.2)'[q1=0.61;q2=0.23;q3=0.16]'/31.6:0.1,CAW98_1_S441:0.5);"
    cleaned=m.strip_2023_internal_node_figure_annotations(raw.decode())
    assert "q1" not in cleaned
    assert ":0.1" in cleaned
    assert m.tip_labels(raw)==["CAW123_1_S454","CAW98_1_S441","cardBL_1_S443"]
