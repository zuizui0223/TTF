"""Synthetic source-only Costus clade holdout design tests."""
import importlib.util
from pathlib import Path
import pytest
dendropy = pytest.importorskip("dendropy")

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_clade_holdout_source_gate_20261115.py"
s=importlib.util.spec_from_file_location("costus_tree_clade",P)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def small_tree():
    return dendropy.Tree.get(
        data="(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);",
        schema="newick",rooting="force-rooted",preserve_underscores=True,
    )

def test_fixed_tree_split_is_disjoint_without_scent_outcome():
    root=small_tree()
    allnames={x.taxon.label.lower() for x in root.leaf_node_iter()}
    for k in [2,4]:
        groups=m.split_largest_subclades(root,k)
        assert len(groups)==k
        assert set().union(*(set(g) for g in groups))==allnames
        assert sum(map(len,groups))==8
        assert not any(set(groups[i])&set(groups[j]) for i in range(len(groups)) for j in range(i+1,len(groups)))

def test_without_two_pollinator_classes_clade_gate_stops():
    labels={k:("Bee" if k in "abcd" else "Hummingbird") for k in "abcdefgh"}
    groups=m.split_largest_subclades(small_tree(),2)
    d=m.clade_folds(groups,labels)
    assert d["status"]=="INSUFFICIENT_POLLINATOR_OR_CLADE_SUPPORT"
    assert d["n_folds"]==2

def test_fixed_k_and_min_support_are_source_declared():
    assert m.CANDIDATE_CLUSTERS==[4,3,2]
    assert m.MIN_TEST_TIPS_PER_CLADE==3
    assert m.MIN_TRAIN_PER_SOURCE_POLLINATOR_CLASS==4
