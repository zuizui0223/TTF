"""Synthetic unrooted Costus distance block and no-outcome source gate."""
import importlib.util
from pathlib import Path
import pytest
dendropy = pytest.importorskip("dendropy")

F=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_costus_unrooted_neighborhood_gate_20261115.py"
spec=importlib.util.spec_from_file_location("costus_unrooted",F)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def tree():
    return dendropy.Tree.get(data="((A:1,B:1):1,(C:1,D:1):1);",schema="newick",rooting="force-unrooted")

def test_symmetric_patristic_pairwise_distance_only():
    data=tree()
    labels=["a","b","c","d"]
    distances=m.exact_patristic_distances(data,labels)
    assert len(distances)==6
    assert distances[("a","b")]==2.0
    assert distances[("a","d")]==4.0

def test_average_linkage_groups_deterministically_without_root():
    labels=["a","b","c","d"]
    distances=m.exact_patristic_distances(tree(),labels)
    a=m.avg_linkage_partitions(labels,distances)
    b=m.avg_linkage_partitions(labels,distances)
    assert a==b
    assert len(a[4])==4 and len(a[3])==3 and len(a[2])==2
    assert sorted(map(tuple,a[2]))==[("a","b"),("c","d")]

def test_fixed_grouping_family():
    assert m.K_ORDER==[4,3,2]
    assert m.DISTANCE_ROUND_DIGITS==12
