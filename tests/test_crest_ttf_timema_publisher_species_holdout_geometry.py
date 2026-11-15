"""Synthetic source-geometry tests; no Timema biology response opened."""
import importlib.util
from pathlib import Path
import pytest

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_timema_publisher_species_holdout_geometry_20261115.py"
spec=importlib.util.spec_from_file_location("ttf_timema_species_gate",FILE)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def make_rows():
    return [
      {"pair_id":"1","species_1":"A","species_2":"A","population_pair":"p1-p2"},
      {"pair_id":"2","species_1":"A","species_2":"B","population_pair":"p3-p4"},
      {"pair_id":"3","species_1":"B","species_2":"B","population_pair":"p5-p6"},
      {"pair_id":"4","species_1":"C","species_2":"D","population_pair":"p7-p8"},
      {"pair_id":"5","species_1":"C","species_2":"C","population_pair":"p9-p10"},
      {"pair_id":"6","species_1":"D","species_2":"D","population_pair":"p11-p12"},
    ]

def test_split_is_species_disjoint_not_dyad_count():
    g=m.split_geometry(make_rows())
    assert g["species_count"]==4
    assert g["within_species_pair_count"]==4
    assert g["between_species_pair_count"]==2
    assert g["partitions_with_nonempty_species_disjoint_train_and_test"]>=1
    p=g["best_balanced_species_split_source_only"]
    assert set(p["train_species"]).isdisjoint(set(p["test_species"]))
    assert p["train_dyads"]>0 and p["test_dyads"]>0
    assert (p["train_dyads"]+p["test_dyads"]+p["cross_boundary_dyads_excluded"])==6
    assert g["comparisons_are_not_independent_evolutionary_replicates"] is True

def test_species_count_fail_closed():
    with pytest.raises(ValueError,match="species"):
        m.split_geometry(make_rows()[:1])
