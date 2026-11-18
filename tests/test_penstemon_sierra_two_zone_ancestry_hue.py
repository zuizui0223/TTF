"""Known-truth source SHA, zone leakage and fixed baseline checks."""
from pathlib import Path
import importlib.util
import math
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_sierra_two_zone_ancestry_hue_20261118.py"
spec=importlib.util.spec_from_file_location("sierra_ttf",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def example_raw():
    h="newID,pop,class,subpopulation,modified_hue,K1,K2,region\n"
    return (h+"1,1,hybrid,hybrids_virginialakes,1.0,0.0,1.0,sierras\n"+
           "2,1,hybrid,hybrids_virginialakes,3.0,1.0,0.0,sierras\n"+
           "3,2,hybrid,hybrids_gemlakes,3.0,0.0,1.0,sierras\n"+
           "4,2,hybrid,hybrids_gemlakes,5.0,1.0,0.0,sierras\n"+
           "5,1,davidsonii,davidsonii_virginialakes,9.0,0.7,0.3,sierras\n").encode()


def synthetic_audit(monkeypatch, raw):
    monkeypatch.setattr(m,"EXPECTED_SOURCE_ROWS",5)
    monkeypatch.setattr(m,"MIN_HYBRIDS_PER_ZONE",2)
    monkeypatch.setattr(m,"SOURCE_BLOB",m.git_blob_sha(raw))
    return m.parse_source(raw)


def test_exact_source_grain_and_two_zones(monkeypatch):
    s=synthetic_audit(monkeypatch,example_raw())
    assert s["source_count"]==5
    assert s["complete_case_hybrid_counts"]=={
       "hybrids_gemlakes":2,"hybrids_virginialakes":2
    }
    assert s["excluded_nonhybrids_and_unassigned"]==1
    assert len(s["rows"])==4


def test_train_only_ols_and_baseline_direction_not_flipped(monkeypatch):
    monkeypatch.setattr(m,"MIN_HYBRIDS_PER_ZONE",2)
    train=[{"x":0.0,"y":1.0},{"x":0.5,"y":2.0},{"x":1.0,"y":3.0}]
    test=[{"x":0.0,"y":3.0},{"x":0.5,"y":4.0},{"x":1.0,"y":5.0}]
    a,b=m.ols_fit([r["x"] for r in train],[r["y"] for r in train])
    assert math.isclose(a,1.0) and math.isclose(b,2.0)
    d=m.one_fold(train,test,"train","heldout")
    assert d["heldout_sse_skill_vs_training_mean"]<1
    assert d["fitted_training_ancestry_slope"]==2.0
    assert d["heldout_n_individuals"]==3


def test_wrong_original_blob_or_duplicate_id_stops(monkeypatch):
    raw=example_raw()
    with pytest.raises(ValueError,match="SHA mismatch"):
        m.parse_source(raw)
    duplicate=raw.replace(b"5,1,davidsonii",b"1,1,davidsonii")
    monkeypatch.setattr(m,"EXPECTED_SOURCE_ROWS",5)
    monkeypatch.setattr(m,"MIN_HYBRIDS_PER_ZONE",2)
    monkeypatch.setattr(m,"SOURCE_BLOB",m.git_blob_sha(duplicate))
    with pytest.raises(ValueError,match="duplicate"):
        m.parse_source(duplicate)


def test_zero_training_variation_is_not_predicted():
    with pytest.raises(ValueError,match="no ancestry K1 variation"):
        m.ols_fit([0.2]*21,[1.0]*21)


def test_training_zone_labels_were_selected_without_hue():
    assert m.ZONES==("hybrids_virginialakes","hybrids_gemlakes")
    assert m.MIN_HYBRIDS_PER_ZONE==20
    assert m.EXPECTED_SOURCE_ROWS==421
