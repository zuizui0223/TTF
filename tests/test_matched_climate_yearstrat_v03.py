"""Outcome-blind tests for all-398 yearstrat matched climate v0.3."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pytest
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/"scripts"))
from host_geometry_yearstrat_v03 import load_host_snapshot,EXPECTED_SUMMARY_SHA256
from matched_climate_yearstrat_v03 import matched_metrics, match_sampling_indices,residual_fraction

def test_frozen_518_plant_taxa_and_v02_checksums():
    root=HERE/"benchmarks/exploratory/lgm_host_reference_yearstrat_v02"
    assert hashlib.sha256((root/"summary.json").read_bytes()).hexdigest()==EXPECTED_SUMMARY_SHA256
    completed,plants,review,failures,plan=load_host_snapshot(root)
    assert (len(completed),len(plants),len(review),len(failures))==(518,506,1,11)
    assert sum(map(len,plants.values()))==296697
    assert sum(bool(v) for v in plants.values())==498
    assert plan["binding"]["not_equivalent_to_v01"] is True

def test_host_v01_cannot_be_loaded_as_v02():
    with pytest.raises(ValueError):
        load_host_snapshot(HERE/"benchmarks/exploratory/lgm_host_reference_v01")

def test_both_cloud_reference_counts_are_equal_and_seed_deterministic():
    a,b=match_sampling_indices(191,83,"Butterfly example",7)
    a2,b2=match_sampling_indices(191,83,"Butterfly example",7)
    assert np.array_equal(a,a2) and np.array_equal(b,b2)
    assert len(a)==len(b)==83 and len(np.unique(a))==83

def test_climate_geographic_controls_reduce_or_preserve_variance():
    rng=np.random.default_rng(404)
    current=rng.normal(size=(90,5,4))
    lgm=current+rng.normal(size=(90,5,4),scale=.3)
    host=rng.normal(size=(110,4))
    self_=rng.normal(size=(100,4))
    distance=np.log1p(np.arange(90,dtype=float)+1)
    v=matched_metrics(current,lgm,host,self_,distance,"Climate comparison",0)
    assert 0<=v["unique_history_fraction"]<=v["unique_before_geography"]+1e-12<=1+1e-12
    assert v["reference_n"]==100

def test_identical_host_and_self_cloud_is_rank_deficient_not_positive_evidence():
    rng=np.random.default_rng(22)
    current=rng.normal(size=(50,5,4))
    lgm=current+rng.normal(size=(50,5,4),scale=.2)
    cloud=rng.normal(size=(70,4))
    out=matched_metrics(current,lgm,cloud,cloud,np.arange(50)/3,"Identical",0)
    assert out["relative_host_minus_self"]==pytest.approx(0.,abs=1e-12)
    assert out["unique_history_fraction"] is None

def test_constant_host_history_has_no_identifiable_information():
    rng=np.random.default_rng(23)
    controls=rng.normal(size=(40,4))
    assert residual_fraction(np.ones(40),controls) is None
