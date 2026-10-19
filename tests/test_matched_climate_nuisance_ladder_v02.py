"""Tests for the post-v0.1 exploratory matched-reference nuisance ladder."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from matched_climate_nuisance_ladder_v02 import matched_metrics,residual_fraction,match_sampling_indices

def test_geographic_control_cannot_increase_retained_variance():
    rng=np.random.default_rng(344)
    current=rng.normal(size=(66,5,4));lgm=current+rng.normal(scale=.4,size=(66,5,4))
    host=rng.normal(size=(89,4)); insect=rng.normal(size=(91,4))
    for draw in range(5):
        q=matched_metrics(current,lgm,host,insect,np.linspace(.1,4.,66),'Species test',draw)
        assert q['unique_before_geography'] is not None
        assert q['unique_history_fraction'] is not None
        assert 0 <= q['unique_history_fraction'] <= q['unique_before_geography']+1e-12 <= 1+1e-12

def test_same_references_produce_zero_host_minus_insect_change():
    rng=np.random.default_rng(12)
    current=rng.normal(size=(50,5,4))
    shift=rng.normal(size=current.shape)*.3
    cloud=rng.normal(size=(70,4))
    q=matched_metrics(current,current+shift,cloud,cloud,np.arange(50)/10,'Identical clouds',0)
    assert q['relative_host_minus_self']==pytest.approx(0.,abs=1e-12)
    assert q['unique_before_geography'] is not None
    assert q['unique_history_fraction']==pytest.approx(0.,abs=1e-12)

def test_added_geo_removes_designed_nuisance():
    rng=np.random.default_rng(44);n=100
    x=rng.normal(size=(n,3))
    geo=rng.normal(size=n)
    y=3*geo+.05*rng.normal(size=n)
    before=residual_fraction(y,x)
    after=residual_fraction(y,np.column_stack([x,geo]))
    assert before>0.90 and after<0.01

def test_matched_cloud_subsampling_is_deterministic_and_preserved():
    hi,si=match_sampling_indices(89,91,'Test taxon',2)
    assert np.array_equal(hi,match_sampling_indices(89,91,'Test taxon',2)[0])
    assert np.array_equal(si,match_sampling_indices(89,91,'Test taxon',2)[1])
    assert len(hi)==len(si)==89 and len(set(hi))==89 and len(set(si))==89
