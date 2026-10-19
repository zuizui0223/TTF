import numpy as np
import pytest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from same_niche_null_v04 import (draw_independent_clouds,
    score_analog_history,score_real_and_placebo)

def test_independent_disjoint_reference_clouds_and_repeatable():
    real,self_,placebo=draw_independent_clouds(131,200,'Epinotia tedella',11)
    again=draw_independent_clouds(131,200,'Epinotia tedella',11)
    assert len(real)==len(self_)==len(placebo)==65
    assert len(set(self_)&set(placebo))==0
    assert len(set(real))==65
    assert all(np.array_equal(a,b) for a,b in zip((real,self_,placebo),again))
    real2,_,_=draw_independent_clouds(131,200,'Epinotia tedella',12)
    assert not np.array_equal(real,real2)

@pytest.mark.parametrize('insect,host',[(59,100),(80,29),(0,100),(59,29)])
def test_insufficient_cloud_fails_closed(insect,host):
    with pytest.raises(ValueError):draw_independent_clouds(insect,host,'X y',0)

def test_identical_actual_and_placebo_has_exact_zero_contrast():
    r=np.random.default_rng(2)
    x=r.normal(size=(30,5,4));z=x+r.normal(scale=.25,size=x.shape)
    host=r.normal(size=(40,4));self_=r.normal(size=(40,4))
    geo=np.log1p(np.arange(30))
    actual,placebo=score_real_and_placebo(x,z,host,self_,host,geo)
    assert actual['unique_fraction']==pytest.approx(placebo['unique_fraction'],abs=1e-12)
    assert actual['residual_variance']==pytest.approx(placebo['residual_variance'],abs=1e-12)
    assert actual['host_current_mean']==pytest.approx(placebo['host_current_mean'],abs=1e-12)

def test_no_glacial_change_is_nonidentifiable_not_effect():
    r=np.random.default_rng(5)
    x=r.normal(size=(20,5,4));h=r.normal(size=(30,4));s=r.normal(size=(30,4))
    out=score_analog_history(x,x,h,s,np.linspace(0,2,20))
    assert out['unique_fraction'] is None and out['residual_variance'] is None

def test_raw_residual_fraction_normalization():
    r=np.random.default_rng(12)
    x=r.normal(size=(40,5,4));l=x+r.normal(size=x.shape,scale=.3)
    h=r.normal(size=(80,4));s=r.normal(size=(80,4))
    metric=score_analog_history(x,l,h,s,np.log1p(np.arange(40)))
    assert metric['unique_fraction'] is not None and metric['history_variance']>0
    assert metric['residual_variance']>=0
    assert metric['unique_fraction']==pytest.approx(
        metric['residual_variance']/metric['history_variance'],rel=1e-10)
