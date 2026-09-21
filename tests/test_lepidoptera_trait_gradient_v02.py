import numpy as np
import pytest

from ttf.lepidoptera_trait_gradient_v02 import (
    EVALUATION_TAG,
    REFERENCE_TAG,
    envelope_pvalues,
    frozen_v02_seed,
    target_equal_mean_correlation,
    target_fe_geometry_residual,
    wilson_interval,
)


def test_v02_seed_namespaces_are_disjoint_and_locked():
    a=frozen_v02_seed(20260920,REFERENCE_TAG,"private",0)
    b=frozen_v02_seed(20260920,EVALUATION_TAG,"private",0)
    assert a!=b
    assert a==frozen_v02_seed(20260920,REFERENCE_TAG,"private",0)
    with pytest.raises(ValueError):
        frozen_v02_seed(20260920,"spent-v01","private",0)


def test_target_fe_residual_is_orthogonal_to_centered_geometry():
    target=np.repeat(np.arange(8),8)
    g=np.arange(len(target),dtype=float)
    x=np.column_stack([g,g**2,np.sin(g),np.cos(g),np.sqrt(g+1)])
    y=0.2*g-0.01*g**2+np.tile(np.linspace(-1,1,8),8)
    r=target_fe_geometry_residual(y,x,target)
    for label in np.unique(target):
        idx=np.flatnonzero(target==label)
        assert abs(float(np.mean(r[idx])))<1e-10
    xc=x.copy()
    for label in np.unique(target):
        idx=np.flatnonzero(target==label)
        xc[idx]-=xc[idx].mean(axis=0)
    assert np.max(np.abs(xc.T@r))<1e-7


def test_target_equal_mean_correlation_recovers_positive_relation():
    target=np.repeat(np.arange(10),6)
    x=np.tile(np.linspace(-1,1,6),10)
    y=2*x
    stat,per=target_equal_mean_correlation(y,x,target)
    assert np.isclose(stat,1.0)
    assert len(per)==10


def test_envelope_uses_least_favourable_null():
    private=np.linspace(-.2,.2,99)
    geometry=np.linspace(0,.4,99)
    p=envelope_pvalues(np.asarray([.1,.3]),private,geometry)
    assert p.shape==(2,)
    assert np.all((p>0)&(p<=1))
    assert p[0]>p[1]


def test_wilson_interval_bounds_probability():
    lo,hi=wilson_interval(25,500)
    assert 0<=lo<.05<hi<=1
