"""Synthetic climate-only reference-size tests, no genetic outcomes."""
import sys
from pathlib import Path
import numpy as np
import pytest
from rasterio.transform import from_origin
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_host_reference_size_rarefaction import (
    sampled_reference_indices, reference_corridor_delta, run_rarefaction,
    verify_baseline,
)


def tiny():
    r,c=8,8
    z=np.zeros((r,c,4),float)
    for i in range(r):
        for j in range(c): z[i,j]=[i/8,j/8,(i+j)/16,(i-j)/16]
    v=z+0.3
    clouds=np.array([[0.,0.,0.,0.],[0.1,0.5,0.1,-0.1],[0.7,0.8,0.6,-0.2],
                     [1,1,1,1],[0.4,0.6,0.5,0.2],[0.2,0.9,0.3,-0.5]])
    pts=np.array([[1,1],[2,4],[6,5]],int)
    edges=np.array([[0,1],[0,2],[1,2]],int)
    tf=from_origin(-86,12,0.1,0.1)
    return z,v,clouds,pts,edges,tf


def test_reference_draws_fixed_and_without_replacement():
    x=sampled_reference_indices(665,26,50,integer_seed=15)
    y=sampled_reference_indices(665,26,50,integer_seed=15)
    np.testing.assert_array_equal(x,y)
    assert x.shape==(50,26)
    assert all(len(set(row))==26 for row in x)
    assert np.all(x>=0) and np.all(x<665)
    assert not np.array_equal(x[0],x[1])


def test_illegal_sample_or_count_rejected():
    with pytest.raises(ValueError):sampled_reference_indices(3,4,5,integer_seed=1)
    with pytest.raises(ValueError):sampled_reference_indices(20,5,0,integer_seed=1)


def test_identical_present_and_lgm_has_zero_delta():
    x,_,h,p,e,t=tiny()
    present,lgm,delta=reference_corridor_delta(x,x,h,t,p,e)
    assert np.allclose(delta,0,atol=1e-12)
    assert np.allclose(present,lgm)


def test_removing_host_points_never_decreases_individual_raster_analogue_cost():
    x,y,h,p,e,t=tiny()
    from scipy.spatial import cKDTree
    original=cKDTree(h).query(x.reshape(-1,4),k=1)[0]
    sub=cKDTree(h[[0,2,4]]).query(x.reshape(-1,4),k=1)[0]
    assert np.all(sub>=original-1e-14)
    o=reference_corridor_delta(x,y,h,t,p,e)
    s=reference_corridor_delta(x,y,h[[0,2,4]],t,p,e)
    assert np.all(s[0]>=o[0]-1e-9)
    assert np.all(s[1]>=o[1]-1e-9)


def test_refit_draws_and_summary_are_deterministic_and_descriptive_only():
    x,y,h,p,e,t=tiny()
    baseline=reference_corridor_delta(x,y,h,t,p,e)
    check=verify_baseline(baseline,baseline,atol=1e-12)
    assert check==0
    run1=run_rarefaction(x,y,h,t,p,e,sample_size=3,repetitions=5,integer_seed=55)
    run2=run_rarefaction(x,y,h,t,p,e,sample_size=3,repetitions=5,integer_seed=55)
    assert run1==run2
    assert len(run1['host_delta_mean_by_draw'])==5
    assert run1['new_empirical_genetic_test'] is False
    assert run1['new_p_values']==0
    with pytest.raises(RuntimeError,match='replay mismatch'):
        verify_baseline(baseline,(baseline[0]+1,baseline[1],baseline[2]))
