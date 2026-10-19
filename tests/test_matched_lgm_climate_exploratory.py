import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from matched_climate_diagnostic import match_sampling_indices,matched_metrics,residual_fraction,analogue,edge_dist_km
from reference_helpers_min import join_hosts,filter_host_points,normalize_id

def test_equal_count_seeded_without_replacement():
    h,s=match_sampling_indices(200,37,'Species test',0)
    assert len(h)==len(s)==37 and len(set(h))==37
    assert np.array_equal(h,match_sampling_indices(200,37,'Species test',0)[0])
    assert np.array_equal(s,np.arange(37))
    assert not np.array_equal(h,match_sampling_indices(200,37,'Species test',1)[0])

def test_identical_climate_references_have_identical_delta():
    rng=np.random.default_rng(5)
    cur=rng.normal(size=(15,5,4));lgm=cur+rng.normal(scale=.2,size=(15,5,4))
    cloud=rng.normal(size=(70,4))
    result=matched_metrics(cur,lgm,cloud,cloud,np.linspace(.2,2.,15),'A b',0)
    assert abs(result['relative_host_minus_self'])<1e-12
    assert result['reference_n']==70

def test_residual_nuisance_projection_is_bounded():
    rng=np.random.default_rng(3)
    ctrl=rng.normal(size=(50,4));y=3*ctrl[:,0]+rng.normal(size=50)
    q=residual_fraction(y,ctrl)
    assert 0<q<1
    assert residual_fraction(np.ones(50),ctrl) is None

def test_missing_reference_is_not_promoted_to_valid_species():
    with pytest.raises(ValueError): match_sampling_indices(29,100,'A b',0)
    with pytest.raises(ValueError): match_sampling_indices(100,29,'A b',0)

def test_five_interior_positions_are_required():
    with pytest.raises(ValueError):
        matched_metrics(np.ones((4,4,4)),np.ones((4,4,4)),np.zeros((50,4)),np.zeros((50,4)),np.ones(4),'Species',0)

def test_nonfinite_climate_fails_closed():
    with pytest.raises(ValueError): analogue(np.ones((3,4)),np.array([[np.nan,0,0,0]]))

def test_geodesic_km_control():
    assert abs(edge_dist_km([[0,0],[0,1]])-111.195)<.01

def test_host_source_deduplication_and_spatial_exclusion():
    raw={'10':[{'gbif_id':1,'latitude':0.,'longitude':0.,'uncertainty_m':0},
               {'gbif_id':2,'latitude':1.,'longitude':0.,'uncertainty_m':0}],
         '20':[{'gbif_id':2,'latitude':1.,'longitude':0.,'uncertainty_m':0}]}
    v=join_hosts(['10','20'],raw)
    assert len(v)==2
    r=filter_host_points(v,np.array([[0.,0.]]),buffer_km=50,thinning_km=10)
    assert [x['gbif_id'] for x in r]==[2]

def test_host_rank_id_exact_not_float_rounded():
    assert normalize_id('4e+05')=='400000'
    with pytest.raises(ValueError): normalize_id('1.2')
