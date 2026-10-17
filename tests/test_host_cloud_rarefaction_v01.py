"""Synthetic, completely genotype-free checks for host cloud rarefaction pilot."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

TARGET=Path(__file__).resolve().parents[1]/'scripts/pilot_host_cloud_size_rarefaction.py'


def module():
    spec=importlib.util.spec_from_file_location('host_rarefaction',TARGET)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_exact_seed_and_reproducible_without_replacement():
    m=module()
    a=m.fixed_subsample_indices(665,26,64)
    b=m.fixed_subsample_indices(665,26,64)
    np.testing.assert_array_equal(a,b)
    assert a.shape==(64,26)
    assert a.min()>=0 and a.max()<665
    assert all(len(set(row))==26 for row in a)
    assert len(set(tuple(sorted(row)) for row in a))==64


def test_sample_n_never_adapts_to_available_cloud():
    m=module()
    with pytest.raises(ValueError,match='665'):
        m.fixed_subsample_indices(25,26,64)
    with pytest.raises(ValueError,match='64'):
        m.fixed_subsample_indices(665,26,63)
    with pytest.raises(ValueError,match='665'):
        m.fixed_subsample_indices(664,26,64)


def test_summary_uses_fixed_self_mean_and_full_count_not_subset():
    m=module()
    host=np.array([[1,2,3],[5,6,7]],float)
    # Each row is a replicate across original 3 edges, self is fixed per edge.
    self_change=np.array([3.,3.,3.])
    result=m.summarize_rarefied(host,self_change)
    assert result['rarefied_host_delta_cost_km']['means_per_draw']==[2.0,6.0]
    assert result['rarefied_host_minus_self_delta_cost_km']['means_per_draw']==[-1.,3.]
    assert result['rarefied_host_minus_self_delta_cost_km']['fraction_positive_draws']==0.5
    assert result['draw_n']==2
    assert result['original_edge_n']==3


def test_invalid_shape_and_nonfinite_costs_fail_closed():
    m=module()
    with pytest.raises(ValueError,match='nonfinite'):
        m.summarize_rarefied(np.array([[1,np.nan]]),np.ones(2))
    with pytest.raises(ValueError,match='shape'):
        m.summarize_rarefied(np.ones((2,3)),np.ones(4))


def test_original_full_cloud_and_self_costs_must_match_preexisting_edges():
    m=module()
    original=[{'edge_index':str(i), 'R_host_0BP_km':str(1+i),
        'R_host_LGM_km':str(4+i),'R_self_0BP_km':str(2+i),
        'R_self_LGM_km':str(5+i)} for i in range(3)]
    scores={('host','0BP'):np.array([1.,2.,3.]), ('host','LGM'):np.array([4.,5.,6.]),
            ('self','0BP'):np.array([2.,3.,4.]), ('self','LGM'):np.array([5.,6.,7.])}
    assert m.verify_original_edges(scores,original,atol=1e-8)==0.0
    scores[('host','LGM')][1]+=0.01
    with pytest.raises(RuntimeError,match='baseline replay'):
        m.verify_original_edges(scores,original,atol=1e-8)


def test_no_claim_about_an_actual_biological_barrier_in_receipt():
    m=module()
    assert m.SPECIES=='Bungalotis erythus'
    assert m.HOST_ORIGINAL_N==665
    assert m.RAREFIED_N==26
    assert m.REPLICATES==64
