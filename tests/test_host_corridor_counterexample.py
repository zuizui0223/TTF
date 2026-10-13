"""Synthetic landscape thought experiment; it contains no TTF empirical genetics."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

S=importlib.util.spec_from_file_location('toy_corridor',Path(__file__).resolve().parents[1]/'scripts/probe_host_corridor_nonidentifiability.py')

def mod():
    m=importlib.util.module_from_spec(S)
    S.loader.exec_module(m)
    return m


def test_same_five_point_host_climate_average_different_2d_least_cost_route():
    m=mod()
    z=m.toy_counterexample()
    assert z['class']=='SYNTHETIC_NONIDENTIFIABILITY_EXAMPLE'
    assert z['five_point_straight_line_cost']['with_alternative_corridor']==20
    assert z['five_point_straight_line_cost']['without_alternative_corridor']==20
    assert z['least_cost_distance']['with_alternative_corridor'] < z['least_cost_distance']['without_alternative_corridor']
    assert z['not_an_empirical_LGM_result'] is True


def test_corridor_path_distance_uses_cell_cost_not_crossing_water():
    m=mod()
    cost=np.ones((5,5),float)
    assert m.illustrative_four_neighbor_least_cost(cost,(2,0),(2,4))==pytest.approx(4.0)
    cost[:,2]=np.inf
    with pytest.raises(ValueError,match='unreachable'):
        m.illustrative_four_neighbor_least_cost(cost,(2,0),(2,4))


def test_no_negative_resistance_is_silently_clamped():
    m=mod();cost=np.ones((3,3));cost[1,1]=-1.0
    with pytest.raises(ValueError,match='positive'):
        m.illustrative_four_neighbor_least_cost(cost,(0,0),(2,2))
