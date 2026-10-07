import numpy as np

from ttf.historical_host_connectivity import (
    RESISTANCE_EPSILON,
    historical_unique_fraction,
    host_union,
    information_gate,
    resistance_from_union,
)


def test_host_union_is_probabilistic_union():
    x=np.array([[0.2,0.0],[0.5,1.0]])
    got=host_union(x)
    assert np.allclose(got,[0.6,1.0])


def test_resistance_floor_is_frozen():
    got=resistance_from_union(np.array([0.0,0.5]))
    assert RESISTANCE_EPSILON==1e-6
    assert np.isclose(got[0],-np.log(1e-6))
    assert got[0]>got[1]


def test_historical_unique_fraction_detects_redundancy_and_independence():
    now=np.linspace(-2,2,101)
    assert historical_unique_fraction(now,now)<1e-12
    lgm=np.sin(np.linspace(0,8*np.pi,101))
    assert historical_unique_fraction(lgm,now)>0.9


def test_information_gate_equal_weights_species_not_edges():
    by={}
    for i in range(100):
        n=20 if i else 1000
        now=np.linspace(-2,2,n)
        lgm=np.sin(np.linspace(0,8*np.pi,n)+i*0.01)
        by[f"s{i:03d}"]=(lgm,now)
    out=information_gate(by)
    assert out["status"]=="PASS_TO_HISTORICAL_HOST_SYNTHETIC_QUALIFICATION"
    assert out["species"]==100
    assert out["max_single_species_signal_share"]<0.02
