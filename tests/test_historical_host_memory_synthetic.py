"""No genetic response is read by the historical host-memory synthetic test family."""
import itertools

import numpy as np
import pytest

from scripts.qualify_historical_host_memory_synthetic import (
    assess, endpoint_safe_ibd, prepare_species, simulate_betas, statistic, wilson
)


def toy():
    rng=np.random.default_rng(123)
    nodes={i: xyz for i,xyz in enumerate(rng.normal(size=(7,3))*100)}
    edges={i: e for i,e in enumerate(itertools.combinations(range(7),2))}
    rows=[]
    for i in edges:
        rows.append({"edge_index":str(i),
                     "M_host":str(rng.normal()),
                     "M_self":str(rng.normal()),
                     "D_host_0BP":str(rng.normal())})
    return prepare_species("Toy species", rows, nodes, edges)


def test_endpoint_safe_ibd_matches_explicit_leave_two_nodes_out():
    p=toy()
    rng=np.random.default_rng(33)
    raw=rng.normal(size=(len(p["x"]),3))
    fast=endpoint_safe_ibd(raw,p)
    for i,(a,b) in enumerate(zip(p["left"],p["right"])):
        use=np.array([(a not in (x,y) and b not in (x,y))
                      for x,y in zip(p["left"],p["right"])])
        assert np.sum(use)>=5
        coeff=np.linalg.lstsq(
            np.column_stack([np.ones(np.sum(use)),p["x"][use]]),
            raw[use],rcond=None)[0]
        slow=raw[i]-(coeff[0]+coeff[1]*p["x"][i])
        assert fast[i]==pytest.approx(slow,abs=1e-10)


def test_fixed_seed_and_frozen_groups_are_reproducible():
    p=toy()
    a=simulate_betas(p,"reference",7,block_size=64)
    b=simulate_betas(p,"reference",7,block_size=64)
    c=simulate_betas(p,"positive",7,block_size=64)
    np.testing.assert_array_equal(a,b)
    assert len(a)==7 and np.all(np.isfinite(a))
    assert not np.array_equal(a,c)


def test_studentized_species_equal_weight_and_error_guard():
    arr=np.array([[0.1,0.3],[0.2,0.2],[0.5,0.1]])
    obs=statistic(arr)
    exp=arr.mean(axis=0)/(arr.std(axis=0,ddof=1)/np.sqrt(3))
    np.testing.assert_allclose(obs,exp)
    with pytest.raises(RuntimeError,match="multiple independent"):
        statistic(arr[:1])


def test_mc_gate_exact_one_sided_tail_and_wilson_bounds():
    reference=np.arange(1999,dtype=float)
    # All evaluation null worlds are extreme; no false-positive rejection.
    null_eval=np.full(500,-10.0)
    positive=np.full(500,3000.0)
    result=assess(reference,null_eval,positive)
    assert result["null_rejections"]==0
    assert result["positive_rejections"]==500
    assert result["type1_gate_pass"] is True
    assert result["positive_power_gate_pass"] is True
    assert result["decision"]=="PASS_TO_CONFIRMATORY_CHARACTER_MASK_ONLY"
    assert wilson(0,500)[1]<0.10
    with pytest.raises(RuntimeError,match="frozen Monte Carlo"):
        assess(reference[:1998],null_eval,positive)


def test_singular_ibd_disjoint_fit_fails_closed():
    p=toy()
    rows=[{"edge_index":str(i),"M_host":str(i),"M_self":str(i*i),
           "D_host_0BP":str(i**3)} for i in range(21)]
    identical={i:np.zeros(3) for i in range(7)}
    with pytest.raises(RuntimeError,match="invariant"):
        prepare_species("Toy bad",rows,identical,
                        {i:e for i,e in enumerate(itertools.combinations(range(7),2))})
