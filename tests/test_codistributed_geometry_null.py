import numpy as np

from ttf.codistributed_geometry_null import (
    GeometryNullCenter,
    center_dyad_scores,
    fit_geometry_null_center,
    frozen_uint64_seed,
    least_favourable_beta_pvalues,
    prepare_fixed_dyad_transfer_cache,
    score_fixed_dyad_transfer_cache,
    source_only_transfer_scores_batch,
)
from ttf.core import SpeciesEdges
from ttf.relational_genetic_empirical import RelationalSpeciesResponse, source_only_transfer_scores


def _edges(name, x):
    x=np.asarray(x,dtype=float)
    nodes=np.asarray([[0,1],[1,2],[2,3],[0,3]],dtype=np.int64)
    start=x[nodes[:,0]]
    end=x[nodes[:,1]]
    return SpeciesEdges(
        species=name,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5*(start+end),
        length=np.linalg.norm(end-start,axis=1),
        turnover=np.zeros(len(nodes),dtype=float),
    )


def test_seed_namespaces_are_deterministic_and_separate():
    a=frozen_uint64_seed(20261003,"center","A2",7)
    b=frozen_uint64_seed(20261003,"center","A2",7)
    c=frozen_uint64_seed(20261003,"reference","A2",7)
    d=frozen_uint64_seed(20261003,"center","A2",8)
    assert a==b
    assert len({a,c,d})==3
    assert 0 <= a < 2**64


def test_geometry_null_center_weights_amplitude_cells_equally():
    cells={
        0.5:np.asarray([[0.0,2.0],[2.0,4.0]]),
        3.0:np.asarray([[10.0],[20.0]]),
    }
    center=fit_geometry_null_center(cells)
    # Cell means are [1,3] and [10,20], then amplitudes receive equal weight.
    np.testing.assert_allclose(center.mu0,[5.5,11.5])
    assert center.amplitudes==(0.5,3.0)
    assert center.worlds_per_amplitude==(2,1)
    x=np.asarray([[6.5,4.5],[12.5,10.5]])
    np.testing.assert_allclose(center_dyad_scores(x,center),[[1,-1],[1,-1]])


def test_batch_single_source_operator_matches_scalar_reference():
    source=_edges("s",[[0,0],[1,0],[2,0],[3,0]])
    target1=_edges("t1",[[0,0.2],[1,0.2],[2,0.2],[3,0.2]])
    target2=_edges("t2",[[0,-0.3],[1,-0.3],[2,-0.3],[3,-0.3]])
    edge_map={"s":source,"t1":target1,"t2":target2}
    pairs=[("s","t1"),("s","t2")]
    train_worlds=np.asarray([
        [0.1,0.8],
        [0.9,0.2],
        [0.3,0.6],
        [0.7,0.4],
    ])
    eval1=np.asarray([
        [0.2,0.7],
        [0.8,0.1],
        [0.4,0.9],
        [0.6,0.3],
    ])
    eval2=np.asarray([
        [0.7,0.4],
        [0.1,0.8],
        [0.9,0.2],
        [0.3,0.6],
    ])
    batch=source_only_transfer_scores_batch(
        edge_map,
        {"s":train_worlds},
        {"t1":eval1,"t2":eval2},
        pairs,
        bandwidth_km=2.0,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )
    expected=np.empty_like(batch)
    for j in range(2):
        responses={
            "s":RelationalSpeciesResponse("s",train_worlds[:,j],train_worlds[:,j]),
            "t1":RelationalSpeciesResponse("t1",eval1[:,j],eval1[:,j]),
            "t2":RelationalSpeciesResponse("t2",eval2[:,j],eval2[:,j]),
        }
        expected[:,j]=source_only_transfer_scores(
            edge_map,responses,pairs,
            bandwidth_km=2.0,
            prior_strength=0.25,
            prior_mean=0.0,
            segment_points=5,
        )
    np.testing.assert_allclose(batch,expected,rtol=0,atol=1e-14)




def test_cached_batch_operator_matches_uncached_batch():
    # The cache also freezes endpoint-safe IBD geometry for later synthetic use,
    # so use a toy graph that genuinely satisfies the >=3 disjoint-edge rule.
    def dense_edges(name, y):
        x=np.column_stack((np.arange(6,dtype=float),np.full(6,float(y))))
        nodes=np.asarray(
            [(i,j) for i in range(6) for j in range(i+1,6)],
            dtype=np.int64,
        )
        start=x[nodes[:,0]]
        end=x[nodes[:,1]]
        return SpeciesEdges(
            species=name,
            nodes=nodes,
            start=start,
            end=end,
            midpoint=0.5*(start+end),
            length=np.linalg.norm(end-start,axis=1),
            turnover=np.zeros(len(nodes),dtype=float),
        )

    source=dense_edges("s",0.0)
    target1=dense_edges("t1",0.2)
    target2=dense_edges("t2",-0.3)
    edge_map={"s":source,"t1":target1,"t2":target2}
    pairs=[("s","t1"),("s","t2")]
    n=source.n_edges
    base=np.linspace(0.05,0.95,n)
    train={"s":np.column_stack((base,base[::-1]))}
    evaluation={
        "t1":np.column_stack((np.roll(base,1),np.roll(base[::-1],2))),
        "t2":np.column_stack((np.roll(base,3),np.roll(base[::-1],1))),
    }
    uncached=source_only_transfer_scores_batch(
        edge_map,train,evaluation,pairs,
        bandwidth_km=2.0,prior_strength=0.25,prior_mean=0.0,segment_points=5,
    )
    cache=prepare_fixed_dyad_transfer_cache(
        edge_map,pairs,
        bandwidth_km=2.0,prior_strength=0.25,prior_mean=0.0,segment_points=5,
        min_training_edges=3,
    )
    cached=score_fixed_dyad_transfer_cache(cache,train,evaluation)
    np.testing.assert_allclose(cached,uncached,rtol=0,atol=1e-14)

def test_least_favourable_beta_pvalue_takes_maximum_component():
    refs={
        "A0p5":np.linspace(-1,1,99),
        "A3":np.linspace(0,2,99),
    }
    p,labels=least_favourable_beta_pvalues(np.asarray([0.5]),refs)
    # A3 has more mass above 0.5 and is therefore least favourable.
    assert labels==("A3",)
    assert 0.0 < p[0] <= 1.0


def test_center_rejects_row_drift():
    center=GeometryNullCenter(np.zeros(3),(1.0,),(2,))
    try:
        center_dyad_scores(np.zeros((2,4)),center)
    except ValueError as exc:
        assert "one row" in str(exc)
    else:
        raise AssertionError("row drift was not rejected")
