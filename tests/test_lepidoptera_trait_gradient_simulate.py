import numpy as np
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_trait_gradient_simulate import frozen_seed,simulate_trait_gradient_world

def geometries():
    out={}
    for j in range(8):
        x=np.linspace(0,10,12)
        xyz=np.column_stack([x,np.sin(x+j),np.cos(x+j)])
        out[f"s{j}"]=prepare_density_scaled_genetic_geometry(xyz,neighbor_fraction=.15)
    return out

def kernels(names):
    n=len(names)
    trait=np.fromfunction(lambda i,j: np.exp(-np.abs(i-j)/2.0),(n,n))
    geom=np.fromfunction(lambda i,j: np.exp(-((i-j)/3.0)**2),(n,n))
    return trait,geom

def test_seed_namespace_is_cell_disjoint():
    assert len({frozen_seed(20260920,c,0) for c in ("private","trait_gradient_positive","geometry_confounded_trap")})==3

def test_world_is_deterministic_and_complete():
    g=geometries(); names=tuple(sorted(g)); t,q=kernels(names)
    a=simulate_trait_gradient_world(g,names,t,q,cell="trait_gradient_positive",seed=3)
    b=simulate_trait_gradient_world(g,names,t,q,cell="trait_gradient_positive",seed=3)
    assert set(a.edge_response)==set(g)
    for n in g: assert np.array_equal(a.edge_response[n],b.edge_response[n])

def test_private_ignores_supplied_similarity_kernels():
    g=geometries(); names=tuple(sorted(g)); n=len(names)
    eye=np.eye(n); ones=np.ones((n,n))
    a=simulate_trait_gradient_world(g,names,eye,eye,cell="private",seed=9)
    b=simulate_trait_gradient_world(g,names,ones,ones,cell="private",seed=9)
    for name in g: assert np.array_equal(a.edge_response[name],b.edge_response[name])

def test_positive_and_trap_use_distinct_kernels():
    g=geometries(); names=tuple(sorted(g)); t,q=kernels(names)
    a=simulate_trait_gradient_world(g,names,t,q,cell="trait_gradient_positive",seed=5)
    b=simulate_trait_gradient_world(g,names,t,q,cell="geometry_confounded_trap",seed=5)
    assert any(not np.array_equal(a.edge_response[n],b.edge_response[n]) for n in names)
