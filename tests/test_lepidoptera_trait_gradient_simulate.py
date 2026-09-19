import numpy as np
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_trait_gradient_simulate import frozen_seed,simulate_trait_gradient_world

def geometries():
    out={}
    for j in range(8):
        x=np.linspace(0,10,12); xyz=np.column_stack([x,np.sin(x+j),np.cos(x+j)])
        out[f"s{j}"]=prepare_density_scaled_genetic_geometry(xyz,neighbor_fraction=.15)
    return out

def test_seed_namespace_is_cell_disjoint():
    assert len({frozen_seed(20260920,c,0) for c in ("private","trait_gradient_positive","geometry_confounded_trap")})==3

def test_world_is_deterministic_and_complete():
    g=geometries(); t={n:i/7 for i,n in enumerate(sorted(g))}; q={n:1-i/7 for i,n in enumerate(sorted(g))}
    a=simulate_trait_gradient_world(g,t,q,cell="trait_gradient_positive",seed=3)
    b=simulate_trait_gradient_world(g,t,q,cell="trait_gradient_positive",seed=3)
    assert set(a.edge_response)==set(g)
    for n in g: assert np.array_equal(a.edge_response[n],b.edge_response[n])

def test_private_ignores_trait_and_geometry_scores():
    g=geometries(); z={n:0. for n in g}; o={n:1. for n in g}
    a=simulate_trait_gradient_world(g,z,z,cell="private",seed=9)
    b=simulate_trait_gradient_world(g,o,o,cell="private",seed=9)
    for n in g: assert np.array_equal(a.edge_response[n],b.edge_response[n])
