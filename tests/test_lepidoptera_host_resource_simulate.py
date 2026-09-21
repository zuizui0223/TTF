from __future__ import annotations

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_host_resource_simulate import (
    prepare_host_resource_simulator,
    simulate_prepared_host_resource_world,
)


def geometries():
    out={}
    for j in range(7):
        x=np.linspace(0,12,14)
        xyz=np.column_stack([x,np.sin(x+j),np.cos(x-j)])
        out[f"s{j}"]=prepare_density_scaled_genetic_geometry(xyz,neighbor_fraction=.15)
    return out


def kernels(n):
    idx=np.arange(n)
    g=np.exp(-np.abs(idx[:,None]-idx[None,:])/2)
    b=np.exp(-((idx[:,None]-idx[None,:])/3)**2)
    h=np.eye(n)*.5+.5*np.exp(-np.abs(idx[:,None]-idx[None,:]))
    np.fill_diagonal(h,1)
    return g,b,h


def test_host_resource_simulator_is_deterministic_and_cell_specific():
    geo=geometries()
    names=tuple(sorted(geo))
    g,b,h=kernels(len(names))
    prepared=prepare_host_resource_simulator(
        geo,names,
        geometry_kernel=g,
        resource_breadth_kernel=b,
        host_resource_positive_kernel=h,
    )
    a=simulate_prepared_host_resource_world(prepared,cell="host_resource_gradient_positive",seed=7)
    c=simulate_prepared_host_resource_world(prepared,cell="host_resource_gradient_positive",seed=7)
    d=simulate_prepared_host_resource_world(prepared,cell="geometry_confounded_trap",seed=7)
    assert all(np.array_equal(a.edge_response[n],c.edge_response[n]) for n in names)
    assert any(not np.array_equal(a.edge_response[n],d.edge_response[n]) for n in names)


def test_private_cell_does_not_depend_on_shared_kernels():
    geo=geometries()
    names=tuple(sorted(geo))
    n=len(names)
    a=prepare_host_resource_simulator(
        geo,names,
        geometry_kernel=np.eye(n),
        resource_breadth_kernel=np.eye(n),
        host_resource_positive_kernel=np.eye(n),
    )
    b=prepare_host_resource_simulator(
        geo,names,
        geometry_kernel=np.ones((n,n)),
        resource_breadth_kernel=np.ones((n,n)),
        host_resource_positive_kernel=np.ones((n,n)),
    )
    wa=simulate_prepared_host_resource_world(a,cell="private",seed=11)
    wb=simulate_prepared_host_resource_world(b,cell="private",seed=11)
    assert all(np.array_equal(wa.edge_response[n],wb.edge_response[n]) for n in names)
