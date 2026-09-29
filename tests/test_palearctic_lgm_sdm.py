import numpy as np

from ttf.palearctic_lgm_sdm import (
    deterministic_background_indices,
    ordered_nonself_pairs,
    palearctic_core_grid,
    schoener_d,
    within_occurrence_buffer,
)


def test_palearctic_grid_and_ordered_pairs_are_deterministic():
    grid=palearctic_core_grid(
        lat_min=30,lat_max=31,lon_min=0,lon_max=1,step_degrees=0.5
    )
    assert grid.shape==(9,2)
    assert ordered_nonself_pairs(["b","a","c"])[0]==("a","b")
    assert len(ordered_nonself_pairs(["a","b","c"]))==6


def test_occurrence_buffer_contains_nearby_not_far_cells():
    grid=np.asarray([[50.0,10.0],[50.1,10.1],[60.0,30.0]])
    occ=np.asarray([[50.0,10.0]])
    mask=within_occurrence_buffer(grid,occ,radius_km=50)
    assert mask.tolist()==[True,True,False]


def test_background_subsampling_is_order_stable():
    grid=np.asarray([[50+i*0.1,10.0] for i in range(20)])
    eligible=np.ones(20,dtype=bool)
    a=deterministic_background_indices("Species x",grid,eligible,maximum=5)
    b=deterministic_background_indices("Species x",grid,eligible,maximum=5)
    assert np.array_equal(a,b)
    assert len(a)==5


def test_schoener_d_bounds_and_identity():
    a=np.asarray([1.0,2.0,1.0])
    b=np.asarray([0.0,1.0,3.0])
    assert schoener_d(a,a)==1.0
    assert 0.0<=schoener_d(a,b)<=1.0
