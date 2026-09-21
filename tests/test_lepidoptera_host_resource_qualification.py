from __future__ import annotations

import numpy as np

from ttf.lepidoptera_host_resource_qualification import (
    frozen_alignment_indices,
    geometry_kernel_from_features,
    resource_breadth_kernel,
    resource_cosine_kernel,
    resource_jaccard_kernel,
)


def test_resource_jaccard_kernel():
    h=np.asarray([
        [1,1,0,0],
        [0,1,1,0],
        [0,0,0,1],
    ],dtype=np.uint8)
    k=resource_jaccard_kernel(h)
    assert np.isclose(k[0,1],1/3)
    assert k[0,2]==0
    assert np.allclose(np.diag(k),1)


def test_breadth_and_geometry_kernels_are_symmetric_unit_diagonal():
    b=resource_breadth_kernel(np.asarray([1,3,9.]),np.asarray([4,12,30.]))
    g=geometry_kernel_from_features(np.asarray([[0,0],[1,0],[2,3]],float))
    for k in (b,g):
        assert np.allclose(k,k.T)
        assert np.allclose(np.diag(k),1)


def test_alignment_is_deterministic_and_capped():
    target=np.column_stack([np.arange(300,dtype=float),np.zeros(300)])
    source=np.column_stack([np.arange(300,dtype=float),np.zeros(300)])
    a=frozen_alignment_indices(
        target,source,target_name="A a",source_name="B b",
        radius=500,maximum_rows=128,
    )
    b=frozen_alignment_indices(
        target,source,target_name="A a",source_name="B b",
        radius=500,maximum_rows=128,
    )
    assert len(a[0])==128
    assert np.array_equal(a[0],b[0])
    assert np.array_equal(a[1],b[1])


def test_resource_cosine_kernel_is_psd_and_unit_diagonal():
    h=np.asarray([
        [1,1,0,0],
        [0,1,1,0],
        [0,0,0,1],
    ],dtype=np.uint8)
    k=resource_cosine_kernel(h)
    assert np.allclose(k,k.T)
    assert np.allclose(np.diag(k),1)
    assert np.min(np.linalg.eigvalsh(k))>-1e-10
    assert np.isclose(k[0,1],0.5)
