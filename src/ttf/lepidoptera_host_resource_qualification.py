from __future__ import annotations

import hashlib
from typing import Mapping

import numpy as np

from .lepidoptera_trait_gradient_v02 import target_fe_geometry_residual
from .private_null_inference import envelope_upper_pvalues


REFERENCE_TAG="lepidoptera-host-resource-v01-reference"
EVALUATION_TAG="lepidoptera-host-resource-v01-evaluation"
ALLOWED_CELLS=("private","geometry_confounded_trap","resource_breadth_trap","host_resource_positive")


def resource_jaccard_kernel(host_presence: np.ndarray) -> np.ndarray:
    h=np.asarray(host_presence,dtype=np.uint8)
    if h.ndim!=2 or len(h)<2:
        raise ValueError("host_presence must be a species x unit matrix")
    intersection=h.astype(np.int64)@h.astype(np.int64).T
    size=np.sum(h,axis=1,dtype=np.int64)
    union=size[:,None]+size[None,:]-intersection
    out=np.zeros_like(intersection,dtype=float)
    valid=union>0
    out[valid]=intersection[valid]/union[valid]
    np.fill_diagonal(out,1.0)
    return out


def resource_breadth_kernel(
    native_host_species_count: np.ndarray,
    footprint_unit_count: np.ndarray,
) -> np.ndarray:
    hosts=np.asarray(native_host_species_count,float)
    units=np.asarray(footprint_unit_count,float)
    if hosts.ndim!=1 or units.ndim!=1 or len(hosts)!=len(units):
        raise ValueError("breadth arrays drift")
    def zlog(x):
        y=np.log1p(np.maximum(x,0.0))
        sd=float(np.std(y))
        return np.zeros_like(y) if sd<=np.sqrt(np.finfo(float).eps) else (y-float(np.mean(y)))/sd
    zh=zlog(hosts); zu=zlog(units)
    return 0.5*(
        np.exp(-np.abs(zh[:,None]-zh[None,:]))
        +np.exp(-np.abs(zu[:,None]-zu[None,:]))
    )


def geometry_kernel_from_features(features: np.ndarray) -> np.ndarray:
    x=np.asarray(features,float)
    if x.ndim!=2 or len(x)<2 or not np.isfinite(x).all():
        raise ValueError("geometry features must be a finite matrix")
    scale=np.std(x,axis=0,ddof=0)
    scale[scale<=np.sqrt(np.finfo(float).eps)]=1.0
    z=(x-np.mean(x,axis=0))/scale
    d2=np.mean((z[:,None,:]-z[None,:,:])**2,axis=2)
    return np.exp(-0.5*d2)


def residualize_host_resource_similarity(
    raw_jaccard: np.ndarray,
    nuisance_covariates: np.ndarray,
    target_index: np.ndarray,
) -> np.ndarray:
    return target_fe_geometry_residual(
        np.asarray(raw_jaccard,float),
        np.asarray(nuisance_covariates,float),
        np.asarray(target_index,np.int64),
    )


def frozen_alignment_indices(
    target_midpoint: np.ndarray,
    source_midpoint: np.ndarray,
    *,
    target_name: str,
    source_name: str,
    radius: float=500.0,
    maximum_rows: int=128,
) -> tuple[np.ndarray,np.ndarray]:
    target=np.asarray(target_midpoint,float)
    source=np.asarray(source_midpoint,float)
    if target.ndim!=2 or source.ndim!=2 or target.shape[1]!=source.shape[1]:
        raise ValueError("midpoint dimensions drift")
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        delta=target[:,None,:]-source[None,:,:]
        d2=np.sum(delta*delta,axis=2)
        nearest=np.argmin(d2,axis=1)
        distance=np.sqrt(d2[np.arange(len(target)),nearest])
    else:
        distance,nearest=cKDTree(source).query(target,k=1)
    valid=np.flatnonzero(distance<=float(radius))
    if len(valid)<3:
        raise ValueError("fewer than three aligned target edges")
    if len(valid)>int(maximum_rows):
        payload=f"host-resource-v03-alignment|{target_name}|{source_name}".encode("utf-8")
        seed=int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")
        rng=np.random.default_rng(seed)
        valid=np.sort(rng.choice(valid,size=int(maximum_rows),replace=False))
    return valid.astype(np.int64),np.asarray(nearest[valid],dtype=np.int64)


def frozen_seed(master_seed: int, tag: str, cell: str, replicate: int) -> int:
    if tag not in {REFERENCE_TAG,EVALUATION_TAG}:
        raise ValueError("unauthorized host-resource seed namespace")
    if cell not in ALLOWED_CELLS:
        raise ValueError("unknown host-resource synthetic cell")
    if replicate<0:
        raise ValueError("replicate must be non-negative")
    payload=f"{int(master_seed)}|{tag}|{cell}|{int(replicate)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")


def nuisance_envelope_pvalues(
    observed: np.ndarray,
    *,
    private_reference: np.ndarray,
    geometry_reference: np.ndarray,
    breadth_reference: np.ndarray,
) -> np.ndarray:
    p,_=envelope_upper_pvalues(
        np.asarray(observed,float),
        {
            "private":np.asarray(private_reference,float),
            "geometry_confounded_trap":np.asarray(geometry_reference,float),
            "resource_breadth_trap":np.asarray(breadth_reference,float),
        },
    )
    return np.asarray(p,float)


__all__=[
    "ALLOWED_CELLS",
    "EVALUATION_TAG",
    "REFERENCE_TAG",
    "frozen_alignment_indices",
    "frozen_seed",
    "geometry_kernel_from_features",
    "nuisance_envelope_pvalues",
    "residualize_host_resource_similarity",
    "resource_breadth_kernel",
    "resource_jaccard_kernel",
]
