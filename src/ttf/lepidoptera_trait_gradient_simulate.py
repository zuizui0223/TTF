from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping
import hashlib
import numpy as np

from .genetic_geometry import GeneticSamplingGeometry

SEED_TAG = "lepidoptera-trait-gradient-v01-formal"

@dataclass(frozen=True)
class TraitGradientSyntheticWorld:
    edge_response: Mapping[str, np.ndarray]
    cell: str
    seed: int

def frozen_seed(master_seed: int, cell: str, replicate: int) -> int:
    allowed={"private","trait_gradient_positive","geometry_confounded_trap"}
    if cell not in allowed: raise ValueError("unknown frozen cell")
    payload=f"{master_seed}|{SEED_TAG}|{cell}|{replicate}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")

def _psd_factor(kernel: np.ndarray) -> np.ndarray:
    k=np.asarray(kernel,float)
    if k.ndim!=2 or k.shape[0]!=k.shape[1] or not np.isfinite(k).all():
        raise ValueError("kernel must be finite square")
    k=0.5*(k+k.T)
    diag=np.sqrt(np.maximum(np.diag(k),np.finfo(float).tiny))
    k=k/diag[:,None]/diag[None,:]
    values,vectors=np.linalg.eigh(k)
    values=np.maximum(values,1e-10)
    return vectors @ np.diag(np.sqrt(values))

def _standardized_coordinates(geometries):
    labels=tuple(sorted(geometries))
    pooled=np.vstack([geometries[n].coordinates for n in labels])
    center=np.median(pooled,axis=0); scale=np.std(pooled,axis=0)
    scale[scale<=np.sqrt(np.finfo(float).eps)]=1.
    return {n:(geometries[n].coordinates-center)/scale for n in labels}

def _edge_field(g,coords,normal,offset,width):
    nodes=np.asarray(g.edge_nodes,int)
    state=np.tanh(((coords@normal)-offset)/float(width))
    return state[nodes[:,0]]-state[nodes[:,1]]

def simulate_trait_gradient_world(
    geometries: Mapping[str, GeneticSamplingGeometry],
    species_order: tuple[str,...],
    trait_kernel: np.ndarray,
    geometry_kernel: np.ndarray,
    *,
    cell: str,
    seed: int,
    shared_fraction: float = 0.85,
    private_amplitude: float = 0.35,
    noise_sd: float = 0.10,
    transition_width: float = 0.20,
    latent_fields: int = 6,
) -> TraitGradientSyntheticWorld:
    labels=tuple(species_order)
    if tuple(sorted(geometries))!=tuple(sorted(labels)) or len(set(labels))!=len(labels):
        raise ValueError("species order does not match geometry labels")
    n=len(labels)
    tk=np.asarray(trait_kernel,float); gk=np.asarray(geometry_kernel,float)
    if tk.shape!=(n,n) or gk.shape!=(n,n): raise ValueError("kernel dimension drift")
    if not 0<=shared_fraction<=1 or latent_fields<1: raise ValueError("invalid simulator setting")
    if cell=="private":
        kernel=np.eye(n)
    elif cell=="trait_gradient_positive":
        kernel=(1-shared_fraction)*np.eye(n)+shared_fraction*tk
    elif cell=="geometry_confounded_trap":
        kernel=(1-shared_fraction)*np.eye(n)+shared_fraction*gk
    else:
        raise ValueError("unknown frozen cell")
    factor=_psd_factor(kernel)
    z=_standardized_coordinates(geometries)
    rng=np.random.default_rng(int(seed)); dim=next(iter(z.values())).shape[1]
    responses={name:np.zeros(geometries[name].n_edges,float) for name in labels}
    for _ in range(int(latent_fields)):
        normal=rng.normal(size=dim); normal/=np.linalg.norm(normal)
        pooled=np.concatenate([z[name]@normal for name in labels]); offset=float(np.median(pooled))
        loadings=factor @ rng.normal(size=n)
        for i,name in enumerate(labels):
            responses[name]+=float(loadings[i])*_edge_field(
                geometries[name],z[name],normal,offset,transition_width
            )
    for name in labels:
        pn=rng.normal(size=dim); pn/=np.linalg.norm(pn)
        po=float(np.median(z[name]@pn))
        responses[name]+=float(private_amplitude)*_edge_field(
            geometries[name],z[name],pn,po,transition_width
        )
        responses[name]+=rng.normal(0,float(noise_sd),size=geometries[name].n_edges)
        sd=float(np.std(responses[name]))
        if sd>np.finfo(float).tiny:
            responses[name]=(responses[name]-float(np.mean(responses[name])))/sd
    return TraitGradientSyntheticWorld(edge_response=responses,cell=cell,seed=int(seed))

def make_worlds(
    geometries,
    species_order,
    trait_kernel,
    geometry_kernel,
    *,
    cell,
    start,
    count,
    master_seed=20260920,
    **kwargs,
):
    if start<0 or count<1: raise ValueError("invalid shard")
    return tuple(
        simulate_trait_gradient_world(
            geometries,species_order,trait_kernel,geometry_kernel,
            cell=cell,seed=frozen_seed(master_seed,cell,start+i),**kwargs
        )
        for i in range(count)
    )
