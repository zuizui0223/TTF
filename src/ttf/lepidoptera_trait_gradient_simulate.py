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

@dataclass(frozen=True)
class PreparedTraitGradientSimulator:
    geometries: Mapping[str, GeneticSamplingGeometry]
    species_order: tuple[str, ...]
    standardized_coordinates: Mapping[str, np.ndarray]
    factors: Mapping[str, np.ndarray]

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

def prepare_trait_gradient_simulator(
    geometries: Mapping[str, GeneticSamplingGeometry],
    species_order: tuple[str,...],
    trait_kernel: np.ndarray,
    geometry_kernel: np.ndarray,
    *,
    shared_fraction: float = 0.85,
) -> PreparedTraitGradientSimulator:
    labels=tuple(species_order)
    if tuple(sorted(geometries))!=tuple(sorted(labels)) or len(set(labels))!=len(labels):
        raise ValueError("species order does not match geometry labels")
    n=len(labels)
    tk=np.asarray(trait_kernel,float); gk=np.asarray(geometry_kernel,float)
    if tk.shape!=(n,n) or gk.shape!=(n,n): raise ValueError("kernel dimension drift")
    if not 0<=shared_fraction<=1: raise ValueError("invalid shared fraction")
    factors={
        "private": np.eye(n),
        "trait_gradient_positive": _psd_factor((1-shared_fraction)*np.eye(n)+shared_fraction*tk),
        "geometry_confounded_trap": _psd_factor((1-shared_fraction)*np.eye(n)+shared_fraction*gk),
    }
    return PreparedTraitGradientSimulator(
        geometries=geometries,
        species_order=labels,
        standardized_coordinates=_standardized_coordinates(geometries),
        factors=factors,
    )

def _edge_field(g,coords,normal,offset,width):
    nodes=np.asarray(g.edge_nodes,int)
    state=np.tanh(((coords@normal)-offset)/float(width))
    return state[nodes[:,0]]-state[nodes[:,1]]

def simulate_prepared_trait_gradient_world(
    prepared: PreparedTraitGradientSimulator,
    *,
    cell: str,
    seed: int,
    private_amplitude: float = 0.35,
    noise_sd: float = 0.10,
    transition_width: float = 0.20,
    latent_fields: int = 6,
) -> TraitGradientSyntheticWorld:
    labels=prepared.species_order
    if cell not in prepared.factors: raise ValueError("unknown frozen cell")
    if latent_fields<1: raise ValueError("latent_fields must be positive")
    factor=prepared.factors[cell]
    z=prepared.standardized_coordinates
    rng=np.random.default_rng(int(seed)); dim=next(iter(z.values())).shape[1]
    responses={name:np.zeros(prepared.geometries[name].n_edges,float) for name in labels}
    pooled_cache={}
    for _ in range(int(latent_fields)):
        normal=rng.normal(size=dim); normal/=np.linalg.norm(normal)
        key=tuple(normal.tolist())
        pooled=np.concatenate([z[name]@normal for name in labels]); offset=float(np.median(pooled))
        loadings=factor @ rng.normal(size=len(labels))
        for i,name in enumerate(labels):
            responses[name]+=float(loadings[i])*_edge_field(
                prepared.geometries[name],z[name],normal,offset,transition_width
            )
    for name in labels:
        pn=rng.normal(size=dim); pn/=np.linalg.norm(pn)
        po=float(np.median(z[name]@pn))
        responses[name]+=float(private_amplitude)*_edge_field(
            prepared.geometries[name],z[name],pn,po,transition_width
        )
        responses[name]+=rng.normal(0,float(noise_sd),size=prepared.geometries[name].n_edges)
        sd=float(np.std(responses[name]))
        if sd>np.finfo(float).tiny:
            responses[name]=(responses[name]-float(np.mean(responses[name])))/sd
    return TraitGradientSyntheticWorld(edge_response=responses,cell=cell,seed=int(seed))

def simulate_trait_gradient_world(
    geometries: Mapping[str, GeneticSamplingGeometry],
    species_order: tuple[str,...],
    trait_kernel: np.ndarray,
    geometry_kernel: np.ndarray,
    *,
    cell: str,
    seed: int,
    shared_fraction: float = 0.85,
    **kwargs,
) -> TraitGradientSyntheticWorld:
    prepared=prepare_trait_gradient_simulator(
        geometries,species_order,trait_kernel,geometry_kernel,shared_fraction=shared_fraction
    )
    return simulate_prepared_trait_gradient_world(prepared,cell=cell,seed=seed,**kwargs)

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
    shared_fraction=0.85,
    **kwargs,
):
    if start<0 or count<1: raise ValueError("invalid shard")
    prepared=prepare_trait_gradient_simulator(
        geometries,species_order,trait_kernel,geometry_kernel,shared_fraction=shared_fraction
    )
    return tuple(
        simulate_prepared_trait_gradient_world(
            prepared,cell=cell,seed=frozen_seed(master_seed,cell,start+i),**kwargs
        )
        for i in range(count)
    )
