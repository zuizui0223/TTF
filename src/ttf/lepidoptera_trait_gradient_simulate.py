from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence
import hashlib
import numpy as np

from .genetic_geometry import GeneticSamplingGeometry


SEED_TAG="lepidoptera-trait-gradient-v01"


@dataclass(frozen=True)
class TraitGradientSyntheticWorld:
    edge_response: Mapping[str,np.ndarray]
    cell: str
    seed: int


def frozen_seed(master_seed:int,cell:str,replicate:int)->int:
    if cell not in {"private","trait_gradient_positive","geometry_confounded_trap"}:
        raise ValueError("unknown frozen cell")
    x=f"{master_seed}|{SEED_TAG}|{cell}|{replicate}".encode()
    return int.from_bytes(hashlib.sha256(x).digest()[:8],"big")


def _standardized_coordinates(geometries):
    labels=tuple(sorted(geometries)); pooled=np.vstack([geometries[n].coordinates for n in labels])
    center=np.median(pooled,axis=0); scale=np.std(pooled,axis=0); scale[scale<=np.sqrt(np.finfo(float).eps)]=1.
    return {n:(geometries[n].coordinates-center)/scale for n in labels}


def _edge_boundary(g:GeneticSamplingGeometry,coords:np.ndarray,normal:np.ndarray,offset:float,width:.2)->np.ndarray:
    nodes=np.asarray(g.edge_nodes,int); state=np.tanh((coords@normal-offset)/float(width))
    return np.abs(state[nodes[:,0]]-state[nodes[:,1]])


def simulate_trait_gradient_world(
    geometries:Mapping[str,GeneticSamplingGeometry],
    trait_score:Mapping[str,float],
    geometry_score:Mapping[str,float],
    *,
    cell:str,
    seed:int,
    private_amplitude:float=.35,
    shared_amplitude:float=1.0,
    noise_sd:float=.10,
    transition_width:float=.20,
)->TraitGradientSyntheticWorld:
    labels=tuple(sorted(geometries))
    if set(labels)!=set(trait_score) or set(labels)!=set(geometry_score): raise ValueError("score labels drift")
    if cell not in {"private","trait_gradient_positive","geometry_confounded_trap"}: raise ValueError("unknown cell")
    z=_standardized_coordinates(geometries); rng=np.random.default_rng(int(seed)); dim=next(iter(z.values())).shape[1]
    normal=rng.normal(size=dim); normal/=np.linalg.norm(normal)
    pooled=np.concatenate([z[n]@normal for n in labels]); offset=float(np.median(pooled))
    out={}
    for n in labels:
        g=geometries[n]; nodes=np.asarray(g.edge_nodes,int)
        shared=_edge_boundary(g,z[n],normal,offset,transition_width)
        pn=rng.normal(size=dim); pn/=np.linalg.norm(pn); po=float(np.median(z[n]@pn))
        private=_edge_boundary(g,z[n],pn,po,transition_width)
        if cell=="private": amp=0.
        elif cell=="trait_gradient_positive": amp=float(shared_amplitude)*float(np.clip(trait_score[n],0,1))
        else: amp=float(shared_amplitude)*float(np.clip(geometry_score[n],0,1))
        noise=np.abs(rng.normal(0,float(noise_sd),size=g.n_edges))
        out[n]=float(private_amplitude)*private+amp*shared+noise
    return TraitGradientSyntheticWorld(edge_response=out,cell=cell,seed=int(seed))


def make_worlds(geometries,trait_score,geometry_score,*,cell,start,count,master_seed=20260920,**kwargs):
    if start<0 or count<1: raise ValueError("invalid shard")
    return tuple(simulate_trait_gradient_world(geometries,trait_score,geometry_score,cell=cell,seed=frozen_seed(master_seed,cell,start+i),**kwargs) for i in range(count))
