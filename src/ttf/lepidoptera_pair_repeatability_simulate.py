from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Mapping

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry


REFERENCE_TAG="lepidoptera-pair-repeatability-v01-reference"
EVALUATION_TAG="lepidoptera-pair-repeatability-v01-evaluation"
ALLOWED_CELLS=(
    "private",
    "geometry_confounded_trap",
    "host_breadth_confounded_trap",
    "host_resource_gradient_trap",
    "relational_latent_positive",
)


@dataclass(frozen=True)
class PairRepeatabilitySyntheticWorld:
    edge_response: Mapping[str,np.ndarray]
    cell: str
    seed: int


@dataclass(frozen=True)
class PreparedPairRepeatabilitySimulator:
    geometries: Mapping[str,GeneticSamplingGeometry]
    species_order: tuple[str,...]
    standardized_coordinates: Mapping[str,np.ndarray]
    factors: Mapping[str,np.ndarray]
    relational_kernel: np.ndarray


def frozen_seed(master_seed:int,tag:str,cell:str,replicate:int)->int:
    if tag not in {REFERENCE_TAG,EVALUATION_TAG}:
        raise ValueError("unauthorized pair-repeatability seed namespace")
    if cell not in ALLOWED_CELLS:
        raise ValueError("unknown pair-repeatability cell")
    if replicate<0:
        raise ValueError("replicate must be non-negative")
    payload=f"{int(master_seed)}|{tag}|{cell}|{int(replicate)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")


def _psd_factor(kernel:np.ndarray)->np.ndarray:
    k=np.asarray(kernel,float)
    if k.ndim!=2 or k.shape[0]!=k.shape[1] or not np.isfinite(k).all():
        raise ValueError("kernel must be finite square")
    k=.5*(k+k.T)
    diag=np.sqrt(np.maximum(np.diag(k),np.finfo(float).tiny))
    k=k/diag[:,None]/diag[None,:]
    values,vectors=np.linalg.eigh(k)
    values=np.maximum(values,1e-10)
    return vectors@np.diag(np.sqrt(values))


def relational_latent_kernel(species_order:tuple[str,...],dimensions:int=8)->np.ndarray:
    if dimensions<2:
        raise ValueError("dimensions must be >=2")
    rows=[]
    for name in species_order:
        seed=int.from_bytes(
            hashlib.sha256(
                f"pair-repeatability-relational-kernel-v01|{name}".encode("utf-8")
            ).digest()[:8],
            "big",
        )
        rng=np.random.default_rng(seed)
        v=rng.normal(size=int(dimensions))
        norm=float(np.linalg.norm(v))
        if norm<=np.finfo(float).tiny:
            raise RuntimeError("zero relational embedding")
        rows.append(v/norm)
    e=np.asarray(rows,float)
    k=e@e.T
    np.fill_diagonal(k,1.0)
    return k


def _standardized_coordinates(geometries,labels):
    pooled=np.vstack([np.asarray(geometries[n].coordinates,float) for n in labels])
    center=np.median(pooled,axis=0)
    scale=np.std(pooled,axis=0)
    scale[scale<=np.sqrt(np.finfo(float).eps)]=1.0
    return {
        n:(np.asarray(geometries[n].coordinates,float)-center)/scale
        for n in labels
    }


def prepare_pair_repeatability_simulator(
    geometries:Mapping[str,GeneticSamplingGeometry],
    species_order:tuple[str,...],
    *,
    geometry_kernel:np.ndarray,
    resource_breadth_kernel:np.ndarray,
    host_resource_kernel:np.ndarray,
    relational_dimensions:int=8,
)->PreparedPairRepeatabilitySimulator:
    labels=tuple(map(str,species_order))
    if tuple(sorted(geometries))!=tuple(sorted(labels)):
        raise ValueError("species/geometry drift")
    n=len(labels)
    eye=np.eye(n)
    g=np.asarray(geometry_kernel,float)
    b=np.asarray(resource_breadth_kernel,float)
    h=np.asarray(host_resource_kernel,float)
    if any(k.shape!=(n,n) for k in (g,b,h)):
        raise ValueError("kernel dimension drift")
    r=relational_latent_kernel(labels,dimensions=relational_dimensions)
    factors={
        "private":eye,
        "geometry_confounded_trap":_psd_factor(.15*eye+.85*g),
        "host_breadth_confounded_trap":_psd_factor(.15*eye+.85*b),
        "host_resource_gradient_trap":_psd_factor(.10*eye+.90*h),
        "relational_latent_positive":_psd_factor(.10*eye+.90*r),
    }
    return PreparedPairRepeatabilitySimulator(
        geometries=geometries,
        species_order=labels,
        standardized_coordinates=_standardized_coordinates(geometries,labels),
        factors=factors,
        relational_kernel=r,
    )


def _edge_field(geometry,coordinates,normal,offset,width):
    nodes=np.asarray(geometry.edge_nodes,np.int64)
    state=np.tanh(((coordinates@normal)-float(offset))/float(width))
    return state[nodes[:,0]]-state[nodes[:,1]]


def simulate_pair_repeatability_world(
    prepared:PreparedPairRepeatabilitySimulator,
    *,
    cell:str,
    seed:int,
    latent_fields:int=6,
    private_amplitude:float=.35,
    noise_sd:float=.10,
    transition_width:float=.20,
)->PairRepeatabilitySyntheticWorld:
    if cell not in prepared.factors:
        raise ValueError("unknown pair-repeatability cell")
    labels=prepared.species_order
    z=prepared.standardized_coordinates
    factor=prepared.factors[cell]
    rng=np.random.default_rng(int(seed))
    dim=next(iter(z.values())).shape[1]
    response={
        name:np.zeros(prepared.geometries[name].n_edges,float)
        for name in labels
    }
    for _ in range(int(latent_fields)):
        normal=rng.normal(size=dim)
        normal/=np.linalg.norm(normal)
        pooled=np.concatenate([z[n]@normal for n in labels])
        offset=float(np.median(pooled))
        loading=factor@rng.normal(size=len(labels))
        for i,name in enumerate(labels):
            response[name]+=float(loading[i])*_edge_field(
                prepared.geometries[name],z[name],normal,offset,transition_width
            )
    for name in labels:
        normal=rng.normal(size=dim)
        normal/=np.linalg.norm(normal)
        offset=float(np.median(z[name]@normal))
        response[name]+=float(private_amplitude)*_edge_field(
            prepared.geometries[name],z[name],normal,offset,transition_width
        )
        response[name]+=rng.normal(
            0,float(noise_sd),size=prepared.geometries[name].n_edges
        )
        sd=float(np.std(response[name]))
        if sd>np.finfo(float).tiny:
            response[name]=(response[name]-float(np.mean(response[name])))/sd
    return PairRepeatabilitySyntheticWorld(
        edge_response=response,cell=str(cell),seed=int(seed)
    )


__all__=[
    "ALLOWED_CELLS",
    "EVALUATION_TAG",
    "REFERENCE_TAG",
    "PairRepeatabilitySyntheticWorld",
    "PreparedPairRepeatabilitySimulator",
    "frozen_seed",
    "prepare_pair_repeatability_simulator",
    "relational_latent_kernel",
    "simulate_pair_repeatability_world",
]
