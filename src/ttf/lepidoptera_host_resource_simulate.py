from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry


@dataclass(frozen=True)
class HostResourceSyntheticWorld:
    edge_response: Mapping[str,np.ndarray]
    cell: str
    seed: int


@dataclass(frozen=True)
class PreparedHostResourceSimulator:
    geometries: Mapping[str,GeneticSamplingGeometry]
    species_order: tuple[str,...]
    standardized_coordinates: Mapping[str,np.ndarray]
    factors: Mapping[str,np.ndarray]


def _psd_factor(kernel: np.ndarray) -> np.ndarray:
    k=np.asarray(kernel,float)
    if k.ndim!=2 or k.shape[0]!=k.shape[1] or not np.isfinite(k).all():
        raise ValueError("kernel must be a finite square matrix")
    k=.5*(k+k.T)
    diag=np.sqrt(np.maximum(np.diag(k),np.finfo(float).tiny))
    k=k/diag[:,None]/diag[None,:]
    values,vectors=np.linalg.eigh(k)
    values=np.maximum(values,1e-10)
    return vectors@np.diag(np.sqrt(values))


def _standardized_coordinates(
    geometries: Mapping[str,GeneticSamplingGeometry],
    labels: tuple[str,...],
) -> dict[str,np.ndarray]:
    pooled=np.vstack([np.asarray(geometries[n].coordinates,float) for n in labels])
    center=np.median(pooled,axis=0)
    scale=np.std(pooled,axis=0)
    scale[scale<=np.sqrt(np.finfo(float).eps)]=1.0
    return {
        n:(np.asarray(geometries[n].coordinates,float)-center)/scale
        for n in labels
    }


def prepare_host_resource_simulator(
    geometries: Mapping[str,GeneticSamplingGeometry],
    species_order: tuple[str,...],
    *,
    geometry_kernel: np.ndarray,
    resource_breadth_kernel: np.ndarray,
    host_resource_positive_kernel: np.ndarray,
) -> PreparedHostResourceSimulator:
    labels=tuple(map(str,species_order))
    if tuple(sorted(geometries))!=tuple(sorted(labels)) or len(set(labels))!=len(labels):
        raise ValueError("species order does not match geometry labels")
    n=len(labels)
    eye=np.eye(n)
    g=np.asarray(geometry_kernel,float)
    b=np.asarray(resource_breadth_kernel,float)
    h=np.asarray(host_resource_positive_kernel,float)
    if g.shape!=(n,n) or b.shape!=(n,n) or h.shape!=(n,n):
        raise ValueError("kernel dimension drift")
    factors={
        "private":eye,
        "geometry_confounded_trap":_psd_factor(.15*eye+.85*g),
        "host_breadth_confounded_trap":_psd_factor(.15*eye+.85*b),
        "host_resource_gradient_positive":_psd_factor(.10*eye+.90*h),
    }
    return PreparedHostResourceSimulator(
        geometries=geometries,
        species_order=labels,
        standardized_coordinates=_standardized_coordinates(geometries,labels),
        factors=factors,
    )


def _edge_field(
    geometry: GeneticSamplingGeometry,
    coordinates: np.ndarray,
    normal: np.ndarray,
    offset: float,
    width: float,
) -> np.ndarray:
    nodes=np.asarray(geometry.edge_nodes,np.int64)
    state=np.tanh(((coordinates@normal)-float(offset))/float(width))
    return state[nodes[:,0]]-state[nodes[:,1]]


def simulate_prepared_host_resource_world(
    prepared: PreparedHostResourceSimulator,
    *,
    cell: str,
    seed: int,
    latent_fields: int=6,
    private_amplitude: float=.35,
    noise_sd: float=.10,
    transition_width: float=.20,
) -> HostResourceSyntheticWorld:
    if cell not in prepared.factors:
        raise ValueError("unknown host-resource synthetic cell")
    if latent_fields<1:
        raise ValueError("latent_fields must be positive")
    labels=prepared.species_order
    z=prepared.standardized_coordinates
    factor=prepared.factors[cell]
    rng=np.random.default_rng(int(seed))
    dim=next(iter(z.values())).shape[1]
    response={
        name:np.zeros(prepared.geometries[name].n_edges,dtype=float)
        for name in labels
    }
    for _ in range(int(latent_fields)):
        normal=rng.normal(size=dim)
        normal/=np.linalg.norm(normal)
        pooled=np.concatenate([z[name]@normal for name in labels])
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
        response[name]+=rng.normal(0,float(noise_sd),size=prepared.geometries[name].n_edges)
        sd=float(np.std(response[name]))
        if sd>np.finfo(float).tiny:
            response[name]=(response[name]-float(np.mean(response[name])))/sd
    return HostResourceSyntheticWorld(
        edge_response=response,
        cell=str(cell),
        seed=int(seed),
    )


__all__=[
    "HostResourceSyntheticWorld",
    "PreparedHostResourceSimulator",
    "prepare_host_resource_simulator",
    "simulate_prepared_host_resource_world",
]
