from __future__ import annotations

import hashlib
import math
from typing import Iterable

import numpy as np

EARTH_RADIUS_KM=6371.0088


def palearctic_core_grid(
    *,
    lat_min: float=30.0,
    lat_max: float=75.0,
    lon_min: float=-25.0,
    lon_max: float=180.0,
    step_degrees: float=0.5,
) -> np.ndarray:
    if step_degrees<=0:
        raise ValueError("step_degrees must be positive")
    lats=np.arange(lat_min,lat_max+0.5*step_degrees,step_degrees,dtype=float)
    lons=np.arange(lon_min,lon_max+0.5*step_degrees,step_degrees,dtype=float)
    lon_grid,lat_grid=np.meshgrid(lons,lats)
    return np.column_stack((lat_grid.ravel(),lon_grid.ravel()))


def haversine_km(a: np.ndarray,b: np.ndarray) -> np.ndarray:
    a=np.asarray(a,dtype=float)
    b=np.asarray(b,dtype=float)
    lat1=np.deg2rad(a[...,0]); lon1=np.deg2rad(a[...,1])
    lat2=np.deg2rad(b[...,0]); lon2=np.deg2rad(b[...,1])
    dlat=lat2-lat1; dlon=lon2-lon1
    h=np.sin(dlat/2.0)**2+np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2.0)**2
    return 2.0*EARTH_RADIUS_KM*np.arcsin(np.minimum(1.0,np.sqrt(h)))


def within_occurrence_buffer(
    grid_latlon: np.ndarray,
    occurrence_latlon: np.ndarray,
    *,
    radius_km: float=500.0,
    chunk_size: int=512,
) -> np.ndarray:
    grid=np.asarray(grid_latlon,dtype=float)
    occ=np.asarray(occurrence_latlon,dtype=float)
    if grid.ndim!=2 or grid.shape[1]!=2 or occ.ndim!=2 or occ.shape[1]!=2:
        raise ValueError("lat/lon arrays must have shape (n,2)")
    if len(occ)==0:
        return np.zeros(len(grid),dtype=bool)
    out=np.zeros(len(grid),dtype=bool)
    for start in range(0,len(grid),chunk_size):
        stop=min(start+chunk_size,len(grid))
        g=grid[start:stop,None,:]
        o=occ[None,:,:]
        dist=haversine_km(g,o)
        out[start:stop]=np.min(dist,axis=1)<=float(radius_km)
    return out


def deterministic_background_indices(
    species: str,
    grid_latlon: np.ndarray,
    eligible: np.ndarray,
    *,
    maximum: int=2000,
) -> np.ndarray:
    grid=np.asarray(grid_latlon,dtype=float)
    eligible=np.asarray(eligible,dtype=bool)
    idx=np.flatnonzero(eligible)
    if len(idx)<=maximum:
        return idx
    keyed=[]
    for i in idx:
        lat,lon=grid[i]
        token=f"{species}|{lat:.6f}|{lon:.6f}".encode()
        keyed.append((hashlib.sha256(token).hexdigest(),int(i)))
    keyed.sort()
    return np.asarray([i for _,i in keyed[:maximum]],dtype=np.int64)


def standardize_from_background(
    background: np.ndarray,
    *arrays: np.ndarray,
) -> tuple[np.ndarray,...]:
    bg=np.asarray(background,dtype=float)
    mu=bg.mean(axis=0)
    sd=bg.std(axis=0)
    if np.any(~np.isfinite(mu)) or np.any(~np.isfinite(sd)) or np.any(sd<=0):
        raise ValueError("degenerate background climate")
    outputs=[(np.asarray(a,dtype=float)-mu)/sd for a in arrays]
    return (*outputs,mu,sd)


def quadratic_features(z: np.ndarray) -> np.ndarray:
    z=np.asarray(z,dtype=float)
    if z.ndim!=2:
        raise ValueError("z must be 2D")
    return np.column_stack((z,z*z))


def normalize_surface(values: np.ndarray) -> np.ndarray:
    x=np.asarray(values,dtype=float)
    if x.ndim!=1 or np.any(~np.isfinite(x)):
        raise ValueError("surface must be a finite 1D vector")
    x=np.clip(x,0.0,None)
    total=float(x.sum())
    if total<=0:
        raise ValueError("surface has zero total suitability")
    return x/total


def schoener_d(a: np.ndarray,b: np.ndarray) -> float:
    pa=normalize_surface(a)
    pb=normalize_surface(b)
    if pa.shape!=pb.shape:
        raise ValueError("surface shapes differ")
    return float(1.0-0.5*np.abs(pa-pb).sum())


def ordered_nonself_pairs(species: Iterable[str]) -> list[tuple[str,str]]:
    names=sorted(map(str,species))
    return [(s,t) for s in names for t in names if s!=t]
