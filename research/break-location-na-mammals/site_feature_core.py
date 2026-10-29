"""Predictor extraction only: explicit nodata and bounded spherical GIS distances.

Polyline vertices are filled at <=0.01 degrees in wrapped lon/lat. Distance to
nearest vertex overestimates continuous spherical-line distance by <=0.79 km.
WGS84 ellipsoid vs mean-radius sphere is a separate approximately 0.6% bound.
No genetic outcomes are fitted here. Restricted coordinates stay ephemeral.
"""
from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree
import rasterio
from pyproj import Transformer

R_KM=6371.0088
STEP_DEG=.01

def sphere_xyz(lonlat):
    a=np.asarray(lonlat,dtype=float)
    if a.ndim!=2 or a.shape[1]!=2 or not np.isfinite(a).all():
        raise ValueError('finite n x 2 lon/lat coordinates required')
    if (np.abs(a[:,0])>180).any() or (np.abs(a[:,1])>90).any():
        raise ValueError('invalid WGS84 coordinates')
    lo,la=np.deg2rad(a).T
    return np.column_stack((np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)))

def coordinate_paths(g):
    if g.is_empty:return
    kind=g.geom_type
    if kind=='Polygon':
        yield np.asarray(g.exterior.coords)[:,:2]
        for r in g.interiors:yield np.asarray(r.coords)[:,:2]
    elif kind in ('LineString','LinearRing'):
        yield np.asarray(g.coords)[:,:2]
    elif hasattr(g,'geoms'):
        for part in g.geoms:yield from coordinate_paths(part)
    else:raise ValueError('expected polygon or linear geometry')

def dense_vertices(geometries):
    chunks=[]
    for g in geometries:
        for p in coordinate_paths(g):
            if len(p)<2:continue
            delta=np.diff(p,axis=0)
            delta[:,0]=(delta[:,0]+180)%360-180
            steps=np.maximum(1,np.ceil(np.linalg.norm(delta,axis=1)/STEP_DEG).astype(int))
            # Original paths are handled independently; no bridges between rings.
            for start in range(0,len(steps),10000):
                n=steps[start:start+10000]
                idx=np.repeat(np.arange(len(n)),n)
                off=np.repeat(np.r_[0,np.cumsum(n)[:-1]],n)
                f=(np.arange(n.sum())-off)/n[idx]
                v=p[start+idx]+delta[start+idx]*f[:,None]
                v[:,0]=(v[:,0]+180)%360-180
                chunks.append(v)
            chunks.append(p[-1:].copy())
    if not chunks:raise ValueError('empty required barrier geometry')
    out=np.vstack(chunks);sphere_xyz(out)
    return out

def vertex_distances(lonlat,vertices):
    tree=vertices if isinstance(vertices,cKDTree) else cKDTree(sphere_xyz(vertices))
    d,_=tree.query(sphere_xyz(lonlat),workers=2)
    return 2*R_KM*np.arcsin(np.clip(d/2,0,1))

def raster_values(path,lonlat,mask_nodata=True):
    xy=np.asarray(lonlat,dtype=float);sphere_xyz(xy)
    out=np.full(len(xy),np.nan)
    with rasterio.open(path) as ds:
        if ds.crs is None or ds.count!=1:raise ValueError('single-band referenced raster required')
        t=Transformer.from_crs(4326,ds.crs,always_xy=True)
        xx,yy=t.transform(xy[:,0],xy[:,1])
        keep=np.isfinite(xx)&np.isfinite(yy)&(xx>=ds.bounds.left)&(xx<ds.bounds.right)&(yy>ds.bounds.bottom)&(yy<=ds.bounds.top)
        ids=np.flatnonzero(keep)
        # Sort for raster block locality without changing row order in returned data.
        order=ids[np.lexsort((xx[ids],yy[ids]))]
        for i,v in zip(order,ds.sample(zip(xx[order],yy[order]),indexes=1,masked=mask_nodata)):
            if np.ma.is_masked(v[0]):continue
            val=float(v[0])
            if mask_nodata and ds.nodata is not None and val==ds.nodata:continue
            if np.isfinite(val):out[i]=val
    return out

def complete_rows(matrix):
    a=np.asarray(matrix,float)
    if a.ndim!=2:raise ValueError('feature matrix required')
    return np.isfinite(a).all(axis=1)
