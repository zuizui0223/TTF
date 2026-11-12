"""GIS distance construction for the Jensen baseline; stand-alone, no genetics.

All inputs are WGS84 source geometries whose publisher SHA has been verified.
Pre-clip to species survey extent before calling for memory efficiency.
This module only defines deterministic geographic covariates, never fits breaks.
"""
from __future__ import annotations
from collections.abc import Iterable
import numpy as np
from pyproj import CRS,Transformer
from shapely.geometry import Point
from shapely.ops import transform,unary_union

def local_transformer(center_lon:float,center_lat:float)->Transformer:
    if not -180<=center_lon<=180 or not -90<=center_lat<=90:
        raise ValueError("invalid local projection anchor")
    crs=CRS.from_proj4("+proj=aeqd +lon_0="+str(center_lon)+" +lat_0="+str(center_lat)+" +datum=WGS84 +units=m +no_defs")
    return Transformer.from_crs("EPSG:4326",crs,always_xy=True)

def boundary_network(polygons):
    """Do NOT take the boundary of the dissolved union: it destroys ecotones."""
    if not polygons:raise ValueError("no ecological regions")
    edges=[g.boundary for g in polygons if g is not None and not g.is_empty]
    if not edges:raise ValueError("empty ecoregion boundaries")
    return unary_union(edges)

def distance_covariates_wgs84(points_lon_lat, mountain_polys, river_lines,
                              large_lake_polys, ecoregion_l1_polys,
                              ecoregion_l2_polys):
    points=list(points_lon_lat)
    if not points:raise ValueError("empty sites")
    if any(not (-180<=p[0]<=180 and -90<=p[1]<=90) for p in points):
        raise ValueError("bad geographic site coordinate")
    mx=float(np.median([p[0] for p in points]))
    my=float(np.median([p[1] for p in points]))
    t=local_transformer(mx,my)
    def projected_set(geom):
        valid=[transform(t.transform,g) for g in geom if g is not None and not g.is_empty]
        if not valid:raise ValueError("missing required GIS barrier input")
        return valid
    mountain=unary_union(projected_set(mountain_polys))
    water=unary_union(projected_set(list(river_lines)+list(large_lake_polys)))
    eco1=boundary_network(projected_set(ecoregion_l1_polys))
    eco2=boundary_network(projected_set(ecoregion_l2_polys))
    result=np.empty((len(points),4),dtype=float)
    for i,(lon,lat) in enumerate(points):
        p=Point(*t.transform(lon,lat))
        result[i,:]=[p.distance(obj)/1000.0 for obj in (mountain,water,eco1,eco2)]
    if not np.isfinite(result).all() or (result<0).any():raise ValueError("nonfinite source distances")
    return result

def original_distance_transform(distances_km):
    if np.any(np.asarray(distances_km)<0):raise ValueError("negative GIS distance")
    return np.sqrt(np.asarray(distances_km,dtype=float))
