"""Frozen 50km break-line positional jitter within surveyed range. No genetic effect fit."""
import numpy as np
import shapely

def jitter_cases(points,polygon,rng,radius_m=50000,max_attempts=20):
    if radius_m!=50000 or max_attempts!=20:
        raise ValueError("Frozen 50km 20-attempt uncertainty rule must not change")
    pts=np.asarray(points,float)
    if pts.ndim!=2 or pts.shape[1]!=2 or not np.isfinite(pts).all():
        raise ValueError("Finite n x 2 projected case-point coordinates needed")
    if polygon.is_empty or not polygon.is_valid:
        raise ValueError("Source surveyed range is invalid")
    if not shapely.covers(polygon,shapely.points(pts)).all():
        raise ValueError("Original break case point outside surveyed range")
    out=pts.copy();pending=np.ones(len(pts),dtype=bool)
    attempted=0
    for trial in range(max_attempts):
        ids=np.flatnonzero(pending)
        if not len(ids):break
        angle=rng.uniform(0,2*np.pi,size=len(ids))
        length=radius_m*np.sqrt(rng.uniform(size=len(ids)))
        candidate=pts[ids]+np.column_stack([length*np.cos(angle),length*np.sin(angle)])
        okay=shapely.contains_xy(polygon,candidate[:,0],candidate[:,1])
        out[ids[okay]]=candidate[okay]
        pending[ids[okay]]=False
        attempted+=len(ids)
    distance=np.linalg.norm(out-pts,axis=1)
    return out,{"attempted_draws":int(attempted),"successful":int((~pending).sum()),
        "unmoved":int(pending.sum()),
        "displacement_km_median":float(np.median(distance)/1000) if len(distance) else None,
        "displacement_km_p95":float(np.quantile(distance,.95)/1000) if len(distance) else None}


def repair_survey_polygon(polygon):
    """Repair projection-induced ring self-intersections, keep polygonal area only."""
    from shapely.validation import make_valid
    from shapely.ops import unary_union
    if polygon.is_empty:raise ValueError("empty projected source surveyed range")
    fixed=make_valid(polygon) if not polygon.is_valid else polygon
    polys=[]
    def collect(g):
        if g.geom_type=="Polygon":polys.append(g)
        elif g.geom_type=="MultiPolygon":polys.extend(g.geoms)
        elif hasattr(g,"geoms"):
            for part in g.geoms:collect(part)
    collect(fixed)
    if not polys:raise ValueError("no polygonal area after projected survey repair")
    output=unary_union(polys)
    if not output.is_valid or output.area<=0:
        raise ValueError("surveyed range remains invalid")
    return output
