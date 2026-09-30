#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from shapely.geometry import shape
from shapely import intersects_xy


VARIABLES=("bio1","bio7","bio12","bio15")


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def load_union(path: Path):
    payload=json.loads(path.read_text())
    if payload.get("type")=="FeatureCollection":
        features=payload.get("features") or []
        if len(features)!=1:
            raise RuntimeError("Palearctic union must contain exactly one feature")
        return shape(features[0]["geometry"])
    if payload.get("type")=="Feature":
        return shape(payload["geometry"])
    return shape(payload)


def global_centres(resolution: float) -> tuple[np.ndarray,np.ndarray]:
    if np.isclose(resolution,0.25):
        lon=-179.875+0.25*np.arange(1440,dtype=float)
        lat=-89.875+0.25*np.arange(720,dtype=float)
        return lon,lat
    if np.isclose(resolution,0.1):
        lon=-179.95+0.1*np.arange(3600,dtype=float)
        lat=-89.95+0.1*np.arange(1800,dtype=float)
        return lon,lat
    raise ValueError("unsupported frozen grid resolution")


def palearctic_grid_points(geometry, resolution: float) -> np.ndarray:
    lon,lat=global_centres(resolution)
    minx,miny,maxx,maxy=geometry.bounds
    lon=lon[(lon>=minx-0.05)&(lon<=maxx+0.05)]
    lat=lat[(lat>=miny-0.05)&(lat<=maxy+0.05)]
    xx,yy=np.meshgrid(lon,lat,indexing="xy")
    flat_x=xx.ravel()
    flat_y=yy.ravel()
    keep=np.asarray(intersects_xy(geometry,flat_x,flat_y),dtype=bool)
    return np.column_stack((flat_y[keep],flat_x[keep]))


def harmonize(values: np.ndarray, variable: str, dataset: str) -> np.ndarray:
    x=np.asarray(values,dtype=float).copy()
    if dataset in {"lgm","trace0"} and variable=="bio1":
        x-=273.15
    return x


def sample_scaled_raster(path: Path, latlon: np.ndarray, *, chunk_size: int=50000):
    import rasterio

    points=np.asarray(latlon,dtype=float)
    if points.ndim!=2 or points.shape[1]!=2:
        raise ValueError("latlon must be n x 2")
    out=np.empty(len(points),dtype=float)
    with rasterio.open(path) as ds:
        if ds.crs is None or not ds.crs.is_geographic:
            raise RuntimeError(f"raster is not geographic: {path}")
        scale=float(ds.scales[0]) if ds.scales else 1.0
        offset=float(ds.offsets[0]) if ds.offsets else 0.0
        for start in range(0,len(points),int(chunk_size)):
            stop=min(len(points),start+int(chunk_size))
            coords=[(float(lon),float(lat)) for lat,lon in points[start:stop]]
            vals=[]
            for sample in ds.sample(coords,masked=True):
                value=sample[0]
                if np.ma.is_masked(value):
                    vals.append(np.nan)
                else:
                    vals.append(float(value)*scale+offset)
            out[start:stop]=vals
        meta={
            "path":str(path),
            "sha256":sha256_path(path),
            "dtype":str(ds.dtypes[0]),
            "nodata":None if ds.nodata is None else float(ds.nodata),
            "gdal_scale":scale,
            "gdal_offset":offset,
            "width":int(ds.width),
            "height":int(ds.height),
            "crs":str(ds.crs),
        }
    return out,meta


def sample_environment(paths: dict[str,Path], points: np.ndarray, *, dataset: str):
    cols=[]; metas={}
    for var in VARIABLES:
        values,meta=sample_scaled_raster(paths[var],points)
        cols.append(harmonize(values,var,dataset))
        meta["logical_variable"]=var
        meta["harmonization"]=(
            "kelvin_to_celsius_minus_273.15"
            if dataset in {"lgm","trace0"} and var=="bio1"
            else "identity"
        )
        metas[var]=meta
    return np.column_stack(cols),metas


def compatibility(current: np.ndarray, trace0: np.ndarray) -> dict:
    out={}
    for i,var in enumerate(VARIABLES):
        a=np.asarray(current[:,i],float)
        b=np.asarray(trace0[:,i],float)
        valid=np.isfinite(a)&np.isfinite(b)
        if valid.sum()<3:
            out[var]={"n":int(valid.sum()),"pearson":None}
            continue
        diff=b[valid]-a[valid]
        out[var]={
            "n":int(valid.sum()),
            "pearson":float(np.corrcoef(a[valid],b[valid])[0,1]),
            "median_difference":float(np.median(diff)),
            "q05_difference":float(np.quantile(diff,0.05)),
            "q95_difference":float(np.quantile(diff,0.95)),
        }
    return out


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--occurrences",type=Path,required=True)
    ap.add_argument("--panel-metadata",type=Path,required=True)
    ap.add_argument("--realm-geojson",type=Path,required=True)
    ap.add_argument("--climate-rule",type=Path,required=True)
    for prefix in ("current","lgm","trace0"):
        for var in VARIABLES:
            ap.add_argument(f"--{prefix}-{var}",dest=f"{prefix}_{var}",type=Path,required=True)
    ap.add_argument("--output-occurrence-climate",type=Path,required=True)
    ap.add_argument("--output-grid",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.climate_rule.read_text())
    if rule.get("schema") not in {
        "ttf_genetic_palearctic_lgm_climate_input_rule_v0.1",
        "ttf_genetic_palearctic_lgm_climate_input_rule_v0.2",
    }:
        raise RuntimeError("unexpected frozen climate input rule")
    panel=json.loads(args.panel_metadata.read_text())
    if panel.get("schema")!="ttf_genetic_palearctic_lgm_panel_metadata_v0.1":
        raise RuntimeError("unexpected panel metadata")
    meta={row["species"]:row for row in panel["rows"]}
    if len(meta)!=44:
        raise RuntimeError("panel metadata is not the frozen 44 species")

    occ=list(csv.DictReader(args.occurrences.open(encoding="utf-8")))
    occ_points=[]
    occ_rows=[]
    for row in occ:
        sp=str(row["species"])
        if sp not in meta:
            continue
        try:
            lat=float(row["latitude"]); lon=float(row["longitude"])
        except (TypeError,ValueError):
            continue
        occ_points.append((lat,lon)); occ_rows.append(row)
    occ_points=np.asarray(occ_points,dtype=float)

    current_paths={var:getattr(args,f"current_{var}") for var in VARIABLES}
    lgm_paths={var:getattr(args,f"lgm_{var}") for var in VARIABLES}
    trace0_paths={var:getattr(args,f"trace0_{var}") for var in VARIABLES}

    occ_env,current_meta=sample_environment(current_paths,occ_points,dataset="current")
    valid_occ=np.isfinite(occ_env).all(axis=1)

    args.output_occurrence_climate.parent.mkdir(parents=True,exist_ok=True)
    fields=[
        "species","role","source_key","latitude","longitude",
        "priority_rank","priority_sha256",*VARIABLES,
    ]
    climate_rows=[]
    for row,env,keep in zip(occ_rows,occ_env,valid_occ):
        if not keep:
            continue
        out={k:row[k] for k in fields if k in row}
        for var,value in zip(VARIABLES,env):
            out[var]=float(value)
        climate_rows.append(out)

    counts={sp:0 for sp in meta}
    for row in climate_rows:
        counts[row["species"]]+=1
    valid_species=sorted(sp for sp,n in counts.items() if n>=30)
    source=[sp for sp in valid_species if meta[sp]["role"]=="source"]
    target=[sp for sp in valid_species if meta[sp]["role"]=="target"]
    req=rule["occurrence_climate_validity"]["after_drop_require"]
    coverage_pass=(
        len(valid_species)>=int(req["minimum_total_species"])
        and len(source)>=int(req["minimum_source_clusters"])
        and len(target)>=int(req["minimum_target_clusters"])
    )
    if not coverage_pass:
        status=str(rule["occurrence_climate_validity"]["failure"])
        payload={
            "schema":"ttf_genetic_palearctic_lgm_climate_inputs_v0.1",
            "status":status,
            "valid_species":len(valid_species),
            "source_clusters":len(source),
            "target_clusters":len(target),
            "climate_valid_counts":counts,
            "response_firewall":{
                "species_level_genetic_scores_used":False,
                "pairwise_subpanel_T_st_computed":False,
                "beta_LGM_computed":False,
            },
        }
        args.output_summary.parent.mkdir(parents=True,exist_ok=True)
        args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
        print(json.dumps({"status":status,"valid_species":len(valid_species)},sort_keys=True))
        return 0

    valid_set=set(valid_species)
    climate_rows=[row for row in climate_rows if row["species"] in valid_set]
    climate_rows.sort(key=lambda r:(r["species"],int(r["priority_rank"]),int(r["source_key"])))
    with args.output_occurrence_climate.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields)
        writer.writeheader(); writer.writerows(climate_rows)

    geom=load_union(args.realm_geojson)
    resolution=float(rule["grid"]["resolution_degrees"])
    grid_points=palearctic_grid_points(geom,resolution)
    current_grid,_=sample_environment(current_paths,grid_points,dataset="current")
    lgm_grid,lgm_meta=sample_environment(lgm_paths,grid_points,dataset="lgm")
    trace0_grid,trace0_meta=sample_environment(trace0_paths,grid_points,dataset="trace0")

    grid_valid=np.isfinite(current_grid).all(axis=1)&np.isfinite(lgm_grid).all(axis=1)
    current_valid=current_grid[grid_valid]
    lgm_valid=lgm_grid[grid_valid]
    point_valid=grid_points[grid_valid]
    trace0_valid=trace0_grid[grid_valid]

    np.savez_compressed(
        args.output_grid,
        latitude=point_valid[:,0],
        longitude=point_valid[:,1],
        current_environment=current_valid,
        lgm_environment=lgm_valid,
    )

    payload={
        "schema":"ttf_genetic_palearctic_lgm_climate_inputs_v0.1",
        "status":"PASS_TO_RESPONSE_BLIND_LGM_RELATION_BUILD",
        "formal_panel_species":44,
        "climate_valid_species":len(valid_species),
        "source_clusters":len(source),
        "target_clusters":len(target),
        "dropped_for_climate_coverage":sorted(set(meta)-valid_set),
        "climate_valid_counts":counts,
        "occurrence_rows_before_current_climate_validity":len(occ_rows),
        "occurrence_rows_after_current_climate_validity_and_species_gate":len(climate_rows),
        "grid":{
            "candidate_palearctic_centroids":int(len(grid_points)),
            "common_current_lgm_finite_cells":int(grid_valid.sum()),
            "resolution_degrees":resolution,
            "paleocoastline_reconstruction":False,
        },
        "compatibility_current_V2_1_vs_TraCE_0BP":compatibility(current_valid,trace0_valid),
        "raster_metadata":{
            "current":current_meta,
            "lgm":lgm_meta,
            "trace0":trace0_meta,
        },
        "inputs":{
            "occurrences_sha256":sha256_path(args.occurrences),
            "panel_metadata_sha256":sha256_path(args.panel_metadata),
            "realm_geojson_sha256":sha256_path(args.realm_geojson),
            "climate_rule_sha256":sha256_path(args.climate_rule),
        },
        "outputs":{
            "occurrence_climate_sha256":sha256_path(args.output_occurrence_climate),
            "grid_npz_sha256":sha256_path(args.output_grid),
        },
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_subpanel_T_st_computed":False,
            "beta_LGM_computed":False,
        },
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "climate_valid_species":len(valid_species),
        "grid_cells":int(grid_valid.sum()),
        "compatibility":payload["compatibility_current_V2_1_vs_TraCE_0BP"],
    },sort_keys=True))


if __name__=="__main__":
    raise SystemExit(main())
