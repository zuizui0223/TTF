#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ttf.palearctic_lgm import (
    climate_features,
    fit_climate_scaling,
    fit_weighted_ridge_logistic,
    palearctic_role_order,
    predict_relative_suitability,
    schoener_d,
    target_group_background_cells,
)


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def sample_rasters(coords_lon_lat: np.ndarray, paths: list[Path]) -> tuple[np.ndarray,np.ndarray]:
    import rasterio

    coords=np.asarray(coords_lon_lat,dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be lon/lat n x 2")
    datasets=[rasterio.open(path) for path in paths]
    try:
        cols=[]
        valid=np.ones(len(coords),dtype=bool)
        tuples=[(float(lon),float(lat)) for lon,lat in coords]
        for ds in datasets:
            if ds.crs is None or not ds.crs.is_geographic:
                raise RuntimeError(f"climate raster is not geographic: {ds.name}")
            values=np.asarray([float(row[0]) for row in ds.sample(tuples)],dtype=float)
            if ds.nodata is not None:
                valid &= ~np.isclose(values,float(ds.nodata),rtol=0.0,atol=0.0)
            valid &= np.isfinite(values)
            cols.append(values)
        return np.column_stack(cols),valid
    finally:
        for ds in datasets:
            ds.close()


def realm_mask(polygon, coords_lon_lat: np.ndarray) -> np.ndarray:
    coords=np.asarray(coords_lon_lat,dtype=float)
    try:
        from shapely import covers_xy
        return np.asarray(covers_xy(polygon,coords[:,0],coords[:,1]),dtype=bool)
    except ImportError:
        from shapely.geometry import Point
        return np.asarray(
            [bool(polygon.covers(Point(float(lon),float(lat)))) for lon,lat in coords],
            dtype=bool,
        )


def realm_grid(polygon, *, resolution: float) -> np.ndarray:
    minx,miny,maxx,maxy=map(float,polygon.bounds)
    i0=max(0,int(math.floor((minx+180.0)/resolution))-1)
    i1=min(int(round(360.0/resolution))-1,int(math.ceil((maxx+180.0)/resolution))+1)
    j0=max(0,int(math.floor((miny+90.0)/resolution))-1)
    j1=min(int(round(180.0/resolution))-1,int(math.ceil((maxy+90.0)/resolution))+1)
    ix=np.arange(i0,i1+1,dtype=np.int64)
    iy=np.arange(j0,j1+1,dtype=np.int64)
    lon=-180.0+(ix+0.5)*resolution
    lat=-90.0+(iy+0.5)*resolution
    xx,yy=np.meshgrid(lon,lat,indexing="xy")
    coords=np.column_stack((xx.ravel(),yy.ravel()))
    return coords[realm_mask(polygon,coords)]


def pairwise_geographic_coop(
    source_coords: np.ndarray,
    target_coords: np.ndarray,
    *,
    radius_km: float=500.0,
) -> float:
    s=np.asarray(source_coords,dtype=float)
    t=np.asarray(target_coords,dtype=float)
    if len(s)==0 or len(t)==0:
        raise ValueError("empty geographic occurrence set")
    # Vectorized haversine in blocks over target points.
    slat=np.radians(s[:,1])
    slon=np.radians(s[:,0])
    radius=6371.0088
    hit=0
    for start in range(0,len(t),128):
        block=t[start:start+128]
        tlat=np.radians(block[:,1])[:,None]
        tlon=np.radians(block[:,0])[:,None]
        dlat=slat[None,:]-tlat
        dlon=slon[None,:]-tlon
        h=np.sin(dlat/2.0)**2 + np.cos(tlat)*np.cos(slat[None,:])*np.sin(dlon/2.0)**2
        d=2.0*radius*np.arcsin(np.minimum(1.0,np.sqrt(h)))
        hit += int(np.count_nonzero(np.min(d,axis=1)<=radius_km))
    return hit/len(t)


def quantiles(values: np.ndarray) -> dict[str,float]:
    q=np.quantile(np.asarray(values,dtype=float),[0,.1,.25,.5,.75,.9,1])
    return {
        key:float(value)
        for key,value in zip(("min","q10","q25","median","q75","q90","max"),q)
    }


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--panel",type=Path,required=True)
    ap.add_argument("--occurrences",type=Path,required=True)
    ap.add_argument("--realm-geojson",type=Path,required=True)
    ap.add_argument("--current-raster",type=Path,action="append",required=True)
    ap.add_argument("--lgm-raster",type=Path,action="append",required=True)
    ap.add_argument("--output-relation",type=Path,required=True)
    ap.add_argument("--output-bootstrap-design",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    args=ap.parse_args()

    contract=json.loads(args.contract.read_text())
    if contract.get("schema")!="ttf_palearctic_insect_lgm_refugia_v0.4":
        raise RuntimeError("unexpected Palearctic insect LGM contract")
    if len(args.current_raster)!=4 or len(args.lgm_raster)!=4:
        raise RuntimeError("exactly four current and four LGM rasters are required")

    panel=json.loads(args.panel.read_text())
    if panel.get("schema")!="ttf_palearctic_insect_lgm_panel_v0.1":
        raise RuntimeError("unexpected formal panel artifact")
    if panel.get("status")!="PASS_RESPONSE_BLIND_PALEARCTIC_INSECT_PANEL":
        raise RuntimeError("formal realm panel did not pass")
    if any(bool(v) for v in panel["response_firewall"].values()):
        raise RuntimeError("response firewall is open")

    from shapely.geometry import shape
    from shapely.ops import unary_union
    realm=json.loads(args.realm_geojson.read_text())
    field=contract["response_blind_panel_eligibility"]["palearctic_realm"]["realm_field"]
    value=contract["response_blind_panel_eligibility"]["palearctic_realm"]["required_value"]
    features=[
        f for f in realm.get("features",[])
        if str(f.get("properties",{}).get(field))==str(value)
    ]
    if not features:
        raise RuntimeError("no formal Palearctic features")
    polygon=unary_union([shape(f["geometry"]) for f in features])

    by_species: dict[str,list[tuple[float,float]]]={}
    with args.occurrences.open(newline="",encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            by_species.setdefault(str(row["species"]),[]).append(
                (float(row["longitude"]),float(row["latitude"]))
            )

    # Target-group background is based on the full frozen 143-species parent
    # occurrence pool, not on the outcome of SDM eligibility.
    pooled=np.asarray(
        [coord for coords in by_species.values() for coord in coords],
        dtype=float,
    )
    pooled=pooled[realm_mask(polygon,pooled)]
    bg_cells=target_group_background_cells(
        pooled,
        resolution_degrees=float(contract["sdm"]["grid_resolution_degrees"]),
        maximum_cells=int(contract["sdm"]["background_cells"]),
    )
    bg_cells=bg_cells[realm_mask(polygon,bg_cells)]
    bg_climate,bg_valid=sample_rasters(bg_cells,args.current_raster)
    bg_cells=bg_cells[bg_valid]
    bg_climate=bg_climate[bg_valid]
    if len(bg_cells)<100:
        raise RuntimeError("too few climate-valid target-group background cells")
    scaling=fit_climate_scaling(bg_climate)
    bg_features=climate_features(bg_climate,scaling)

    grid=realm_grid(
        polygon,
        resolution=float(contract["sdm"]["grid_resolution_degrees"]),
    )
    current_grid,current_valid=sample_rasters(grid,args.current_raster)
    lgm_grid,lgm_valid=sample_rasters(grid,args.lgm_raster)
    grid_valid=current_valid & lgm_valid
    grid=grid[grid_valid]
    current_grid=current_grid[grid_valid]
    lgm_grid=lgm_grid[grid_valid]
    if len(grid)<1000:
        raise RuntimeError("too few climate-valid Palearctic prediction cells")
    current_grid_features=climate_features(current_grid,scaling)
    lgm_grid_features=climate_features(lgm_grid,scaling)

    formal_species=[
        str(row["species"])
        for row in panel["species_ledger"]
        if bool(row["eligible"])
    ]
    ledger_lookup={str(row["species"]):row for row in panel["species_ledger"]}
    min_presence=int(
        contract["sdm_admissibility"]["minimum_finite_current_climate_presence_points"]
    )

    presence_features: dict[str,np.ndarray]={}
    inrealm_coords: dict[str,np.ndarray]={}
    sdm_ledger=[]
    final_species=[]
    for species in formal_species:
        coords=np.asarray(by_species.get(species,[]),dtype=float)
        coords=coords[realm_mask(polygon,coords)]
        climate,valid=sample_rasters(coords,args.current_raster)
        valid_coords=coords[valid]
        valid_climate=climate[valid]
        pass_support=len(valid_coords)>=min_presence
        sdm_ledger.append({
            "species":species,
            "formal_realm_occurrences":int(len(coords)),
            "finite_current_climate_presence_points":int(len(valid_coords)),
            "minimum_required":min_presence,
            "pass_sdm_support":bool(pass_support),
        })
        if pass_support:
            final_species.append(species)
            inrealm_coords[species]=valid_coords
            presence_features[species]=climate_features(valid_climate,scaling)

    minimum_final=int(contract["sdm_admissibility"]["minimum_final_sdm_species"])
    if len(final_species)<minimum_final:
        summary={
            "schema":"ttf_palearctic_insect_lgm_relation_v0.1",
            "status":"NOT_EVALUABLE_PALEARCTIC_INSECT_LGM_SDM_SUPPORT",
            "formal_realm_species":len(formal_species),
            "final_sdm_species":len(final_species),
            "minimum_required":minimum_final,
            "sdm_ledger":sdm_ledger,
            "genetic_response_used":False,
        }
        args.output_summary.parent.mkdir(parents=True,exist_ok=True)
        args.output_summary.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
        print(json.dumps(summary,sort_keys=True))
        return 0

    ordered=palearctic_role_order(final_species)
    cut=len(ordered)//2
    sources=ordered[:cut]
    targets=ordered[cut:]
    species_index={name:i for i,name in enumerate(ordered)}

    fits={}
    current_surface={}
    lgm_surface={}
    for species in ordered:
        fit=fit_weighted_ridge_logistic(
            presence_features[species],
            bg_features,
            ridge_lambda=1.0,
            presence_total_weight=.5,
            background_total_weight=.5,
            max_iter=2000,
        )
        if not fit.converged:
            summary={
                "schema":"ttf_palearctic_insect_lgm_relation_v0.1",
                "status":"NOT_EVALUABLE_PALEARCTIC_INSECT_LGM_SDM_FIT",
                "failed_species":species,
                "final_sdm_species":len(ordered),
                "genetic_response_used":False,
            }
            args.output_summary.parent.mkdir(parents=True,exist_ok=True)
            args.output_summary.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
            print(json.dumps(summary,sort_keys=True))
            return 0
        fits[species]=fit
        current_surface[species]=predict_relative_suitability(current_grid_features,fit)
        lgm_surface[species]=predict_relative_suitability(lgm_grid_features,fit)

    src_idx=[]
    tgt_idx=[]
    r_lgm=[]
    r_current=[]
    coverage=[]
    same_order=[]
    same_family=[]
    count_ratio=[]
    for source in sources:
        for target in targets:
            src_idx.append(species_index[source])
            tgt_idx.append(species_index[target])
            r_lgm.append(schoener_d(lgm_surface[source],lgm_surface[target]))
            r_current.append(schoener_d(current_surface[source],current_surface[target]))
            coverage.append(pairwise_geographic_coop(inrealm_coords[source],inrealm_coords[target]))
            srow=ledger_lookup[source]
            trow=ledger_lookup[target]
            same_order.append(float(
                bool(srow.get("gbif_order"))
                and srow.get("gbif_order")==trow.get("gbif_order")
            ))
            same_family.append(float(
                bool(srow.get("gbif_family"))
                and srow.get("gbif_family")==trow.get("gbif_family")
            ))
            count_ratio.append(len(inrealm_coords[source])/len(inrealm_coords[target]))

    src_idx=np.asarray(src_idx,dtype=np.int64)
    tgt_idx=np.asarray(tgt_idx,dtype=np.int64)
    r_lgm=np.asarray(r_lgm,dtype=float)
    r_current=np.asarray(r_current,dtype=float)
    coverage=np.asarray(coverage,dtype=float)
    same_order=np.asarray(same_order,dtype=float)
    same_family=np.asarray(same_family,dtype=float)
    count_ratio=np.asarray(count_ratio,dtype=float)

    args.output_relation.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_relation,
        species_order=np.asarray(ordered,dtype="U160"),
        source_species=np.asarray(sources,dtype="U160"),
        target_species=np.asarray(targets,dtype="U160"),
        source_index=src_idx,
        target_index=tgt_idx,
        R_LGM=r_lgm,
        R_current=r_current,
        geographic_coop=coverage,
        same_order=same_order,
        same_family=same_family,
        presence_count_ratio=count_ratio,
    )

    concatenated=np.vstack([presence_features[name] for name in ordered])
    offsets=[0]
    for name in ordered:
        offsets.append(offsets[-1]+len(presence_features[name]))
    coefficients=np.vstack([fits[name].coefficients for name in ordered])
    np.savez_compressed(
        args.output_bootstrap_design,
        species_order=np.asarray(ordered,dtype="U160"),
        source_index=src_idx,
        target_index=tgt_idx,
        source_species=np.asarray(sources,dtype="U160"),
        target_species=np.asarray(targets,dtype="U160"),
        presence_features=concatenated,
        presence_offsets=np.asarray(offsets,dtype=np.int64),
        background_features=bg_features,
        current_grid_features=current_grid_features,
        lgm_grid_features=lgm_grid_features,
        fitted_coefficients=coefficients,
    )

    summary={
        "schema":"ttf_palearctic_insect_lgm_relation_v0.1",
        "status":"PASS_RESPONSE_BLIND_PALEARCTIC_INSECT_LGM_RELATION",
        "formal_realm_species":len(formal_species),
        "final_sdm_species":len(ordered),
        "source_count":len(sources),
        "target_count":len(targets),
        "directed_dyads":len(r_lgm),
        "target_group_background_cells":len(bg_cells),
        "prediction_grid_cells":len(grid),
        "R_LGM_quantiles":quantiles(r_lgm),
        "R_current_quantiles":quantiles(r_current),
        "R_LGM_current_correlation":float(np.corrcoef(r_lgm,r_current)[0,1]),
        "sdm_ledger":sdm_ledger,
        "inputs":{
            "contract_sha256":sha256_path(args.contract),
            "panel_sha256":sha256_path(args.panel),
            "occurrences_sha256":sha256_path(args.occurrences),
            "realm_geojson_sha256":sha256_path(args.realm_geojson),
            "current_rasters_sha256":{p.name:sha256_path(p) for p in args.current_raster},
            "lgm_rasters_sha256":{p.name:sha256_path(p) for p in args.lgm_raster},
        },
        "outputs":{
            "relation_npz_sha256":sha256_path(args.output_relation),
            "bootstrap_design_npz_sha256":sha256_path(args.output_bootstrap_design),
        },
        "response_firewall":{
            "species_level_Phase4_scores_used":False,
            "pairwise_genetic_distances_used":False,
            "subpanel_T_st_computed":False,
            "subpanel_beta_LGM_computed":False,
        },
        "next_step":"run 100 frozen occurrence-bootstrap relation replicates and TTF-Q qualification before any subpanel genetic response",
    }
    args.output_summary.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":summary["status"],
        "final_sdm_species":len(ordered),
        "sources":len(sources),
        "targets":len(targets),
        "dyads":len(r_lgm),
        "R_LGM_current_correlation":summary["R_LGM_current_correlation"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
