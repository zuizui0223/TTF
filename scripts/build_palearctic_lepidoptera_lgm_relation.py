#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from ttf.external_geographic_opportunity import directed_occurrence_coverage
from ttf.palearctic_lgm_sdm import (
    deterministic_background_indices,
    ordered_nonself_pairs,
    palearctic_core_grid,
    quadratic_features,
    schoener_d,
    standardize_from_background,
    within_occurrence_buffer,
)


def sample_rasters(latlon: np.ndarray, paths: list[Path]) -> tuple[np.ndarray,np.ndarray]:
    import rasterio

    latlon=np.asarray(latlon,dtype=float)
    coords=[(float(lon),float(lat)) for lat,lon in latlon]
    datasets=[rasterio.open(p) for p in paths]
    try:
        valid=np.ones(len(coords),dtype=bool)
        cols=[]
        for ds in datasets:
            if ds.crs is None or not ds.crs.is_geographic:
                raise RuntimeError(f"CHELSA raster is not geographic: {ds.name}")
            values=np.asarray([float(x[0]) for x in ds.sample(coords)],dtype=float)
            if ds.nodata is not None:
                valid &= ~np.isclose(values,float(ds.nodata),rtol=0.0,atol=0.0)
            valid &= np.isfinite(values)
            cols.append(values)
        return np.column_stack(cols),valid
    finally:
        for ds in datasets:
            ds.close()


def quantiles(x: np.ndarray) -> dict[str,float]:
    q=np.quantile(np.asarray(x,dtype=float),[0,.1,.25,.5,.75,.9,1])
    keys=("min","q10","q25","median","q75","q90","max")
    return {k:float(v) for k,v in zip(keys,q)}


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--occurrences",type=Path,required=True)
    ap.add_argument("--occurrence-execution",type=Path,required=True)
    ap.add_argument("--screen",type=Path,required=True)
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--current-raster",type=Path,action="append",required=True)
    ap.add_argument("--lgm-raster",type=Path,action="append",required=True)
    ap.add_argument("--output-npz",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()

    if len(args.current_raster)!=4 or len(args.lgm_raster)!=4:
        raise RuntimeError("exactly four current and four LGM rasters are required")
    contract=json.loads(args.contract.read_text())
    screen=json.loads(args.screen.read_text())
    execution=json.loads(args.occurrence_execution.read_text())
    if contract.get("schema")!="ttf_palearctic_lepidoptera_lgm_sdm_rule_v0.1":
        raise RuntimeError("unexpected LGM SDM rule")
    if execution.get("status")!="PASS_TO_LGM_CLIMATE_EXTRACTION":
        raise RuntimeError("occurrence execution did not authorize climate extraction")
    if any(bool(v) for v in execution["response_firewall"].values()):
        raise RuntimeError("occurrence response firewall is open")

    meta={str(x["species"]):dict(x) for x in screen["primary_species"]}
    passing=list(map(str,execution["predictor_admissible_species_names"]))
    if not set(passing).issubset(meta):
        raise RuntimeError("occurrence-admissible species outside frozen primary panel")

    rows=list(csv.DictReader(args.occurrences.open(newline="",encoding="utf-8")))
    by_species=defaultdict(list)
    for row in rows:
        if row["species"] in passing:
            by_species[row["species"]].append(row)
    for name in by_species:
        by_species[name].sort(key=lambda r:(int(r["priority_rank"]),int(r["source_key"])))

    bbox=screen["screen"]["sampling_domain"]
    step=float(contract["pairwise_primary_relation"]["grid_resolution_degrees"])
    full_grid=palearctic_core_grid(
        lat_min=float(bbox["latitude_min"]),
        lat_max=float(bbox["latitude_max"]),
        lon_min=float(bbox["longitude_min"]),
        lon_max=float(bbox["longitude_max"]),
        step_degrees=step,
    )
    current_grid_raw,current_valid=sample_rasters(full_grid,list(args.current_raster))
    lgm_grid_raw,lgm_valid=sample_rasters(full_grid,list(args.lgm_raster))
    grid_valid=current_valid & lgm_valid
    grid=full_grid[grid_valid]
    current_grid=current_grid_raw[grid_valid]
    lgm_grid=lgm_grid_raw[grid_valid]
    if len(grid)<1000:
        raise RuntimeError("too few finite Palearctic grid cells")

    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    suitability_current={}
    suitability_lgm={}
    occurrence_n={}
    background_n={}
    auc={}
    failed={}
    retained_occ={}
    max_bg=int(contract["sdm"]["maximum_background_cells"])
    min_bg=int(contract["sdm"]["minimum_background_cells"])
    C=float(contract["sdm"]["ridge_logistic_C"])

    for species in passing:
        srows=by_species.get(species,[])
        latlon=np.asarray(
            [[float(r["latitude"]),float(r["longitude"])] for r in srows],
            dtype=float,
        )
        occurrence_n[species]=int(len(latlon))
        if len(latlon)<int(contract["ecological_occurrences"]["minimum_retained_per_species"]):
            failed[species]="too_few_occurrences_after_transport"
            continue

        positive_raw,positive_valid=sample_rasters(latlon,list(args.current_raster))
        positive=positive_raw[positive_valid]
        valid_latlon=latlon[positive_valid]
        if len(positive)<int(contract["ecological_occurrences"]["minimum_retained_per_species"]):
            failed[species]="too_few_occurrences_after_current_climate_sampling"
            continue

        accessible=within_occurrence_buffer(
            grid,
            valid_latlon,
            radius_km=500.0,
        )
        bg_idx=deterministic_background_indices(
            species,grid,accessible,maximum=max_bg
        )
        if len(bg_idx)<min_bg:
            failed[species]="too_few_accessible_background_cells"
            continue
        background=current_grid[bg_idx]
        try:
            z_pos,z_bg,z_cur,z_lgm,mu,sd=standardize_from_background(
                background,positive,background,current_grid,lgm_grid
            )
        except ValueError:
            failed[species]="degenerate_background_climate"
            continue

        X_pos=quadratic_features(z_pos)
        X_bg=quadratic_features(z_bg)
        X=np.vstack((X_pos,X_bg))
        y=np.concatenate((np.ones(len(X_pos),dtype=int),np.zeros(len(X_bg),dtype=int)))
        model=LogisticRegression(
            C=C,
            penalty="l2",
            solver="lbfgs",
            class_weight="balanced",
            max_iter=2000,
        )
        try:
            model.fit(X,y)
            cur=model.predict_proba(quadratic_features(z_cur))[:,1]
            lgm=model.predict_proba(quadratic_features(z_lgm))[:,1]
            train_prob=model.predict_proba(X)[:,1]
        except Exception as exc:
            failed[species]=f"model_fit_error:{type(exc).__name__}"
            continue
        if (
            np.any(~np.isfinite(cur)) or np.any(~np.isfinite(lgm))
            or float(np.std(cur))<=1e-12 or float(np.std(lgm))<=1e-12
        ):
            failed[species]="nonfinite_or_constant_prediction"
            continue

        suitability_current[species]=cur.astype(np.float64)
        suitability_lgm[species]=lgm.astype(np.float64)
        background_n[species]=int(len(bg_idx))
        auc[species]=float(roc_auc_score(y,train_prob))
        retained_occ[species]=valid_latlon

    modeled=sorted(suitability_lgm)
    minimum=int(contract["ecological_occurrences"]["minimum_predictor_admissible_species"])
    if len(modeled)<minimum:
        payload={
            "schema":"ttf_palearctic_lepidoptera_lgm_relation_v0.1",
            "status":"STOP_LGM_SDM_MODEL_PANEL_TOO_SMALL",
            "occurrence_admissible_species":len(passing),
            "modeled_species":len(modeled),
            "minimum_required":minimum,
            "failed_species":failed,
            "response_firewall":{
                "subpanel_genetic_response_used":False,
            },
        }
        args.output_json.parent.mkdir(parents=True,exist_ok=True)
        args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
        print(json.dumps(payload,sort_keys=True))
        return 0

    pairs=ordered_nonself_pairs(modeled)
    index={name:i for i,name in enumerate(modeled)}
    source_idx=[]; target_idx=[]
    R_lgm=[]; R_current=[]; geo=[]; same_family=[]; log_ratio=[]
    for source,target in pairs:
        source_idx.append(index[source]); target_idx.append(index[target])
        R_lgm.append(schoener_d(suitability_lgm[source],suitability_lgm[target]))
        R_current.append(schoener_d(suitability_current[source],suitability_current[target]))
        geo.append(directed_occurrence_coverage(
            retained_occ[target],retained_occ[source],radius_km=500.0
        ))
        same_family.append(float(meta[source]["family"]==meta[target]["family"]))
        log_ratio.append(abs(float(np.log(occurrence_n[source]/occurrence_n[target]))))

    current_matrix=np.vstack([suitability_current[name] for name in modeled])
    lgm_matrix=np.vstack([suitability_lgm[name] for name in modeled])

    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        species_order=np.asarray(modeled,dtype="U160"),
        family=np.asarray([meta[x]["family"] for x in modeled],dtype="U96"),
        occurrence_n=np.asarray([occurrence_n[x] for x in modeled],dtype=np.int64),
        grid_latlon=grid.astype(np.float32),
        current_suitability=current_matrix.astype(np.float32),
        lgm_suitability=lgm_matrix.astype(np.float32),
        source_index=np.asarray(source_idx,dtype=np.int64),
        target_index=np.asarray(target_idx,dtype=np.int64),
        R_LGM=np.asarray(R_lgm,dtype=float),
        R_current=np.asarray(R_current,dtype=float),
        geographic_coverage=np.asarray(geo,dtype=float),
        same_family=np.asarray(same_family,dtype=float),
        abs_log_occurrence_ratio=np.asarray(log_ratio,dtype=float),
    )

    payload={
        "schema":"ttf_palearctic_lepidoptera_lgm_relation_v0.1",
        "status":"PASS_RESPONSE_BLIND_LGM_RELATION",
        "occurrence_admissible_species":len(passing),
        "modeled_species":len(modeled),
        "modeled_species_names":modeled,
        "failed_species":failed,
        "grid_cells":int(len(grid)),
        "grid_resolution_degrees":step,
        "ordered_nonself_dyads":len(pairs),
        "per_species":{
            name:{
                "occurrences":occurrence_n[name],
                "background_cells":background_n[name],
                "descriptive_resubstitution_auc":auc[name],
            }
            for name in modeled
        },
        "R_LGM_quantiles":quantiles(np.asarray(R_lgm)),
        "R_current_quantiles":quantiles(np.asarray(R_current)),
        "geographic_coverage_quantiles":quantiles(np.asarray(geo)),
        "response_firewall":{
            "subpanel_nucleotide_identity_read":False,
            "subpanel_pairwise_genetic_distance_read":False,
            "subpanel_transfer_response_constructed":False,
        },
        "next_step":"Run TTF-Q information decomposition and calibrated detectability on this exact predictor design before any subpanel genetic response.",
    }
    args.output_json.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "modeled_species":len(modeled),
        "dyads":len(pairs),
        "R_LGM":payload["R_LGM_quantiles"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
