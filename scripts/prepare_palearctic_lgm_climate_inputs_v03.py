#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

try:
    from scripts.prepare_palearctic_lgm_climate_inputs import (
        VARIABLES, compatibility, load_union, palearctic_grid_points,
        sample_environment, sha256_path,
    )
except ModuleNotFoundError:
    from prepare_palearctic_lgm_climate_inputs import (
        VARIABLES, compatibility, load_union, palearctic_grid_points,
        sample_environment, sha256_path,
    )


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
    if rule.get("schema")!="ttf_genetic_palearctic_lgm_climate_input_rule_v0.3":
        raise RuntimeError("unexpected v0.3 climate rule")
    panel=json.loads(args.panel_metadata.read_text())
    if panel.get("schema")!="ttf_genetic_palearctic_lgm_panel_metadata_v0.3":
        raise RuntimeError("unexpected v0.3 panel metadata")
    if any(bool(v) for v in panel["response_firewall"].values()):
        raise RuntimeError("v0.3 panel response firewall is open")
    meta={str(row["species"]):dict(row) for row in panel["rows"]}
    if len(meta)!=23:
        raise RuntimeError("v0.3 panel metadata species drift")

    occ=list(csv.DictReader(args.occurrences.open(encoding="utf-8")))
    occ_points=[]; occ_rows=[]
    for row in occ:
        sp=str(row.get("species",""))
        if sp not in meta:
            continue
        try:
            lat=float(row["latitude"]); lon=float(row["longitude"])
        except (TypeError,ValueError):
            continue
        row=dict(row); row["role"]=str(meta[sp]["role"])
        occ_points.append((lat,lon)); occ_rows.append(row)
    if not occ_points:
        raise RuntimeError("no v0.3 occurrences available for climate sampling")
    occ_points=np.asarray(occ_points,dtype=float)

    current_paths={var:getattr(args,f"current_{var}") for var in VARIABLES}
    lgm_paths={var:getattr(args,f"lgm_{var}") for var in VARIABLES}
    trace0_paths={var:getattr(args,f"trace0_{var}") for var in VARIABLES}
    occ_env,current_meta=sample_environment(current_paths,occ_points,dataset="current")
    valid_occ=np.isfinite(occ_env).all(axis=1)

    fields=["species","role","source_key","latitude","longitude","priority_rank","priority_sha256",*VARIABLES]
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
    min_rows=int(rule["occurrence_climate_validity"]["minimum_rows_after_current_climate_validity_per_species"])
    valid_species=sorted(sp for sp,n in counts.items() if n>=min_rows)
    source=[sp for sp in valid_species if meta[sp]["role"]=="source"]
    target=[sp for sp in valid_species if meta[sp]["role"]=="target"]
    req=rule["occurrence_climate_validity"]["after_drop_require"]
    coverage_pass=(
        len(valid_species)>=int(req["minimum_total_species"])
        and len(source)>=int(req["minimum_source_clusters"])
        and len(target)>=int(req["minimum_target_clusters"])
    )
    if not coverage_pass:
        payload={
            "schema":"ttf_genetic_palearctic_lgm_climate_inputs_v0.3",
            "status":str(rule["occurrence_climate_validity"]["failure"]),
            "formal_panel_species":23,
            "climate_valid_species":len(valid_species),
            "source_clusters":len(source),
            "target_clusters":len(target),
            "climate_valid_counts":counts,
            "response_firewall":{
                "species_level_genetic_scores_used":False,
                "pairwise_v03_T_st_computed":False,
                "beta_LGM_computed":False,
            },
        }
        args.output_summary.parent.mkdir(parents=True,exist_ok=True)
        args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
        print(json.dumps({"status":payload["status"],"climate_valid_species":len(valid_species)},sort_keys=True))
        return 0

    valid_set=set(valid_species)
    climate_rows=[row for row in climate_rows if row["species"] in valid_set]
    climate_rows.sort(key=lambda r:(r["species"],int(r["priority_rank"]),int(r["source_key"])))
    args.output_occurrence_climate.parent.mkdir(parents=True,exist_ok=True)
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
        "schema":"ttf_genetic_palearctic_lgm_climate_inputs_v0.3",
        "status":"PASS_TO_V03_RESPONSE_BLIND_LGM_RELATION_BUILD",
        "formal_panel_species":23,
        "climate_valid_species":len(valid_species),
        "source_clusters":len(source),
        "target_clusters":len(target),
        "directed_dyads":len(source)*len(target),
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
        "raster_metadata":{"current":current_meta,"lgm":lgm_meta,"trace0":trace0_meta},
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
            "pairwise_v03_T_st_computed":False,
            "beta_LGM_computed":False,
        },
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "climate_valid_species":len(valid_species),
        "source_clusters":len(source),
        "target_clusters":len(target),
        "grid_cells":int(grid_valid.sum()),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
