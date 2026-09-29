#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.external_geographic_opportunity import directed_occurrence_coverage
from ttf.lgm_climate_envelope import (
    climatic_relation,
    current_lgm_suitability,
    extrapolation_fraction,
    fit_common_climate_space,
)


CLIMATE_FIELDS=("bio1","bio7","bio12","bio15")


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20),b""):
            h.update(chunk)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    values=np.asarray(x,dtype=float)
    sd=float(values.std())
    if not np.isfinite(sd) or sd<=np.finfo(float).eps:
        raise ValueError("cannot z-score constant/non-finite predictor")
    return (values-float(values.mean()))/sd


def read_panel_metadata(path: Path):
    payload=json.loads(path.read_text())
    if payload.get("schema")!="ttf_genetic_palearctic_lgm_panel_metadata_v0.1":
        raise RuntimeError("unexpected frozen panel metadata")
    if any(bool(v) for v in payload["response_firewall"].values()):
        raise RuntimeError("panel metadata response firewall is open")
    rows=list(payload["rows"])
    meta={str(row["species"]):dict(row) for row in rows}
    sources=sorted(sp for sp,row in meta.items() if row["role"]=="source")
    targets=sorted(sp for sp,row in meta.items() if row["role"]=="target")
    if len(meta)!=44 or len(sources)!=22 or len(targets)!=22:
        raise RuntimeError("frozen 44-species role geometry drift")
    if set(sources)&set(targets):
        raise RuntimeError("source/target role overlap")
    return sources,targets,meta


def read_occurrence_climate(path: Path, retained: set[str]):
    rows=list(csv.DictReader(path.open(encoding="utf-8")))
    by_species={name:[] for name in retained}
    for row in rows:
        sp=str(row.get("species",""))
        if sp not in by_species:
            continue
        try:
            lat=float(row["latitude"]); lon=float(row["longitude"])
            env=np.asarray([float(row[k]) for k in CLIMATE_FIELDS],dtype=float)
        except (KeyError,TypeError,ValueError):
            continue
        if not np.isfinite(env).all() or not np.isfinite(lat) or not np.isfinite(lon):
            continue
        by_species[sp].append((lat,lon,env))
    return by_species


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel-metadata",type=Path,required=True)
    ap.add_argument("--occurrence-climate-csv",type=Path,required=True)
    ap.add_argument("--grid-npz",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--external-rule",type=Path,required=True)
    ap.add_argument("--climate-rule",type=Path,required=True)
    ap.add_argument("--output-npz",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    ap.add_argument("--chunk-size",type=int,default=4096)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!="ttf_genetic_palearctic_lgm_subpanel_rule_v0.1":
        raise RuntimeError("unexpected Palearctic-LGM rule")
    external=json.loads(args.external_rule.read_text())
    if external.get("schema")!="ttf_genetic_palearctic_lgm_external_data_rule_v0.1":
        raise RuntimeError("unexpected external-data rule")
    climate=json.loads(args.climate_rule.read_text())
    if climate.get("schema")!="ttf_genetic_palearctic_lgm_climate_input_rule_v0.1":
        raise RuntimeError("unexpected climate-input rule")

    formal_sources,formal_targets,meta=read_panel_metadata(args.panel_metadata)
    formal_species=sorted(set(formal_sources)|set(formal_targets))

    occurrences=read_occurrence_climate(args.occurrence_climate_csv,set(formal_species))
    min_occ=int(climate["occurrence_climate_validity"]["minimum_rows_after_climate_validity_per_species"])
    species=sorted(
        sp for sp in formal_species
        if len(occurrences.get(sp,[]))>=min_occ
    )
    sources=[sp for sp in formal_sources if sp in set(species)]
    targets=[sp for sp in formal_targets if sp in set(species)]
    req=climate["occurrence_climate_validity"]["after_drop_require"]
    passed=(
        len(species)>=int(req["minimum_total_species"])
        and len(sources)>=int(req["minimum_source_clusters"])
        and len(targets)>=int(req["minimum_target_clusters"])
    )
    if not passed:
        payload={
            "schema":"ttf_genetic_palearctic_lgm_relation_design_v0.1",
            "status":str(climate["occurrence_climate_validity"]["failure"]),
            "formal_species":44,
            "climate_valid_species":len(species),
            "source_clusters":len(sources),
            "target_clusters":len(targets),
            "failed_species":{
                sp:len(occurrences.get(sp,[]))
                for sp in formal_species if sp not in set(species)
            },
            "response_firewall":{
                "species_level_genetic_scores_used":False,
                "pairwise_subpanel_T_st_computed":False,
                "beta_LGM_computed":False,
            },
        }
        args.output_summary.parent.mkdir(parents=True,exist_ok=True)
        args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
        print(json.dumps({"status":payload["status"],"climate_valid_species":len(species)},sort_keys=True))
        return 0

    pooled=np.vstack([
        np.vstack([item[2] for item in occurrences[sp]])
        for sp in species
    ])
    space=fit_common_climate_space(pooled,axes=4)

    grid=np.load(args.grid_npz,allow_pickle=False)
    current_grid=np.asarray(grid["current_environment"],dtype=float)
    lgm_grid=np.asarray(grid["lgm_environment"],dtype=float)
    if current_grid.shape!=lgm_grid.shape or current_grid.ndim!=2 or current_grid.shape[1]!=4:
        raise RuntimeError("current/LGM grid environments must be equal-shape n x 4")
    if not np.isfinite(current_grid).all() or not np.isfinite(lgm_grid).all():
        raise RuntimeError("climate grid contains non-finite values")

    current_density={}
    lgm_density={}
    species_extrapolation={}
    for sp in species:
        env=np.vstack([item[2] for item in occurrences[sp]])
        now,past=current_lgm_suitability(
            env,current_grid,lgm_grid,space,chunk_size=args.chunk_size
        )
        current_density[sp]=now.astype(np.float32)
        lgm_density[sp]=past.astype(np.float32)
        species_extrapolation[sp]={
            "current_grid_outside_pooled_current_range":extrapolation_fraction(
                current_grid,pooled,space
            ),
            "lgm_grid_outside_pooled_current_range":extrapolation_fraction(
                lgm_grid,pooled,space
            ),
        }

    src_out=[]; tgt_out=[]; r_lgm=[]; r_present=[]; coverage=[]
    same_order=[]; same_family=[]; count_ratio=[]
    occurrence_latlon={
        sp:np.asarray([[x[0],x[1]] for x in occurrences[sp]],dtype=float)
        for sp in species
    }
    counts={sp:len(occurrences[sp]) for sp in species}

    for s in sources:
        for t in targets:
            src_out.append(species.index(s))
            tgt_out.append(species.index(t))
            r_lgm.append(climatic_relation(lgm_density[s],lgm_density[t]))
            r_present.append(climatic_relation(current_density[s],current_density[t]))
            coverage.append(directed_occurrence_coverage(
                occurrence_latlon[t],occurrence_latlon[s],radius_km=500.0
            ))
            same_order.append(float(meta[s]["order"]==meta[t]["order"]))
            same_family.append(float(meta[s]["family"]==meta[t]["family"]))
            count_ratio.append(float(counts[s]/counts[t]))

    src_out=np.asarray(src_out,np.int64)
    tgt_out=np.asarray(tgt_out,np.int64)
    r_lgm=np.asarray(r_lgm,float)
    r_present=np.asarray(r_present,float)
    coverage=np.asarray(coverage,float)
    same_order=np.asarray(same_order,float)
    same_family=np.asarray(same_family,float)
    count_ratio=np.asarray(count_ratio,float)

    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        species_order=np.asarray(species,dtype="U160"),
        source_species=np.asarray(sources,dtype="U160"),
        target_species=np.asarray(targets,dtype="U160"),
        source_index=src_out,
        target_index=tgt_out,
        R_LGM=r_lgm,
        R_present=r_present,
        geographic_coverage=coverage,
        same_order=same_order,
        same_family=same_family,
        retained_occurrence_count_ratio=count_ratio,
        climate_mean=space.mean,
        climate_sd=space.sd,
        climate_eigenvalues=space.eigval,
        climate_eigenvectors=space.eigvec,
    )

    corr=float(np.corrcoef(r_lgm,r_present)[0,1])
    payload={
        "schema":"ttf_genetic_palearctic_lgm_relation_design_v0.1",
        "status":"PASS_TO_TTF_Q_RESPONSE_BLIND_CHARACTERIZATION",
        "formal_species":44,
        "species":len(species),
        "source_clusters":len(sources),
        "target_clusters":len(targets),
        "dropped_before_relation":sorted(set(formal_species)-set(species)),
        "dyads":len(r_lgm),
        "R_LGM":{
            "min":float(r_lgm.min()),
            "median":float(np.median(r_lgm)),
            "max":float(r_lgm.max()),
        },
        "R_present":{
            "min":float(r_present.min()),
            "median":float(np.median(r_present)),
            "max":float(r_present.max()),
        },
        "pearson_R_LGM_vs_R_present":corr,
        "climate_grid_cells":int(len(current_grid)),
        "pooled_current_occurrences":int(len(pooled)),
        "maximum_lgm_grid_extrapolation_fraction":float(max(
            row["lgm_grid_outside_pooled_current_range"]
            for row in species_extrapolation.values()
        )),
        "inputs":{
            "panel_metadata_sha256":sha256_path(args.panel_metadata),
            "occurrence_climate_csv_sha256":sha256_path(args.occurrence_climate_csv),
            "grid_npz_sha256":sha256_path(args.grid_npz),
            "rule_sha256":sha256_path(args.rule),
            "external_rule_sha256":sha256_path(args.external_rule),
            "climate_rule_sha256":sha256_path(args.climate_rule),
        },
        "design_npz_sha256":sha256_path(args.output_npz),
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_subpanel_T_st_computed":False,
            "beta_LGM_computed":False,
        },
        "next_step":"Run TTF-Q information survival, concentration, null qualification and detectability. Do not compute subgroup genetic response unless the frozen opening rule passes.",
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "species":payload["species"],
        "dyads":payload["dyads"],
        "pearson_R_LGM_vs_R_present":corr,
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
