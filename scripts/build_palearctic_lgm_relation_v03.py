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
        for chunk in iter(lambda: fh.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


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
    ap.add_argument("--climate-rule",type=Path,required=True)
    ap.add_argument("--output-npz",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    ap.add_argument("--chunk-size",type=int,default=4096)
    args=ap.parse_args()

    rule=json.loads(args.rule.read_text())
    if rule.get("schema")!="ttf_genetic_palearctic_lgm_subpanel_rule_v0.3":
        raise RuntimeError("unexpected v0.3 parent rule")
    climate=json.loads(args.climate_rule.read_text())
    if climate.get("schema")!="ttf_genetic_palearctic_lgm_climate_input_rule_v0.4":
        raise RuntimeError("unexpected v0.3 climate rule")
    panel=json.loads(args.panel_metadata.read_text())
    if panel.get("schema")!="ttf_genetic_palearctic_lgm_panel_metadata_v0.3":
        raise RuntimeError("unexpected v0.3 panel metadata")
    if any(bool(v) for v in panel["response_firewall"].values()):
        raise RuntimeError("v0.3 panel response firewall is open")

    meta={str(row["species"]):dict(row) for row in panel["rows"]}
    formal_species=sorted(meta)
    formal_sources=sorted(sp for sp,row in meta.items() if row["role"]=="source")
    formal_targets=sorted(sp for sp,row in meta.items() if row["role"]=="target")
    if len(formal_species)!=23 or len(formal_sources)!=11 or len(formal_targets)!=12:
        raise RuntimeError("v0.3 formal source-target geometry drift")
    if set(formal_sources)&set(formal_targets):
        raise RuntimeError("v0.3 source-target roles must be disjoint")

    occurrences=read_occurrence_climate(args.occurrence_climate_csv,set(formal_species))
    min_occ=int(climate["occurrence_climate_validity"]["minimum_rows_after_current_climate_validity_per_species"])
    species=sorted(sp for sp in formal_species if len(occurrences.get(sp,[]))>=min_occ)
    species_set=set(species)
    sources=[sp for sp in formal_sources if sp in species_set]
    targets=[sp for sp in formal_targets if sp in species_set]
    req=climate["occurrence_climate_validity"]["after_drop_require"]
    passed=(
        len(species)>=int(req["minimum_total_species"])
        and len(sources)>=int(req["minimum_source_clusters"])
        and len(targets)>=int(req["minimum_target_clusters"])
    )
    if not passed:
        payload={
            "schema":"ttf_genetic_palearctic_lgm_relation_design_v0.3",
            "status":str(climate["occurrence_climate_validity"]["failure"]),
            "formal_species":23,
            "climate_valid_species":len(species),
            "source_clusters":len(sources),
            "target_clusters":len(targets),
            "response_firewall":{
                "species_level_genetic_scores_used":False,
                "pairwise_v03_T_st_computed":False,
                "beta_LGM_computed":False,
            },
        }
        args.output_summary.parent.mkdir(parents=True,exist_ok=True)
        args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
        print(json.dumps({"status":payload["status"],"climate_valid_species":len(species)},sort_keys=True))
        return 0

    pooled=np.vstack([np.vstack([item[2] for item in occurrences[sp]]) for sp in species])
    space=fit_common_climate_space(pooled,axes=4)

    grid=np.load(args.grid_npz,allow_pickle=False)
    current_grid=np.asarray(grid["current_environment"],dtype=float)
    lgm_grid=np.asarray(grid["lgm_environment"],dtype=float)
    if current_grid.shape!=lgm_grid.shape or current_grid.ndim!=2 or current_grid.shape[1]!=4:
        raise RuntimeError("current/LGM grid environments must be equal-shape n x 4")
    if not np.isfinite(current_grid).all() or not np.isfinite(lgm_grid).all():
        raise RuntimeError("climate grid contains non-finite values")

    current_density={}; lgm_density={}
    current_extrap=extrapolation_fraction(current_grid,pooled,space)
    lgm_extrap=extrapolation_fraction(lgm_grid,pooled,space)
    for sp in species:
        env=np.vstack([item[2] for item in occurrences[sp]])
        now,past=current_lgm_suitability(
            env,current_grid,lgm_grid,space,chunk_size=args.chunk_size
        )
        current_density[sp]=now.astype(np.float32)
        lgm_density[sp]=past.astype(np.float32)

    occurrence_latlon={
        sp:np.asarray([[x[0],x[1]] for x in occurrences[sp]],dtype=float)
        for sp in species
    }
    counts={sp:len(occurrences[sp]) for sp in species}

    src_out=[]; tgt_out=[]; r_lgm=[]; r_present=[]; coverage=[]
    same_order=[]; same_family=[]; count_ratio=[]
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

    payload={
        "schema":"ttf_genetic_palearctic_lgm_relation_design_v0.3",
        "status":"PASS_TO_V03_TTF_Q_RESPONSE_BLIND_CHARACTERIZATION",
        "formal_species":23,
        "species":len(species),
        "source_clusters":len(sources),
        "target_clusters":len(targets),
        "dyads":len(r_lgm),
        "dropped_before_relation":sorted(set(formal_species)-species_set),
        "R_LGM":{"min":float(r_lgm.min()),"median":float(np.median(r_lgm)),"max":float(r_lgm.max())},
        "R_present":{"min":float(r_present.min()),"median":float(np.median(r_present)),"max":float(r_present.max())},
        "pearson_R_LGM_vs_R_present":float(np.corrcoef(r_lgm,r_present)[0,1]),
        "climate_grid_cells":int(len(current_grid)),
        "pooled_current_occurrences":int(len(pooled)),
        "current_grid_extrapolation_fraction":current_extrap,
        "lgm_grid_extrapolation_fraction":lgm_extrap,
        "inputs":{
            "panel_metadata_sha256":sha256_path(args.panel_metadata),
            "occurrence_climate_csv_sha256":sha256_path(args.occurrence_climate_csv),
            "grid_npz_sha256":sha256_path(args.grid_npz),
            "rule_sha256":sha256_path(args.rule),
            "climate_rule_sha256":sha256_path(args.climate_rule),
        },
        "design_npz_sha256":sha256_path(args.output_npz),
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_v03_T_st_computed":False,
            "beta_LGM_computed":False,
        },
        "next_step":"Run frozen v0.3 TTF-Q. Do not compute subgroup genetic response unless its opening rule passes.",
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "species":len(species),
        "source_clusters":len(sources),
        "target_clusters":len(targets),
        "dyads":len(r_lgm),
        "pearson_R_LGM_vs_R_present":payload["pearson_R_LGM_vs_R_present"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
