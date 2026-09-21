#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from ttf.relational_dyadic import prepare_dyadic_regression


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    x=np.asarray(x,float)
    sd=float(x.std())
    if not np.isfinite(sd) or sd<=np.finfo(float).eps:
        raise ValueError("constant/non-finite control")
    return (x-float(x.mean()))/sd


def nearest_coverage(target: np.ndarray, source: np.ndarray, radius: float=500.0, chunk: int=256) -> float:
    t=np.asarray(target,float); s=np.asarray(source,float)
    if t.ndim!=2 or s.ndim!=2 or t.shape[1]!=3 or s.shape[1]!=3 or not len(t) or not len(s):
        raise ValueError("edge midpoints must be non-empty n x 3")
    hit=0
    r2=float(radius)**2
    for start in range(0,len(t),chunk):
        block=t[start:start+chunk]
        delta=block[:,None,:]-s[None,:,:]
        nearest=np.min(np.sum(delta*delta,axis=2),axis=1)
        hit += int(np.count_nonzero(nearest <= r2))
    return hit/len(t)


def load_edges(path: Path) -> dict[str,np.ndarray]:
    grouped=defaultdict(list)
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            grouped[row["species"]].append([
                float(row["mid_x_km"]),float(row["mid_y_km"]),float(row["mid_z_km"])
            ])
    return {name:np.asarray(values,float) for name,values in grouped.items()}


def load_metadata(path: Path) -> dict[str,dict[str,str]]:
    rows=list(csv.DictReader(path.open(encoding="utf-8")))
    return {str(r["species"]):r for r in rows}


def remap(values: np.ndarray) -> np.ndarray:
    unique=sorted(map(int,np.unique(values)))
    lookup={v:i for i,v in enumerate(unique)}
    return np.asarray([lookup[int(v)] for v in values],dtype=np.int64)


def spearman(x: np.ndarray,y: np.ndarray) -> float:
    def ranks(a):
        a=np.asarray(a,float); order=np.argsort(a,kind="stable"); out=np.empty(len(a),float)
        sx=a[order]; starts=np.r_[0,1+np.flatnonzero(sx[1:]!=sx[:-1])]; stops=np.r_[starts[1:],len(a)]
        for st,sp in zip(starts,stops): out[order[st:sp]]=0.5*((st+1)+sp)
        return out
    a=ranks(x);b=ranks(y);a-=a.mean();b-=b.mean();den=np.sqrt(np.dot(a,a)*np.dot(b,b))
    return 0.0 if den<=np.finfo(float).eps else float(np.dot(a,b)/den)


def process_panel(
    label: str,
    species: np.ndarray,
    metadata: dict[str,dict[str,str]],
    edge_map: dict[str,np.ndarray],
    source_idx: np.ndarray,
    target_idx: np.ndarray,
    relation: np.ndarray,
    locality_n: dict[str,int],
    *,
    radius: float,
    minimum_coverage: float,
    minimum_sources: int,
) -> tuple[dict[str,np.ndarray],dict]:
    coverage=np.empty(len(relation),float)
    same_class=np.empty(len(relation),bool)
    same_order=np.empty(len(relation),bool)
    same_family=np.empty(len(relation),bool)
    locality_ratio=np.empty(len(relation),float)
    for i,(s,t) in enumerate(zip(source_idx,target_idx)):
        sn=str(species[int(s)]); tn=str(species[int(t)])
        coverage[i]=nearest_coverage(edge_map[tn],edge_map[sn],radius=radius)
        same_class[i]=metadata[sn]["class"]==metadata[tn]["class"]
        same_order[i]=metadata[sn]["order"]==metadata[tn]["order"]
        same_family[i]=metadata[sn]["family"]==metadata[tn]["family"]
        locality_ratio[i]=locality_n[sn]/locality_n[tn]
        if (i+1)%5000==0:
            print(json.dumps({"panel":label,"coverage_dyads_done":i+1,"total":len(relation)}))

    candidate=coverage>=minimum_coverage
    counts=Counter(map(int,target_idx[candidate]))
    supported_targets={t for t,n in counts.items() if n>=minimum_sources}
    eligible=candidate & np.asarray([int(t) in supported_targets for t in target_idx],bool)

    s=source_idx[eligible];t=target_idx[eligible];r=relation[eligible];g=coverage[eligible]
    sc=same_class[eligible].astype(float);so=same_order[eligible].astype(float);sf=same_family[eligible].astype(float)
    lr=locality_ratio[eligible]
    if not len(r):
        raise RuntimeError(f"{label}: no eligible dyads")
    predictors=np.column_stack((
        zscore(r),zscore(g),
        sc-sc.mean(),so-so.mean(),sf-sf.mean(),
        zscore(np.abs(np.log(lr))),
    ))
    prepared=prepare_dyadic_regression(remap(s),remap(t),predictors,primary_index=0)
    quant=np.quantile(r,[0,.1,.25,.5,.75,.9,1])
    source_names={str(species[int(x)]) for x in s}
    target_names={str(species[int(x)]) for x in t}
    summary={
      "panel":label,"candidate_dyads":int(len(relation)),"eligible_dyads":int(len(r)),
      "supported_sources":len(source_names),"supported_targets":len(target_names),
      "R_env_quantiles":{"min":float(quant[0]),"q10":float(quant[1]),"q25":float(quant[2]),"median":float(quant[3]),"q75":float(quant[4]),"q90":float(quant[5]),"max":float(quant[6])},
      "R_env_effective_unique_rounded_1e6":int(len(np.unique(np.round(r,6)))),
      "coverage_quantiles":{k:float(v) for k,v in zip(["min","q10","q25","median","q75","q90","max"],np.quantile(g,[0,.1,.25,.5,.75,.9,1]))},
      "spearman_R_vs_coverage":spearman(r,g),
      "spearman_R_vs_abs_log_locality_ratio":spearman(r,np.abs(np.log(lr))),
      "same_class_fraction":float(sc.mean()),"same_order_fraction":float(so.mean()),"same_family_fraction":float(sf.mean()),
      "predictor_condition_number_after_two_way_FE":float(prepared.condition_number),
      "source_order_counts":dict(Counter(metadata[x]["order"] for x in source_names).most_common()),
      "target_order_counts":dict(Counter(metadata[x]["order"] for x in target_names).most_common()),
    }
    arrays={
      f"{label}_source_index":s,f"{label}_target_index":t,f"{label}_R_env":r,
      f"{label}_coverage":g,f"{label}_same_class":sc.astype(np.int8),f"{label}_same_order":so.astype(np.int8),
      f"{label}_same_family":sf.astype(np.int8),f"{label}_locality_count_ratio":lr,
    }
    return arrays,summary


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--environment-design",type=Path,required=True)
    ap.add_argument("--edges",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--opportunity-rule",type=Path,required=True)
    ap.add_argument("--output-npz",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.opportunity_rule.read_text())
    if rule.get("schema")!="ttf_relational_environment_opportunity_rule_v0.1": raise RuntimeError("bad opportunity rule")
    if any(bool(v) for v in rule["response_firewall"].values()): raise RuntimeError("Study B response firewall open")
    data=np.load(args.environment_design,allow_pickle=False)
    species=np.asarray(data["species_order"]).astype(str)
    metadata=load_metadata(args.candidates)
    edge_map=load_edges(args.edges)
    missing=set(map(str,species))-set(edge_map)
    if missing: raise RuntimeError(f"missing frozen genetic geometry: {len(missing)} species")
    locality_n={name:int(metadata[name]["n_localities"]) for name in metadata}
    geo=rule["directed_geographic_opportunity"]
    radius=float(geo["support_radius_km"]); mincov=float(geo["minimum_target_coverage"]); mins=int(geo["minimum_source_species_per_target"])

    arrays={}
    summaries={}
    for label,prefix in [("development","development"),("confirmatory","confirmatory")]:
        a,s=process_panel(
          label,species,metadata,edge_map,
          np.asarray(data[f"{prefix}_source_index"],np.int64),np.asarray(data[f"{prefix}_target_index"],np.int64),np.asarray(data[f"{prefix}_R_env"],float),
          locality_n,radius=radius,minimum_coverage=mincov,minimum_sources=mins)
        arrays.update(a);summaries[label]=s

    gates=rule["structural_gates_before_synthetic_qualification"];dev=summaries["development"]
    structural=(dev["supported_targets"]>=int(gates["minimum_supported_target_species"]) and dev["supported_sources"]>=int(gates["minimum_supported_source_species"]) and dev["eligible_dyads"]>=int(gates["minimum_supported_directed_dyads"]))
    status="PASS_TO_DEVELOPMENT_SYNTHETIC_QUALIFICATION" if structural else "NOT_EVALUABLE_ENVIRONMENT_OPPORTUNITY_GEOMETRY"
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output_npz,species_order=species,**arrays)
    payload={
      "schema":"ttf_relational_environment_opportunity_design_v0.1","status":status,
      "environment_design_sha256":sha256_path(args.environment_design),"edge_geometry_sha256":sha256_path(args.edges),
      "candidate_csv_sha256":sha256_path(args.candidates),"opportunity_rule_sha256":sha256_path(args.opportunity_rule),
      "panels":summaries,"structural_gate_pass":bool(structural),"design_npz_sha256":sha256_path(args.output_npz),
      "response_firewall":{"Study_B_sequence_identity_opened":False,"Study_B_pairwise_genetic_distances_opened":False,"Study_B_T_st_computed":False,"Study_B_beta_R_computed":False},
      "next_step":"Freeze and run development-only synthetic Type-I/power qualification on this exact eligible dyad geometry. Confirmatory genetic response remains closed."
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":status,"development":dev,"confirmatory":summaries["confirmatory"]},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
