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
from ttf.relational_environment import nearest_coverage


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


def load_edges(path: Path) -> dict[str,np.ndarray]:
    grouped=defaultdict(list)
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            grouped[row["species"]].append([
                float(row["mid_x_km"]),float(row["mid_y_km"]),float(row["mid_z_km"])
            ])
    return {name:np.asarray(values,float) for name,values in grouped.items()}


def load_compact_edges(path: Path) -> tuple[dict[str,np.ndarray], dict[str,int]]:
    data=np.load(path,allow_pickle=False)
    required={
        "species_order","edge_offsets","edge_midpoints_ecef_km",
        "n_localities","graph_k","n_edges","min_endpoint_disjoint_training_edges",
    }
    missing=required-set(data.files)
    if missing:
        raise RuntimeError(f"compact geometry missing arrays: {sorted(missing)}")
    species=np.asarray(data["species_order"]).astype(str)
    offsets=np.asarray(data["edge_offsets"],dtype=np.int64)
    midpoints=np.asarray(data["edge_midpoints_ecef_km"],dtype=float)
    n_localities=np.asarray(data["n_localities"],dtype=np.int64)
    n_edges=np.asarray(data["n_edges"],dtype=np.int64)
    if len(species)!=1000 or len(set(map(str,species)))!=1000:
        raise RuntimeError("compact geometry must contain exact 1000 unique species")
    if offsets.shape!=(1001,) or offsets[0]!=0 or np.any(np.diff(offsets)<0):
        raise RuntimeError("compact geometry edge offsets are invalid")
    if midpoints.ndim!=2 or midpoints.shape[1]!=3 or int(offsets[-1])!=len(midpoints):
        raise RuntimeError("compact geometry midpoint/offset shape drift")
    if n_localities.shape!=(1000,) or n_edges.shape!=(1000,):
        raise RuntimeError("compact geometry per-species vector shape drift")
    if not np.array_equal(np.diff(offsets),n_edges):
        raise RuntimeError("compact geometry edge-count/offset drift")
    if not np.isfinite(midpoints).all():
        raise RuntimeError("compact geometry contains non-finite edge midpoints")
    edge_map={
        str(name):np.asarray(midpoints[int(offsets[i]):int(offsets[i+1])],dtype=float)
        for i,name in enumerate(species)
    }
    locality_n={str(name):int(n_localities[i]) for i,name in enumerate(species)}
    return edge_map,locality_n


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


def taxonomic_breadth_summary(
    names: set[str],
    metadata: dict[str,dict[str,str]],
    gate: dict[str,object],
) -> dict[str,object]:
    counts=Counter(
        (metadata[name].get("order") or "UNKNOWN").strip() or "UNKNOWN"
        for name in names
    )
    total=len(names)
    threshold=float(gate["order_fraction_threshold"])
    largest_max=float(gate["largest_single_order_fraction_max"])
    minimum_orders=int(gate["minimum_orders_at_or_above_fraction_threshold"])
    if total == 0:
        return {
            "species":0,
            "orders":0,
            "top_orders":{},
            "largest_order_fraction":None,
            "largest_single_order_fraction_max":largest_max,
            "order_fraction_threshold":threshold,
            "orders_at_or_above_fraction_threshold":[],
            "orders_at_or_above_fraction_threshold_count":0,
            "minimum_orders_at_or_above_fraction_threshold":minimum_orders,
            "pass":False,
        }
    fractions={name:count/total for name,count in counts.items()}
    qualifying=sorted(name for name,value in fractions.items() if value>=threshold)
    largest=max(fractions.values())
    return {
        "species":total,
        "orders":len(counts),
        "top_orders":dict(counts.most_common()),
        "largest_order_fraction":largest,
        "largest_single_order_fraction_max":largest_max,
        "order_fraction_threshold":threshold,
        "orders_at_or_above_fraction_threshold":qualifying,
        "orders_at_or_above_fraction_threshold_count":len(qualifying),
        "minimum_orders_at_or_above_fraction_threshold":minimum_orders,
        "pass":bool(largest<=largest_max and len(qualifying)>=minimum_orders),
    }


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
    breadth_gate: dict[str,object],
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
    source_names={str(species[int(x)]) for x in s}
    target_names={str(species[int(x)]) for x in t}
    source_breadth=taxonomic_breadth_summary(source_names,metadata,breadth_gate)
    target_breadth=taxonomic_breadth_summary(target_names,metadata,breadth_gate)

    predictor_preparation_pass=False
    predictor_condition=None
    predictor_error=None
    if len(r):
        try:
            predictors=np.column_stack((
                zscore(r),zscore(g),
                sc-sc.mean(),so-so.mean(),sf-sf.mean(),
                zscore(np.abs(np.log(lr))),
            ))
            prepared=prepare_dyadic_regression(remap(s),remap(t),predictors,primary_index=0)
            predictor_condition=float(prepared.condition_number)
            predictor_preparation_pass=bool(np.isfinite(predictor_condition))
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            predictor_error=f"{type(exc).__name__}: {exc}"

    if len(r):
        quant=np.quantile(r,[0,.1,.25,.5,.75,.9,1])
        coverage_quant=np.quantile(g,[0,.1,.25,.5,.75,.9,1])
        r_quant={"min":float(quant[0]),"q10":float(quant[1]),"q25":float(quant[2]),"median":float(quant[3]),"q75":float(quant[4]),"q90":float(quant[5]),"max":float(quant[6])}
        g_quant={k:float(v) for k,v in zip(["min","q10","q25","median","q75","q90","max"],coverage_quant)}
        r_unique=int(len(np.unique(np.round(r,6))))
        rho_coverage=spearman(r,g)
        rho_locality=spearman(r,np.abs(np.log(lr)))
        same_class_fraction=float(sc.mean());same_order_fraction=float(so.mean());same_family_fraction=float(sf.mean())
    else:
        r_quant=None;g_quant=None;r_unique=0;rho_coverage=None;rho_locality=None
        same_class_fraction=None;same_order_fraction=None;same_family_fraction=None

    summary={
      "panel":label,"candidate_dyads":int(len(relation)),"eligible_dyads":int(len(r)),
      "supported_sources":len(source_names),"supported_targets":len(target_names),
      "R_env_quantiles":r_quant,
      "R_env_effective_unique_rounded_1e6":r_unique,
      "coverage_quantiles":g_quant,
      "spearman_R_vs_coverage":rho_coverage,
      "spearman_R_vs_abs_log_locality_ratio":rho_locality,
      "same_class_fraction":same_class_fraction,"same_order_fraction":same_order_fraction,"same_family_fraction":same_family_fraction,
      "predictor_preparation_pass":predictor_preparation_pass,
      "predictor_preparation_error":predictor_error,
      "predictor_condition_number_after_two_way_FE":predictor_condition,
      "source_order_counts":dict(Counter((metadata[x].get("order") or "UNKNOWN") for x in source_names).most_common()),
      "target_order_counts":dict(Counter((metadata[x].get("order") or "UNKNOWN") for x in target_names).most_common()),
      "source_taxonomic_breadth":source_breadth,
      "target_taxonomic_breadth":target_breadth,
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
    geometry=ap.add_mutually_exclusive_group(required=True)
    geometry.add_argument("--edges",type=Path)
    geometry.add_argument("--compact-geometry",type=Path)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--opportunity-rule",type=Path,required=True)
    ap.add_argument("--output-npz",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.opportunity_rule.read_text())
    if rule.get("schema")!="ttf_relational_environment_opportunity_rule_v0.2": raise RuntimeError("bad opportunity rule")
    if any(bool(v) for v in rule["response_firewall"].values()): raise RuntimeError("Study B response firewall open")
    data=np.load(args.environment_design,allow_pickle=False)
    species=np.asarray(data["species_order"]).astype(str)
    metadata=load_metadata(args.candidates)
    if args.compact_geometry is not None:
        edge_map,compact_locality_n=load_compact_edges(args.compact_geometry)
        candidate_species=set(metadata)
        if set(edge_map)!=candidate_species:
            raise RuntimeError("compact geometry species set does not match frozen candidate 1000")
        for name in candidate_species:
            if int(metadata[name]["n_localities"])!=int(compact_locality_n[name]):
                raise RuntimeError(f"compact geometry locality-count drift: {name}")
        geometry_path=args.compact_geometry
        geometry_source="compact_npz"
    else:
        edge_map=load_edges(args.edges)
        compact_locality_n=None
        geometry_path=args.edges
        geometry_source="canonical_edge_csv"
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
          locality_n,radius=radius,minimum_coverage=mincov,minimum_sources=mins,
          breadth_gate=rule["structural_gates_before_synthetic_qualification"]["taxonomic_breadth"])
        arrays.update(a);summaries[label]=s

    gates=rule["structural_gates_before_synthetic_qualification"];dev=summaries["development"];con=summaries["confirmatory"]
    numeric_support=(
        dev["supported_targets"]>=int(gates["minimum_supported_target_species"])
        and dev["supported_sources"]>=int(gates["minimum_supported_source_species"])
        and dev["eligible_dyads"]>=int(gates["minimum_supported_directed_dyads"])
    )
    taxonomic_support=all([
        dev["source_taxonomic_breadth"]["pass"],
        dev["target_taxonomic_breadth"]["pass"],
        con["source_taxonomic_breadth"]["pass"],
        con["target_taxonomic_breadth"]["pass"],
    ])
    predictor_preparation=bool(dev["predictor_preparation_pass"] and con["predictor_preparation_pass"])
    structural=bool(numeric_support and taxonomic_support and predictor_preparation)
    status="PASS_TO_DEVELOPMENT_SYNTHETIC_QUALIFICATION" if structural else "NOT_EVALUABLE_ENVIRONMENT_OPPORTUNITY_GEOMETRY"
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output_npz,species_order=species,**arrays)
    payload={
      "schema":"ttf_relational_environment_opportunity_design_v0.2","status":status,
      "environment_design_sha256":sha256_path(args.environment_design),
      "edge_geometry_sha256":sha256_path(geometry_path),
      "edge_geometry_source":geometry_source,
      "candidate_csv_sha256":sha256_path(args.candidates),"opportunity_rule_sha256":sha256_path(args.opportunity_rule),
      "panels":summaries,
      "structural_gates":{
        "numeric_support_pass":bool(numeric_support),
        "taxonomic_breadth_pass":bool(taxonomic_support),
        "predictor_preparation_pass":bool(predictor_preparation),
        "overall_pass":bool(structural),
      },
      "structural_gate_pass":bool(structural),"design_npz_sha256":sha256_path(args.output_npz),
      "response_firewall":{"Study_B_sequence_identity_opened":False,"Study_B_pairwise_genetic_distances_opened":False,"Study_B_T_st_computed":False,"Study_B_beta_R_computed":False},
      "next_step":"Freeze and run development-only synthetic Type-I/power qualification on this exact eligible dyad geometry. Confirmatory genetic response remains closed."
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":status,"development":dev,"confirmatory":summaries["confirmatory"]},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
