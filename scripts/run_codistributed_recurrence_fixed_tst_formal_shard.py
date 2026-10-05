#!/usr/bin/env python3
"""Execution-only sharding for the frozen codistributed-recurrence formal gate.

This file changes no scientific quantity. It implements the shard boundaries
frozen before the development-screen result in
genetic_codistributed_recurrence_formal_sharding_v0.1.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from run_codistributed_recurrence_fixed_tst_screen import (
    EXPECTED_DESIGN_SCHEMA,
    EXPECTED_RULE_SCHEMA,
    beta_for_chunks,
    build_design,
    label_A,
    read_csv,
    reconstruct_geometries,
    simulate_chunks,
    ranked_role_digest,
)
from ttf.codistributed_formal_shards import frozen_seed_range, validate_exact_ranges
from ttf.codistributed_geometry_null import (
    GeometryNullCenter,
    least_favourable_beta_pvalues,
    prepare_fixed_dyad_transfer_cache,
    template_edges_from_geometries,
)
from ttf.relational_dyadic import wilson_interval


EXPECTED_SHARD_SCHEMA="ttf_genetic_codistributed_recurrence_formal_sharding_v0.1"
EXPECTED_SCREEN_SCHEMA="ttf_genetic_codistributed_recurrence_fixed_tst_screen_v0.3"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def expected_ranges(sharding: dict, component: str) -> list[tuple[int,int]]:
    if component in {"reference","evaluation"}:
        raw=sharding["shards"][component]["per_amplitude"]
    elif component=="positive":
        raw=sharding["shards"]["positive"]["ranges"]
    elif component=="center":
        raw=sharding["shards"]["center"]["per_amplitude"]
    else:
        raise ValueError(f"unknown component {component}")
    return [(int(a),int(b)) for a,b in raw]


def require_range(sharding: dict, component: str, start: int, stop: int) -> None:
    pair=(int(start),int(stop))
    if pair not in expected_ranges(sharding,component):
        raise RuntimeError(f"{component} range {pair} is not frozen in sharding contract")


def load_screen(path: Path) -> dict:
    p=json.loads(path.read_text())
    if p.get("schema")!=EXPECTED_SCREEN_SCHEMA:
        raise RuntimeError("unexpected screen schema")
    if p.get("status")!="PASS_TO_FORMAL_FIXED_TST_QUALIFICATION":
        raise RuntimeError("formal sharding is not authorized by screen result")
    if not bool(p.get("gate",{}).get("overall_pass",False)):
        raise RuntimeError("screen result lacks overall PASS")
    if any(bool(v) for v in p["response_firewall"].values()):
        raise RuntimeError("screen response firewall is open")
    return p


def load_contracts(design_path: Path, rule_path: Path, sharding_path: Path):
    design=json.loads(design_path.read_text())
    rule=json.loads(rule_path.read_text())
    sharding=json.loads(sharding_path.read_text())
    if design.get("schema")!=EXPECTED_DESIGN_SCHEMA:
        raise RuntimeError("unexpected design schema")
    if rule.get("schema")!=EXPECTED_RULE_SCHEMA:
        raise RuntimeError("unexpected geometry-null rule schema")
    if sharding.get("schema")!=EXPECTED_SHARD_SCHEMA:
        raise RuntimeError("unexpected formal sharding schema")
    if sharding.get("status")!="FROZEN_BEFORE_DEVELOPMENT_SCREEN_RESULT":
        raise RuntimeError("formal sharding contract is not frozen prospectively")
    if int(sharding["invariants"]["world_batch_size"])!=200:
        raise RuntimeError("formal world batch size drift")
    return design,rule,sharding


def load_development_context(args, design: dict, rule: dict):
    loc=read_csv(args.localities)
    ed=read_csv(args.edges)
    cand=read_csv(args.candidates)
    geometries,midpoint,role,meta=reconstruct_geometries(loc,ed,cand)
    sources=sorted(n for n,r in role.items() if r=="source")
    targets=sorted(n for n,r in role.items() if r=="target")
    if len(sources)!=190 or len(targets)!=191 or len(geometries)!=381:
        raise RuntimeError("development role counts drift")
    sel=design["selection"]["development"]
    if ranked_role_digest(sources,"development")!=sel["source_digest_sha256"]:
        raise RuntimeError("source digest drift")
    if ranked_role_digest(targets,"development")!=sel["target_digest_sha256"]:
        raise RuntimeError("target digest drift")
    sources,targets,pairs,G,prepared=build_design(
        midpoint,geometries,role,meta,
        radius=float(design["focal_relation"]["radius_km"]),
        expected_condition=float(
            json.loads(
                Path("benchmarks/frozen/genetic_codistributed_recurrence_dyadic_qualification_result_v0.1.json").read_text()
            )["geometry"]["predictor_condition_number"]
        ),
    )
    if len(pairs)!=36290:
        raise RuntimeError("dyad count drift")
    op=rule["actual_operator"]
    fixed_cache=prepare_fixed_dyad_transfer_cache(
        template_edges_from_geometries(geometries),
        pairs,
        bandwidth_km=float(op["bandwidth_km"]),
        prior_strength=float(op["prior_strength"]),
        prior_mean=float(op["prior_mean"]),
        segment_points=int(op["segment_points"]),
        min_training_edges=5,
    )
    return geometries,pairs,G,prepared,fixed_cache


def load_center(path: Path, rule: dict) -> GeometryNullCenter:
    p=np.load(path,allow_pickle=False)
    mu0=np.asarray(p["mu0"],dtype=float)
    amplitudes=tuple(map(float,np.asarray(p["amplitudes"],dtype=float)))
    counts=tuple(map(int,np.asarray(p["worlds_per_amplitude"],dtype=np.int64)))
    expected=tuple(float(a) for a in rule["geometry_null_center"]["private_amplitudes"])
    n=int(rule["formal_qualification"]["center_worlds_per_private_amplitude"])
    if amplitudes!=expected or counts!=tuple([n]*len(expected)):
        raise RuntimeError("formal center metadata drift")
    if mu0.shape!=(36290,) or not np.isfinite(mu0).all():
        raise RuntimeError("formal center dyad vector drift")
    return GeometryNullCenter(mu0=mu0,amplitudes=amplitudes,worlds_per_amplitude=counts)


def run_center(args,design,rule,sharding):
    require_range(sharding,"center",0,int(rule["formal_qualification"]["center_worlds_per_private_amplitude"]))
    geometries,pairs,G,prepared,fixed_cache=load_development_context(args,design,rule)
    cfg=rule["formal_qualification"]
    master=int(rule["synthetic_parameters"]["master_seed"])
    amplitudes=[float(a) for a in rule["geometry_null_center"]["private_amplitudes"]]
    n=int(cfg["center_worlds_per_private_amplitude"])
    means=[]; counts=[]
    for A in amplitudes:
        cell=label_A(A)
        seed_list=frozen_seed_range(master,cfg["namespaces"]["center"],cell,0,n)
        chunks=list(simulate_chunks(
            geometries,pairs,seed_list,0.0,A,rule,fixed_cache,
            batch_size=int(sharding["invariants"]["world_batch_size"]),
        ))
        if len(chunks)!=1 or chunks[0].shape!=(len(pairs),n):
            raise RuntimeError("center cell no longer matches frozen single-batch operation")
        means.append(np.sum(chunks[0],axis=1)/n)
        counts.append(n)
    center=GeometryNullCenter(
        mu0=np.mean(np.vstack(means),axis=0),
        amplitudes=tuple(amplitudes),
        worlds_per_amplitude=tuple(counts),
    )
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        mu0=center.mu0,
        amplitudes=np.asarray(center.amplitudes,dtype=float),
        worlds_per_amplitude=np.asarray(center.worlds_per_amplitude,dtype=np.int64),
    )
    summary={
        "schema":"ttf_genetic_codistributed_recurrence_formal_center_v0.1",
        "status":"PASS_FROZEN_FORMAL_CENTER",
        "screen_result_sha256":sha256_path(args.screen_result),
        "sharding_contract_sha256":sha256_path(args.sharding),
        "worlds_per_amplitude":n,
        "amplitudes":amplitudes,
        "mu0_summary":{
            "mean":float(np.mean(center.mu0)),
            "sd":float(np.std(center.mu0,ddof=0)),
            "min":float(np.min(center.mu0)),
            "max":float(np.max(center.mu0)),
        },
        "response_firewall":{
            "development_nucleotide_identity_opened":False,
            "confirmatory_nucleotide_identity_opened":False,
            "empirical_T_st_opened":False,
            "empirical_beta_G_opened":False,
        },
    }
    args.output_json.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":summary["status"],"mu0":summary["mu0_summary"]},sort_keys=True))


def run_beta_shard(args,design,rule,sharding,component):
    require_range(sharding,component,args.start,args.stop)
    geometries,pairs,G,prepared,fixed_cache=load_development_context(args,design,rule)
    center=load_center(args.center_npz,rule)
    cfg=rule["formal_qualification"]
    master=int(rule["synthetic_parameters"]["master_seed"])
    payload={
        "component":np.asarray(component),
        "start":np.asarray(int(args.start),dtype=np.int64),
        "stop":np.asarray(int(args.stop),dtype=np.int64),
    }
    if component in {"reference","evaluation"}:
        namespace=cfg["namespaces"][component]
        for A in map(float,rule["geometry_null_center"]["private_amplitudes"]):
            cell=label_A(A)
            beta=beta_for_chunks(
                prepared,center,
                simulate_chunks(
                    geometries,pairs,
                    frozen_seed_range(master,namespace,cell,args.start,args.stop),
                    0.0,A,rule,fixed_cache,
                    batch_size=int(sharding["invariants"]["world_batch_size"]),
                ),
            )
            if len(beta)!=(args.stop-args.start):
                raise RuntimeError("formal beta shard world-count drift")
            payload[cell]=np.asarray(beta,dtype=float)
    elif component=="positive":
        pos=rule["synthetic_parameters"]["positive_world"]
        beta=beta_for_chunks(
            prepared,center,
            simulate_chunks(
                geometries,pairs,
                frozen_seed_range(master,cfg["namespaces"]["positive"],pos["name"],args.start,args.stop),
                float(pos["shared_fraction"]),float(pos["residual_amplitude"]),
                rule,fixed_cache,
                batch_size=int(sharding["invariants"]["world_batch_size"]),
            ),
        )
        if len(beta)!=(args.stop-args.start):
            raise RuntimeError("formal positive shard world-count drift")
        payload["positive"]=np.asarray(beta,dtype=float)
    else:
        raise RuntimeError("unsupported beta shard component")
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output_npz,**payload)
    print(json.dumps({"component":component,"start":args.start,"stop":args.stop},sort_keys=True))


def read_shard_dir(path: Path,component: str,expected: list[tuple[int,int]],keys: tuple[str,...]):
    found={}
    for file in sorted(path.glob("*.npz")):
        p=np.load(file,allow_pickle=False)
        comp=str(np.asarray(p["component"]).item())
        if comp!=component:
            continue
        start=int(np.asarray(p["start"]).item())
        stop=int(np.asarray(p["stop"]).item())
        key=(start,stop)
        if key in found:
            raise RuntimeError(f"duplicate {component} shard {key}")
        data={}
        for k in keys:
            if k not in p.files:
                raise RuntimeError(f"{component} shard {key} missing {k}")
            x=np.asarray(p[k],dtype=float)
            if x.shape!=(stop-start,) or not np.isfinite(x).all():
                raise RuntimeError(f"{component} shard {key} shape drift for {k}")
            data[k]=x
        found[key]=data
    validate_exact_ranges(found.keys(), expected)
    return found


def concat_key(parts,expected,key):
    return np.concatenate([parts[r][key] for r in expected])


def run_aggregate(args,design,rule,sharding):
    screen=load_screen(args.screen_result)
    center=load_center(args.center_npz,rule)
    labels=tuple(label_A(float(a)) for a in rule["geometry_null_center"]["private_amplitudes"])
    ref_ranges=expected_ranges(sharding,"reference")
    eval_ranges=expected_ranges(sharding,"evaluation")
    pos_ranges=expected_ranges(sharding,"positive")
    refs_parts=read_shard_dir(args.reference_dir,"reference",ref_ranges,labels)
    eval_parts=read_shard_dir(args.evaluation_dir,"evaluation",eval_ranges,labels)
    pos_parts=read_shard_dir(args.positive_dir,"positive",pos_ranges,("positive",))

    refs={k:concat_key(refs_parts,ref_ranges,k) for k in labels}
    evaluations={k:concat_key(eval_parts,eval_ranges,k) for k in labels}
    positive=concat_key(pos_parts,pos_ranges,"positive")
    cfg=rule["formal_qualification"]
    if any(len(v)!=int(cfg["private_reference_worlds_per_amplitude"]) for v in refs.values()):
        raise RuntimeError("formal reference count drift")
    if any(len(v)!=int(cfg["private_evaluation_worlds_per_amplitude"]) for v in evaluations.values()):
        raise RuntimeError("formal evaluation count drift")
    if len(positive)!=int(cfg["positive_evaluation_worlds"]):
        raise RuntimeError("formal positive count drift")

    reference_summary={
        k:{
            "worlds":int(len(v)),
            "mean_beta_G":float(np.mean(v)),
            "sd_beta_G":float(np.std(v,ddof=0)),
            "q95_beta_G":float(np.quantile(v,0.95)),
        } for k,v in refs.items()
    }

    private_eval={}
    private_pass=True
    alpha=float(rule["inference"]["alpha"])
    for k in labels:
        p,least=least_favourable_beta_pvalues(evaluations[k],refs)
        reject=int(np.count_nonzero(p<=alpha))
        n=int(len(p)); lo,hi=wilson_interval(reject,n)
        passed=bool(hi<=0.10)
        private_pass &= passed
        private_eval[k]={
            "worlds":n,
            "rejections":reject,
            "rejection_rate":reject/n,
            "wilson95_low":lo,
            "wilson95_high":hi,
            "stage_pass":passed,
            "least_favourable_counts":dict(Counter(least)),
        }

    pos_p,pos_least=least_favourable_beta_pvalues(positive,refs)
    pos_reject=int(np.count_nonzero(pos_p<=alpha))
    pos_n=int(len(pos_p)); pos_lo,pos_hi=wilson_interval(pos_reject,pos_n)
    positive_pass=bool(pos_lo>=0.80)
    overall=bool(private_pass and positive_pass)
    pos=rule["synthetic_parameters"]["positive_world"]

    payload={
        "schema":"ttf_genetic_codistributed_recurrence_fixed_tst_qualification_v0.3",
        "status":"PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION" if overall else "NOT_EVALUABLE_CODISTRIBUTED_RECURRENCE_FIXED_TST",
        "stage":"formal",
        "rule":"docs/supporting/genetic_codistributed_recurrence_geometry_null_rule_v0.3.json",
        "screen_result":"benchmarks/frozen/genetic_codistributed_recurrence_fixed_tst_screen_result_v0.3.json",
        "screen_result_sha256":sha256_path(args.screen_result),
        "formal_sharding_contract":"docs/supporting/genetic_codistributed_recurrence_formal_sharding_v0.1.json",
        "formal_sharding_contract_sha256":sha256_path(args.sharding),
        "execution":{
            "mode":"predeclared_sharded_formal_v0.1",
            "world_batch_size":int(sharding["invariants"]["world_batch_size"]),
            "reference_ranges":[list(x) for x in ref_ranges],
            "evaluation_ranges":[list(x) for x in eval_ranges],
            "positive_ranges":[list(x) for x in pos_ranges],
        },
        "development_geometry":screen["development_geometry"],
        "geometry_null_center":{
            "amplitudes":list(center.amplitudes),
            "worlds_per_amplitude":list(center.worlds_per_amplitude),
            "mu0_summary":{
                "mean":float(np.mean(center.mu0)),
                "sd":float(np.std(center.mu0,ddof=0)),
                "min":float(np.min(center.mu0)),
                "max":float(np.max(center.mu0)),
            },
        },
        "private_references":reference_summary,
        "private_evaluation":private_eval,
        "positive_evaluation":{
            "cell":pos["name"],
            "worlds":pos_n,
            "rejections":pos_reject,
            "power":pos_reject/pos_n,
            "wilson95_low":pos_lo,
            "wilson95_high":pos_hi,
            "stage_pass":positive_pass,
            "least_favourable_counts":dict(Counter(pos_least)),
            "mean_beta_G":float(np.mean(positive)),
        },
        "gate":{
            "private_pass":private_pass,
            "positive_pass":positive_pass,
            "overall_pass":overall,
            "scientific_authority":"method qualification only; no empirical genetic conclusion",
            "private_rule":"Wilson95 upper <=0.10 in every private cell",
            "positive_rule":"Wilson95 lower >=0.80 for globally_shared_place_A2",
            "if_pass":"proceed only to confirmatory character-mask preparation and survivor-geometry requalification",
            "if_fail":"close v0.3 fixed-T_st centering estimator without tuning on these worlds",
        },
        "response_firewall":{
            "development_nucleotide_identity_opened":False,
            "confirmatory_nucleotide_identity_opened":False,
            "empirical_T_st_opened":False,
            "empirical_beta_G_opened":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":payload["status"],"gate":payload["gate"]},sort_keys=True))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--component",choices=("center","reference","evaluation","positive","aggregate"),required=True)
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--sharding",type=Path,required=True)
    ap.add_argument("--screen-result",type=Path,required=True)
    ap.add_argument("--localities",type=Path)
    ap.add_argument("--edges",type=Path)
    ap.add_argument("--candidates",type=Path)
    ap.add_argument("--center-npz",type=Path)
    ap.add_argument("--reference-dir",type=Path)
    ap.add_argument("--evaluation-dir",type=Path)
    ap.add_argument("--positive-dir",type=Path)
    ap.add_argument("--start",type=int)
    ap.add_argument("--stop",type=int)
    ap.add_argument("--output-npz",type=Path)
    ap.add_argument("--output-json",type=Path)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()

    load_screen(args.screen_result)
    design,rule,sharding=load_contracts(args.design,args.rule,args.sharding)

    if args.component=="center":
        if not all((args.localities,args.edges,args.candidates,args.output_npz,args.output_json)):
            raise RuntimeError("center component missing required paths")
        run_center(args,design,rule,sharding)
    elif args.component in {"reference","evaluation","positive"}:
        if None in (args.start,args.stop) or not all((args.localities,args.edges,args.candidates,args.center_npz,args.output_npz)):
            raise RuntimeError("beta shard component missing required paths")
        run_beta_shard(args,design,rule,sharding,args.component)
    else:
        if not all((args.center_npz,args.reference_dir,args.evaluation_dir,args.positive_dir,args.output)):
            raise RuntimeError("aggregate component missing required paths")
        run_aggregate(args,design,rule,sharding)


if __name__=="__main__":
    main()
