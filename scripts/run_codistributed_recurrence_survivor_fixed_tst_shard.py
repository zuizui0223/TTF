#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from run_codistributed_recurrence_fixed_tst_screen import (
    beta_for_chunks,
    fast_symmetric_coverage,
    label_A,
    read_csv,
    reconstruct_geometries,
    simulate_chunks,
    zscore,
)
from ttf.codistributed_formal_shards import (
    aggregate_equal_amplitude_center,
    frozen_seed_range,
    validate_exact_ranges,
)
from ttf.codistributed_geometry_null import (
    GeometryNullCenter,
    least_favourable_beta_pvalues,
    prepare_fixed_dyad_transfer_cache,
    template_edges_from_geometries,
)
from ttf.relational_dyadic import prepare_dyadic_regression, wilson_interval

EXPECTED_MASK="ttf_genetic_codistributed_recurrence_confirmatory_mask_result_v0.1"
EXPECTED_INFO="ttf_genetic_codistributed_recurrence_survivor_information_v0.1"
EXPECTED_RULE="ttf_genetic_codistributed_recurrence_geometry_null_rule_v0.3"
EXPECTED_MASK_RULE="ttf_genetic_codistributed_recurrence_confirmatory_mask_rule_v0.1"
EXPECTED_SHARD="ttf_genetic_codistributed_recurrence_formal_sharding_v0.1"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def digest_labels(names) -> str:
    return hashlib.sha256("\n".join(sorted(map(str,names))).encode()).hexdigest()


def expected_ranges(sharding: dict, component: str) -> list[tuple[int,int]]:
    raw=(
        sharding["shards"][component]["per_amplitude"]
        if component in {"center","reference","evaluation"}
        else sharding["shards"]["positive"]["ranges"]
    )
    return [(int(a),int(b)) for a,b in raw]


def require_range(sharding: dict, component: str, start: int, stop: int) -> None:
    if (int(start),int(stop)) not in expected_ranges(sharding,component):
        raise RuntimeError(f"{component} range is not frozen")


def load_contracts(args):
    mask=json.loads(args.mask_result.read_text())
    info=json.loads(args.survivor_info.read_text())
    rule=json.loads(args.rule.read_text())
    mask_rule=json.loads(args.mask_rule.read_text())
    sharding=json.loads(args.sharding.read_text())
    if mask.get("schema")!=EXPECTED_MASK or mask.get("status")!="PASS_TO_SURVIVOR_INFORMATION_GATE":
        raise RuntimeError("confirmatory mask result not authorized")
    if (
        info.get("schema")!=EXPECTED_INFO
        or info.get("status")!="PASS_TO_SURVIVOR_FIXED_TST_REQUALIFICATION"
        or info.get("information_gate",{}).get("overall_pass") is not True
    ):
        raise RuntimeError("survivor information gate not authorized")
    if rule.get("schema")!=EXPECTED_RULE or mask_rule.get("schema")!=EXPECTED_MASK_RULE:
        raise RuntimeError("survivor rule schema drift")
    if sharding.get("schema")!=EXPECTED_SHARD:
        raise RuntimeError("sharding schema drift")
    if any(bool(v) for v in mask["response_firewall"].values()):
        raise RuntimeError("mask response firewall open")
    if any(bool(v) for v in info["response_firewall"].values()):
        raise RuntimeError("survivor information response firewall open")
    q=mask_rule["exact_survivor_fixed_tst_requalification"]
    cfg=rule["formal_qualification"]
    key_pairs=(
        ("center_worlds_per_amplitude","center_worlds_per_private_amplitude"),
        ("private_reference_worlds_per_amplitude","private_reference_worlds_per_amplitude"),
        ("private_evaluation_worlds_per_amplitude","private_evaluation_worlds_per_amplitude"),
    )
    for survivor_key,formal_key in key_pairs:
        if int(q[survivor_key])!=int(cfg[formal_key]):
            raise RuntimeError(f"survivor world-count drift: {survivor_key}")
    if int(q["positive_worlds"])!=int(cfg["positive_evaluation_worlds"]):
        raise RuntimeError("survivor positive world-count drift")
    return mask,info,rule,mask_rule,sharding


def load_context(args,mask,info,rule):
    candidate_rows=read_csv(args.candidates)
    geometries,midpoint,role,meta=reconstruct_geometries(
        read_csv(args.localities),read_csv(args.edges),candidate_rows
    )
    sources=tuple(sorted(n for n,r in role.items() if r=="source"))
    targets=tuple(sorted(n for n,r in role.items() if r=="target"))
    if (len(geometries),len(sources),len(targets))!=(326,166,160):
        raise RuntimeError("survivor role counts drift")
    if digest_labels(sources)!=mask["surviving_source_digest_sha256"]:
        raise RuntimeError("survivor source digest drift")
    if digest_labels(targets)!=mask["surviving_target_digest_sha256"]:
        raise RuntimeError("survivor target digest drift")

    pairs=[(s,t) for s in sources for t in targets]
    G=fast_symmetric_coverage(midpoint,sources,targets,500.0)
    center={n:np.mean(midpoint[n],axis=0) for n in geometries}
    centroid=np.empty(len(pairs),dtype=float)
    edge_ratio=np.empty(len(pairs),dtype=float)
    locality_ratio=np.empty(len(pairs),dtype=float)
    same_order=np.empty(len(pairs),dtype=float)
    for i,(s,t) in enumerate(pairs):
        centroid[i]=np.log1p(float(np.linalg.norm(center[s]-center[t]))/500.0)
        edge_ratio[i]=abs(np.log(geometries[s].n_edges/geometries[t].n_edges))
        locality_ratio[i]=abs(np.log(geometries[s].n_localities/geometries[t].n_localities))
        same_order[i]=float(meta[s]["order"]==meta[t]["order"])
    X=np.column_stack((
        zscore(G),zscore(centroid),zscore(edge_ratio),zscore(locality_ratio),zscore(same_order)
    ))
    src=np.repeat(np.arange(len(sources),dtype=np.int64),len(targets))
    tgt=np.tile(np.arange(len(targets),dtype=np.int64),len(sources))
    prepared=prepare_dyadic_regression(src,tgt,X,primary_index=0)
    gs=info["G_summary"]
    observed=(float(np.mean(G)),float(np.median(G)),float(np.quantile(G,0.9)),float(np.max(G)))
    expected=(float(gs["mean"]),float(gs["median"]),float(gs["q90"]),float(gs["max"]))
    if not np.allclose(observed,expected,rtol=0,atol=1e-14):
        raise RuntimeError("survivor G summary drift")
    op=rule["actual_operator"]
    cache=prepare_fixed_dyad_transfer_cache(
        template_edges_from_geometries(geometries),pairs,
        bandwidth_km=float(op["bandwidth_km"]),
        prior_strength=float(op["prior_strength"]),
        prior_mean=float(op["prior_mean"]),
        segment_points=int(op["segment_points"]),
        min_training_edges=5,
    )
    return geometries,pairs,prepared,cache


def load_center(path: Path,rule: dict,mask_rule: dict,n_dyads: int) -> GeometryNullCenter:
    p=np.load(path,allow_pickle=False)
    mu0=np.asarray(p["mu0"],dtype=float)
    amps=tuple(map(float,np.asarray(p["amplitudes"],dtype=float)))
    counts=tuple(map(int,np.asarray(p["worlds_per_amplitude"],dtype=np.int64)))
    expected=tuple(map(float,rule["geometry_null_center"]["private_amplitudes"]))
    n=int(mask_rule["exact_survivor_fixed_tst_requalification"]["center_worlds_per_amplitude"])
    if amps!=expected or counts!=tuple([n]*len(expected)) or mu0.shape!=(n_dyads,):
        raise RuntimeError("survivor center drift")
    return GeometryNullCenter(mu0=mu0,amplitudes=amps,worlds_per_amplitude=counts)


def run_center(args,mask,info,rule,mask_rule,sharding):
    q=mask_rule["exact_survivor_fixed_tst_requalification"]
    n=int(q["center_worlds_per_amplitude"])
    require_range(sharding,"center",0,n)
    geometries,pairs,prepared,cache=load_context(args,mask,info,rule)
    master=int(rule["synthetic_parameters"]["master_seed"])
    namespaces=q["survivor_seed_namespaces"]
    amplitudes=[float(a) for a in rule["geometry_null_center"]["private_amplitudes"]]
    means=[]
    for A in amplitudes:
        cell=label_A(A)
        seeds=frozen_seed_range(master,namespaces["center"],cell,0,n)
        chunks=list(simulate_chunks(
            geometries,pairs,seeds,0.0,A,rule,cache,batch_size=200
        ))
        if len(chunks)!=1 or chunks[0].shape!=(len(pairs),n):
            raise RuntimeError("survivor center batch drift")
        means.append(np.sum(chunks[0],axis=1)/n)
    center=GeometryNullCenter(
        mu0=np.mean(np.vstack(means),axis=0),
        amplitudes=tuple(amplitudes),
        worlds_per_amplitude=tuple([n]*len(amplitudes)),
    )
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_npz,mu0=center.mu0,
        amplitudes=np.asarray(center.amplitudes,dtype=float),
        worlds_per_amplitude=np.asarray(center.worlds_per_amplitude,dtype=np.int64),
    )
    args.output_json.write_text(json.dumps({
        "schema":"ttf_genetic_codistributed_recurrence_survivor_formal_center_v0.1",
        "status":"PASS_SURVIVOR_CENTER",
        "dyads":len(pairs),
        "predictor_condition_number":prepared.condition_number,
        "mu0_summary":{
            "mean":float(np.mean(center.mu0)),
            "sd":float(np.std(center.mu0,ddof=0)),
            "min":float(np.min(center.mu0)),
            "max":float(np.max(center.mu0)),
        },
        "response_firewall":info["response_firewall"],
    },indent=2,sort_keys=True)+"\n")


def run_center_cell(args,mask,info,rule,mask_rule,sharding):
    q=mask_rule["exact_survivor_fixed_tst_requalification"]
    n=int(q["center_worlds_per_amplitude"])
    require_range(sharding,"center",0,n)
    amplitudes=tuple(float(a) for a in rule["geometry_null_center"]["private_amplitudes"])
    A=float(args.amplitude)
    if A not in amplitudes:
        raise RuntimeError("center amplitude is not frozen")
    geometries,pairs,prepared,cache=load_context(args,mask,info,rule)
    cell=label_A(A)
    seeds=frozen_seed_range(
        int(rule["synthetic_parameters"]["master_seed"]),
        q["survivor_seed_namespaces"]["center"],
        cell,0,n,
    )
    chunks=list(simulate_chunks(
        geometries,pairs,seeds,0.0,A,rule,cache,batch_size=200
    ))
    if len(chunks)!=1 or chunks[0].shape!=(len(pairs),n):
        raise RuntimeError("survivor center cell no longer matches frozen single-batch operation")
    mean=np.sum(chunks[0],axis=1)/n
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        component=np.asarray("center-cell"),
        amplitude=np.asarray(A,dtype=float),
        cell=np.asarray(cell),
        worlds=np.asarray(n,dtype=np.int64),
        mean=np.asarray(mean,dtype=float),
        predictor_condition_number=np.asarray(prepared.condition_number,dtype=float),
    )
    print(json.dumps({
        "component":"center-cell","amplitude":A,"cell":cell,
        "worlds":n,"dyads":len(pairs),
        "mean_summary":{
            "mean":float(np.mean(mean)),
            "sd":float(np.std(mean,ddof=0)),
            "min":float(np.min(mean)),
            "max":float(np.max(mean)),
        },
    },sort_keys=True))


def run_center_aggregate(args,mask,info,rule,mask_rule,sharding):
    q=mask_rule["exact_survivor_fixed_tst_requalification"]
    n=int(q["center_worlds_per_amplitude"])
    amplitudes=tuple(float(a) for a in rule["geometry_null_center"]["private_amplitudes"])
    expected_cells={label_A(a):a for a in amplitudes}
    found={}
    conditions=[]
    for file in sorted(args.center_dir.glob("*.npz")):
        p=np.load(file,allow_pickle=False)
        if str(np.asarray(p["component"]).item())!="center-cell":
            continue
        A=float(np.asarray(p["amplitude"]).item())
        cell=str(np.asarray(p["cell"]).item())
        worlds=int(np.asarray(p["worlds"]).item())
        if cell not in expected_cells or A!=expected_cells[cell] or worlds!=n:
            raise RuntimeError(f"center cell metadata drift in {file}")
        if A in found:
            raise RuntimeError(f"duplicate center amplitude {A:g}")
        mean=np.asarray(p["mean"],dtype=float)
        if mean.shape!=(26560,) or not np.isfinite(mean).all():
            raise RuntimeError(f"center cell dyad payload drift for A={A:g}")
        found[A]=mean
        conditions.append(float(np.asarray(p["predictor_condition_number"]).item()))
    if set(found)!=set(amplitudes):
        raise RuntimeError("center amplitude shard set incomplete")
    if not np.allclose(conditions,conditions[0],rtol=0,atol=1e-12):
        raise RuntimeError("center cell predictor condition-number drift")
    center=aggregate_equal_amplitude_center(
        found,amplitudes,worlds_per_amplitude=n
    )
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_npz,mu0=center.mu0,
        amplitudes=np.asarray(center.amplitudes,dtype=float),
        worlds_per_amplitude=np.asarray(center.worlds_per_amplitude,dtype=np.int64),
    )
    args.output_json.write_text(json.dumps({
        "schema":"ttf_genetic_codistributed_recurrence_survivor_formal_center_v0.1",
        "status":"PASS_SURVIVOR_CENTER",
        "execution_mode":"predeclared_parallel_amplitude_cells",
        "dyads":26560,
        "predictor_condition_number":conditions[0],
        "mu0_summary":{
            "mean":float(np.mean(center.mu0)),
            "sd":float(np.std(center.mu0,ddof=0)),
            "min":float(np.min(center.mu0)),
            "max":float(np.max(center.mu0)),
        },
        "response_firewall":info["response_firewall"],
    },indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":"PASS_SURVIVOR_CENTER",
        "execution_mode":"predeclared_parallel_amplitude_cells",
        "mu0_summary":{
            "mean":float(np.mean(center.mu0)),
            "sd":float(np.std(center.mu0,ddof=0)),
            "min":float(np.min(center.mu0)),
            "max":float(np.max(center.mu0)),
        },
    },sort_keys=True))


def run_beta(args,mask,info,rule,mask_rule,sharding,component):
    require_range(sharding,component,args.start,args.stop)
    geometries,pairs,prepared,cache=load_context(args,mask,info,rule)
    center=load_center(args.center_npz,rule,mask_rule,len(pairs))
    q=mask_rule["exact_survivor_fixed_tst_requalification"]
    namespaces=q["survivor_seed_namespaces"]
    master=int(rule["synthetic_parameters"]["master_seed"])
    payload={
        "component":np.asarray(component),
        "start":np.asarray(args.start,dtype=np.int64),
        "stop":np.asarray(args.stop,dtype=np.int64),
    }
    if component in {"reference","evaluation"}:
        for A in map(float,rule["geometry_null_center"]["private_amplitudes"]):
            cell=label_A(A)
            payload[cell]=np.asarray(beta_for_chunks(
                prepared,center,
                simulate_chunks(
                    geometries,pairs,
                    frozen_seed_range(master,namespaces[component],cell,args.start,args.stop),
                    0.0,A,rule,cache,batch_size=200,
                ),
            ),dtype=float)
    else:
        pos=rule["synthetic_parameters"]["positive_world"]
        payload["positive"]=np.asarray(beta_for_chunks(
            prepared,center,
            simulate_chunks(
                geometries,pairs,
                frozen_seed_range(master,namespaces["positive"],pos["name"],args.start,args.stop),
                float(pos["shared_fraction"]),float(pos["residual_amplitude"]),
                rule,cache,batch_size=200,
            ),
        ),dtype=float)
    args.output_npz.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output_npz,**payload)


def read_parts(path,component,ranges,keys):
    found={}
    for file in sorted(path.glob("*.npz")):
        p=np.load(file,allow_pickle=False)
        if str(np.asarray(p["component"]).item())!=component:
            continue
        key=(int(np.asarray(p["start"]).item()),int(np.asarray(p["stop"]).item()))
        if key in found:
            raise RuntimeError(f"duplicate {component} shard {key}")
        data={}
        for k in keys:
            x=np.asarray(p[k],dtype=float)
            if x.shape!=(key[1]-key[0],) or not np.isfinite(x).all():
                raise RuntimeError(f"{component} shard {key} payload drift for {k}")
            data[k]=x
        found[key]=data
    validate_exact_ranges(found.keys(),ranges)
    return found


def run_aggregate(args,mask,info,rule,mask_rule,sharding):
    labels=tuple(label_A(float(a)) for a in rule["geometry_null_center"]["private_amplitudes"])
    rr=expected_ranges(sharding,"reference")
    er=expected_ranges(sharding,"evaluation")
    pr=expected_ranges(sharding,"positive")
    rp=read_parts(args.reference_dir,"reference",rr,labels)
    ep=read_parts(args.evaluation_dir,"evaluation",er,labels)
    pp=read_parts(args.positive_dir,"positive",pr,("positive",))
    refs={k:np.concatenate([rp[r][k] for r in rr]) for k in labels}
    evaluations={k:np.concatenate([ep[r][k] for r in er]) for k in labels}
    positive=np.concatenate([pp[r]["positive"] for r in pr])
    q=mask_rule["exact_survivor_fixed_tst_requalification"]
    if any(len(v)!=int(q["private_reference_worlds_per_amplitude"]) for v in refs.values()):
        raise RuntimeError("survivor reference world-count drift")
    if any(len(v)!=int(q["private_evaluation_worlds_per_amplitude"]) for v in evaluations.values()):
        raise RuntimeError("survivor evaluation world-count drift")
    if len(positive)!=int(q["positive_worlds"]):
        raise RuntimeError("survivor positive world-count drift")

    alpha=float(rule["inference"]["alpha"])
    private={}
    private_pass=True
    for cell in labels:
        p,least=least_favourable_beta_pvalues(evaluations[cell],refs)
        rejected=int(np.count_nonzero(p<=alpha))
        lo,hi=wilson_interval(rejected,len(p))
        passed=bool(hi<=0.10)
        private_pass &= passed
        private[cell]={
            "worlds":int(len(p)),
            "rejections":rejected,
            "rejection_rate":rejected/len(p),
            "wilson95_low":lo,
            "wilson95_high":hi,
            "stage_pass":passed,
            "least_favourable_counts":dict(Counter(least)),
        }

    p,least=least_favourable_beta_pvalues(positive,refs)
    rejected=int(np.count_nonzero(p<=alpha))
    lo,hi=wilson_interval(rejected,len(p))
    positive_pass=bool(lo>=0.80)
    overall=bool(private_pass and positive_pass)
    payload={
        "schema":"ttf_genetic_codistributed_recurrence_survivor_fixed_tst_qualification_v0.1",
        "status":(
            "PASS_TO_ONE_SHOT_CONFIRMATORY_EMPIRICAL_OPENING"
            if overall else
            "NOT_EVALUABLE_CODISTRIBUTED_RECURRENCE_SURVIVOR_GEOMETRY"
        ),
        "private_evaluation":private,
        "positive_evaluation":{
            "worlds":int(len(p)),
            "rejections":rejected,
            "power":rejected/len(p),
            "wilson95_low":lo,
            "wilson95_high":hi,
            "stage_pass":positive_pass,
            "least_favourable_counts":dict(Counter(least)),
            "mean_beta_G":float(np.mean(positive)),
        },
        "gate":{
            "private_pass":bool(private_pass),
            "positive_pass":bool(positive_pass),
            "overall_pass":overall,
            "private_rule":"Wilson95 upper <=0.10 in every private cell",
            "positive_rule":"Wilson95 lower >=0.80 for globally_shared_place_A2",
        },
        "response_firewall":info["response_firewall"],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":payload["status"],"gate":payload["gate"]},sort_keys=True))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--component",choices=("center","center-cell","center-aggregate","reference","evaluation","positive","aggregate"),required=True)
    for name in (
        "localities","edges","candidates","mask-result","survivor-info",
        "rule","mask-rule","sharding","center-npz","center-dir","reference-dir",
        "evaluation-dir","positive-dir","output-npz","output-json","output",
    ):
        ap.add_argument("--"+name,type=Path)
    ap.add_argument("--start",type=int)
    ap.add_argument("--stop",type=int)
    ap.add_argument("--amplitude",type=float)
    args=ap.parse_args()
    mask,info,rule,mask_rule,sharding=load_contracts(args)
    if args.component=="center":
        run_center(args,mask,info,rule,mask_rule,sharding)
    elif args.component=="center-cell":
        if args.amplitude is None or not all((args.localities,args.edges,args.candidates,args.output_npz)):
            raise RuntimeError("center-cell component missing required inputs")
        run_center_cell(args,mask,info,rule,mask_rule,sharding)
    elif args.component=="center-aggregate":
        if not all((args.center_dir,args.output_npz,args.output_json)):
            raise RuntimeError("center-aggregate component missing required inputs")
        run_center_aggregate(args,mask,info,rule,mask_rule,sharding)
    elif args.component in {"reference","evaluation","positive"}:
        run_beta(args,mask,info,rule,mask_rule,sharding,args.component)
    else:
        run_aggregate(args,mask,info,rule,mask_rule,sharding)


if __name__=="__main__":
    main()
