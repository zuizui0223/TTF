#!/usr/bin/env python3
"""Development-only actual-T_st qualification screen for codistributed recurrence.

This runner is response-blind. It accepts exact reconstructed development
locality/edge geometry and frozen taxonomy metadata, generates synthetic genetic
worlds only, and evaluates the v0.3 dyad-specific geometry-null-centering rule.
It never reads aligned nucleotide sequence characters.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from numba import njit

from ttf.codistributed_geometry_null import (
    GeometryNullCenter,
    center_dyad_scores,
    frozen_uint64_seed,
    least_favourable_beta_pvalues,
    prepare_fixed_dyad_transfer_cache,
    simulate_fixed_dyad_tst_batch,
    template_edges_from_geometries,
)
from ttf.genetic_geometry import GeneticSamplingGeometry, endpoint_disjoint_training_counts
from ttf.relational_dyadic import (
    batch_primary_test,
    prepare_dyadic_regression,
    wilson_interval,
)

EXPECTED_DESIGN_SCHEMA="ttf_genetic_codistributed_recurrence_response_blind_design_v0.1"
EXPECTED_RULE_SCHEMA="ttf_genetic_codistributed_recurrence_geometry_null_rule_v0.3"


SOURCE_ARCHIVE_SHA256="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
PANEL_TAG="place-recurrence-insecta-fresh-v0.1"


def digest_labels(names):
    return hashlib.sha256("\n".join(sorted(map(str,names))).encode()).hexdigest()


def ranked_role_digest(names, panel_label):
    tag=f"{PANEL_TAG}|{panel_label}-role"
    ordered=sorted(
        map(str,names),
        key=lambda name:(
            hashlib.sha256(
                f"{tag}|{SOURCE_ARCHIVE_SHA256}|{name}".encode()
            ).hexdigest(),
            name,
        ),
    )
    return hashlib.sha256("\n".join(ordered).encode()).hexdigest()


def read_csv(path: Path):
    with path.open(newline="",encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def zscore(x):
    a=np.asarray(x,dtype=float)
    sd=float(np.std(a,ddof=0))
    if not np.isfinite(sd) or sd <= np.finfo(float).eps:
        raise ValueError("cannot z-score constant/non-finite predictor")
    return (a-float(np.mean(a)))/sd


def coverage_fraction(target,source,radius=500.0,chunk=64):
    t=np.asarray(target,dtype=float); s=np.asarray(source,dtype=float)
    if len(t)==0 or len(s)==0:
        raise ValueError("empty edge midpoint set")
    threshold2=float(radius)**2
    hit=0
    for i in range(0,len(t),int(chunk)):
        x=t[i:i+int(chunk)]
        d=x[:,None,:]-s[None,:,:]
        nearest=np.min(np.sum(d*d,axis=2),axis=1)
        hit += int(np.count_nonzero(nearest <= threshold2))
    return hit/len(t)


@njit(cache=False)
def _directed_coverage_flat(points, offsets, query_index, reference_index, threshold2):
    q0=offsets[query_index]
    q1=offsets[query_index+1]
    r0=offsets[reference_index]
    r1=offsets[reference_index+1]
    covered=0
    for i in range(q0,q1):
        for j in range(r0,r1):
            dx=points[i,0]-points[j,0]
            dy=points[i,1]-points[j,1]
            dz=points[i,2]-points[j,2]
            if dx*dx+dy*dy+dz*dz <= threshold2:
                covered += 1
                break
    return covered / (q1-q0)


@njit(cache=False)
def _symmetric_coverage_grid(
    points, offsets, centers, radii, source_index, target_index, radius
):
    out=np.empty(len(source_index)*len(target_index),dtype=np.float64)
    threshold2=radius*radius
    row=0
    for si in range(len(source_index)):
        s=source_index[si]
        for ti in range(len(target_index)):
            t=target_index[ti]
            dx=centers[s,0]-centers[t,0]
            dy=centers[s,1]-centers[t,1]
            dz=centers[s,2]-centers[t,2]
            cd=(dx*dx+dy*dy+dz*dz)**0.5
            if cd > radii[s]+radii[t]+radius:
                out[row]=0.0
            else:
                target_to_source=_directed_coverage_flat(
                    points,offsets,t,s,threshold2
                )
                source_to_target=_directed_coverage_flat(
                    points,offsets,s,t,threshold2
                )
                out[row]=min(target_to_source,source_to_target)
            row += 1
    return out


def fast_symmetric_coverage(midpoint, sources, targets, radius):
    labels=tuple(sorted(midpoint))
    label_index={name:i for i,name in enumerate(labels)}
    offsets=np.zeros(len(labels)+1,dtype=np.int64)
    arrays=[]
    centers=np.empty((len(labels),3),dtype=float)
    radii=np.empty(len(labels),dtype=float)
    for i,name in enumerate(labels):
        pts=np.asarray(midpoint[name],dtype=float)
        if pts.ndim!=2 or pts.shape[1]!=3 or len(pts)==0:
            raise ValueError(f"invalid midpoint geometry for {name}")
        arrays.append(pts)
        offsets[i+1]=offsets[i]+len(pts)
        center=np.mean(pts,axis=0)
        centers[i]=center
        radii[i]=float(np.max(np.linalg.norm(pts-center,axis=1)))
    points=np.vstack(arrays)
    source_index=np.asarray([label_index[name] for name in sources],dtype=np.int64)
    target_index=np.asarray([label_index[name] for name in targets],dtype=np.int64)
    return _symmetric_coverage_grid(
        points,offsets,centers,radii,source_index,target_index,float(radius)
    )


def reconstruct_geometries(locality_rows,edge_rows,candidate_rows):
    meta={r["species"]:r for r in candidate_rows}
    loc_by={}
    role={}
    for r in locality_rows:
        name=r["species"]
        loc_by.setdefault(name,[]).append(r)
        role.setdefault(name,r["role"])
        if role[name] != r["role"]:
            raise RuntimeError(f"role drift for {name}")
    edge_by={}
    for r in edge_rows:
        name=r["species"]
        edge_by.setdefault(name,[]).append(r)
        if name in role and role[name] != r["role"]:
            raise RuntimeError(f"edge/locality role drift for {name}")
        role.setdefault(name,r["role"])

    if set(loc_by)!=set(edge_by):
        raise RuntimeError("locality/edge species mismatch")
    geometries={}
    midpoint={}
    for name in sorted(loc_by):
        if name not in meta:
            raise RuntimeError(f"missing candidate metadata for {name}")
        loc=sorted(loc_by[name],key=lambda r:int(r["locality_index"]))
        idx=[int(r["locality_index"]) for r in loc]
        if idx != list(range(len(loc))):
            raise RuntimeError(f"locality index drift for {name}")
        xyz=np.asarray([[float(r["x_km"]),float(r["y_km"]),float(r["z_km"])] for r in loc])
        ed=sorted(edge_by[name],key=lambda r:int(r["edge_index"]))
        eidx=[int(r["edge_index"]) for r in ed]
        if eidx != list(range(len(ed))):
            raise RuntimeError(f"edge index drift for {name}")
        nodes=np.asarray([[int(r["node_left"]),int(r["node_right"])] for r in ed],dtype=np.int64)
        mids=np.asarray([[float(r["mid_x_km"]),float(r["mid_y_km"]),float(r["mid_z_km"])] for r in ed])
        direct=0.5*(xyz[nodes[:,0]]+xyz[nodes[:,1]])
        if not np.allclose(mids,direct,rtol=0,atol=1e-9):
            raise RuntimeError(f"edge midpoint drift for {name}")
        counts=endpoint_disjoint_training_counts(nodes)
        g=GeneticSamplingGeometry(
            coordinates=xyz,
            record_to_locality=np.arange(len(xyz),dtype=np.int64),
            records_per_locality=np.ones(len(xyz),dtype=np.int64),
            edge_nodes=nodes,
            endpoint_disjoint_training_edges=counts,
            graph_k=int(meta[name]["graph_k"]),
        )
        expected=(int(meta[name]["n_localities"]),int(meta[name]["edges"]),int(meta[name]["min_endpoint_disjoint_training_edges"]))
        observed=(g.n_localities,g.n_edges,g.min_endpoint_disjoint_training_edges)
        if observed != expected:
            raise RuntimeError(f"candidate geometry invariant drift for {name}: {observed} != {expected}")
        geometries[name]=g
        midpoint[name]=direct
    return geometries,midpoint,role,meta


def build_design(midpoint,geometries,role,meta,radius,expected_condition):
    sources=tuple(sorted(name for name in geometries if role[name]=="source"))
    targets=tuple(sorted(name for name in geometries if role[name]=="target"))
    if set(sources)&set(targets):
        raise RuntimeError("source/target overlap")
    pairs=[(s,t) for s in sources for t in targets]

    center={name:np.mean(midpoint[name],axis=0) for name in geometries}

    G=fast_symmetric_coverage(midpoint,sources,targets,float(radius))
    centroid=np.empty(len(pairs),dtype=float)
    edge_ratio=np.empty(len(pairs),dtype=float)
    locality_ratio=np.empty(len(pairs),dtype=float)
    same_order=np.empty(len(pairs),dtype=float)
    for i,(s,t) in enumerate(pairs):
        cd=float(np.linalg.norm(center[s]-center[t]))
        centroid[i]=np.log1p(cd/float(radius))
        edge_ratio[i]=abs(np.log(geometries[s].n_edges/geometries[t].n_edges))
        locality_ratio[i]=abs(np.log(geometries[s].n_localities/geometries[t].n_localities))
        same_order[i]=float(meta[s]["order"]==meta[t]["order"])

    X=np.column_stack((
        zscore(G),zscore(centroid),zscore(edge_ratio),zscore(locality_ratio),zscore(same_order)
    ))
    src_index=np.repeat(np.arange(len(sources),dtype=np.int64),len(targets))
    tgt_index=np.tile(np.arange(len(targets),dtype=np.int64),len(sources))
    prepared=prepare_dyadic_regression(src_index,tgt_index,X,primary_index=0)
    if not np.isclose(prepared.condition_number,float(expected_condition),rtol=0,atol=1e-7):
        raise RuntimeError(
            f"predictor condition number drift: {prepared.condition_number} != {expected_condition}"
        )
    return sources,targets,pairs,G,prepared


def seeds(master,namespace,cell,n):
    return [
        frozen_uint64_seed(master,namespace,cell,i)
        for i in range(int(n))
    ]


def label_A(a):
    return "A"+str(a).replace(".","p").replace("p0","")


def simulate_chunks(
    geometries,pairs,seed_list,shared_fraction,A,rule,prepared_cache,
    *,batch_size=25,
):
    syn=rule["synthetic_parameters"]
    op=rule["actual_operator"]
    all_seeds=tuple(int(x) for x in seed_list)
    for start in range(0,len(all_seeds),int(batch_size)):
        chunk=all_seeds[start:start+int(batch_size)]
        yield simulate_fixed_dyad_tst_batch(
            geometries,pairs,chunk,
            shared_fraction=float(shared_fraction),
            residual_amplitude=float(A),
            ibd_strength=float(syn["ibd_strength"]),
            noise_sd=float(syn["noise_sd"]),
            transition_width=float(syn["transition_width"]),
            noise_dimensions=int(syn["noise_dimensions"]),
            min_training_edges=5,
            bandwidth_km=float(op["bandwidth_km"]),
            prior_strength=float(op["prior_strength"]),
            prior_mean=float(op["prior_mean"]),
            segment_points=int(op["segment_points"]),
            prepared_cache=prepared_cache,
        )


def beta_for_chunks(prepared,center,chunks):
    pieces=[]
    for tst in chunks:
        beta=batch_primary_test(prepared,center_dyad_scores(tst,center)).coefficient
        pieces.append(np.asarray(beta,dtype=float))
    if not pieces:
        raise RuntimeError("synthetic batch produced no beta_G values")
    return np.concatenate(pieces)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--localities",type=Path,required=True)
    ap.add_argument("--edges",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--stage",choices=("screen","formal"),default="screen")
    ap.add_argument("--screen-result",type=Path)
    ap.add_argument(
        "--world-batch-size",
        type=int,
        default=200,
        help="Execution-only number of synthetic worlds scored per geometry pass; does not alter seeds or estimand.",
    )
    args=ap.parse_args()
    if args.world_batch_size < 1:
        raise RuntimeError("--world-batch-size must be positive")

    design=json.loads(args.design.read_text())
    rule=json.loads(args.rule.read_text())
    if design.get("schema")!=EXPECTED_DESIGN_SCHEMA: raise RuntimeError("unexpected design schema")
    if rule.get("schema")!=EXPECTED_RULE_SCHEMA: raise RuntimeError("unexpected rule schema")
    if rule.get("status")!="FROZEN_BEFORE_FIXED_DYAD_T_ST_SYNTHETIC_RESULTS":
        raise RuntimeError("fixed-T_st screen rule is not frozen")

    if args.stage=="formal":
        if args.screen_result is None:
            raise RuntimeError("--screen-result is required for formal qualification")
        prior_screen=json.loads(args.screen_result.read_text())
        if prior_screen.get("schema")!="ttf_genetic_codistributed_recurrence_fixed_tst_screen_v0.3":
            raise RuntimeError("unexpected screen-result schema")
        if prior_screen.get("status")!="PASS_TO_FORMAL_FIXED_TST_QUALIFICATION":
            raise RuntimeError("formal qualification is not authorized by the screen result")
        if not bool(prior_screen.get("gate",{}).get("overall_pass",False)):
            raise RuntimeError("screen result does not record overall PASS")

    loc=read_csv(args.localities); ed=read_csv(args.edges); cand=read_csv(args.candidates)
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
            json.loads(Path("benchmarks/frozen/genetic_codistributed_recurrence_dyadic_qualification_result_v0.1.json").read_text())
            ["geometry"]["predictor_condition_number"]
        ),
    )
    if len(pairs)!=36290: raise RuntimeError("dyad count drift")

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

    config=(
        rule["development_screen"]
        if args.stage=="screen"
        else rule["formal_qualification"]
    )
    master=int(rule["synthetic_parameters"]["master_seed"])
    amplitudes=[float(a) for a in rule["geometry_null_center"]["private_amplitudes"]]

    cell_means=[]
    center_counts=[]
    for A in amplitudes:
        cell=label_A(A)
        total=np.zeros(len(pairs),dtype=float)
        count=0
        for tst in simulate_chunks(
            geometries,pairs,
            seeds(master,config["namespaces"]["center"],cell,int(config["center_worlds_per_private_amplitude"])),
            0.0,A,rule,fixed_cache,batch_size=args.world_batch_size,
        ):
            total += np.sum(tst,axis=1)
            count += int(tst.shape[1])
        if count != int(config["center_worlds_per_private_amplitude"]):
            raise RuntimeError("center world-count drift")
        cell_means.append(total/count)
        center_counts.append(count)
    center=GeometryNullCenter(
        mu0=np.mean(np.vstack(cell_means),axis=0),
        amplitudes=tuple(amplitudes),
        worlds_per_amplitude=tuple(center_counts),
    )

    private_refs={}
    reference_summary={}
    for A in amplitudes:
        cell=label_A(A)
        beta=beta_for_chunks(
            prepared,center,
            simulate_chunks(
                geometries,pairs,
                seeds(master,config["namespaces"]["reference"],cell,int(config["private_reference_worlds_per_amplitude"])),
                0.0,A,rule,fixed_cache,batch_size=args.world_batch_size,
            ),
        )
        private_refs[cell]=np.asarray(beta,dtype=float)
        reference_summary[cell]={
            "worlds":int(len(beta)),
            "mean_beta_G":float(np.mean(beta)),
            "sd_beta_G":float(np.std(beta,ddof=0)),
            "q95_beta_G":float(np.quantile(beta,0.95)),
        }

    private_eval={}
    all_private_pass=True
    for A in amplitudes:
        cell=label_A(A)
        beta=beta_for_chunks(
            prepared,center,
            simulate_chunks(
                geometries,pairs,
                seeds(master,config["namespaces"]["evaluation"],cell,int(config["private_evaluation_worlds_per_amplitude"])),
                0.0,A,rule,fixed_cache,batch_size=args.world_batch_size,
            ),
        )
        p,least=least_favourable_beta_pvalues(beta,private_refs)
        reject=int(np.count_nonzero(p<=float(rule["inference"]["alpha"])))
        n=len(beta); lo,hi=wilson_interval(reject,n)
        if args.stage=="screen":
            cell_pass=reject<=5
        else:
            cell_pass=hi<=0.10
        all_private_pass &= cell_pass
        private_eval[cell]={
            "worlds":int(n),"rejections":reject,"rejection_rate":reject/n,
            "wilson95_low":lo,"wilson95_high":hi,
            "stage_pass":bool(cell_pass),
            "least_favourable_counts":dict(Counter(least)),
        }

    pos=rule["synthetic_parameters"]["positive_world"]
    pos_beta=beta_for_chunks(
        prepared,center,
        simulate_chunks(
            geometries,pairs,
            seeds(master,config["namespaces"]["positive"],pos["name"],int(config["positive_evaluation_worlds"])),
            float(pos["shared_fraction"]),float(pos["residual_amplitude"]),rule,fixed_cache,batch_size=args.world_batch_size,
        ),
    )
    pos_p,pos_least=least_favourable_beta_pvalues(pos_beta,private_refs)
    pos_reject=int(np.count_nonzero(pos_p<=float(rule["inference"]["alpha"])))
    pos_n=len(pos_beta); pos_lo,pos_hi=wilson_interval(pos_reject,pos_n)
    positive_pass=(
        pos_reject>=35
        if args.stage=="screen"
        else pos_lo>=0.80
    )

    overall=bool(all_private_pass and positive_pass)
    if args.stage=="screen":
        schema="ttf_genetic_codistributed_recurrence_fixed_tst_screen_v0.3"
        status="PASS_TO_FORMAL_FIXED_TST_QUALIFICATION" if overall else "STOP_FIXED_TST_V03_SCREEN"
    else:
        schema="ttf_genetic_codistributed_recurrence_fixed_tst_qualification_v0.3"
        status="PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION" if overall else "NOT_EVALUABLE_CODISTRIBUTED_RECURRENCE_FIXED_TST"
    payload={
        "schema":schema,
        "status":status,
        "stage":args.stage,
        "rule":str(args.rule),
        "execution":{"world_batch_size":int(args.world_batch_size)},
        "development_geometry":{
            "species":len(geometries),"sources":len(sources),"targets":len(targets),"dyads":len(pairs),
            "predictor_condition_number":prepared.condition_number,
            "G_summary":{
                "mean":float(np.mean(G)),"median":float(np.median(G)),
                "q90":float(np.quantile(G,0.9)),"max":float(np.max(G)),
            }
        },
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
            "cell":pos["name"],"worlds":pos_n,"rejections":pos_reject,
            "power":pos_reject/pos_n,"wilson95_low":pos_lo,"wilson95_high":pos_hi,
            "stage_pass":bool(positive_pass),
            "least_favourable_counts":dict(Counter(pos_least)),
            "mean_beta_G":float(np.mean(pos_beta)),
        },
        "gate":{
            "private_pass":bool(all_private_pass),
            "positive_pass":bool(positive_pass),
            "overall_pass":overall,
            "scientific_authority":"development screen only" if args.stage=="screen" else "method qualification only; no empirical genetic conclusion",
            "private_rule":"<=5/50 rejections" if args.stage=="screen" else "Wilson95 upper <=0.10 in every private cell",
            "positive_rule":">=35/50 rejections" if args.stage=="screen" else "Wilson95 lower >=0.80 for globally_shared_place_A2",
            "if_pass":(
                "execute already-predeclared formal qualification with disjoint formal namespaces"
                if args.stage=="screen"
                else "proceed only to confirmatory character-mask preparation and survivor-geometry requalification"
            ),
            "if_fail":"close v0.3 fixed-T_st centering estimator without tuning on these worlds",
        },
        "response_firewall":{
            "development_nucleotide_identity_opened":False,
            "confirmatory_nucleotide_identity_opened":False,
            "empirical_T_st_opened":False,
            "empirical_beta_G_opened":False,
        }
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"stage":args.stage,"status":payload["status"],"private_pass":all_private_pass,"positive_pass":positive_pass},sort_keys=True))


if __name__=="__main__":
    main()
