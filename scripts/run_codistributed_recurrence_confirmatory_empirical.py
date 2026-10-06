#!/usr/bin/env python3
"""One-shot empirical opening for fresh codistributed-recurrence Insecta.

This runner is inert until a separate authorization file is created after the
survivor fixed-Tst synthetic requalification passes. It serializes no nucleotide
identity, edge genetic-distance vectors, or dyad-level empirical T_st values.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from run_codistributed_recurrence_fixed_tst_screen import (
    fast_symmetric_coverage,
    label_A,
    read_csv,
    reconstruct_geometries,
    zscore,
)
from ttf.codistributed_formal_shards import validate_exact_ranges
from ttf.codistributed_geometry_null import (
    GeometryNullCenter,
    center_dyad_scores,
    least_favourable_beta_pvalues,
    prepare_fixed_dyad_transfer_cache,
    score_fixed_dyad_transfer_cache,
    template_edges_from_geometries,
)
from ttf.phylogatr_confirmatory import read_occurrence_rows
from ttf.phylogatr_empirical import (
    canonical_mask_sha256_from_identity,
    frozen_edge_mean_p_distances,
    read_aligned_nucleotide_identity,
    sequences_by_frozen_locality,
)
from ttf.relational_dyadic import batch_primary_test, prepare_dyadic_regression
from ttf.relational_genetic_empirical import post_ibd_responses

EXPECTED_AUTH="ttf_genetic_codistributed_recurrence_identity_opening_authorization_v0.1"
EXPECTED_OPENING_RULE="ttf_genetic_codistributed_recurrence_empirical_opening_rule_v0.1"
EXPECTED_QUAL="ttf_genetic_codistributed_recurrence_survivor_fixed_tst_qualification_v0.1"
EXPECTED_MASK="ttf_genetic_codistributed_recurrence_confirmatory_mask_result_v0.1"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def digest_labels(names) -> str:
    return hashlib.sha256("\n".join(sorted(map(str,names))).encode()).hexdigest()


def npz_hashes(path: Path) -> dict[str,str]:
    files=sorted(path.glob("*.npz"))
    if not files:
        raise RuntimeError("no NPZ artifacts found")
    return {file.name:sha256_path(file) for file in files}


def read_metadata(path: Path) -> dict[str,dict[str,str]]:
    rows=list(csv.DictReader(path.open(newline="",encoding="utf-8")))
    out={str(r["species"]):r for r in rows}
    if len(out)!=len(rows):
        raise RuntimeError("duplicate candidate species")
    return out


def read_full_mask(path: Path,expected_sha: str) -> dict:
    if sha256_path(path)!=expected_sha:
        raise RuntimeError("full confirmatory mask ledger SHA drift")
    p=json.loads(path.read_text())
    if p.get("schema")!=EXPECTED_MASK or p.get("status")!="PASS_TO_SURVIVOR_INFORMATION_GATE":
        raise RuntimeError("unexpected full mask ledger")
    if any(bool(v) for v in p["response_firewall"].values()):
        raise RuntimeError("full mask ledger firewall open")
    return p


def expected_ranges(sharding: dict) -> list[tuple[int,int]]:
    return [(int(a),int(b)) for a,b in sharding["shards"]["reference"]["per_amplitude"]]


def load_private_references(path: Path,sharding: dict,labels: tuple[str,...]) -> dict[str,np.ndarray]:
    ranges=expected_ranges(sharding)
    found={}
    for file in sorted(path.glob("*.npz")):
        p=np.load(file,allow_pickle=False)
        if str(np.asarray(p["component"]).item())!="reference":
            continue
        key=(int(np.asarray(p["start"]).item()),int(np.asarray(p["stop"]).item()))
        if key in found:
            raise RuntimeError(f"duplicate reference shard {key}")
        data={}
        for cell in labels:
            x=np.asarray(p[cell],dtype=float)
            if x.shape!=(key[1]-key[0],) or not np.isfinite(x).all():
                raise RuntimeError(f"reference shard payload drift {key} {cell}")
            data[cell]=x
        found[key]=data
    validate_exact_ranges(found.keys(),ranges)
    refs={cell:np.concatenate([found[r][cell] for r in ranges]) for cell in labels}
    if any(len(x)!=999 for x in refs.values()):
        raise RuntimeError("private reference world-count drift")
    return refs


def component_pvalues(beta: float,refs: dict[str,np.ndarray]) -> dict[str,float]:
    out={}
    for cell,x in refs.items():
        out[cell]=float((1+np.count_nonzero(np.asarray(x)>=float(beta)))/(len(x)+1))
    return out


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,required=True)
    ap.add_argument("--source-archive-sha256",required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--localities",type=Path,required=True)
    ap.add_argument("--edges",type=Path,required=True)
    ap.add_argument("--full-mask-ledger",type=Path,required=True)
    ap.add_argument("--survivor-qualification",type=Path,required=True)
    ap.add_argument("--center-npz",type=Path,required=True)
    ap.add_argument("--reference-dir",type=Path,required=True)
    ap.add_argument("--geometry-rule",type=Path,required=True)
    ap.add_argument("--sharding",type=Path,required=True)
    ap.add_argument("--opening-rule",type=Path,required=True)
    ap.add_argument("--authorization",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    # No nucleotide identity is read before every authorization check below passes.
    auth=json.loads(args.authorization.read_text())
    if auth.get("schema")!=EXPECTED_AUTH:
        raise RuntimeError("unexpected identity-opening authorization schema")
    if auth.get("status")!="AUTHORIZE_ONE_SHOT_CODISTRIBUTED_NUCLEOTIDE_IDENTITY_OPENING":
        raise RuntimeError("one-shot codistributed identity opening is not authorized")

    opening=json.loads(args.opening_rule.read_text())
    if opening.get("schema")!=EXPECTED_OPENING_RULE:
        raise RuntimeError("unexpected empirical opening rule")
    if opening.get("status")!="FROZEN_BEFORE_SURVIVOR_FIXED_TST_REQUALIFICATION_RESULT":
        raise RuntimeError("empirical opening rule was not prospectively frozen")
    if args.source_archive_sha256!=opening["source_binding"]["archive_sha256"]:
        raise RuntimeError("source archive SHA drift")

    qual=json.loads(args.survivor_qualification.read_text())
    if qual.get("schema")!=EXPECTED_QUAL:
        raise RuntimeError("unexpected survivor fixed-Tst qualification schema")
    if qual.get("status")!="PASS_TO_ONE_SHOT_CONFIRMATORY_EMPIRICAL_OPENING":
        raise RuntimeError("survivor synthetic requalification did not authorize identity opening")
    if qual.get("gate",{}).get("overall_pass") is not True:
        raise RuntimeError("survivor synthetic gate lacks overall PASS")
    if any(bool(v) for v in qual["response_firewall"].values()):
        raise RuntimeError("survivor qualification response firewall open")
    binding=qual.get("artifact_binding")
    if not isinstance(binding,dict):
        raise RuntimeError("survivor qualification lacks synthetic artifact binding")
    current_reference_hashes=npz_hashes(args.reference_dir)
    if sha256_path(args.center_npz)!=binding.get("center_npz_sha256"):
        raise RuntimeError("survivor center differs from qualification artifact")
    if current_reference_hashes!=binding.get("reference_shards_sha256"):
        raise RuntimeError("survivor private references differ from qualification artifacts")

    actual_inputs={
        "candidates":sha256_path(args.candidates),
        "localities":sha256_path(args.localities),
        "edges":sha256_path(args.edges),
        "full_mask_ledger":sha256_path(args.full_mask_ledger),
        "survivor_qualification":sha256_path(args.survivor_qualification),
        "center_npz":sha256_path(args.center_npz),
        "geometry_rule":sha256_path(args.geometry_rule),
        "sharding":sha256_path(args.sharding),
        "opening_rule":sha256_path(args.opening_rule),
    }
    if actual_inputs!=auth["inputs_sha256"]:
        raise RuntimeError("authorized empirical input hash drift")
    if current_reference_hashes!=auth.get("reference_shards_sha256"):
        raise RuntimeError("authorized private-reference artifact hash drift")
    for path,expected in auth.get("code_sha256",{}).items():
        p=Path(path)
        if not p.is_file() or sha256_path(p)!=expected:
            raise RuntimeError(f"authorized code drift: {path}")

    root=args.root.resolve()
    if sha256_path(root/"genes.txt")!=auth["source"]["genes_txt_sha256"]:
        raise RuntimeError("genes.txt drift")
    if sha256_path(root/"cite.txt")!=auth["source"]["cite_txt_sha256"]:
        raise RuntimeError("cite.txt drift")

    full_mask=read_full_mask(
        args.full_mask_ledger,
        opening["frozen_empirical_domain"]["full_mask_ledger_sha256"],
    )
    mask_by_species={str(row["species"]):row for row in full_mask["species"]}

    candidate_rows=read_csv(args.candidates)
    metadata={str(row["species"]):row for row in candidate_rows}
    locality_rows=read_csv(args.localities)
    edge_rows=read_csv(args.edges)
    geometries,midpoint,role,meta=reconstruct_geometries(
        locality_rows,edge_rows,candidate_rows
    )
    sources=tuple(sorted(n for n,r in role.items() if r=="source"))
    targets=tuple(sorted(n for n,r in role.items() if r=="target"))
    frozen=opening["frozen_empirical_domain"]
    if (len(geometries),len(sources),len(targets),len(sources)*len(targets))!=(
        frozen["species"],frozen["sources"],frozen["targets"],frozen["directed_dyads"]
    ):
        raise RuntimeError("frozen empirical domain count drift")
    if digest_labels(sources)!=frozen["source_digest_sha256"]:
        raise RuntimeError("empirical source digest drift")
    if digest_labels(targets)!=frozen["target_digest_sha256"]:
        raise RuntimeError("empirical target digest drift")

    latlon_by={}
    for row in locality_rows:
        latlon_by.setdefault(row["species"],[]).append(row)
    canonical_latlon={}
    for name,rows in latlon_by.items():
        rr=sorted(rows,key=lambda r:int(r["locality_index"]))
        canonical_latlon[name]=np.asarray([
            [float(r["latitude"]),float(r["longitude"])] for r in rr
        ],dtype=float)

    pairs=[(s,t) for s in sources for t in targets]
    G=fast_symmetric_coverage(midpoint,sources,targets,500.0)
    centers={name:np.mean(midpoint[name],axis=0) for name in geometries}
    centroid=np.empty(len(pairs),dtype=float)
    edge_ratio=np.empty(len(pairs),dtype=float)
    locality_ratio=np.empty(len(pairs),dtype=float)
    same_order=np.empty(len(pairs),dtype=float)
    for i,(s,t) in enumerate(pairs):
        centroid[i]=np.log1p(float(np.linalg.norm(centers[s]-centers[t]))/500.0)
        edge_ratio[i]=abs(np.log(geometries[s].n_edges/geometries[t].n_edges))
        locality_ratio[i]=abs(np.log(geometries[s].n_localities/geometries[t].n_localities))
        same_order[i]=float(meta[s]["order"]==meta[t]["order"])
    X=np.column_stack((
        zscore(G),zscore(centroid),zscore(edge_ratio),zscore(locality_ratio),zscore(same_order)
    ))
    src=np.repeat(np.arange(len(sources),dtype=np.int64),len(targets))
    tgt=np.tile(np.arange(len(targets),dtype=np.int64),len(sources))
    prepared=prepare_dyadic_regression(src,tgt,X,primary_index=0)

    rule=json.loads(args.geometry_rule.read_text())
    sharding=json.loads(args.sharding.read_text())
    amps=tuple(map(float,rule["geometry_null_center"]["private_amplitudes"]))
    center_file=np.load(args.center_npz,allow_pickle=False)
    center=GeometryNullCenter(
        mu0=np.asarray(center_file["mu0"],dtype=float),
        amplitudes=tuple(map(float,center_file["amplitudes"])),
        worlds_per_amplitude=tuple(map(int,center_file["worlds_per_amplitude"])),
    )
    if center.amplitudes!=amps or center.worlds_per_amplitude!=(199,199,199,199):
        raise RuntimeError("authorized survivor center metadata drift")
    if center.mu0.shape!=(len(pairs),):
        raise RuntimeError("authorized survivor center dyad drift")
    labels=tuple(label_A(a) for a in amps)
    refs=load_private_references(args.reference_dir,sharding,labels)

    op=opening["response_operator"]
    edge_map=template_edges_from_geometries(geometries)
    cache=prepare_fixed_dyad_transfer_cache(
        edge_map,pairs,
        bandwidth_km=float(op["bandwidth_km"]),
        prior_strength=float(op["prior_strength"]),
        prior_mean=float(op["prior_mean"]),
        segment_points=int(op["segment_points"]),
        min_training_edges=int(op["min_endpoint_disjoint_training_edges"]),
    )

    # Single authorized nucleotide-identity opening starts here.
    response_map={}
    pair_counts={}
    empirical_species=tuple(sorted(set(sources)|set(targets)))
    for index,name in enumerate(empirical_species,1):
        mask_row=mask_by_species.get(name)
        if not mask_row or mask_row.get("status")!="PASS_MASK":
            raise RuntimeError(f"empirical species lacks frozen PASS mask: {name}")
        row=metadata[name]
        fasta=root/str(row["raw_dir"])/f"{row['raw_gene']}.afa"
        occurrence=root/str(row["raw_dir"])/"occurrences.txt"
        alignment=read_aligned_nucleotide_identity(fasta)
        if canonical_mask_sha256_from_identity(alignment)!=str(mask_row["mask_sha256"]):
            raise RuntimeError(f"canonical mask digest drift after identity opening: {name}")
        grouped=sequences_by_frozen_locality(
            alignment,read_occurrence_rows(occurrence),canonical_latlon[name]
        )
        distance=frozen_edge_mean_p_distances(
            grouped,geometries[name].edge_nodes,
            alignment_length=alignment.alignment_length,
            minimum_comparable_fraction=0.50,
        )
        response_map[name]=post_ibd_responses(
            edge_map[name],distance.genetic_distance,
            min_training_edges=int(op["min_endpoint_disjoint_training_edges"]),
        )
        pair_counts[name]={
            "edges":int(len(distance.genetic_distance)),
            "valid_sequence_pairs_total":int(np.sum(distance.valid_pair_counts)),
            "valid_sequence_pairs_min_per_edge":int(np.min(distance.valid_pair_counts)),
        }
        if index%50==0:
            print(json.dumps({"identity_species_done":index,"total":len(empirical_species)}))

    train={name:response_map[name].train_turnover[:,None] for name in sources}
    evaluation={name:response_map[name].eval_turnover[:,None] for name in targets}
    t_st=score_fixed_dyad_transfer_cache(cache,train,evaluation)[:,0]
    centered=center_dyad_scores(t_st,center)
    fit=batch_primary_test(prepared,centered[:,None])
    beta=float(fit.coefficient[0])
    se=float(fit.standard_error[0])
    z=float(fit.z_score[0])

    component=component_pvalues(beta,refs)
    envelope,least=least_favourable_beta_pvalues(np.asarray([beta]),refs)
    p=float(envelope[0])
    least_cell=str(least[0])
    alpha=float(opening["primary_inference"]["alpha"])
    positive=bool(beta>0 and p<=alpha)
    decision=(
        "CONFIRMATORY_POSITIVE_GEOGRAPHIC_RECURRENCE_CONDITIONAL_ON_SHARED_OPPORTUNITY"
        if positive else
        "CONFIRMATORY_QUALIFIED_NULL_SHARED_GEOGRAPHIC_OPPORTUNITY_INSUFFICIENT"
    )

    tq=np.quantile(t_st,[0,.1,.25,.5,.75,.9,1])
    eq=np.quantile(centered,[0,.1,.25,.5,.75,.9,1])
    payload={
        "schema":"ttf_genetic_codistributed_recurrence_empirical_result_v0.1",
        "status":"ONE_SHOT_CONFIRMATORY_EMPIRICAL_RESULT_OPENED",
        "authorization_sha256":sha256_path(args.authorization),
        "inputs_sha256":actual_inputs,
        "empirical_domain":{
            "species":len(empirical_species),
            "sources":len(sources),
            "targets":len(targets),
            "directed_dyads":len(pairs),
        },
        "primary":{
            "estimand":"beta_G for z(G_st) on centered E_st",
            "coefficient":beta,
            "two_way_cluster_standard_error_descriptive":se,
            "two_way_cluster_z_descriptive":z,
            "component_monte_carlo_p":component,
            "least_favourable_cell":least_cell,
            "primary_envelope_p_value":p,
            "alpha":alpha,
            "positive":positive,
        },
        "T_st_distribution":{
            "mean":float(np.mean(t_st)),"sd":float(np.std(t_st,ddof=0)),
            "min":float(tq[0]),"q10":float(tq[1]),"q25":float(tq[2]),
            "median":float(tq[3]),"q75":float(tq[4]),"q90":float(tq[5]),"max":float(tq[6]),
        },
        "centered_E_distribution":{
            "mean":float(np.mean(centered)),"sd":float(np.std(centered,ddof=0)),
            "min":float(eq[0]),"q10":float(eq[1]),"q25":float(eq[2]),
            "median":float(eq[3]),"q75":float(eq[4]),"q90":float(eq[5]),"max":float(eq[6]),
        },
        "predictor_condition_number":float(prepared.condition_number),
        "decision":decision,
        "pair_count_summary_by_species":pair_counts,
        "outcome_state":{
            "confirmatory_nucleotide_identity_opened":True,
            "confirmatory_pairwise_genetic_distances_computed_in_memory":True,
            "serialized_nucleotide_identity":False,
            "serialized_edge_genetic_distance_vectors":False,
            "serialized_dyad_T_st":False,
            "confirmatory_T_st_computed":True,
            "confirmatory_beta_G_computed":True,
        },
        "one_shot":{
            "post_result_retuning_allowed":False,
            "alternate_subgroup_result_allowed":False,
            "result_selection_rerun_allowed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "decision":decision,"beta_G":beta,"primary_p":p,
        "least_favourable_cell":least_cell,"dyads":len(pairs),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
