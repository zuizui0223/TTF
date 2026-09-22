#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

try:
    from scripts.reconstruct_relational_environment_geometry import reconstruct_candidate, verify_candidate
except ModuleNotFoundError:
    from reconstruct_relational_environment_geometry import reconstruct_candidate, verify_candidate

from ttf.core import SpeciesEdges
from ttf.phylogatr_confirmatory import read_occurrence_rows
from ttf.phylogatr_empirical import (
    canonical_mask_sha256_from_identity,
    frozen_edge_mean_p_distances,
    read_aligned_nucleotide_identity,
    sequences_by_frozen_locality,
)
from ttf.relational_dyadic import batch_primary_test, prepare_dyadic_regression
from ttf.relational_genetic_empirical import post_ibd_responses, source_only_transfer_scores


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def digest_lines(values) -> str:
    return hashlib.sha256(("\n".join(map(str,values))+"\n").encode()).hexdigest()


def zscore(x):
    x=np.asarray(x,float)
    sd=float(x.std())
    if not np.isfinite(sd) or sd<=np.finfo(float).eps:
        raise ValueError("constant/non-finite predictor")
    return (x-float(x.mean()))/sd


def remap(values):
    vals=sorted(map(int,np.unique(values)))
    lookup={v:i for i,v in enumerate(vals)}
    return np.asarray([lookup[int(v)] for v in values],np.int64)


def read_candidates(path: Path):
    rows=list(csv.DictReader(path.open(encoding="utf-8")))
    out={str(row["species"]):row for row in rows}
    if len(out)!=len(rows):
        raise RuntimeError("duplicate Study-C candidate species")
    return out


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,required=True)
    ap.add_argument("--source-archive-sha256",required=True)
    ap.add_argument("--survivor-design",type=Path,required=True)
    ap.add_argument("--survivor-summary",type=Path,required=True)
    ap.add_argument("--mask-result",type=Path,required=True)
    ap.add_argument("--survivor-qualification",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--opening-rule",type=Path,required=True)
    ap.add_argument("--authorization",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    auth=json.loads(args.authorization.read_text())
    if auth.get("schema")!="ttf_relational_historical_climate_identity_opening_authorization_v0.1":
        raise RuntimeError("unexpected Study-C identity-opening authorization")
    if auth.get("status")!="AUTHORIZE_ONE_SHOT_STUDY_C_NUCLEOTIDE_IDENTITY_OPENING":
        raise RuntimeError("Study-C identity opening is not authorized")

    actual_inputs={
        "survivor_design":sha256_path(args.survivor_design),
        "survivor_summary":sha256_path(args.survivor_summary),
        "mask_result":sha256_path(args.mask_result),
        "survivor_qualification":sha256_path(args.survivor_qualification),
        "candidates":sha256_path(args.candidates),
        "opening_rule":sha256_path(args.opening_rule),
    }
    if actual_inputs!=auth["inputs_sha256"]:
        raise RuntimeError("Study-C authorized input hash drift")
    if args.source_archive_sha256!=auth["source"]["archive_sha256"]:
        raise RuntimeError("Study-C source archive SHA drift")

    root=args.root.resolve()
    if sha256_path(root/"genes.txt")!=auth["source"]["genes_txt_sha256"]:
        raise RuntimeError("genes.txt drift since Study-C authorization")
    if sha256_path(root/"cite.txt")!=auth["source"]["cite_txt_sha256"]:
        raise RuntimeError("cite.txt drift since Study-C authorization")
    for path,expected in auth["code_sha256"].items():
        p=Path(path)
        if not p.is_file() or sha256_path(p)!=expected:
            raise RuntimeError(f"Study-C authorized code drift: {path}")

    rule=json.loads(args.opening_rule.read_text())
    if rule.get("schema")!="ttf_relational_historical_climate_empirical_opening_rule_v0.1":
        raise RuntimeError("unexpected Study-C empirical opening rule")
    mask=json.loads(args.mask_result.read_text())
    summary=json.loads(args.survivor_summary.read_text())
    qual=json.loads(args.survivor_qualification.read_text())
    for payload in (mask,summary,qual):
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("Study-C pre-identity response firewall is open")
    if qual.get("status")!="PASS_TO_HISTORICAL_EMPIRICAL_IDENTITY_OPENING_PREPARATION":
        raise RuntimeError("Study-C survivor requalification did not pass")

    data=np.load(args.survivor_design,allow_pickle=False)
    species=np.asarray(data["species_order"]).astype(str)
    s=np.asarray(data["survivor_source_index"],np.int64)
    t=np.asarray(data["survivor_target_index"],np.int64)
    r_hist=np.asarray(data["survivor_R_hist"],float)
    r_current=np.asarray(data["survivor_R_current"],float)
    coverage=np.asarray(data["survivor_coverage"],float)
    same_class=np.asarray(data["survivor_same_class"],float)
    same_order=np.asarray(data["survivor_same_order"],float)
    same_family=np.asarray(data["survivor_same_family"],float)
    locality_ratio=np.asarray(data["survivor_locality_count_ratio"],float)
    n=len(s)
    if any(len(x)!=n for x in (t,r_hist,r_current,coverage,same_class,same_order,same_family,locality_ratio)):
        raise RuntimeError("Study-C survivor dyad-array drift")

    source_names=sorted(set(map(str,species[s])))
    target_names=sorted(set(map(str,species[t])))
    empirical_species=sorted(set(source_names)|set(target_names))
    dyad_labels=[f"{species[int(a)]}\t{species[int(b)]}" for a,b in zip(s,t)]
    frozen=auth["empirical_design"]
    checks={
        "source_species":len(source_names),
        "target_species":len(target_names),
        "species_union":len(empirical_species),
        "directed_dyads":len(dyad_labels),
        "source_digest_sha256":digest_lines(source_names),
        "target_digest_sha256":digest_lines(target_names),
        "union_digest_sha256":digest_lines(empirical_species),
        "dyad_digest_sha256":digest_lines(dyad_labels),
    }
    if checks!=frozen:
        raise RuntimeError("Study-C authorized empirical design drift")
    if set(source_names)&set(target_names):
        raise RuntimeError("Study-C source/target role overlap")

    metadata=read_candidates(args.candidates)
    mask_by_species={str(row["species"]):row for row in mask["species"]}
    if not set(empirical_species).issubset(metadata) or not set(empirical_species).issubset(mask_by_species):
        raise RuntimeError("Study-C empirical species missing frozen metadata/mask")

    edge_map={}
    canonical_latlon={}
    for name in empirical_species:
        geometry,latlon=reconstruct_candidate(root,metadata[name],neighbor_fraction=0.15)
        verify_candidate(metadata[name],geometry)
        nodes=np.asarray(geometry.edge_nodes,np.int64)
        coords=np.asarray(geometry.coordinates,float)
        start=coords[nodes[:,0]]
        end=coords[nodes[:,1]]
        edge_map[name]=SpeciesEdges(
            species=name,nodes=nodes,start=start,end=end,
            midpoint=0.5*(start+end),
            length=np.linalg.norm(end-start,axis=1),
            turnover=np.zeros(len(nodes),float),
        )
        canonical_latlon[name]=latlon

    # Single authorized Study-C nucleotide identity opening.
    response_map={}
    pair_count_summary={}
    for i,name in enumerate(empirical_species,1):
        row=metadata[name]
        mask_row=mask_by_species[name]
        if mask_row.get("status")!="PASS_MASK":
            raise RuntimeError(f"Study-C empirical species lacks PASS mask: {name}")
        fasta=root/str(row["raw_dir"])/f"{row['raw_gene']}.afa"
        occurrence=root/str(row["raw_dir"])/"occurrences.txt"
        alignment=read_aligned_nucleotide_identity(fasta)
        if canonical_mask_sha256_from_identity(alignment)!=str(mask_row["mask_sha256"]):
            raise RuntimeError(f"Study-C canonical mask digest drift: {name}")
        grouped=sequences_by_frozen_locality(
            alignment,read_occurrence_rows(occurrence),canonical_latlon[name]
        )
        distance=frozen_edge_mean_p_distances(
            grouped,
            edge_map[name].nodes,
            alignment_length=alignment.alignment_length,
            minimum_comparable_fraction=0.50,
        )
        response_map[name]=post_ibd_responses(
            edge_map[name],distance.genetic_distance,min_training_edges=5
        )
        pair_count_summary[name]={
            "edges":int(len(distance.genetic_distance)),
            "valid_sequence_pairs_total":int(np.sum(distance.valid_pair_counts)),
            "valid_sequence_pairs_min_per_edge":int(np.min(distance.valid_pair_counts)),
        }
        if i%50==0:
            print(json.dumps({"identity_species_done":i,"total":len(empirical_species)}))

    pairs=[(str(species[int(a)]),str(species[int(b)])) for a,b in zip(s,t)]
    kernel=auth["genetic_response"]["kernel"]
    t_st=source_only_transfer_scores(
        edge_map,response_map,pairs,
        bandwidth_km=float(kernel["bandwidth_km"]),
        prior_strength=float(kernel["prior_strength"]),
        prior_mean=float(kernel["prior_mean"]),
        segment_points=int(kernel["segment_points"]),
    )

    predictors=np.column_stack((
        zscore(r_hist),
        zscore(r_current),
        zscore(coverage),
        same_class-same_class.mean(),
        same_order-same_order.mean(),
        same_family-same_family.mean(),
        zscore(np.abs(np.log(locality_ratio))),
    ))
    prepared=prepare_dyadic_regression(remap(s),remap(t),predictors,primary_index=0)
    expected_condition=float(qual["geometry"]["predictor_condition_number"])
    if not np.isclose(prepared.condition_number,expected_condition,rtol=0.0,atol=1e-9):
        raise RuntimeError("Study-C authorized predictor condition-number drift")
    fit=batch_primary_test(prepared,t_st[:,None])
    beta=float(fit.coefficient[0])
    se=float(fit.standard_error[0])
    z=float(fit.z_score[0])
    p=float(fit.p_value_one_sided[0])
    alpha=float(auth["relational_model"]["alpha_one_sided"])
    positive=bool(beta>0 and p<=alpha)
    decision=(
        "STUDY_C_HISTORICAL_RELATIONAL_POSITIVE"
        if positive else
        "STUDY_C_HISTORICAL_RELATIONAL_NULL_WITH_QUALIFIED_POWER"
    )

    q=np.quantile(t_st,[0,.1,.25,.5,.75,.9,1])
    dyads=[
        {
            "source":str(species[int(a)]),
            "target":str(species[int(b)]),
            "T_st":float(score),
            "R_hist":float(rh),
            "R_current":float(rc),
            "geographic_coverage":float(gg),
            "same_class":bool(sc),
            "same_order":bool(so),
            "same_family":bool(sf),
            "locality_count_ratio":float(lr),
        }
        for a,b,score,rh,rc,gg,sc,so,sf,lr in zip(
            s,t,t_st,r_hist,r_current,coverage,same_class,same_order,same_family,locality_ratio
        )
    ]

    payload={
        "schema":"ttf_relational_historical_climate_empirical_result_v0.1",
        "status":"EMPIRICAL_STUDY_C_RESULT_OPENED_UNDER_FROZEN_AUTHORIZATION",
        "authorization_sha256":sha256_path(args.authorization),
        "input_sha256":actual_inputs,
        "empirical_design":checks,
        "primary":{
            "estimand":"beta_hist for z_R_hist controlling z_R_current",
            "coefficient":beta,
            "standard_error_two_way_cluster":se,
            "z_score":z,
            "p_value_one_sided":p,
            "alpha":alpha,
            "positive":positive,
        },
        "T_st_distribution":{
            "mean":float(np.mean(t_st)),
            "sd":float(np.std(t_st)),
            "min":float(q[0]),"q10":float(q[1]),"q25":float(q[2]),
            "median":float(q[3]),"q75":float(q[4]),"q90":float(q[5]),"max":float(q[6]),
        },
        "decision":decision,
        "pair_count_summary_by_species":pair_count_summary,
        "dyads":dyads,
        "outcome_state":{
            "Study_C_sequence_identity_opened":True,
            "Study_C_pairwise_genetic_distances_computed_in_memory":True,
            "serialized_sequence_identity":False,
            "serialized_edge_genetic_distance_vectors":False,
            "Study_C_T_st_computed":True,
            "Study_C_beta_hist_computed":True,
        },
        "one_shot":{
            "post_result_retuning_allowed":False,
            "result_selection_rerun_allowed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "decision":decision,"beta_hist":beta,"se":se,"z":z,"p":p,"dyads":len(t_st)
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
