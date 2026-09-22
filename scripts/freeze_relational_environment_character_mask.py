#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

try:
    from scripts.reconstruct_relational_environment_geometry import reconstruct_candidate
except ModuleNotFoundError:
    from reconstruct_relational_environment_geometry import reconstruct_candidate

from ttf.phylogatr_character_mask import (
    CharacterMaskError,
    canonical_mask_sha256,
    edge_mask_support,
    masks_by_frozen_locality,
    read_canonical_mask_alignment,
)
from ttf.phylogatr_confirmatory import read_occurrence_rows


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def read_candidates(path: Path) -> dict[str,dict[str,str]]:
    rows=list(csv.DictReader(path.open(encoding="utf-8")))
    out={str(row["species"]):row for row in rows}
    if len(out)!=len(rows):
        raise RuntimeError("duplicate candidate species")
    return out


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,required=True)
    ap.add_argument("--opportunity-design",type=Path,required=True)
    ap.add_argument("--opportunity-summary",type=Path,required=True)
    ap.add_argument("--development-qualification",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--mask-rule",type=Path,required=True)
    ap.add_argument("--output-mask-result",type=Path,required=True)
    ap.add_argument("--output-survivor-design",type=Path,required=True)
    ap.add_argument("--output-survivor-summary",type=Path,required=True)
    args=ap.parse_args()

    mask_rule=json.loads(args.mask_rule.read_text())
    if mask_rule.get("schema")!="ttf_relational_environment_character_mask_rule_v0.1":
        raise RuntimeError("unexpected Study-B mask rule")
    if any(bool(v) for v in mask_rule["response_firewall"].values()):
        raise RuntimeError("Study-B mask-rule response firewall is open")

    opportunity=json.loads(args.opportunity_summary.read_text())
    if opportunity.get("status")!="PASS_TO_DEVELOPMENT_SYNTHETIC_QUALIFICATION":
        raise RuntimeError("Study-B opportunity design did not pass")
    if any(bool(v) for v in opportunity["response_firewall"].values()):
        raise RuntimeError("Study-B opportunity response firewall is open")

    qualification=json.loads(args.development_qualification.read_text())
    if qualification.get("schema")!="ttf_relational_environment_qualification_result_v0.2":
        raise RuntimeError("unexpected Study-B development qualification result")
    if qualification.get("status")!="PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION":
        raise RuntimeError("Study-B development qualification did not authorize mask opening")
    if qualification.get("gates",{}).get("overall_pass") is not True:
        raise RuntimeError("Study-B development qualification gate is not PASS")
    if any(bool(v) for v in qualification["response_firewall"].values()):
        raise RuntimeError("Study-B development qualification response firewall is open")
    if qualification.get("opportunity_design_sha256")!=sha256_path(args.opportunity_design):
        raise RuntimeError("Study-B opportunity-design hash drift at mask opening")

    data=np.load(args.opportunity_design,allow_pickle=False)
    species=np.asarray(data["species_order"]).astype(str)
    source=np.asarray(data["confirmatory_source_index"],np.int64)
    target=np.asarray(data["confirmatory_target_index"],np.int64)
    arrays={
        "R_env":np.asarray(data["confirmatory_R_env"],float),
        "coverage":np.asarray(data["confirmatory_coverage"],float),
        "same_class":np.asarray(data["confirmatory_same_class"],np.int8),
        "same_order":np.asarray(data["confirmatory_same_order"],np.int8),
        "same_family":np.asarray(data["confirmatory_same_family"],np.int8),
        "locality_count_ratio":np.asarray(data["confirmatory_locality_count_ratio"],float),
    }
    n=len(source)
    if any(len(v)!=n for v in arrays.values()) or len(target)!=n:
        raise RuntimeError("Study-B confirmatory opportunity-array width drift")

    empirical_union=sorted(set(map(str,species[source]))|set(map(str,species[target])))
    metadata=read_candidates(args.candidates)
    if not set(empirical_union).issubset(metadata):
        raise RuntimeError("Study-B confirmatory species absent from frozen candidate metadata")

    fraction=float(mask_rule["inherited_parser_contract"]["minimum_comparable_fraction"])
    if fraction!=0.50:
        raise RuntimeError("Study-B mask comparable-site fraction drift")

    root=args.root.resolve()
    ledger={}
    survivors=[]
    for i,name in enumerate(empirical_union,1):
        row=metadata[name]
        geometry,canonical_latlon=reconstruct_candidate(
            root,row,neighbor_fraction=0.15
        )
        fasta=root/str(row["raw_dir"])/f"{row['raw_gene']}.afa"
        occurrence=root/str(row["raw_dir"])/"occurrences.txt"
        if not fasta.is_file() or not occurrence.is_file():
            raise RuntimeError(f"missing frozen source file for {name}")
        try:
            alignment=read_canonical_mask_alignment(fasta)
        except CharacterMaskError as exc:
            ledger[name]={
                "species":name,
                "status":"NOT_EVALUABLE_CHARACTER_SUPPORT",
                "reason":str(exc),
                "survives":False,
                "nucleotide_identity_persisted":False,
            }
            continue
        grouped=masks_by_frozen_locality(
            alignment,read_occurrence_rows(occurrence),canonical_latlon
        )
        support=edge_mask_support(
            grouped,
            geometry.edge_nodes,
            alignment_length=alignment.alignment_length,
            minimum_comparable_fraction=fraction,
        )
        survives=bool(support.all_edges_valid)
        ledger[name]={
            "species":name,
            "status":"PASS_MASK" if survives else "NOT_EVALUABLE_CHARACTER_SUPPORT",
            "reason":None if survives else "one_or_more_frozen_edges_lack_valid_cross_locality_pair",
            "survives":survives,
            "alignment_length":int(alignment.alignment_length),
            "aligned_records":int(len(alignment.headers)),
            "frozen_edges":int(len(support.valid_edges)),
            "valid_edges":int(np.count_nonzero(support.valid_edges)),
            "invalid_edges":int(np.count_nonzero(~support.valid_edges)),
            "minimum_comparable_columns":int(support.minimum_comparable_columns),
            "mask_sha256":canonical_mask_sha256(alignment),
            "nucleotide_identity_persisted":False,
        }
        if survives:
            survivors.append(name)
        if i%50==0:
            print(json.dumps({"mask_species_done":i,"total":len(empirical_union)}))

    survivor_set=set(survivors)
    pair_keep=np.asarray([
        str(species[int(s)]) in survivor_set and str(species[int(t)]) in survivor_set
        for s,t in zip(source,target)
    ],bool)
    target_counts=Counter(map(int,target[pair_keep]))
    minimum_sources=int(mask_rule["post_mask_dyad_support"]["minimum_surviving_sources_per_target"])
    supported_targets={t for t,count in target_counts.items() if count>=minimum_sources}
    pair_keep &= np.asarray([int(t) in supported_targets for t in target],bool)

    ss=source[pair_keep]
    tt=target[pair_keep]
    survivor_arrays={key:value[pair_keep] for key,value in arrays.items()}
    source_species=sorted(set(map(str,species[ss])))
    target_species=sorted(set(map(str,species[tt])))
    empirical_species=sorted(set(source_species)|set(target_species))

    args.output_survivor_design.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(
        args.output_survivor_design,
        species_order=species,
        survivor_source_index=ss,
        survivor_target_index=tt,
        survivor_R_env=survivor_arrays["R_env"],
        survivor_coverage=survivor_arrays["coverage"],
        survivor_same_class=survivor_arrays["same_class"],
        survivor_same_order=survivor_arrays["same_order"],
        survivor_same_family=survivor_arrays["same_family"],
        survivor_locality_count_ratio=survivor_arrays["locality_count_ratio"],
    )

    mask_payload={
        "schema":"ttf_relational_environment_character_mask_result_v0.1",
        "status":"PASS_TO_CONFIRMATORY_SURVIVOR_REQUALIFICATION",
        "opportunity_design_sha256":sha256_path(args.opportunity_design),
        "development_qualification_sha256":sha256_path(args.development_qualification),
        "mask_rule_sha256":sha256_path(args.mask_rule),
        "species_before_mask":len(empirical_union),
        "surviving_species_before_dyad_support":len(survivors),
        "surviving_empirical_union_after_dyad_support":len(empirical_species),
        "source_clusters_after_mask":len(source_species),
        "target_clusters_after_mask":len(target_species),
        "directed_dyads_after_mask":int(len(ss)),
        "species": [ledger[name] for name in sorted(ledger)],
        "response_firewall":{
            "Study_B_sequence_identity_opened":False,
            "Study_B_pairwise_genetic_distances_opened":False,
            "Study_B_T_st_computed":False,
            "Study_B_beta_R_computed":False,
        },
    }
    args.output_mask_result.parent.mkdir(parents=True,exist_ok=True)
    args.output_mask_result.write_text(json.dumps(mask_payload,indent=2,sort_keys=True)+"\n")

    summary={
        "schema":"ttf_relational_environment_survivor_opportunity_design_v0.1",
        "status":"PASS_TO_CONFIRMATORY_SURVIVOR_SYNTHETIC_REQUALIFICATION",
        "opportunity_design_sha256":sha256_path(args.opportunity_design),
        "mask_result_sha256":sha256_path(args.output_mask_result),
        "design_npz_sha256":sha256_path(args.output_survivor_design),
        "survivor":{
            "mask_survivors":len(survivors),
            "empirical_union":len(empirical_species),
            "source_clusters":len(source_species),
            "target_clusters":len(target_species),
            "dyads":int(len(ss)),
        },
        "response_firewall":{
            "Study_B_sequence_identity_opened":False,
            "Study_B_pairwise_genetic_distances_opened":False,
            "Study_B_T_st_computed":False,
            "Study_B_beta_R_computed":False,
        },
    }
    args.output_survivor_summary.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":summary["status"],**summary["survivor"]},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
