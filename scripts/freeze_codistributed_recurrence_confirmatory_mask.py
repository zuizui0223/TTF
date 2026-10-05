#!/usr/bin/env python3
"""Freeze the confirmatory canonical-character mask for codistributed recurrence.

Authorized only after the development actual-T_st formal qualification passes.
The script observes canonical validity (A/C/G/T versus everything else) only.
It never computes or persists nucleotide identity, pairwise genetic distance,
T_st, or beta_G.
"""
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

from ttf.phylogatr_character_mask import (
    CharacterMaskError,
    canonical_mask_sha256,
    edge_mask_support,
    masks_by_frozen_locality,
    read_canonical_mask_alignment,
)
from ttf.phylogatr_confirmatory import read_occurrence_rows


SOURCE_SHA="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
TAG="place-recurrence-insecta-fresh-v0.1"
EXPECTED_DESIGN="ttf_genetic_codistributed_recurrence_response_blind_design_v0.1"
EXPECTED_RULE="ttf_genetic_codistributed_recurrence_confirmatory_mask_rule_v0.1"
EXPECTED_FORMAL="ttf_genetic_codistributed_recurrence_fixed_tst_qualification_v0.3"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def digest_labels(names) -> str:
    return hashlib.sha256("\n".join(sorted(map(str,names))).encode()).hexdigest()


def rank_species(name: str, tag: str) -> str:
    return hashlib.sha256(f"{tag}|{SOURCE_SHA}|{name}".encode()).hexdigest()


def read_candidates(path: Path) -> list[dict[str,str]]:
    rows=list(csv.DictReader(path.open(encoding="utf-8")))
    names=[str(r["species"]) for r in rows]
    if len(names)!=len(set(names)):
        raise RuntimeError("duplicate candidate species")
    return rows


def read_excluded(path: Path) -> set[str]:
    p=json.loads(path.read_text())
    names=p.get("species")
    if not isinstance(names,list):
        raise RuntimeError(f"exclusion file lacks species list: {path}")
    return set(map(str,names))


def frozen_confirmatory_roles(rows, s1: set[str], s2: set[str], design: dict):
    eligible=[
        str(r["species"]) for r in rows
        if str(r["species"]) not in (s1|s2) and str(r["class"])=="Insecta"
    ]
    if len(eligible)!=762:
        raise RuntimeError(f"fresh Insecta domain drift: {len(eligible)} != 762")
    ordered=sorted(eligible,key=lambda s:(rank_species(s,TAG),s))
    if digest_labels(ordered)!=design["domain"]["species_digest_sha256"]:
        raise RuntimeError("fresh Insecta domain digest drift")
    confirmatory=ordered[381:]
    if digest_labels(confirmatory)!=design["selection"]["confirmatory_species_digest_sha256"]:
        raise RuntimeError("confirmatory panel digest drift")
    role_order=sorted(
        confirmatory,
        key=lambda s:(rank_species(s,f"{TAG}|confirmatory-role"),s),
    )
    source=tuple(role_order[:190])
    target=tuple(role_order[190:])
    if digest_labels(source)!=design["selection"]["confirmatory"]["source_digest_sha256"]:
        raise RuntimeError("confirmatory source digest drift")
    if digest_labels(target)!=design["selection"]["confirmatory"]["target_digest_sha256"]:
        raise RuntimeError("confirmatory target digest drift")
    return source,target


def authorize_formal(path: Path) -> dict:
    p=json.loads(path.read_text())
    if p.get("schema")!=EXPECTED_FORMAL:
        raise RuntimeError("unexpected formal qualification schema")
    if p.get("status")!="PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION":
        raise RuntimeError("formal qualification does not authorize character-mask opening")
    if p.get("gate",{}).get("overall_pass") is not True:
        raise RuntimeError("formal qualification overall gate is not PASS")
    if any(bool(v) for v in p["response_firewall"].values()):
        raise RuntimeError("formal qualification response firewall is open")
    return p


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--formal-qualification",type=Path,required=True)
    ap.add_argument("--mask-rule",type=Path,required=True)
    ap.add_argument("--s1-exclusion",type=Path,required=True)
    ap.add_argument("--s2-exclusion",type=Path,required=True)
    ap.add_argument("--output-mask-result",type=Path,required=True)
    ap.add_argument("--output-survivors",type=Path,required=True)
    args=ap.parse_args()

    formal=authorize_formal(args.formal_qualification)
    design=json.loads(args.design.read_text())
    rule=json.loads(args.mask_rule.read_text())
    if design.get("schema")!=EXPECTED_DESIGN:
        raise RuntimeError("unexpected codistributed design schema")
    if rule.get("schema")!=EXPECTED_RULE:
        raise RuntimeError("unexpected confirmatory mask rule schema")
    if rule.get("status")!="FROZEN_BEFORE_CONFIRMATORY_CHARACTER_MASK_AND_BEFORE_ANY_CONFIRMATORY_GENETIC_RESPONSE":
        raise RuntimeError("confirmatory mask rule is not prospectively frozen")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("confirmatory mask-rule response firewall is open")
    if design["source"]["archive_sha256"]!=SOURCE_SHA:
        raise RuntimeError("source archive binding drift")

    rows=read_candidates(args.candidates)
    meta={str(r["species"]):r for r in rows}
    source,target=frozen_confirmatory_roles(
        rows,
        read_excluded(args.s1_exclusion),
        read_excluded(args.s2_exclusion),
        design,
    )
    roles={name:"source" for name in source}
    roles.update({name:"target" for name in target})
    if len(roles)!=381:
        raise RuntimeError("confirmatory role union drift")

    fraction=float(rule["inherited_parser_contract"]["minimum_comparable_fraction"])
    if fraction!=0.50:
        raise RuntimeError("minimum comparable fraction drift")

    root=args.root.resolve()
    ledger=[]
    survivors=[]
    for index,name in enumerate(sorted(roles),1):
        row=meta[name]
        geometry,canonical_latlon=reconstruct_candidate(root,row,neighbor_fraction=0.15)
        verify_candidate(row,geometry)
        fasta=root/str(row["raw_dir"])/f"{row['raw_gene']}.afa"
        occurrence=root/str(row["raw_dir"])/"occurrences.txt"
        if not fasta.is_file() or not occurrence.is_file():
            raise RuntimeError(f"missing exact frozen source file for {name}")
        try:
            alignment=read_canonical_mask_alignment(fasta)
            grouped=masks_by_frozen_locality(
                alignment,
                read_occurrence_rows(occurrence),
                canonical_latlon,
            )
            support=edge_mask_support(
                grouped,
                geometry.edge_nodes,
                alignment_length=alignment.alignment_length,
                minimum_comparable_fraction=fraction,
            )
            survives=bool(support.all_edges_valid)
            item={
                "species":name,
                "role":roles[name],
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
        except CharacterMaskError as exc:
            survives=False
            item={
                "species":name,
                "role":roles[name],
                "status":"NOT_EVALUABLE_CHARACTER_SUPPORT",
                "reason":str(exc),
                "survives":False,
                "nucleotide_identity_persisted":False,
            }
        ledger.append(item)
        if survives:
            survivors.append((name,roles[name]))
        if index%50==0:
            print(json.dumps({"mask_species_done":index,"total":381},sort_keys=True))

    surviving_sources=sorted(n for n,r in survivors if r=="source")
    surviving_targets=sorted(n for n,r in survivors if r=="target")
    args.output_survivors.parent.mkdir(parents=True,exist_ok=True)
    with args.output_survivors.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["species","role"],lineterminator="\n")
        w.writeheader()
        for name,role in sorted(survivors):
            w.writerow({"species":name,"role":role})

    payload={
        "schema":"ttf_genetic_codistributed_recurrence_confirmatory_mask_result_v0.1",
        "status":"PASS_TO_SURVIVOR_INFORMATION_GATE",
        "formal_qualification_sha256":sha256_path(args.formal_qualification),
        "mask_rule_sha256":sha256_path(args.mask_rule),
        "design_sha256":sha256_path(args.design),
        "source_archive_sha256":SOURCE_SHA,
        "species_before_mask":381,
        "sources_before_mask":190,
        "targets_before_mask":191,
        "surviving_species":len(survivors),
        "surviving_sources":len(surviving_sources),
        "surviving_targets":len(surviving_targets),
        "surviving_source_digest_sha256":digest_labels(surviving_sources),
        "surviving_target_digest_sha256":digest_labels(surviving_targets),
        "survivor_csv_sha256":sha256_path(args.output_survivors),
        "species":ledger,
        "formal_gate":{
            "status":formal["status"],
            "overall_pass":formal["gate"]["overall_pass"],
        },
        "response_firewall":{
            "confirmatory_nucleotide_identity_opened":False,
            "confirmatory_pairwise_genetic_distances_opened":False,
            "confirmatory_T_st_opened":False,
            "confirmatory_beta_G_opened":False,
        },
    }
    args.output_mask_result.parent.mkdir(parents=True,exist_ok=True)
    args.output_mask_result.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "surviving_species":len(survivors),
        "surviving_sources":len(surviving_sources),
        "surviving_targets":len(surviving_targets),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
