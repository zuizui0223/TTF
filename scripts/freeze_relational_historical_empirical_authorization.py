#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def digest_lines(values) -> str:
    return hashlib.sha256(("\n".join(map(str,values))+"\n").encode()).hexdigest()


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
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rule=json.loads(args.opening_rule.read_text())
    if rule.get("schema")!="ttf_relational_historical_climate_empirical_opening_rule_v0.1":
        raise RuntimeError("unexpected Study-C empirical opening rule")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("Study-C opening-rule firewall is open")

    summary=json.loads(args.survivor_summary.read_text())
    mask=json.loads(args.mask_result.read_text())
    qual=json.loads(args.survivor_qualification.read_text())
    if summary.get("status")!="PASS_TO_HISTORICAL_CONFIRMATORY_SURVIVOR_SYNTHETIC_REQUALIFICATION":
        raise RuntimeError("Study-C survivor summary is not eligible for requalification")
    if mask.get("status")!="PASS_TO_HISTORICAL_CONFIRMATORY_SURVIVOR_REQUALIFICATION":
        raise RuntimeError("Study-C mask result did not pass")
    if qual.get("schema")!="ttf_relational_historical_climate_qualification_result_v0.1":
        raise RuntimeError("unexpected Study-C survivor qualification schema")
    if qual.get("qualification_stage")!="confirmatory_survivor":
        raise RuntimeError("Study-C qualification is not survivor-stage")
    if qual.get("status")!="PASS_TO_HISTORICAL_EMPIRICAL_IDENTITY_OPENING_PREPARATION":
        raise RuntimeError("Study-C survivor requalification did not pass")
    if qual.get("gates",{}).get("overall_pass") is not True:
        raise RuntimeError("Study-C survivor requalification overall gate is not PASS")
    for payload in (summary,mask,qual):
        if any(bool(v) for v in payload["response_firewall"].values()):
            raise RuntimeError("Study-C pre-identity response firewall is open")

    design_sha=sha256_path(args.survivor_design)
    if summary.get("design_npz_sha256")!=design_sha:
        raise RuntimeError("Study-C survivor design hash drift")
    if qual.get("opportunity_design_sha256")!=design_sha:
        raise RuntimeError("Study-C survivor qualification design hash drift")

    data=np.load(args.survivor_design,allow_pickle=False)
    species=np.asarray(data["species_order"]).astype(str)
    s=np.asarray(data["survivor_source_index"],np.int64)
    t=np.asarray(data["survivor_target_index"],np.int64)
    if len(s)!=len(t) or len(s)==0:
        raise RuntimeError("invalid Study-C survivor dyads")
    source_names=sorted(set(map(str,species[s])))
    target_names=sorted(set(map(str,species[t])))
    if set(source_names)&set(target_names):
        raise RuntimeError("Study-C source and target roles overlap")
    union=sorted(set(source_names)|set(target_names))
    dyad_labels=[f"{species[int(a)]}\t{species[int(b)]}" for a,b in zip(s,t)]

    root=args.root.resolve()
    genes=root/"genes.txt"; cite=root/"cite.txt"
    if not genes.is_file() or not cite.is_file():
        raise FileNotFoundError("source root lacks genes.txt/cite.txt")

    code_paths=[
        Path("scripts/run_relational_historical_empirical.py"),
        Path("src/ttf/relational_genetic_empirical.py"),
        Path("src/ttf/relational_dyadic.py"),
        Path("src/ttf/phylogatr_empirical.py"),
        Path("src/ttf/phylogatr_compact_ibd.py"),
        Path("src/ttf/geometry_control.py"),
        Path("src/ttf/transfer.py"),
    ]
    missing=[str(p) for p in code_paths if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"missing frozen Study-C empirical code: {missing}")

    payload={
        "schema":"ttf_relational_historical_climate_identity_opening_authorization_v0.1",
        "status":"AUTHORIZE_ONE_SHOT_STUDY_C_NUCLEOTIDE_IDENTITY_OPENING",
        "source":{
            "archive_sha256":args.source_archive_sha256,
            "genes_txt_sha256":sha256_path(genes),
            "cite_txt_sha256":sha256_path(cite),
        },
        "inputs_sha256":{
            "survivor_design":design_sha,
            "survivor_summary":sha256_path(args.survivor_summary),
            "mask_result":sha256_path(args.mask_result),
            "survivor_qualification":sha256_path(args.survivor_qualification),
            "candidates":sha256_path(args.candidates),
            "opening_rule":sha256_path(args.opening_rule),
        },
        "empirical_design":{
            "source_species":len(source_names),
            "target_species":len(target_names),
            "species_union":len(union),
            "directed_dyads":len(dyad_labels),
            "source_digest_sha256":digest_lines(source_names),
            "target_digest_sha256":digest_lines(target_names),
            "union_digest_sha256":digest_lines(union),
            "dyad_digest_sha256":digest_lines(dyad_labels),
        },
        "code_sha256":{str(p):sha256_path(p) for p in code_paths},
        "relational_model":rule["relational_model"],
        "genetic_response":rule["genetic_response"],
        "one_shot_rules":rule["one_shot_rules"],
        "opening_state":{
            "Study_C_sequence_identity_opened":False,
            "Study_C_pairwise_genetic_distances_opened":False,
            "Study_C_T_st_computed":False,
            "Study_C_beta_hist_computed":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":payload["status"],**payload["empirical_design"]},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
