#!/usr/bin/env python3
"""Create the one-shot codistributed-recurrence identity-opening authorization.

This script never reads aligned FASTA sequence identity. It binds the qualified
survivor synthetic artifacts, frozen inputs, exact source provenance, and exact
empirical code bytes into a single authorization receipt.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.codistributed_formal_shards import validate_exact_ranges


EXPECTED_QUAL="ttf_genetic_codistributed_recurrence_survivor_fixed_tst_qualification_v0.1"
EXPECTED_OPENING="ttf_genetic_codistributed_recurrence_empirical_opening_rule_v0.1"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def reference_ranges(sharding: dict) -> list[tuple[int,int]]:
    return [
        (int(a),int(b))
        for a,b in sharding["shards"]["reference"]["per_amplitude"]
    ]


def validate_reference_artifacts(path: Path,sharding: dict) -> dict[str,str]:
    expected=reference_ranges(sharding)
    found={}
    labels=("A0p5","A1","A2","A3")
    for file in sorted(path.glob("*.npz")):
        p=np.load(file,allow_pickle=False)
        if str(np.asarray(p["component"]).item())!="reference":
            continue
        key=(int(np.asarray(p["start"]).item()),int(np.asarray(p["stop"]).item()))
        if key in found:
            raise RuntimeError(f"duplicate reference shard {key}")
        for label in labels:
            if label not in p.files:
                raise RuntimeError(f"reference shard {key} missing {label}")
            x=np.asarray(p[label],dtype=float)
            if x.shape!=(key[1]-key[0],) or not np.isfinite(x).all():
                raise RuntimeError(f"reference shard {key} payload drift for {label}")
        found[key]=file
    validate_exact_ranges(found.keys(),expected)
    return {
        found[key].name:sha256_path(found[key])
        for key in expected
    }


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
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    opening=json.loads(args.opening_rule.read_text())
    if opening.get("schema")!=EXPECTED_OPENING:
        raise RuntimeError("unexpected empirical opening rule")
    if opening.get("status")!="FROZEN_BEFORE_SURVIVOR_FIXED_TST_REQUALIFICATION_RESULT":
        raise RuntimeError("empirical opening rule is not prospectively frozen")
    if args.source_archive_sha256!=opening["source_binding"]["archive_sha256"]:
        raise RuntimeError("source archive SHA drift")
    if sha256_path(args.localities)!=opening["source_binding"]["geometry_localities_sha256"]:
        raise RuntimeError("survivor locality SHA drift")
    if sha256_path(args.edges)!=opening["source_binding"]["geometry_edges_sha256"]:
        raise RuntimeError("survivor edge SHA drift")
    if sha256_path(args.full_mask_ledger)!=opening["frozen_empirical_domain"]["full_mask_ledger_sha256"]:
        raise RuntimeError("full confirmatory mask ledger SHA drift")

    qual=json.loads(args.survivor_qualification.read_text())
    if qual.get("schema")!=EXPECTED_QUAL:
        raise RuntimeError("unexpected survivor qualification schema")
    if qual.get("status")!="PASS_TO_ONE_SHOT_CONFIRMATORY_EMPIRICAL_OPENING":
        raise RuntimeError("survivor qualification does not authorize identity opening")
    if qual.get("gate",{}).get("overall_pass") is not True:
        raise RuntimeError("survivor qualification lacks overall PASS")
    if any(bool(v) for v in qual["response_firewall"].values()):
        raise RuntimeError("survivor qualification response firewall is open")

    sharding=json.loads(args.sharding.read_text())
    reference_hashes=validate_reference_artifacts(args.reference_dir,sharding)
    binding=qual.get("artifact_binding")
    if not isinstance(binding,dict):
        raise RuntimeError("qualification lacks synthetic artifact binding")
    if sha256_path(args.center_npz)!=binding.get("center_npz_sha256"):
        raise RuntimeError("center NPZ is not the qualified center")
    if reference_hashes!=binding.get("reference_shards_sha256"):
        raise RuntimeError("reference shards are not the qualified references")

    root=args.root.resolve()
    genes=root/"genes.txt"
    cite=root/"cite.txt"
    if not genes.is_file() or not cite.is_file():
        raise RuntimeError("exact phylogatR source root is incomplete")

    inputs={
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
    code_paths=(
        "scripts/run_codistributed_recurrence_confirmatory_empirical.py",
        "src/ttf/phylogatr_empirical.py",
        "src/ttf/phylogatr_character_mask.py",
        "src/ttf/relational_genetic_empirical.py",
        "src/ttf/codistributed_geometry_null.py",
        "src/ttf/relational_dyadic.py",
    )
    code_hashes={}
    for name in code_paths:
        path=Path(name)
        if not path.is_file():
            raise RuntimeError(f"authorized code file missing: {name}")
        code_hashes[name]=sha256_path(path)

    center=np.load(args.center_npz,allow_pickle=False)
    if np.asarray(center["mu0"]).shape!=(26560,):
        raise RuntimeError("qualified center dyad-count drift")
    if tuple(map(float,np.asarray(center["amplitudes"])))!=(0.5,1.0,2.0,3.0):
        raise RuntimeError("qualified center amplitude drift")
    if tuple(map(int,np.asarray(center["worlds_per_amplitude"])))!=(199,199,199,199):
        raise RuntimeError("qualified center world-count drift")

    payload={
        "schema":"ttf_genetic_codistributed_recurrence_identity_opening_authorization_v0.1",
        "status":"AUTHORIZE_ONE_SHOT_CODISTRIBUTED_NUCLEOTIDE_IDENTITY_OPENING",
        "source":{
            "archive_sha256":args.source_archive_sha256,
            "genes_txt_sha256":sha256_path(genes),
            "cite_txt_sha256":sha256_path(cite),
        },
        "inputs_sha256":inputs,
        "reference_shards_sha256":reference_hashes,
        "qualification":{
            "sha256":sha256_path(args.survivor_qualification),
            "status":qual["status"],
            "overall_pass":qual["gate"]["overall_pass"],
        },
        "code_sha256":code_hashes,
        "opening_rule_sha256":sha256_path(args.opening_rule),
        "identity_opening":{
            "authorized_once":True,
            "nucleotide_identity_opened_at_authorization_time":False,
            "pairwise_genetic_distances_opened_at_authorization_time":False,
            "empirical_T_st_opened_at_authorization_time":False,
            "empirical_beta_G_opened_at_authorization_time":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "qualification_sha256":payload["qualification"]["sha256"],
        "reference_shards":len(reference_hashes),
        "code_files":len(code_hashes),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
