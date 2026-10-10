#!/usr/bin/env python3
"""Confirmatory-only Boolean canonical-mask support under frozen TTF host-memory rules.

Reads FASTA into A/C/G/T-validity masks; never stores a nucleotide identity,
genetic distance, post-IBD response or historical host genetic coefficient.
Must run only after exact archive/source and synthetic PASS authorization.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from ttf.phylogatr_character_mask import (
    CharacterMaskError, canonical_mask_sha256, edge_mask_support,
    masks_by_frozen_locality, read_canonical_mask_alignment
)
from ttf.phylogatr_confirmatory import read_occurrence_rows

CAN_SHA="c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b"
ROLE_SHA="4168e7572d3378e7f1595ae562d53c81bb2ef165cbe50eb2ed8ec4e7406fc6ae"
LOC_SHA="037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9"
EDGE_SHA="ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9"
ARCH_SHA="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce61a5"
SYNTHETIC_RECEIPT_SHA="0498f71d636274a7f2e654f42591545cf26cbd0b97e603902377ea2ef702495e"
MASK_RULE_BLOB_SHA1="1d5bec837157955cd2f989b2dd127d502d57318c"
SYNTHETIC_WORLD_SHA256="2aa61a0fd3e0fbf9fdf44f1ef7c4cffe4a6137cc91543d9971902d1349ebbf71"


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob_sha1(path: Path) -> str:
    """Pin even the mask-only rule bytes, not merely its JSON schema name."""
    raw = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\\0" + raw).hexdigest()


def mask_source_bindings(mask_rule_sha256: str) -> dict[str, str]:
    """Immutable provenance passed forward to both survivor gates."""
    return {
        "candidate_table_sha256": CAN_SHA,
        "panel_roles_sha256": ROLE_SHA,
        "locality_geometry_sha256": LOC_SHA,
        "edge_geometry_sha256": EDGE_SHA,
        "synthetic_pass_receipt_sha256": SYNTHETIC_RECEIPT_SHA,
        "mask_rule_git_blob_sha1": MASK_RULE_BLOB_SHA1,
        "mask_rule_sha256": mask_rule_sha256,
    }


def safe_raw_path(root: Path, raw_dir: str, filename: str) -> Path:
    raw = Path(raw_dir)
    if raw.is_absolute() or ".." in raw.parts or Path(filename).name!=filename:
        raise RuntimeError("frozen raw archive path traversal")
    base = root.resolve()
    result = (base/raw/filename).resolve()
    if base not in result.parents:
        raise RuntimeError("frozen raw path escapes exact archive")
    return result


def read_frozen_geo(local_csv: Path, edge_csv: Path, target: set[str]):
    ll=defaultdict(dict)
    ee=defaultdict(dict)
    with local_csv.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            sp=row["species"]
            if sp in target:
                idx=int(row["locality_index"])
                if idx in ll[sp]:raise RuntimeError("duplicate frozen locality")
                ll[sp][idx]=(float(row["latitude"]),float(row["longitude"]))
    with edge_csv.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            sp=row["species"]
            if sp in target:
                idx=int(row["edge_index"])
                if idx in ee[sp]:raise RuntimeError("duplicate frozen edge")
                ee[sp][idx]=(int(row["node_left"]),int(row["node_right"]))
    if set(ll)!=target or set(ee)!=target:
        raise RuntimeError("confirmatory geometry source does not cover exact panel")
    return ll, ee


def score_species(row: dict, root: Path, loc: dict, edges: dict, comparable_fraction=0.5):
    name=row["species"]
    if sorted(loc)!=list(range(len(loc))) or sorted(edges)!=list(range(len(edges))):
        raise RuntimeError(f"noncontiguous frozen graph/locality IDs for {name}")
    if int(row["n_localities"])!=len(loc) or int(row["edges"])!=len(edges):
        raise RuntimeError(f"frozen candidate graph geometry changed for {name}")
    ordered_local=np.asarray([loc[k] for k in sorted(loc)],dtype=float)
    nodes=np.asarray([edges[k] for k in sorted(edges)],dtype=int)
    if np.any(nodes<0) or np.any(nodes>=len(loc)):
        raise RuntimeError(f"frozen graph edge references nonexistent locality in {name}")
    fasta=safe_raw_path(root,row["raw_dir"],row["raw_gene"]+".afa")
    occurrence=safe_raw_path(root,row["raw_dir"],"occurrences.txt")
    if not fasta.is_file() or not occurrence.is_file():
        raise FileNotFoundError(f"exact raw source absent for {name}")
    try:
        alignment=read_canonical_mask_alignment(fasta)
    except CharacterMaskError as exc:
        return {"species":name,"status":"NOT_EVALUABLE_CHARACTER_SUPPORT",
                "reason_code":"malformed_or_duplicate_frozen_alignment",
                "survives":False,"nucleotide_identity_persisted":False}
    if int(row["n_headers"])!=len(alignment.headers):
        raise RuntimeError(f"frozen FASTA header count changed for {name}")
    masks=masks_by_frozen_locality(
        alignment,read_occurrence_rows(occurrence),ordered_local
    )
    score=edge_mask_support(masks,nodes,alignment_length=alignment.alignment_length,
                            minimum_comparable_fraction=comparable_fraction)
    passes=bool(score.all_edges_valid)
    return {
        "species":name,
        "status":"PASS_MASK" if passes else "NOT_EVALUABLE_CHARACTER_SUPPORT",
        "reason_code":None if passes else "one_or_more_frozen_edges_have_no_valid_cross_locality_pair",
        "survives":passes,
        "frozen_edges":len(edges),
        "valid_edges":int(np.count_nonzero(score.valid_edges)),
        "invalid_edges":int(np.count_nonzero(~score.valid_edges)),
        "alignment_length":int(alignment.alignment_length),
        "aligned_records":len(alignment.headers),
        "minimum_comparable_columns":int(score.minimum_comparable_columns),
        "mask_sha256":canonical_mask_sha256(alignment),
        "nucleotide_identity_persisted":False,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-archive",type=Path,required=True)
    ap.add_argument("--extracted-root",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--roles",type=Path,required=True)
    ap.add_argument("--localities",type=Path,required=True)
    ap.add_argument("--edges",type=Path,required=True)
    ap.add_argument("--synthetic-receipt",type=Path,required=True)
    ap.add_argument("--rule",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    pins={"source_archive":ARCH_SHA,"candidates":CAN_SHA,"roles":ROLE_SHA,
          "localities":LOC_SHA,"edges":EDGE_SHA,
          "synthetic_receipt":SYNTHETIC_RECEIPT_SHA}
    for k,d in pins.items():
        p=getattr(args,k)
        if k=="source_archive" and p.stat().st_size!=274_988_692:
            raise RuntimeError("exact raw archive size mismatch")
        if sha256(p)!=d:
            raise RuntimeError(f"frozen {k} SHA256 mismatch")
    rule=json.loads(args.rule.read_text())
    synth=json.loads(args.synthetic_receipt.read_text())
    if rule.get("schema")!="ttf_historical_host_memory_character_mask_rule_v0.1":
        raise RuntimeError("mask opening rule not frozen")
    if git_blob_sha1(args.rule)!=MASK_RULE_BLOB_SHA1:
        raise RuntimeError("exact frozen mask-rule byte identity mismatch")
    if (synth.get("schema")!="ttf_historical_host_memory_synthetic_qualification_v0.1"
            or synth.get("contract_sha256")!=SYNTHETIC_WORLD_SHA256
            or synth.get("confirmatory_species")!=321
            or synth.get("no_empirical_genetic_response_read") is not True):
        raise RuntimeError("wrong original synthetic authorization source")
    if synth["calibration"]["decision"]!="PASS_TO_CONFIRMATORY_CHARACTER_MASK_ONLY":
        raise RuntimeError("synthetic result does not authorize mask exposure")
    if synth["nucleotide_identity_opened"] or rule["firewall"]["confirmatory_nucleotide_identity_opened"]:
        raise RuntimeError("genetic identity firewall already open")
    if float(rule["minimum_fraction_comparable"])!=0.5:
        raise RuntimeError("frozen canonical mask overlap threshold drift")
    with args.candidates.open(newline="",encoding="utf-8") as f:
        candidates={r["species"]:r for r in csv.DictReader(f)}
    if len(candidates)!=642:
        raise RuntimeError("frozen candidate count drift")
    with args.roles.open(newline="",encoding="utf-8") as f:
        role=list(csv.DictReader(f))
    confirm=sorted(r["species"] for r in role if r["panel"]=="confirmatory")
    if len(confirm)!=321 or len(set(confirm))!=321 or not set(confirm).issubset(candidates):
        raise RuntimeError("confirmatory panel role/source mismatch")
    ll,ee=read_frozen_geo(args.localities,args.edges,set(confirm))
    ledger=[score_species(candidates[sp],args.extracted_root,ll[sp],ee[sp]) for sp in confirm]
    survivors=[x["species"] for x in ledger if x["survives"]]
    minimum=int(rule["minimum_confirmatory_survivors"])
    if minimum!=200:
        raise RuntimeError("frozen survivor floor changed")
    output={
        "schema":"ttf_historical_host_memory_confirmatory_mask_result_v0.1",
        "status":("PASS_TO_EXACT_SURVIVOR_INFORMATION_AND_SYNTHETIC_REQUALIFICATION"
                  if len(survivors)>=minimum else
                  "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_CHARACTER_SUPPORT"),
        "source_archive_sha256":ARCH_SHA,
        "source_bindings":mask_source_bindings(sha256(args.rule)),
        "confirmatory_before_mask":321,
        "surviving_species":len(survivors),
        "minimum_survivors":minimum,
        "survivor_names":survivors,
        "ledger":ledger,
        "no_graph_repair":True,"no_backfill":True,"no_resplitting":True,
        "nucleotide_identity_persisted":False,
        "genetic_distances_opened":False,"empirical_beta_host_opened":False,
        "direct_genetic_opening_authorized":False
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":output["status"],
                      "surviving_species":len(survivors)},sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
