#!/usr/bin/env python3
"""Response-blind survivor-information gate for codistributed recurrence.

This step is authorized only by a frozen confirmatory character-mask result.
It uses unchanged survivor sampling geometry and taxonomy metadata only. It
never reads nucleotide identity, genetic distance, empirical T_st, or beta_G.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from run_codistributed_recurrence_fixed_tst_screen import (
    fast_symmetric_coverage,
    read_csv,
    reconstruct_geometries,
)
from ttf.codistributed_survivor import survivor_information_gate_dict


EXPECTED_MASK_SCHEMA="ttf_genetic_codistributed_recurrence_confirmatory_mask_result_v0.1"
EXPECTED_RULE_SCHEMA="ttf_genetic_codistributed_recurrence_confirmatory_mask_rule_v0.1"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def digest_labels(names) -> str:
    return hashlib.sha256("\n".join(sorted(map(str,names))).encode()).hexdigest()


def read_survivors(path: Path) -> dict[str,str]:
    rows=list(csv.DictReader(path.open(newline="",encoding="utf-8")))
    if not rows or set(rows[0])!={"species","role"}:
        raise RuntimeError("survivor CSV schema drift")
    out={}
    for row in rows:
        name=str(row["species"])
        role=str(row["role"])
        if role not in {"source","target"} or name in out:
            raise RuntimeError("invalid or duplicate survivor role row")
        out[name]=role
    return out


def authorize_mask(mask_path: Path, survivor_path: Path, rule_path: Path):
    mask=json.loads(mask_path.read_text())
    rule=json.loads(rule_path.read_text())
    if mask.get("schema")!=EXPECTED_MASK_SCHEMA:
        raise RuntimeError("unexpected confirmatory mask schema")
    if mask.get("status")!="PASS_TO_SURVIVOR_INFORMATION_GATE":
        raise RuntimeError("confirmatory mask does not authorize survivor information gate")
    if rule.get("schema")!=EXPECTED_RULE_SCHEMA:
        raise RuntimeError("unexpected confirmatory mask rule schema")
    if rule.get("status")!="FROZEN_BEFORE_CONFIRMATORY_CHARACTER_MASK_AND_BEFORE_ANY_CONFIRMATORY_GENETIC_RESPONSE":
        raise RuntimeError("confirmatory mask rule is not frozen")
    if any(bool(v) for v in mask["response_firewall"].values()):
        raise RuntimeError("confirmatory mask response firewall is open")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("confirmatory mask-rule response firewall is open")
    if sha256_path(survivor_path)!=mask["survivor_csv_sha256"]:
        raise RuntimeError("survivor CSV hash drift")
    survivors=read_survivors(survivor_path)
    sources=sorted(n for n,r in survivors.items() if r=="source")
    targets=sorted(n for n,r in survivors.items() if r=="target")
    if len(survivors)!=int(mask["surviving_species"]):
        raise RuntimeError("surviving species count drift")
    if len(sources)!=int(mask["surviving_sources"]) or len(targets)!=int(mask["surviving_targets"]):
        raise RuntimeError("survivor role count drift")
    if digest_labels(sources)!=mask["surviving_source_digest_sha256"]:
        raise RuntimeError("surviving source digest drift")
    if digest_labels(targets)!=mask["surviving_target_digest_sha256"]:
        raise RuntimeError("surviving target digest drift")
    return mask,rule,survivors


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--localities",type=Path,required=True)
    ap.add_argument("--edges",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--mask-result",type=Path,required=True)
    ap.add_argument("--survivors",type=Path,required=True)
    ap.add_argument("--mask-rule",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    mask,rule,survivor_roles=authorize_mask(
        args.mask_result,args.survivors,args.mask_rule
    )
    candidate_rows=read_csv(args.candidates)
    geometries,midpoint,roles,meta=reconstruct_geometries(
        read_csv(args.localities),
        read_csv(args.edges),
        candidate_rows,
    )
    missing=set(survivor_roles)-set(geometries)
    if missing:
        raise RuntimeError(f"survivors missing frozen confirmatory geometry: {sorted(missing)}")
    for name,role in survivor_roles.items():
        if roles.get(name)!=role:
            raise RuntimeError(f"survivor role drift for {name}")

    geometries={n:geometries[n] for n in survivor_roles}
    midpoint={n:midpoint[n] for n in survivor_roles}
    sources=tuple(sorted(n for n,r in survivor_roles.items() if r=="source"))
    targets=tuple(sorted(n for n,r in survivor_roles.items() if r=="target"))

    status="NOT_EVALUABLE_CODISTRIBUTED_RECURRENCE_SURVIVOR_INFORMATION"
    reason=None
    metrics=None
    g_summary=None
    try:
        if len(sources)<2 or len(targets)<2:
            raise ValueError("fewer than two surviving source or target species")
        pairs=[(s,t) for s in sources for t in targets]
        radius=500.0
        G=fast_symmetric_coverage(midpoint,sources,targets,radius)
        center={name:np.mean(midpoint[name],axis=0) for name in survivor_roles}
        controls=np.empty((len(pairs),4),dtype=float)
        for i,(s,t) in enumerate(pairs):
            controls[i,0]=np.log1p(float(np.linalg.norm(center[s]-center[t]))/radius)
            controls[i,1]=abs(np.log(geometries[s].n_edges/geometries[t].n_edges))
            controls[i,2]=abs(np.log(geometries[s].n_localities/geometries[t].n_localities))
            controls[i,3]=float(meta[s]["order"]==meta[t]["order"])
        src_index=np.repeat(np.arange(len(sources),dtype=np.int64),len(targets))
        tgt_index=np.tile(np.arange(len(targets),dtype=np.int64),len(sources))

        gate=rule["response_blind_survivor_information_gate"]
        metrics=survivor_information_gate_dict(
            src_index,tgt_index,G,controls,
            total_unique_fraction_minimum=float(gate["total_unique_fraction_minimum"]),
            maximum_endpoint_signal_share_ceiling=float(gate["maximum_source_signal_share_ceiling"]),
        )
        if float(gate["maximum_target_signal_share_ceiling"]) != float(gate["maximum_source_signal_share_ceiling"]):
            raise RuntimeError("source/target concentration ceilings diverged")
        status=(
            "PASS_TO_SURVIVOR_FIXED_TST_REQUALIFICATION"
            if metrics["overall_pass"]
            else "NOT_EVALUABLE_CODISTRIBUTED_RECURRENCE_SURVIVOR_INFORMATION"
        )
        g_summary={
            "mean":float(np.mean(G)),
            "median":float(np.median(G)),
            "q90":float(np.quantile(G,0.9)),
            "max":float(np.max(G)),
        }
    except ValueError as exc:
        reason=str(exc)

    payload={
        "schema":"ttf_genetic_codistributed_recurrence_survivor_information_v0.1",
        "status":status,
        "mask_result_sha256":sha256_path(args.mask_result),
        "survivor_csv_sha256":sha256_path(args.survivors),
        "mask_rule_sha256":sha256_path(args.mask_rule),
        "surviving_species":len(survivor_roles),
        "surviving_sources":len(sources),
        "surviving_targets":len(targets),
        "survivor_dyads":len(sources)*len(targets),
        "G_summary":g_summary,
        "information_gate":metrics,
        "not_evaluable_reason":reason,
        "scientific_interpretation":"response-blind evaluability only; failure is not a biological null",
        "response_firewall":{
            "confirmatory_nucleotide_identity_opened":False,
            "confirmatory_pairwise_genetic_distances_opened":False,
            "confirmatory_T_st_opened":False,
            "confirmatory_beta_G_opened":False,
        },
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "surviving_sources":len(sources),
        "surviving_targets":len(targets),
        "information_gate":metrics,
        "reason":reason,
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
