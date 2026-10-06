#!/usr/bin/env python3
"""Freeze the response-blind narrow-host panel for historical host connectivity."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

EXPECTED_CANDIDATE_SHA="1c8b35a7a005af639bbabc4d730fbe723df134c8371d04019a564740b095dbef"
PANEL_NAMESPACE="historical-host-connectivity-panel-v0.1"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def digest_names(names: list[str]) -> str:
    return hashlib.sha256(("\n".join(names)+"\n").encode()).hexdigest()


def rank_key(species: str, source_sha256: str) -> tuple[str,str]:
    h=hashlib.sha256(
        f"{PANEL_NAMESPACE}|{source_sha256}|{species}".encode()
    ).hexdigest()
    return h,species


def read_csv(path: Path) -> list[dict[str,str]]:
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def build_panel(
    candidates: list[dict[str,str]],
    interactions: list[dict[str,str]],
    native_rows: list[dict[str,str]],
    *,
    source_sha256: str,
    minimum_species: int=100,
) -> tuple[list[dict[str,object]], dict]:
    candidate_species=[str(r["species"]).strip() for r in candidates]
    if len(candidate_species)!=600 or len(set(candidate_species))!=600:
        raise RuntimeError("candidate panel must contain exactly 600 unique species")

    host_by_insect: dict[str,dict[str,str]]=defaultdict(dict)
    for row in interactions:
        insect=str(row.get("insect_species","")).strip()
        hid=str(row.get("accepted_plant_name_id","")).strip()
        hname=str(row.get("accepted_name","")).strip()
        if insect in set(candidate_species) and hid and hname:
            host_by_insect[insect][hid]=hname

    native_ids={
        str(r.get("accepted_plant_name_id","")).strip()
        for r in native_rows
        if str(r.get("accepted_plant_name_id","")).strip()
    }

    retained=[]
    excluded_counts=defaultdict(int)
    for species in candidate_species:
        hosts=host_by_insect.get(species,{})
        n=len(hosts)
        if n<1:
            excluded_counts["no_wcvp_resolved_species_host"]+=1
            continue
        if n>5:
            excluded_counts["host_breadth_gt_5"]+=1
            continue
        if not set(hosts).issubset(native_ids):
            excluded_counts["host_without_primary_native_wgsrpd3"]+=1
            continue
        retained.append({
            "species":species,
            "n_hosts":n,
            "accepted_host_ids":";".join(sorted(hosts)),
            "accepted_host_names":";".join(hosts[k] for k in sorted(hosts)),
        })

    retained.sort(key=lambda r:rank_key(str(r["species"]),source_sha256))
    if len(retained)<minimum_species:
        return [],{
            "status":"NOT_EVALUABLE_HISTORICAL_HOST_CONNECTIVITY_HOST_DOMAIN_TOO_SMALL",
            "retained_species":len(retained),
            "minimum_species":minimum_species,
            "excluded_counts":dict(excluded_counts),
        }

    for i,row in enumerate(retained):
        row["host_panel_rank"]=i+1

    summary={
        "status":"PASS_TO_HOST_OCCURRENCE_FEASIBILITY",
        "retained_species":len(retained),
        "species_digest_sha256":digest_names(sorted(str(r["species"]) for r in retained)),
        "host_breadth_counts":{
            str(k):sum(int(r["n_hosts"])==k for r in retained)
            for k in range(1,6)
        },
        "excluded_counts":dict(excluded_counts),
    }
    return retained,summary


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--interactions",type=Path,required=True)
    ap.add_argument("--native-distributions",type=Path,required=True)
    ap.add_argument("--source-sha256",required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-receipt",type=Path,required=True)
    args=ap.parse_args()

    if sha256_path(args.candidates)!=EXPECTED_CANDIDATE_SHA:
        raise RuntimeError("historical-host candidate CSV SHA drift")

    panel,summary=build_panel(
        read_csv(args.candidates),
        read_csv(args.interactions),
        read_csv(args.native_distributions),
        source_sha256=args.source_sha256,
    )
    payload={
        "schema":"ttf_genetic_historical_host_connectivity_panel_v0.1",
        **summary,
        "inputs_sha256":{
            "candidates":sha256_path(args.candidates),
            "interactions":sha256_path(args.interactions),
            "native_distributions":sha256_path(args.native_distributions),
        },
        "selection":{
            "host_breadth_min":1,
            "host_breadth_max":5,
            "all_hosts_require_primary_native_wgsrpd3":True,
            "host_panel_namespace":PANEL_NAMESPACE,
            "development_confirmatory_split_deferred_until_after_occurrence_gate":True,
            "no_backfill":True,
        },
        "response_firewall":{
            "fresh_sequence_identity_opened":False,
            "fresh_pairwise_genetic_distances_opened":False,
            "fresh_post_ibd_turnover_computed":False,
            "fresh_beta_host_LGM_computed":False,
        },
    }
    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_csv.open("w",newline="",encoding="utf-8") as f:
        fields=["host_panel_rank","species","n_hosts","accepted_host_ids","accepted_host_names"]
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        for row in panel:
            w.writerow({k:row[k] for k in fields})
    payload["panel_csv_sha256"]=sha256_path(args.output_csv)
    args.output_receipt.parent.mkdir(parents=True,exist_ok=True)
    args.output_receipt.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "retained_species":payload.get("retained_species"),
        "species_digest_sha256":payload.get("species_digest_sha256"),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
