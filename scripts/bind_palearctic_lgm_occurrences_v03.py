#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-csv",type=Path,required=True)
    ap.add_argument("--input-summary",type=Path,required=True)
    ap.add_argument("--panel-metadata",type=Path,required=True)
    ap.add_argument("--external-rule",type=Path,required=True)
    ap.add_argument("--binding",type=Path,required=True)
    ap.add_argument("--output-csv",type=Path,required=True)
    ap.add_argument("--output-summary",type=Path,required=True)
    args=ap.parse_args()

    inp=json.loads(args.input_summary.read_text())
    if inp.get("schema")!="ttf_genetic_palearctic_lgm_occurrence_result_v0.2":
        raise RuntimeError("unexpected bound GBIF summary schema")
    panel=json.loads(args.panel_metadata.read_text())
    if panel.get("schema")!="ttf_genetic_palearctic_lgm_panel_metadata_v0.3":
        raise RuntimeError("unexpected v0.3 panel metadata")
    rule=json.loads(args.external_rule.read_text())
    if rule.get("schema")!="ttf_genetic_palearctic_lgm_external_data_rule_v0.3":
        raise RuntimeError("unexpected v0.3 external rule")
    binding=json.loads(args.binding.read_text())
    if binding.get("schema")!="ttf_genetic_palearctic_lgm_v03_occurrence_producer_binding_v0.1":
        raise RuntimeError("unexpected occurrence producer binding")
    if any(bool(v) for v in binding["response_firewall"].values()):
        raise RuntimeError("v0.3 occurrence binding response firewall is open")

    meta={str(row["species"]):dict(row) for row in panel["rows"]}
    if len(meta)!=23:
        raise RuntimeError("v0.3 panel metadata must contain 23 species")
    sources={sp for sp,row in meta.items() if row["role"]=="source"}
    targets={sp for sp,row in meta.items() if row["role"]=="target"}
    if len(sources)!=11 or len(targets)!=12 or sources & targets:
        raise RuntimeError("v0.3 disjoint role geometry drift")

    request_errors=list(inp.get("request_error_species") or [])
    admissible=sorted(map(str,inp.get("admissible_species_list") or []))
    if not set(admissible).issubset(meta):
        raise RuntimeError("bound GBIF admissible species are outside v0.3 membership")

    rows=list(csv.DictReader(args.input_csv.open(encoding="utf-8")))
    for row in rows:
        sp=str(row["species"])
        if sp not in meta:
            raise RuntimeError(f"bound occurrence row outside v0.3 membership: {sp}")
        row["role"]=str(meta[sp]["role"])

    if request_errors:
        status="INCOMPLETE_TECHNICAL_V03_GBIF_REQUEST_ERRORS"
        retained=[]
    else:
        source_kept=sorted(sp for sp in admissible if sp in sources)
        target_kept=sorted(sp for sp in admissible if sp in targets)
        req=rule["missingness"]
        passed=(
            len(admissible)>=int(req["minimum_total_species"])
            and len(source_kept)>=int(req["minimum_source_clusters"])
            and len(target_kept)>=int(req["minimum_target_clusters"])
        )
        status=(
            "PASS_TO_V03_CURRENT_LGM_CLIMATE_EXTRACTION"
            if passed else str(req["failure"])
        )
        retained=rows if passed else []
    if request_errors:
        source_kept=[]; target_kept=[]

    retained.sort(key=lambda x:(x["species"],int(x["priority_rank"]),int(x["source_key"])))
    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    fields=["species","role","source_key","latitude","longitude","priority_rank","priority_sha256"]
    with args.output_csv.open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=fields)
        writer.writeheader(); writer.writerows(retained)

    payload={
        "schema":"ttf_genetic_palearctic_lgm_occurrence_result_v0.3",
        "status":status,
        "formal_species":23,
        "admissible_species":len(admissible) if not request_errors else 0,
        "source_clusters":len(source_kept),
        "target_clusters":len(target_kept),
        "directed_dyads":len(source_kept)*len(target_kept),
        "admissible_species_list":admissible if not request_errors else [],
        "source_species":source_kept,
        "target_species":target_kept,
        "failed_occurrence_species":inp.get("failed_occurrence_species",[]),
        "rejected_taxon_match_species":inp.get("rejected_taxon_match_species",[]),
        "request_error_species":request_errors,
        "retained_counts":inp.get("retained_counts",{}),
        "retained_occurrence_rows":len(retained),
        "inputs":{
            "bound_v02_summary_sha256":sha256_path(args.input_summary),
            "bound_v02_csv_sha256":sha256_path(args.input_csv),
            "panel_metadata_sha256":sha256_path(args.panel_metadata),
            "external_rule_sha256":sha256_path(args.external_rule),
            "producer_binding_sha256":sha256_path(args.binding),
        },
        "output_csv_sha256":sha256_path(args.output_csv),
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_v03_T_st_computed":False,
            "beta_LGM_computed":False,
        },
        "next_step":(
            "Build v0.3 current/LGM climate inputs."
            if status=="PASS_TO_V03_CURRENT_LGM_CLIMATE_EXTRACTION"
            else (
                "Technical retry only; do not change membership or ecological thresholds."
                if request_errors else
                "STOP. Do not build R_LGM or subgroup genetic response."
            )
        ),
    }
    args.output_summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "admissible_species":payload["admissible_species"],
        "source_clusters":len(source_kept),
        "target_clusters":len(target_kept),
        "directed_dyads":payload["directed_dyads"],
        "request_errors":len(request_errors),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
