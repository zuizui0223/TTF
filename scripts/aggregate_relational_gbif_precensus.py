#!/usr/bin/env python3
"""Aggregate the Relational-TTF GBIF availability pre-census.

This is a response-blind availability gate only. It cannot authorize the final
niche metric, environmental PCA, genetic response opening, or beta_R test.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from collections import Counter
from pathlib import Path


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()


def truthy(value: object) -> bool:
    return str(value).strip().casefold() in {"1","true","yes","y"}


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--counts", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args=ap.parse_args()

    candidates=list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    counts=list(csv.DictReader(args.counts.open(encoding="utf-8")))
    if not candidates or not counts:
        raise ValueError("candidate and count tables must be non-empty")
    if len({r["species"] for r in candidates}) != len(candidates):
        raise ValueError("duplicate candidate species")
    if len({r["species"] for r in counts}) != len(counts):
        raise ValueError("duplicate count species")

    cmap={r["species"]:r for r in candidates}
    qmap={r["species"]:r for r in counts}
    if set(cmap) != set(qmap):
        missing_counts=sorted(set(cmap)-set(qmap))
        extra_counts=sorted(set(qmap)-set(cmap))
        raise RuntimeError(
            f"species-set mismatch: missing_counts={len(missing_counts)}, "
            f"extra_counts={len(extra_counts)}"
        )

    rule=json.loads(args.rule.read_text(encoding="utf-8"))
    minimum=int(rule["external_occurrence_contract"]["minimum_usable_occurrences_candidate"])
    status_counts=Counter(r.get("status","") for r in counts)
    request_errors=int(status_counts.get("REQUEST_ERROR",0))
    accepted=[
        name for name,row in qmap.items()
        if truthy(row.get("match_accepted")) and int(row.get("occurrence_count_2010_endyear") or 0)>=minimum
    ]

    order_counts=Counter(cmap[name]["order"] for name in accepted)
    n=len(accepted)
    largest_order_count=max(order_counts.values(), default=0)
    largest_order_fraction=(largest_order_count/n) if n else None
    orders_ge_5=sum(1 for count in order_counts.values() if n and count/n>=0.05)

    if request_errors:
        decision="INCOMPLETE_PROVIDER_PRECENSUS"
    elif n < 500:
        decision="NOT_EVALUABLE_GBIF_COUNT_FLOOR"
    else:
        decision="PASS_TO_FULL_OCCURRENCE_ACQUISITION"

    payload={
        "schema":"ttf_relational_gbif_precensus_aggregate_v0.1",
        "status":decision,
        "candidate_species":len(candidates),
        "gbif_rows":len(counts),
        "minimum_provider_count_candidate":minimum,
        "provider_status_counts":dict(sorted(status_counts.items())),
        "species_passing_provider_count_gate":n,
        "fraction_passing_provider_count_gate":n/len(candidates),
        "taxonomic_summary_among_count_passers":{
            "orders":len(order_counts),
            "top_orders":dict(order_counts.most_common(12)),
            "largest_order_fraction":largest_order_fraction,
            "orders_with_fraction_ge_0_05":orders_ge_5,
            "breadth_guardrail_evaluated":False,
            "reason":"The frozen breadth gate applies after full occurrence filtering/thinning and environmental admissibility, not at the raw GBIF count triage stage."
        },
        "inputs":{
            "candidate_csv_sha256":sha256_path(args.candidates),
            "counts_csv_sha256":sha256_path(args.counts),
            "rule_sha256":sha256_path(args.rule)
        },
        "claim_boundary":(
            "PASS_TO_FULL_OCCURRENCE_ACQUISITION means only that >=500 fresh species have "
            "at least the candidate number of recent georeferenced GBIF records before exact "
            "duplicate removal, spatial thinning, raster extraction, environmental missingness, "
            "or niche-overlap construction. It does not authorize Schoener's D or genetic opening."
        ),
        "response_firewall":{
            "sequence_identity_opened":False,
            "pairwise_genetic_distances_opened":False,
            "T_st_computed":False,
            "environmental_overlap_computed":False,
            "ecology_genetics_association_computed":False
        }
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":decision,
        "count_passers":n,
        "request_errors":request_errors,
        "largest_order_fraction_raw_count_gate":largest_order_fraction
    },sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
