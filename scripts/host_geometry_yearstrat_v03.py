"""Strict, checkpoint-verified loader for separate v0.2 host year-bin observation data.

Modern occurrence support only: no genetic response, no historical occupancy.
Never merge with the earlier deep-offset v0.1 host observation sample.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
from reference_helpers_min import normalize_id

EXPECTED_SUMMARY_SHA256="cffb1144a61151c9da78737b1c7a84034398fb50d03ea5ecb46642b317d93634"
EXPECTED_MANIFEST="a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044"

def canonical(value):
    return json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(",",":"),allow_nan=False).encode("utf-8")

def load_host_snapshot(directory):
    directory=Path(directory)
    summary_bytes=(directory/"summary.json").read_bytes()
    if hashlib.sha256(summary_bytes).hexdigest()!=EXPECTED_SUMMARY_SHA256:
        raise ValueError("v0.2 plant observation source summary checksum mismatch")
    summary=json.loads(summary_bytes)
    if (summary["status"]!="ALL_TAXA_ATTEMPTED"
        or summary["attempted_taxa"]!=518 or summary["pending_taxa"]!=0
        or summary["quality_records"]!=296697 or summary["old_v01_result_mixed"] is not False
        or summary["native_range_filter_applied"] is not False
        or summary["climate_or_genetic_response_used"] is not False):
        raise ValueError("v0.2 host observation source semantic drift")
    plan=json.loads((directory/"source_bound_plan.json").read_text())
    binding=plan["binding"]
    if (binding!=summary["binding"] or not binding["not_equivalent_to_v01"]
        or binding["manifest_sha256"]!=EXPECTED_MANIFEST
        or len(plan["taxa"])!=518):
        raise ValueError("host v0.2 acquisition binding mismatch")
    expected={normalize_id(p["wcvp_id"]) for p in plan["taxa"]}
    if len(expected)!=518:raise ValueError("duplicate WCVP IDs")
    completed={};plants={};review=set();failures=set()
    checkpoints=list((directory/"checkpoints").glob("*.json"))
    if len(checkpoints)!=518:raise ValueError("not all v0.2 checkpoints present")
    for file in sorted(checkpoints):
        q=json.loads(file.read_text())
        result=q.get("result")
        if (q.get("schema")!="ttf_host_observation_checkpoint_v0.1"
            or q.get("binding")!=binding
            or not isinstance(result,dict)):
            raise ValueError("v0.2 checkpoint source or schema drift")
        taxon=normalize_id(result["wcvp_id"])
        if taxon not in expected or taxon!=file.stem or taxon in completed:
            raise ValueError("v0.2 accepted-plant taxonomic identity drift")
        if hashlib.sha256(canonical(result)).hexdigest()!=q.get("result_sha256"):
            raise ValueError("v0.2 checkpoint digest mismatch")
        state=result["status"]
        pts=result.get("points")
        if not isinstance(pts,list):raise ValueError("missing point list")
        if state=="OBSERVATIONS_ACQUIRED":
            if result["quality_records"]!=len(pts):raise ValueError("observation count mismatch")
            ids=[str(p["gbif_id"]) for p in pts]
            if len(ids)!=len(set(ids)):raise ValueError("duplicate plant GBIF IDs within plant taxon")
            plants[taxon]=pts
        elif state=="TAXONOMIC_RANK_REVIEW":
            if pts or result.get("quality_records") is not None:
                raise ValueError("unresolved rank is not zero observations")
            review.add(taxon)
        elif state=="TAXONOMY_NOT_EXACT":
            if pts or result.get("quality_records") is not None:
                raise ValueError("non-exact plant taxonomy is not observation absence")
            failures.add(taxon)
        else:raise ValueError("unexpected host retrieval state")
        completed[taxon]=result
    if (len(plants),len(review),len(failures),sum(map(len,plants.values())))!=(506,1,11,296697):
        raise ValueError("v0.2 final host source count drift")
    return completed,plants,review,failures,plan
