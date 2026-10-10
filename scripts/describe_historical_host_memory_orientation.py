#!/usr/bin/env python3
"""Descriptive response-blind climate orientation, not a genetic effect test.

This is a retrospective ecological-predictor-only audit after the information
PASS; it cannot modify frozen estimands, predictors, panels, gates, or alpha.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import statistics
from collections import defaultdict
from pathlib import Path
import numpy as np

PRED_SHA="70b748d588d9ab10309ace080f85b8053c6e14aa96d7a84288724c1000272a5b"
ROLE_SHA="4168e7572d3378e7f1595ae562d53c81bb2ef165cbe50eb2ed8ec4e7406fc6ae"
INFO_RULE_SHA="2ec2594510021da653598e22a3d0f4b9e4cb226398aadbf342967ee532909415"


def sha256(path: Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()


def describe(rows, roles):
    source=defaultdict(lambda:defaultdict(list))
    for r in rows:
        species=str(r["species"])
        if species not in roles:
            raise RuntimeError("predictor species missing frozen role")
        for k in ("M_host","M_self"):
            value=float(r[k])
            if not np.isfinite(value):
                raise RuntimeError("nonfinite predictor edge")
            source[species][k].append(value)
    if set(source)!=set(roles):
        raise RuntimeError("frozen role species do not equal complete edge predictors")
    result={}
    for panel in ("development","confirmatory","all"):
        names=sorted(sp for sp in source if panel=="all" or roles[sp]==panel)
        summaries=[]
        for sp in names:
            host=np.asarray(source[sp]["M_host"],dtype=float)
            own=np.asarray(source[sp]["M_self"],dtype=float)
            if len(host)<2 or len(own)!=len(host):
                raise RuntimeError("bad species edge support")
            hstd=float(host.std(ddof=0))
            ostd=float(own.std(ddof=0))
            if hstd<=1e-12 or ostd<=1e-12:
                raise RuntimeError("species history invariant: not qualified")
            summaries.append((
                float(host.mean()),
                float(np.mean(host>0)),
                float(np.corrcoef(host,own)[0,1])
            ))
        arr=np.asarray(summaries,dtype=float)
        if np.any(~np.isfinite(arr)):
            raise RuntimeError("nonfinite predictor-only summary")
        result[panel]={
            "species":len(names),
            "species_mean_host_memory_positive":int(np.sum(arr[:,0]>0)),
            "species_mean_host_memory_negative":int(np.sum(arr[:,0]<0)),
            "median_species_mean_M_host":float(np.median(arr[:,0])),
            "median_fraction_positive_edges_per_species":float(np.median(arr[:,1])),
            "median_within_species_correlation_M_host_M_self":float(np.median(arr[:,2]))
        }
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--edge-predictors",type=Path,required=True)
    p.add_argument("--roles",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if sha256(a.edge_predictors)!=PRED_SHA or sha256(a.roles)!=ROLE_SHA:
        raise RuntimeError("frozen source SHA drift")
    with a.roles.open(newline="",encoding="utf-8") as f:
        roles={row["species"]:row["panel"] for row in csv.DictReader(f)}
    if len(roles)!=641 or sum(x=="development" for x in roles.values())!=320 or sum(x=="confirmatory" for x in roles.values())!=321:
        raise RuntimeError("frozen role panel sizes changed")
    with a.edge_predictors.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=71520:
        raise RuntimeError("frozen edge predictor support changed")
    results=describe(rows,roles)
    out={
        "schema":"ttf_historical_host_memory_response_blind_orientation_v0.1",
        "status":"EXPLORATORY_DESCRIPTIVE_PREDICTOR_ONLY",
        "source_sha256":{"edge_predictors":PRED_SHA,"roles":ROLE_SHA},
        "data":results,
        "n_genetic_edges":len(rows),
        "semantic_note":"Positive M_host means LGM geography is less climatically compatible with the *present* host-resource climate cloud than 0-BP geography. This does not prove ancient occupancy or host use.",
        "scope_note":"Within-species correlation here is ecological predictor co-alignment only. The median of species correlations is not a partial regression beta and is not evidence of genetic memory.",
        "no_change_to_frozen_hypothesis":True,
        "genetic_response_opened":False
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"data":results,"status":out["status"]},sort_keys=True))


if __name__=="__main__":
    main()
