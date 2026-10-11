#!/usr/bin/env python3
"""Audit post-mask sampling representativeness WITHOUT genetic response opening.

Input includes a private canonical-validity ledger, but output contains only
aggregate counts and source digests; no species names, nucleotide masks,
sequence identities, distances or genetic effect sizes are serialized.

This post-mask audit is exploratory and has NO authority to change original
259-species eligibility, effect estimand or stopping thresholds.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

MASK_SHA="b8b80e3edeae07d0a61606ec800d4a0c61296ebb06b379df3a1140ca3ecfa3d5"
CANDIDATE_SHA="c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b"
DIAGNOSTIC_SHA="bee840bf6ae207cab5aef8da4154eae7743f5ec3bbbaf23d7fe5a3dedc50a3e6"
ROLES_SHA="4168e7572d3378e7f1595ae562d53c81bb2ef165cbe50eb2ed8ec4e7406fc6ae"
ARCHIVE_SHA="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
EDGE_BINS=((0,50),(50,100),(100,200),(200,500),(500,None))


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda:file.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def finite_diagnostic(row,key)->float:
    x=float(row[key])
    if not math.isfinite(x):
        raise RuntimeError("nonfinite response-blind predictor diagnostic")
    return x


def aggregate(mask, candidates, diagnostics, roles):
    by_mask={r["species"]:r for r in mask["rows"]}
    by_candidates={r["species"]:r for r in candidates}
    by_diagnostics={r["species"]:r for r in diagnostics}
    by_roles={r["species"]:r for r in roles}

    if (mask.get("schema")!="ttf_historical_host_memory_independent_mask_replay_v0.1"
        or mask.get("source_archive_sha256")!=ARCHIVE_SHA
        or mask.get("genetic_response_not_opened") is not True
        or mask.get("genetic_distance_not_computed") is not True
        or mask.get("nucleotide_base_identity_not_persisted") is not True):
        raise RuntimeError("private mask source not proven response blind")
    if (len(mask["rows"])!=321 or len(by_mask)!=321
        or len(candidates)!=642 or len(by_candidates)!=642
        or len(diagnostics)!=642 or len(by_diagnostics)!=642
        or len(roles)!=641 or len(by_roles)!=641):
        raise RuntimeError("original source count/duplicate drift")
    confirm={sp for sp,r in by_roles.items() if r["panel"]=="confirmatory"}
    if len(confirm)!=321 or set(by_mask)!=confirm:
        raise RuntimeError("private mask does not equal original confirmatory panel")
    for sp,m in by_mask.items():
        if (sp not in by_candidates or by_diagnostics[sp]["status"]!="complete"
            or m["survives"] is not (m["status"]=="PASS_MASK")):
            raise RuntimeError("unqualified or inconsistent original mask species")

    passed=[sp for sp in by_mask if by_mask[sp]["survives"]]
    failed=[sp for sp in by_mask if not by_mask[sp]["survives"]]
    n_all=len(by_mask)
    n_pass=len(passed)
    n_fail=len(failed)
    original_edges=sum(int(by_candidates[sp]["edges"]) for sp in confirm)
    retained_edges=sum(int(by_candidates[sp]["edges"]) for sp in passed)
    invalid_edges=sum(int(m.get("invalid_edges",0)) for m in by_mask.values())
    if (n_pass,n_fail,original_edges,retained_edges,invalid_edges)!=(259,62,41592,20362,521):
        raise RuntimeError("mask population or frozen graph-edge accounting drift")

    out={
      "schema":"ttf_historical_host_memory_mask_attrition_representativeness_v0.1",
      "status":"EXPLORATORY_POST_MASK_REPRESENTATIVENESS_AUDIT",
      "mask_pass_species":n_pass,
      "mask_fail_species":n_fail,
      "confirmatory_species":n_all,
      "mask_pass_species_fraction":n_pass/n_all,
      "original_confirmatory_edges":original_edges,
      "survivor_confirmatory_edges":retained_edges,
      "surviving_original_edge_fraction":retained_edges/original_edges,
      "failed_species_original_edge_count":original_edges-retained_edges,
      "invalid_original_edges_with_finite_masks":invalid_edges,
      "mask_fail_reasons":dict(sorted(Counter(m["reason_code"] for m in mask["rows"] if not m["survives"]).items())),
    }
    def med_for(species,func):
        return float(statistics.median(func(sp) for sp in species))
    out["group_medians"]={
      "pass":{
        "frozen_edge_count":med_for(passed,lambda sp:int(by_candidates[sp]["edges"])),
        "frozen_locality_count":med_for(passed,lambda sp:int(by_candidates[sp]["n_localities"])),
        "host_history_unique_fraction":med_for(passed,lambda sp:finite_diagnostic(by_diagnostics[sp],"unique_fraction_M_host")),
        "predictor_condition_number":med_for(passed,lambda sp:finite_diagnostic(by_diagnostics[sp],"predictor_condition_number"))
      },
      "fail":{
        "frozen_edge_count":med_for(failed,lambda sp:int(by_candidates[sp]["edges"])),
        "frozen_locality_count":med_for(failed,lambda sp:int(by_candidates[sp]["n_localities"])),
        "host_history_unique_fraction":med_for(failed,lambda sp:finite_diagnostic(by_diagnostics[sp],"unique_fraction_M_host")),
        "predictor_condition_number":med_for(failed,lambda sp:finite_diagnostic(by_diagnostics[sp],"predictor_condition_number"))
      }
    }
    out["original_graph_edge_count_bins"]=[]
    for lo,hi in EDGE_BINS:
        elig=[sp for sp in confirm if int(by_candidates[sp]["edges"])>=lo and
              (hi is None or int(by_candidates[sp]["edges"])<hi)]
        n=sum(by_mask[sp]["survives"] for sp in elig)
        out["original_graph_edge_count_bins"].append({
          "minimum_edges_inclusive":lo,
          "maximum_edges_exclusive":hi,
          "original_species":len(elig),
          "mask_survivor_species":n,
          "mask_survival_rate":n/len(elig) if elig else None,
        })
    family=defaultdict(lambda:[0,0])
    for sp in confirm:
        category=by_candidates[sp]["family"]
        if not category:
            raise RuntimeError("missing frozen taxonomic family")
        family[category][0 if by_mask[sp]["survives"] else 1]+=1
    # Families with >=10 original panel species only, no per-species mask leakage.
    rows=[]
    for name,(n_pass,n_fail) in family.items():
        n=n_pass+n_fail
        if n>=10:
            rows.append({"family":name,"original_species":n,
                         "mask_survivor_species":n_pass,"mask_fail_species":n_fail,
                         "mask_survival_rate":n_pass/n})
    out["families_with_at_least_10_original_species"]=sorted(rows,key=lambda r:(-r["original_species"],r["family"]))
    out["families_total_original"]=len(family)
    out["source_sha256"]={"mask_ledger":MASK_SHA,"candidate_table":CANDIDATE_SHA,
                          "species_diagnostics":DIAGNOSTIC_SHA,"roles":ROLES_SHA}
    out["interpretation"]={
      "post_mask_descriptive_only":True,
      "null_significance_or_comparative_genetic_effect_computed":False,
      "mask_failure_is_not_an_ecological_absence":True,
      "families_are_not_an_authorized_subgroup_analysis":True,
      "original_confirmatory_genetic_estimand_is_frozen_to_259_survivors":True,
      "no_extra_exclusions_or_backfill":True,
      "observed_genetic_distance_or_effect_unopened":True,
      "warning":"Longer denser original graphs require all edges to have joint canonical support. Surviving species and edges are not an exchangeable random sample of the original 321-species panel."
    }
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--private-mask-ledger",type=Path,required=True)
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--species-diagnostics",type=Path,required=True)
    ap.add_argument("--roles",type=Path,required=True)
    ap.add_argument("--output-aggregate",type=Path,required=True)
    a=ap.parse_args()
    inputs=((a.private_mask_ledger,MASK_SHA),(a.candidates,CANDIDATE_SHA),
            (a.species_diagnostics,DIAGNOSTIC_SHA),(a.roles,ROLES_SHA))
    for p,h in inputs:
        if sha256(p)!=h:
            raise RuntimeError("frozen input SHA mismatch: "+p.name)
    with a.private_mask_ledger.open(encoding="utf-8") as f:
        mask=json.load(f)
    def load(path):
        with path.open(newline="",encoding="utf-8") as f:
            return list(csv.DictReader(f))
    result=aggregate(mask,load(a.candidates),load(a.species_diagnostics),load(a.roles))
    a.output_aggregate.parent.mkdir(parents=True,exist_ok=True)
    a.output_aggregate.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"pass_species":result["mask_pass_species"],
                      "pass_edge_fraction":result["surviving_original_edge_fraction"]},sort_keys=True))


if __name__=="__main__":
    main()
