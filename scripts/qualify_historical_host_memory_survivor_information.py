#!/usr/bin/env python3
"""Exact-confirmatory-survivor information requalification after mask-only opening.

Cannot open genetic nucleotide identity, genetic distances, empirical turnover,
or beta_host. No replacement species, edge rewiring, or adaptive thresholds.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

PRED_SHA="70b748d588d9ab10309ace080f85b8053c6e14aa96d7a84288724c1000272a5b"
ROLE_SHA="4168e7572d3378e7f1595ae562d53c81bb2ef165cbe50eb2ed8ec4e7406fc6ae"
DIAG_SHA="bee840bf6ae207cab5aef8da4154eae7743f5ec3bbbaf23d7fe5a3dedc50a3e6"
ARCH_SHA="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce61a5"
RULE_SHA="2ec2594510021da653598e22a3d0f4b9e4cb226398aadbf342967ee532909415"
EXPECTED_MASK_BINDINGS={
    "candidate_table_sha256":"c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b",
    "panel_roles_sha256":ROLE_SHA,
    "locality_geometry_sha256":"037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9",
    "edge_geometry_sha256":"ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9",
    "synthetic_pass_receipt_sha256":"0498f71d636274a7f2e654f42591545cf26cbd0b97e603902377ea2ef702495e",
    "mask_rule_git_blob_sha1":"1d5bec837157955cd2f989b2dd127d502d57318c",
}


def file_sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()


def qualify_survivors(mask: dict, role: list[dict], diagnostic: list[dict], rule: dict):
    confirm={r["species"] for r in role if r["panel"]=="confirmatory"}
    if len(confirm)!=321 or len({r["species"] for r in role})!=641:
        raise RuntimeError("frozen pre-mask panel assignments changed")
    if mask.get("schema")!="ttf_historical_host_memory_confirmatory_mask_result_v0.1":
        raise RuntimeError("unexpected confirmatory mask result schema")
    if mask.get("confirmatory_before_mask")!=321 or mask.get("source_archive_sha256")!=ARCH_SHA:
        raise RuntimeError("mask provenance or exact confirmatory sample differs")
    bindings=mask.get("source_bindings")
    if (not isinstance(bindings,dict)
            or any(bindings.get(k)!=v for k,v in EXPECTED_MASK_BINDINGS.items())
            or not isinstance(bindings.get("mask_rule_sha256"),str)
            or len(bindings["mask_rule_sha256"])!=64):
        raise RuntimeError("frozen mask rule or original provenance bindings drift")
    if (mask.get("nucleotide_identity_persisted") is not False
            or mask.get("genetic_distances_opened") is not False
            or mask.get("empirical_beta_host_opened") is not False
            or mask.get("no_graph_repair") is not True
            or mask.get("no_resplitting") is not True
            or mask.get("no_backfill") is not True
            or mask.get("direct_genetic_opening_authorized") is not False):
        raise RuntimeError("unexpected genetic opening or mask design modification")
    ledger=mask.get("ledger",[])
    if len(ledger)!=321 or {r["species"] for r in ledger}!=confirm:
        raise RuntimeError("mask ledger does not cover exact 321 confirmatory species")
    kept=set(mask.get("survivor_names",[]))
    if (kept!={r["species"] for r in ledger if r["survives"]}
            or len(kept)!=int(mask.get("surviving_species",-1))):
        raise RuntimeError("mask survivor identities/counts inconsistent")
    expected_mask_state = (
        "PASS_TO_EXACT_SURVIVOR_INFORMATION_AND_SYNTHETIC_REQUALIFICATION"
        if len(kept) >= 200 else
        "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_CHARACTER_SUPPORT"
    )
    if mask.get("status") != expected_mask_state:
        raise RuntimeError("frozen mask decision/count mismatch")
    for r in ledger:
        if r["status"] != ("PASS_MASK" if r["survives"] else "NOT_EVALUABLE_CHARACTER_SUPPORT"):
            raise RuntimeError("illegal masked species survival state")
        if r["survives"] and (r["invalid_edges"]!=0 or r["valid_edges"]!=r["frozen_edges"]):
            raise RuntimeError("mask survivor has an invalid frozen edge")
    d={r["species"]:r for r in diagnostic}
    if len(d)!=642 or len({r["species"] for r in diagnostic})!=642:
        raise RuntimeError("pre-mask species diagnostics source drift")
    if any(d[s]["status"]!="complete" for s in kept):
        raise RuntimeError("mask survivor must have previously complete predictor")
    if len(kept)<200:
        return {"decision":"NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_CHARACTER_SUPPORT",
                "confirmatory_before":321,"confirmatory_survivors":len(kept),
                "genetic_response_opened":False},kept

    conf=rule["predictor_information"]
    if (float(conf["minimum_species_fraction_with_unique_fraction_at_least_0.05"])!=0.70
        or float(conf["minimum_panel_median_unique_fraction"])!=0.10
        or float(conf["maximum_standardized_predictor_condition_number"])!=30
        or float(conf["minimum_species_fraction_passing_condition_number"])!=0.90
        or int(rule["post_mask_requalification"]["minimum_confirmatory_survivor_species"])!=200):
        raise RuntimeError("frozen post-mask information thresholds changed")
    vals=[float(d[s]["unique_fraction_M_host"]) for s in sorted(kept)]
    cond=[float(d[s]["predictor_condition_number"]) for s in sorted(kept)]
    if any(not math.isfinite(x) for x in vals+cond):
        raise RuntimeError("non-finite survivor information diagnostics")
    fraction=sum(v>=0.05 for v in vals)/len(vals)
    condition=sum(v<=30 for v in cond)/len(cond)
    median=statistics.median(vals)
    checks={
       "minimum_200_species":len(kept)>=200,
       "fraction_nonredundant_ge_0_70":fraction>=0.70,
       "fraction_condition_good_ge_0_90":condition>=0.90,
       "median_unique_ge_0_10":median>=0.10,
    }
    return {"decision":("PASS_TO_EXACT_SURVIVOR_SYNTHETIC_REQUALIFICATION" if all(checks.values())
                         else "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_SURVIVOR_INFORMATION"),
            "confirmatory_before":321,
            "confirmatory_survivors":len(kept),
            "fraction_unique_at_least_0_05":fraction,
            "fraction_condition_at_most_30":condition,
            "survivor_median_unique_fraction":median,
            "checks":checks,
            "survivor_synthetic_qualification_passed":False,
            "genetic_response_opened":False},kept


def main():
    p=argparse.ArgumentParser()
    for k in ("mask_result","roles","diagnostics","predictors","qualification_rule"):
        p.add_argument("--"+k.replace("_","-"),type=Path,required=True)
    p.add_argument("--output-summary",type=Path,required=True)
    p.add_argument("--output-survivor-predictors",type=Path,required=True)
    p.add_argument("--output-survivor-roles",type=Path,required=True)
    a=p.parse_args()
    for key,digest in {"roles":ROLE_SHA,"diagnostics":DIAG_SHA,
                       "predictors":PRED_SHA,"qualification_rule":RULE_SHA}.items():
        if file_sha(getattr(a,key))!=digest:
            raise RuntimeError("frozen predictor/role/rule source digest drift: "+key)
    mask=json.loads(a.mask_result.read_text())
    rule=json.loads(a.qualification_rule.read_text())
    with a.roles.open(newline="",encoding="utf-8") as f:roles=list(csv.DictReader(f))
    with a.diagnostics.open(newline="",encoding="utf-8") as f:diag=list(csv.DictReader(f))
    decision, kept=qualify_survivors(mask,roles,diag,rule)
    decision.update({
        "schema":"ttf_historical_host_memory_survivor_information_v0.1",
        "mask_result_sha256":file_sha(a.mask_result),
        "roles_sha256":ROLE_SHA,
        "diagnostics_sha256":DIAG_SHA,
        "predictor_sha256":PRED_SHA,
        "surviving_species_names":sorted(kept),
        "nucleotide_identity_opened":False,
        "post_IBD_turnover_opened":False,
        "empirical_beta_host_opened":False
    })
    a.output_summary.parent.mkdir(parents=True,exist_ok=True)
    a.output_summary.write_text(json.dumps(decision,indent=2,sort_keys=True)+"\n")
    if decision["decision"]!="PASS_TO_EXACT_SURVIVOR_SYNTHETIC_REQUALIFICATION":
        print(json.dumps({"decision":decision["decision"],"species":len(kept)}))
        return 0

    a.output_survivor_predictors.parent.mkdir(parents=True,exist_ok=True)
    edge_count=0
    with a.predictors.open(newline="",encoding="utf-8") as f, a.output_survivor_predictors.open("w",newline="",encoding="utf-8") as g:
        reader=csv.DictReader(f)
        writer=csv.DictWriter(g,fieldnames=reader.fieldnames,lineterminator="\n")
        writer.writeheader()
        for row in reader:
            if row["species"] in kept:
                writer.writerow(row);edge_count+=1
    with a.output_survivor_roles.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=["species","panel","split_key_sha256"],lineterminator="\n")
        writer.writeheader()
        writer.writerows(r for r in roles if r["species"] in kept)
    decision.update({
        "surviving_frozen_edges":edge_count,
        "survivor_predictor_sha256":file_sha(a.output_survivor_predictors),
        "survivor_roles_sha256":file_sha(a.output_survivor_roles)
    })
    a.output_summary.write_text(json.dumps(decision,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"decision":decision["decision"],"species":len(kept),"edges":edge_count}))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
