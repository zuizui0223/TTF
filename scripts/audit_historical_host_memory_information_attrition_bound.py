#!/usr/bin/env python3
"""Sharp worst-case host-predictor information bound under arbitrary mask attrition.

This response-blind combinatorial audit is valid for every *possible* subset
of at least the frozen 200 of 321 confirmatory species, provided every
survivor retains all its originally frozen graph edges. It never reads
nucleotide characters, sequence masks, empirical genetic distances or effects.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

EXPECTED_DIAGNOSTICS_SHA = "bee840bf6ae207cab5aef8da4154eae7743f5ec3bbbaf23d7fe5a3dedc50a3e6"
EXPECTED_ROLE_SHA = "4168e7572d3378e7f1595ae562d53c81bb2ef165cbe50eb2ed8ec4e7406fc6ae"
EXPECTED_RULE_SHA = "2ec2594510021da653598e22a3d0f4b9e4cb226398aadbf342967ee532909415"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def worst_case_bound(
    diagnostic: list[dict[str,str]], roles: list[dict[str,str]], rule: dict
) -> dict:
    """Exact extremal subset bound, not randomized mask attrition modelling."""
    assignments = {r["species"]:r["panel"] for r in roles}
    if (len(assignments) != 641 or len(roles) != 641
            or sum(r["panel"] == "confirmatory" for r in roles) != 321
            or sum(r["panel"] == "development" for r in roles) != 320):
        raise RuntimeError("frozen 320 / 321 species assignment drift")
    by_species = {r["species"]:r for r in diagnostic}
    if len(diagnostic)!=642 or len(by_species)!=642:
        raise RuntimeError("frozen 642-species diagnostics source drift")
    if not set(assignments).issubset(by_species):
        raise RuntimeError("frozen confirmatory species missing information diagnostics")
    retained = sorted(s for s,panel in assignments.items() if panel == "confirmatory")
    fractions = []
    condition = []
    for s in retained:
        row = by_species[s]
        if row["status"] != "complete":
            raise RuntimeError(f"frozen confirmatory species has noncomplete diagnostics: {s}")
        try:
            a=float(row["unique_fraction_M_host"])
            b=float(row["predictor_condition_number"])
        except (TypeError, ValueError) as exc:
            raise RuntimeError(f"nonparseable original frozen diagnostic: {s}") from exc
        if (not math.isfinite(a) or not 0<=a<=1
                or not math.isfinite(b) or b<1):
            raise RuntimeError(f"invalid frozen predictor diagnostic: {s}")
        fractions.append(a)
        condition.append(b)
    conf=rule["predictor_information"]
    k=int(rule["post_mask_requalification"]["minimum_confirmatory_survivor_species"])
    floor_fraction=float(conf["minimum_species_fraction_with_unique_fraction_at_least_0.05"])
    floor_median=float(conf["minimum_panel_median_unique_fraction"])
    max_condition=float(conf["maximum_standardized_predictor_condition_number"])
    floor_cond_fraction=float(conf["minimum_species_fraction_passing_condition_number"])
    if (k!=200 or floor_fraction!=0.70 or floor_median!=0.10
            or max_condition!=30.0 or floor_cond_fraction!=0.90):
        raise RuntimeError("frozen mask-survivor information thresholds drift")
    # Among all k-species subsets, the minimum possible median is the
    # median of the globally smallest k values, by order-statistic monotonicity.
    ordered_unique=sorted(fractions)
    worst_median=float(statistics.median(ordered_unique[:k]))
    # The largest attainable number of failing species in any k-subset is
    # min(k, total number of bad species), yielding the sharp lower fraction.
    bad_unique=sum(v<0.05 for v in fractions)
    bad_condition=sum(v>max_condition for v in condition)
    worst_unique_fraction=(k-min(k,bad_unique))/k
    worst_condition_fraction=(k-min(k,bad_condition))/k
    # For n>=k these three extremal bounds can only improve, so k suffices.
    decisions={
        "median_unique_ge_0_10":worst_median>=floor_median,
        "fraction_unique_ge_0_70":worst_unique_fraction>=floor_fraction,
        "fraction_condition_ge_0_90":worst_condition_fraction>=floor_cond_fraction,
    }
    return {
        "confirmatory_species":len(retained),
        "minimum_mask_survivors":k,
        "confirmed_original_incomplete":0,
        "unique_fraction_below_0_05_in_full_confirmatory":bad_unique,
        "condition_number_above_30_in_full_confirmatory":bad_condition,
        "worst_case_at_200":{
            "minimum_possible_median_unique_fraction":worst_median,
            "minimum_possible_fraction_unique_ge_0_05":worst_unique_fraction,
            "minimum_possible_fraction_condition_le_30":worst_condition_fraction,
        },
        "threshold_checks":decisions,
        "all_possible_confirmatory_survivor_subsets_ge_200_pass_information":
            all(decisions.values()),
        "argument":"For every subset cardinality n>=200: its median is bounded below by the median of the n smallest original values (itself >= the 200-smallest median); for either fraction, no subset can contain more than the total original count of failing species, and (n-bad)/n increases in n.",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--diagnostics",type=Path,required=True)
    p.add_argument("--roles",type=Path,required=True)
    p.add_argument("--qualification-rule",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    checks=((a.diagnostics,EXPECTED_DIAGNOSTICS_SHA),
            (a.roles,EXPECTED_ROLE_SHA),
            (a.qualification_rule,EXPECTED_RULE_SHA))
    for path,digest in checks:
        if sha256(path)!=digest:
            raise RuntimeError(f"exact frozen information source SHA mismatch: {path}")
    with a.diagnostics.open(newline="",encoding="utf-8") as f:
        diag=list(csv.DictReader(f))
    with a.roles.open(newline="",encoding="utf-8") as f:
        roles=list(csv.DictReader(f))
    rule=json.loads(a.qualification_rule.read_text())
    bound=worst_case_bound(diag,roles,rule)
    result={
        "schema":"ttf_historical_host_memory_arbitrary_mask_attrition_information_bound_v0.1",
        "status":("PASS_INFORMATION_GUARANTEED_FOR_ANY_SURVIVOR_SET_GE_200"
                  if bound["all_possible_confirmatory_survivor_subsets_ge_200_pass_information"]
                  else "NO_WORST_CASE_GUARANTEE"),
        "input_sha256":{"diagnostics":EXPECTED_DIAGNOSTICS_SHA,
                        "roles":EXPECTED_ROLE_SHA,
                        "rule":EXPECTED_RULE_SHA},
        "bound":bound,
        "scope":"Combinatorial response-blind result conditional on >=200 exact original confirmatory species surviving all original frozen edges. Does not certify sequence-mask survival or synthetic power after mask.",
        "changes_to_frozen_mask_rule":False,
        "genetic_response_opened":False
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],
                      "worst_case_at_200":bound["worst_case_at_200"]},sort_keys=True))


if __name__=="__main__":
    main()
