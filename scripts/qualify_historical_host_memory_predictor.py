#!/usr/bin/env python3
"""Response-blind v0.2 information qualification for frozen 642 Lepidoptera.

No sequence identity, pairwise genetic distance, post-IBD response, or
historical-host genetic effect is read here. Does not perform the separately
frozen synthetic null/power qualification.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

SOURCE_SHA256 = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce61a5"
CANDIDATE_SHA256 = "c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b"
SPLIT_NAMESPACE = "historical-host-memory-split-v0.1"
ALLOWED_STATUSES = {
    "complete", "insufficient_host_cloud", "insufficient_self_cloud",
    "nonfinite_edge_climate", "rank_deficient_predictor",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_and_qualify(rows: list[dict[str, str]], candidates: list[str], rule: dict) -> tuple[dict, list[dict]]:
    assert len(candidates) == 642 and len(set(candidates)) == 642
    if len(rows) != len(candidates):
        raise RuntimeError("the frozen 642-row predictor census must contain all candidate species")
    if len({r["species"] for r in rows}) != len(rows) or {r["species"] for r in rows} != set(candidates):
        raise RuntimeError("the predictor census differs from frozen candidate species identity")
    if any(r.get("status") not in ALLOWED_STATUSES for r in rows):
        raise RuntimeError("unknown, non-frozen predictor status")

    by_species = {r["species"]: r for r in rows}
    complete = sorted((r for r in rows if r["status"] == "complete"), key=lambda r:r["species"])
    for r in complete:
        try:
            fraction = float(r["unique_fraction_M_host"])
            condition = float(r["predictor_condition_number"])
            edges = int(r["edges"])
        except (TypeError, ValueError) as exc:
            raise RuntimeError("missing predictor diagnostics for a complete species") from exc
        if (not math.isfinite(fraction) or not math.isfinite(condition)
                or not 0 <= fraction <= 1 or condition < 1 or edges < 1):
            raise RuntimeError("non-finite or invalid complete-species diagnostics")

    config = rule["predictor_information"]
    minimum = int(rule["external_admissibility"]["minimum_species_with_complete_frozen_predictor"])
    min_fraction = float(config["minimum_species_fraction_with_unique_fraction_at_least_0.05"])
    min_median = float(config["minimum_panel_median_unique_fraction"])
    max_cond = float(config["maximum_standardized_predictor_condition_number"])
    min_cond_fraction = float(config["minimum_species_fraction_passing_condition_number"])
    if minimum != 500 or min_fraction != 0.70 or min_median != 0.10 or max_cond != 30 or min_cond_fraction != 0.90:
        raise RuntimeError("frozen information-gate thresholds drift")

    digest_key = lambda sp: (
        hashlib.sha256(f"{SPLIT_NAMESPACE}|{SOURCE_SHA256}|{sp}".encode("utf-8")).hexdigest(), sp
    )
    selected = sorted((r["species"] for r in complete), key=digest_key)
    n_development = len(selected) // 2
    assignments = [
        {"species": sp, "panel": "development" if i < n_development else "confirmatory",
         "split_key_sha256": digest_key(sp)[0]}
        for i, sp in enumerate(selected)
    ]
    panel_by_species = {r["species"]: r["panel"] for r in assignments}

    fractions = [float(r["unique_fraction_M_host"]) for r in complete]
    conds = [float(r["predictor_condition_number"]) for r in complete]
    frac_pass = (sum(x >= 0.05 for x in fractions) / len(fractions)) if fractions else 0.0
    cond_pass = (sum(x <= max_cond for x in conds) / len(conds)) if conds else 0.0
    medians = {
        panel: (statistics.median(
            float(by_species[sp]["unique_fraction_M_host"])
            for sp in selected if panel_by_species[sp] == panel
        ) if any(panel_by_species[sp] == panel for sp in selected) else None)
        for panel in ("development", "confirmatory")
    }
    reasons = []
    if len(complete) < minimum:
        reasons.append("INSUFFICIENT_COMPLETE_EXTERNAL_PREDICTOR")
    if frac_pass < min_fraction:
        reasons.append("INSUFFICIENT_NONREDUNDANT_HOST_VARIANCE")
    if cond_pass < min_cond_fraction:
        reasons.append("PREDICTOR_COLLINEARITY_TOO_HIGH")
    if any(x is None or x < min_median for x in medians.values()):
        reasons.append("PANEL_UNIQUE_HOST_VARIANCE_TOO_SMALL")
    q = {
        "candidate_species": len(candidates),
        "complete_species": len(complete),
        "status_counts": {v: sum(r["status"] == v for r in rows) for v in sorted(ALLOWED_STATUSES)},
        "split_counts": {"development": n_development, "confirmatory": len(selected)-n_development},
        "fraction_unique_at_least_0_05": frac_pass,
        "fraction_condition_at_most_30": cond_pass,
        "panel_median_unique_fraction": medians,
        "min_complete_species": minimum,
        "min_fraction_unique_at_least_0_05": min_fraction,
        "min_panel_median_unique_fraction": min_median,
        "max_condition_number": max_cond,
        "min_fraction_condition_at_most_30": min_cond_fraction,
        "qualification_failures": reasons,
        "decision": ("PASS_TO_SYNTHETIC_QUALIFICATION" if not reasons else
                     "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_PREDICTOR_INFORMATION"),
    }
    return q, assignments


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--species", type=Path, required=True)
    ap.add_argument("--predictor-receipt", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--qualification", type=Path, required=True)
    ap.add_argument("--assignments", type=Path, required=True)
    args = ap.parse_args()

    if sha256(args.candidates) != CANDIDATE_SHA256:
        raise RuntimeError("candidate source SHA mismatch")
    frozen = json.loads(args.rule.read_text())
    if frozen.get("schema") != "ttf_historical_host_memory_qualification_rule_v0.2":
        raise RuntimeError("wrong frozen v0.2 qualification rule")
    if frozen["panel_split"]["ranking"] != "SHA256('historical-host-memory-split-v0.1|<source_sha256>|<species>'), species tie-break":
        raise RuntimeError("frozen split rule changed")
    if float(frozen["panel_split"]["development_fraction"]) != 0.5:
        raise RuntimeError("frozen development fraction changed")
    if any(v is True for v in frozen["retuning"].values()):
        raise RuntimeError("v0.2 permits forbidden scientific retuning")

    receipt = json.loads(args.predictor_receipt.read_text())
    if receipt.get("schema") != "ttf_historical_host_memory_predictor_result_v0.1":
        raise RuntimeError("wrong predictor receipt schema")
    if receipt["inputs_sha256"]["candidates"] != CANDIDATE_SHA256:
        raise RuntimeError("predictor was not built from frozen candidate source")
    if receipt["outputs_sha256"]["species_diagnostics"] != sha256(args.species):
        raise RuntimeError("predictor diagnostics SHA mismatch")

    with args.candidates.open(newline="", encoding="utf-8") as f:
        candidates = [r["species"].strip() for r in csv.DictReader(f)]
    with args.species.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    q, assignments = validate_and_qualify(rows, candidates, frozen)

    args.assignments.parent.mkdir(parents=True, exist_ok=True)
    with args.assignments.open("w", newline="", encoding="utf-8") as f:
        fields = ["species","panel","split_key_sha256"]
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(assignments)
    q.update({
        "schema": "ttf_historical_host_memory_information_qualification_v0.2",
        "genetic_response_opened": False,
        "synthetic_type1_power_qualified": False,
        "source_sha256": SOURCE_SHA256,
        "inputs_sha256": {
            "candidate_csv": CANDIDATE_SHA256,
            "predictor_species_csv": sha256(args.species),
            "predictor_receipt": sha256(args.predictor_receipt),
            "qualification_rule": sha256(args.rule),
        },
        "assignments_sha256": sha256(args.assignments),
    })
    args.qualification.parent.mkdir(parents=True, exist_ok=True)
    args.qualification.write_text(json.dumps(q, indent=2, sort_keys=True) + "\n")
    print(json.dumps(q, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
