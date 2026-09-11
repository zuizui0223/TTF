from __future__ import annotations

import argparse
import json
from pathlib import Path


def aggregate(rule: dict, batch_paths: list[Path]) -> dict:
    if rule["status"] != "frozen_before_development_diagnostic_outcomes":
        raise ValueError("diagnostic rule is not frozen")
    n_expected = int(rule["diagnostic"]["replicates_per_cell"])
    alpha = float(rule["diagnostic"]["alpha_for_descriptive_rejection_fraction"])
    profiles = list(rule["geometry_profiles"])
    responses = [x["id"] for x in rule["response_families"]]
    rows: dict[tuple[str, str], dict[int, dict]] = {(p, r): {} for p in profiles for r in responses}
    for path in batch_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "ttfm_q5_factor_diagnostic_batch_v0.1":
            raise ValueError(f"unexpected batch schema: {path}")
        key = (payload["profile"], payload["response_id"])
        if key not in rows:
            raise ValueError(f"unexpected diagnostic cell: {key}")
        for row in payload["rows"]:
            rep = int(row["replicate"])
            if rep in rows[key]:
                raise ValueError(f"duplicate replicate {rep} in {key}")
            rows[key][rep] = row

    cells = []
    for profile in profiles:
        for response_id in responses:
            key = (profile, response_id)
            rep_map = rows[key]
            if set(rep_map) != set(range(n_expected)):
                raise ValueError(f"coverage mismatch for {key}")
            rr = [rep_map[i] for i in range(n_expected)]
            m_rej = sum(float(x["ttf_m_p_value"]) <= alpha for x in rr)
            c_rej = sum(float(x["ttf_c_p_value"]) <= alpha for x in rr)
            cells.append({
                "profile": profile,
                "response_id": response_id,
                "n_replicates": n_expected,
                "ttf_m_rejection_fraction": m_rej / n_expected,
                "ttf_c_rejection_fraction": c_rej / n_expected,
                "ttf_m_mean_statistic": sum(float(x["ttf_m_statistic"]) for x in rr) / n_expected,
                "ttf_c_mean_statistic": sum(float(x["ttf_c_statistic"]) for x in rr) / n_expected,
                "mean_train_effective_n": sum(float(x["train_effective_n_mean"]) for x in rr) / n_expected,
                "mean_eval_effective_n": sum(float(x["eval_effective_n_mean"]) for x in rr) / n_expected,
                "mean_train_k": sum(float(x["train_k_mean"]) for x in rr) / n_expected,
                "mean_eval_k": sum(float(x["eval_k_mean"]) for x in rr) / n_expected,
            })

    baseline = {(x["response_id"]): x for x in cells if x["profile"] == "balanced"}
    for cell in cells:
        b = baseline[cell["response_id"]]
        cell["ttf_c_rejection_delta_vs_balanced"] = (
            cell["ttf_c_rejection_fraction"] - b["ttf_c_rejection_fraction"]
        )
        cell["ttf_c_statistic_delta_vs_balanced"] = (
            cell["ttf_c_mean_statistic"] - b["ttf_c_mean_statistic"]
        )

    by_response = {}
    for response_id in responses:
        subset = [x for x in cells if x["response_id"] == response_id]
        ranked = sorted(
            subset,
            key=lambda x: abs(float(x["ttf_c_rejection_delta_vs_balanced"])),
            reverse=True,
        )
        by_response[response_id] = {
            "balanced_ttf_c_rejection_fraction": baseline[response_id]["ttf_c_rejection_fraction"],
            "profiles_ranked_by_absolute_rejection_delta": [
                {
                    "profile": x["profile"],
                    "ttf_c_rejection_fraction": x["ttf_c_rejection_fraction"],
                    "delta_vs_balanced": x["ttf_c_rejection_delta_vs_balanced"],
                }
                for x in ranked
            ],
        }

    return {
        "schema": "ttfm_q5_factor_diagnostic_result_v0.1",
        "status": "development_factor_isolation_complete",
        "formal_pass_fail": None,
        "candidate_selected": None,
        "q5_result_remains_immutable_fail": True,
        "cells": cells,
        "by_response": by_response,
        "claim_firewall": rule["claim_firewall"],
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rule", required=True)
    p.add_argument("--batch-dir", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    rule = json.loads(Path(args.rule).read_text(encoding="utf-8"))
    paths = sorted(Path(args.batch_dir).rglob("*.json"))
    if not paths:
        raise ValueError("no diagnostic batches found")
    result = aggregate(rule, paths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"formal_pass_fail": None, "candidate_selected": None, "cells": len(result["cells"])}, sort_keys=True))


if __name__ == "__main__":
    main()
