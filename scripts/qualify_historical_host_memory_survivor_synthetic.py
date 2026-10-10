#!/usr/bin/env python3
"""Independent, response-blind synthetic qualification of exact post-mask survivors.

Uses frozen private-spatial/IBD synthetic world parameters and unchanged test,
but a disjoint predeclared survivor seed namespace. Never reads nucleotide
identity, observed genetic distance, empirical turnover, or empirical beta.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from scripts.qualify_historical_host_memory_synthetic import (
        SOURCE_SHA, GEOMETRY_EDGE_SHA, GEOMETRY_LOCAL_SHA, ROLE_SHA,
        assess, endpoint_safe_ibd, file_sha, load_geometry, prepare_species, statistic, z,
    )
except ModuleNotFoundError:
    from qualify_historical_host_memory_synthetic import (
        SOURCE_SHA, GEOMETRY_EDGE_SHA, GEOMETRY_LOCAL_SHA, ROLE_SHA,
        assess, endpoint_safe_ibd, file_sha, load_geometry, prepare_species, statistic, z,
    )

ORIGINAL_PREDICTOR_SHA = "70b748d588d9ab10309ace080f85b8053c6e14aa96d7a84288724c1000272a5b"
PARENT_WORLD_SHA = "2aa61a0fd3e0fbf9fdf44f1ef7c4cffe4a6137cc91543d9971902d1349ebbf71"
QUAL_RULE_SHA = "2ec2594510021da653598e22a3d0f4b9e4cb226398aadbf342967ee532909415"
SURVIVOR_NAMESPACE = "historical-host-memory-survivor-v0.1"
GROUPS = {"reference": 1999, "evaluation": 500, "positive": 500}


def survivor_seed(category: str, species: str) -> int:
    raw = f"{SURVIVOR_NAMESPACE}|{category}|{SOURCE_SHA}|{species}"
    return int(hashlib.sha256(raw.encode()).hexdigest()[:16], 16)


def survivor_beta_worlds(pre: dict, group: str, worlds: int, block_size: int = 64):
    """Same equations as frozen pre-mask simulator, disjoint seed namespaces."""
    if group not in GROUPS or worlds != GROUPS[group]:
        raise RuntimeError("frozen survivor synthetic group/world size changed")
    rng = np.random.default_rng(survivor_seed(group, pre["name"]))
    out = np.empty(worlds, dtype=float)
    m = len(pre["x"])
    # Rebuild exactly the same private Fourier-basis type with independent seed.
    basis_rng = np.random.default_rng(survivor_seed("spatial-basis", pre["name"]))
    xyz = pre["xyz"]
    frequencies = basis_rng.normal(size=(3, 8)) / 200.0
    phases = xyz @ frequencies
    basis = np.column_stack((np.sin(phases), np.cos(phases)))
    for start in range(0, worlds, block_size):
        stop = min(start + block_size, worlds)
        width = stop - start
        mean_effect = 0.10 if group == "positive" else 0.0
        beta = rng.normal(loc=mean_effect, scale=0.10, size=width)
        weights = rng.normal(size=(16, width))
        field_nodes = basis @ weights / np.sqrt(16.0)
        private = z(np.abs(field_nodes[pre["left"]] - field_nodes[pre["right"]]), axis=0)
        noise = rng.normal(size=(m, width))
        response = (
            0.60 * pre["x"][:, None] +
            0.22 * pre["ms"][:, None] -
            0.18 * pre["dc"][:, None] +
            pre["mh"][:, None] * beta[None, :] +
            0.45 * private + 0.65 * noise
        )
        # Inherited exact leave-two-locality-out IBD correction.
        post = endpoint_safe_ibd(response, pre)
        sd = post.std(axis=0, ddof=0)
        if np.any(~np.isfinite(sd)) or np.any(sd <= 1e-12):
            raise RuntimeError("synthetic survivor response not finite")
        out[start:stop] = (pre["resid"] @ post) / (pre["norm"] * sd)
    if not np.all(np.isfinite(out)):
        raise RuntimeError("nonfinite survivor synthetic beta")
    return out


def validate_source_set(
    summary: dict, survivor_roles: list[dict], predictor_rows: list[dict],
    original_roles: list[dict], mask: dict
) -> tuple[list[str], dict[str, list[dict]]]:
    if summary.get("schema") != "ttf_historical_host_memory_survivor_information_v0.1":
        raise RuntimeError("wrong exact survivor information receipt")
    if summary.get("decision") != "PASS_TO_EXACT_SURVIVOR_SYNTHETIC_REQUALIFICATION":
        raise RuntimeError("survivor information did not qualify; do not simulate")
    if any(summary.get(k) is not False for k in
           ("genetic_response_opened", "nucleotide_identity_opened",
            "post_IBD_turnover_opened", "empirical_beta_host_opened",
            "survivor_synthetic_qualification_passed")):
        raise RuntimeError("survivor information response firewall is open")
    if mask.get("schema") != "ttf_historical_host_memory_confirmatory_mask_result_v0.1":
        raise RuntimeError("wrong confirmatory character-mask receipt")
    if (mask.get("source_bindings",{}).get("candidate_table_sha256")
            != "c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b"
            or mask.get("source_bindings",{}).get("panel_roles_sha256") != ROLE_SHA
            or mask.get("source_bindings",{}).get("locality_geometry_sha256") != GEOMETRY_LOCAL_SHA
            or mask.get("source_bindings",{}).get("edge_geometry_sha256") != GEOMETRY_EDGE_SHA
            or mask.get("source_bindings",{}).get("synthetic_pass_receipt_sha256")
                != "0498f71d636274a7f2e654f42591545cf26cbd0b97e603902377ea2ef702495e"
            or mask.get("source_bindings",{}).get("mask_rule_git_blob_sha1")
                != "1d5bec837157955cd2f989b2dd127d502d57318c"):
        raise RuntimeError("original mask decision provenance binding mismatch")
    if mask.get("status") != "PASS_TO_EXACT_SURVIVOR_INFORMATION_AND_SYNTHETIC_REQUALIFICATION":
        raise RuntimeError("frozen character-mask decision is not PASS")
    if summary.get("mask_result_sha256") is None:
        raise RuntimeError("mask source not bound")
    allowed = {r["species"] for r in original_roles if r["panel"] == "confirmatory"}
    if len(allowed) != 321:
        raise RuntimeError("pre-mask confirmatory set drift")
    kept = sorted(summary["surviving_species_names"])
    if not 200 <= len(kept) <= 321 or len(kept) != len(set(kept)):
        raise RuntimeError("survivor count violates frozen 200-species floor")
    if set(kept) != set(mask["survivor_names"]) or not set(kept).issubset(allowed):
        raise RuntimeError("survivor sample is not the original confirmatory mask subset")
    by_role = {r["species"]: r for r in original_roles}
    if len(survivor_roles) != len(kept) or {r["species"] for r in survivor_roles} != set(kept):
        raise RuntimeError("survivor roles do not cover exact mask survivors")
    if any(r != by_role[r["species"]] for r in survivor_roles):
        raise RuntimeError("forbidden role resplitting or replacement")
    pred = defaultdict(list)
    for r in predictor_rows:
        if r["species"] not in kept:
            raise RuntimeError("non-survivor predictor edge present")
        pred[r["species"]].append(r)
    if set(pred) != set(kept):
        raise RuntimeError("survivor predictor file lacks exact species")
    if len(predictor_rows) != int(summary.get("surviving_frozen_edges", -1)):
        raise RuntimeError("survivor predictor edge count drift")
    return kept, pred


def main() -> int:
    p = argparse.ArgumentParser()
    for field in (
        "localities", "edges", "original_roles", "original_predictors",
        "mask_receipt", "survivor_summary", "survivor_roles",
        "survivor_predictors", "qualification_rule", "original_world", "survivor_world"
    ):
        p.add_argument("--"+field.replace("_","-"), required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()

    pinned = {
        "localities": GEOMETRY_LOCAL_SHA, "edges": GEOMETRY_EDGE_SHA,
        "original_roles": ROLE_SHA, "original_predictors": ORIGINAL_PREDICTOR_SHA,
        "qualification_rule": QUAL_RULE_SHA, "original_world": PARENT_WORLD_SHA
    }
    for key, expected in pinned.items():
        if file_sha(getattr(a, key)) != expected:
            raise RuntimeError(f"frozen survivor simulation source SHA drift: {key}")
    spec = json.loads(a.survivor_world.read_text())
    rule = json.loads(a.qualification_rule.read_text())
    if (spec.get("status") != "FROZEN_BEFORE_ANY_REAL_CONFIRMATORY_MASK_RESULT"
            or spec.get("schema") != "ttf_historical_host_memory_survivor_synthetic_requalification_contract_v0.1"):
        raise RuntimeError("survivor synthetic protocol was not frozen")
    original = json.loads(a.original_world.read_text())
    frozen = rule["synthetic_qualification"]
    if (spec["randomization"]["independent_seed_namespace"] != SURVIVOR_NAMESPACE
            or spec["unmodified_science"]["null_reference_worlds"] != 1999
            or spec["unmodified_science"]["independent_null_evaluation_worlds"] != 500
            or spec["unmodified_science"]["positive_worlds"] != 500
            or spec["unmodified_science"]["alpha"] != 0.05
            or float(spec["unmodified_science"]["maximum_type1_wilson_95_upper"]) != 0.10
            or float(spec["unmodified_science"]["minimum_power_wilson_95_lower"]) != 0.80
            or int(frozen["null_reference_worlds"]) != 1999
            or int(frozen["evaluation_null_worlds"]) != 500
            or int(frozen["positive_worlds"]) != 500
            or original["response_world"]["iid_edge_noise"] != "0.65 times independent standard normal per edge"):
        raise RuntimeError("frozen simulation parameters drift")
    with a.original_roles.open(newline="", encoding="utf-8") as f:
        original_roles = list(csv.DictReader(f))
    with a.survivor_roles.open(newline="", encoding="utf-8") as f:
        roles = list(csv.DictReader(f))
    with a.survivor_predictors.open(newline="", encoding="utf-8") as f:
        predictors = list(csv.DictReader(f))
    mask = json.loads(a.mask_receipt.read_text())
    summary = json.loads(a.survivor_summary.read_text())
    if (file_sha(a.mask_receipt) != summary.get("mask_result_sha256")
            or file_sha(a.survivor_roles) != summary.get("survivor_roles_sha256")
            or file_sha(a.survivor_predictors) != summary.get("survivor_predictor_sha256")):
        raise RuntimeError("mask/survivor input hash mismatch")
    kept, pred = validate_source_set(summary, roles, predictors, original_roles, mask)
    loc, edges = load_geometry(a.localities, a.edges, set(kept))
    betas = {key: np.empty((len(kept), count), dtype=float) for key, count in GROUPS.items()}
    for i, name in enumerate(kept):
        prepared = prepare_species(name, pred[name], loc[name], edges[name])
        prepared["xyz"] = np.stack([loc[name][k] for k in sorted(loc[name])])
        for group, n in GROUPS.items():
            betas[group][i] = survivor_beta_worlds(prepared, group, n)
        if (i+1) % 25 == 0:
            print(json.dumps({"survivor_species_done":i+1,"total":len(kept)}), flush=True)
    test = {group: statistic(b) for group,b in betas.items()}
    outcome = assess(test["reference"],test["evaluation"],test["positive"])
    qualified = (outcome["type1_gate_pass"] and outcome["positive_power_gate_pass"])
    receipt = {
        "schema":"ttf_historical_host_memory_survivor_synthetic_result_v0.1",
        "status":("PASS_TO_SEPARATE_CONFIRMATORY_IDENTITY_AUTHORIZATION_REVIEW"
                  if qualified else "NOT_EVALUABLE_HISTORICAL_HOST_MEMORY_SURVIVOR_SYNTHETIC_QUALIFICATION"),
        "survivor_species":len(kept),
        "frozen_worlds":GROUPS,
        "mask_sha256":file_sha(a.mask_receipt),
        "survivor_summary_sha256":file_sha(a.survivor_summary),
        "survivor_predictor_sha256":file_sha(a.survivor_predictors),
        "survivor_roles_sha256":file_sha(a.survivor_roles),
        "survivor_contract_sha256":file_sha(a.survivor_world),
        "calibration":outcome,
        "nucleotide_identity_opened":False,
        "genetic_distance_opened":False,
        "post_IBD_empirical_turnover_opened":False,
        "empirical_beta_host_opened":False,
        "genetic_opening_authorized_by_this_result":False
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(receipt, sort_keys=True, indent=2)+"\n")
    print(json.dumps({"decision":receipt["status"], "survivors":len(kept),
                      "null_rejections":outcome["null_rejections"],
                      "positive_rejections":outcome["positive_rejections"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
