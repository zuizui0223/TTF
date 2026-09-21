#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


SCHEMA = "ttf_relational_environment_relation_artifact_binding_v0.3"
ARTIFACT = "relational-environment-relation-v0.3"


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workflow-run-id", type=int, required=True)
    ap.add_argument("--artifact-id", type=int, required=True)
    ap.add_argument("--head-sha", required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--design-npz", type=Path, required=True)
    ap.add_argument("--occurrence-ledger", type=Path, required=True)
    ap.add_argument("--relation-rule", type=Path, required=True)
    ap.add_argument("--freshness-rule", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    if args.workflow_run_id <= 0 or args.artifact_id <= 0:
        raise ValueError("workflow-run-id and artifact-id must be positive")
    if len(args.head_sha) != 40 or any(c not in "0123456789abcdef" for c in args.head_sha.lower()):
        raise ValueError("head-sha must be a 40-character hexadecimal commit SHA")

    summary = json.loads(args.summary.read_text())
    if summary.get("schema") != "ttf_relational_environment_relation_design_v0.3":
        raise RuntimeError("unexpected relation summary schema")
    if summary.get("status") not in {
        "PASS_RESPONSE_BLIND_ENVIRONMENT_RELATION_DESIGN",
        "NOT_EVALUABLE_ENVIRONMENT_FRESHNESS",
    }:
        raise RuntimeError("unexpected relation summary status")
    if any(bool(v) for v in summary.get("response_firewall", {}).values()):
        raise RuntimeError("relation summary response firewall is open")

    rule_sha = sha256_path(args.relation_rule)
    freshness_sha = sha256_path(args.freshness_rule)
    if summary.get("status") == "PASS_RESPONSE_BLIND_ENVIRONMENT_RELATION_DESIGN":
        inputs = summary.get("inputs", {})
        if inputs.get("rule_sha256") != rule_sha:
            raise RuntimeError("relation-rule hash drift in relation summary")
        if inputs.get("freshness_rule_sha256") != freshness_sha:
            raise RuntimeError("freshness-rule hash drift in relation summary")
        if summary.get("design_npz_sha256") != sha256_path(args.design_npz):
            raise RuntimeError("design NPZ hash drift relative to relation summary")

    payload = {
        "schema": SCHEMA,
        "status": "FROZEN_RESPONSE_BLIND_RELATION_ARTIFACT",
        "artifact_name": ARTIFACT,
        "workflow_run_id": int(args.workflow_run_id),
        "artifact_id": int(args.artifact_id),
        "head_sha": args.head_sha.lower(),
        "relation_status": summary["status"],
        "files_sha256": {
            "relational_environment_design_v0.3.json": sha256_path(args.summary),
            "relational_environment_design_v0.3.npz": sha256_path(args.design_npz),
            "occurrence_ledger_v0.2.json": sha256_path(args.occurrence_ledger),
        },
        "rules_sha256": {
            "relational_environment_relation_rule_v0.3.json": rule_sha,
            "relational_future_family_freshness_amendment_v0.1.json": freshness_sha,
        },
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "workflow_run_id": payload["workflow_run_id"],
        "artifact_id": payload["artifact_id"],
        "head_sha": payload["head_sha"],
        "relation_status": payload["relation_status"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
