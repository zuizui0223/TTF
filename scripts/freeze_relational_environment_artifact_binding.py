#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


SCHEMA = "ttf_relational_environment_relation_artifact_binding_v0.3"
ARTIFACT = "relational-environment-relation-v0.3"
TRANSPORT_EXECUTION_DEFAULT = Path("benchmarks/frozen/relational_environment_transport_execution_v0.6.json")
RELATION_PRODUCER_DEFAULT = Path("benchmarks/frozen/relational_environment_relation_producer_v0.1.json")


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
    ap.add_argument("--design-npz", type=Path)
    ap.add_argument("--occurrence-ledger", type=Path, required=True)
    ap.add_argument("--relation-rule", type=Path, required=True)
    ap.add_argument("--freshness-rule", type=Path, required=True)
    ap.add_argument("--transport-execution", type=Path, default=TRANSPORT_EXECUTION_DEFAULT)
    ap.add_argument("--relation-producer", type=Path, default=RELATION_PRODUCER_DEFAULT)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    if args.workflow_run_id <= 0 or args.artifact_id <= 0:
        raise ValueError("workflow-run-id and artifact-id must be positive")
    if len(args.head_sha) != 40 or any(c not in "0123456789abcdef" for c in args.head_sha.lower()):
        raise ValueError("head-sha must be a 40-character hexadecimal commit SHA")

    transport = json.loads(args.transport_execution.read_text())
    if transport.get("schema") != "ttf_relational_environment_transport_execution_v0.6":
        raise RuntimeError("unexpected Study-B transport execution schema")
    if transport.get("status") != "FROZEN_AUTHORITATIVE_COMPOSITE_TRANSPORT_BEFORE_RELATION_RESULT":
        raise RuntimeError("Study-B transport execution is not authoritative")
    if transport.get("relation_producer_receipt") != str(args.relation_producer):
        raise RuntimeError("Study-B relation-producer receipt pointer drift")
    execution = json.loads(args.relation_producer.read_text())
    if execution.get("schema") != "ttf_relational_environment_relation_producer_v0.1":
        raise RuntimeError("unexpected Study-B relation producer schema")
    if execution.get("status") != "FROZEN_RELATION_PRODUCER_BEFORE_RESULT":
        raise RuntimeError("Study-B relation producer was not frozen pre-result")
    if execution.get("transport_execution") != str(args.transport_execution):
        raise RuntimeError("Study-B relation producer transport pointer drift")
    if int(execution["workflow_run_id"]) != int(args.workflow_run_id):
        raise RuntimeError("relation artifact does not come from frozen relation-producer run")
    if str(execution["workflow_head_sha"]).lower() != args.head_sha.lower():
        raise RuntimeError("relation artifact head SHA drift from frozen relation-producer run")
    if execution["relation_artifact_name"] != ARTIFACT:
        raise RuntimeError("relation-producer artifact-name drift")

    occurrence = json.loads(args.occurrence_ledger.read_text())
    if occurrence.get("schema") != "ttf_relational_environment_occurrence_acquisition_v0.2":
        raise RuntimeError("unexpected Study-B occurrence acquisition schema")
    if int(occurrence.get("species", -1)) != int(transport["acceptance_gate"]["exact_species_ledgers"]):
        raise RuntimeError("Study-B occurrence ledger does not cover exact frozen species universe")
    request_errors = int((occurrence.get("status_counts") or {}).get("REQUEST_ERROR", 0))
    if request_errors != int(transport["acceptance_gate"]["request_error_count_must_equal"]):
        raise RuntimeError(
            f"Study-B transport is technically incomplete: REQUEST_ERROR={request_errors}"
        )
    if any(bool(v) for v in occurrence.get("response_firewall", {}).values()):
        raise RuntimeError("Study-B occurrence acquisition response firewall is open")

    summary = json.loads(args.summary.read_text())
    if summary.get("schema") != "ttf_relational_environment_relation_design_v0.3":
        raise RuntimeError("unexpected relation summary schema")
    if summary.get("status") not in {
        "PASS_RESPONSE_BLIND_ENVIRONMENT_RELATION_DESIGN",
        "NOT_EVALUABLE_ENVIRONMENT_FRESHNESS",
        "NOT_EVALUABLE_GENERAL_CROSS_TAXON_RELATIONAL_TTF",
    }:
        raise RuntimeError("unexpected relation summary status")
    if any(bool(v) for v in summary.get("response_firewall", {}).values()):
        raise RuntimeError("relation summary response firewall is open")

    rule_sha = sha256_path(args.relation_rule)
    freshness_sha = sha256_path(args.freshness_rule)
    inputs = summary.get("inputs", {})
    if inputs.get("rule_sha256") != rule_sha:
        raise RuntimeError("relation-rule hash drift in relation summary")
    if inputs.get("freshness_rule_sha256") != freshness_sha:
        raise RuntimeError("freshness-rule hash drift in relation summary")
    candidate_sha = sha256_path(Path("benchmarks/frozen/relational_fresh_candidate_1000_v0.1.csv"))
    if inputs.get("candidate_csv_sha256") != candidate_sha:
        raise RuntimeError("candidate CSV hash drift in relation summary")
    if summary.get("status") == "PASS_RESPONSE_BLIND_ENVIRONMENT_RELATION_DESIGN":
        if args.design_npz is None or not args.design_npz.is_file():
            raise RuntimeError("passing relation design requires the bound design NPZ")
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
            "occurrence_ledger_v0.2.json": sha256_path(args.occurrence_ledger),
            **({
                "relational_environment_design_v0.3.npz": sha256_path(args.design_npz)
            } if args.design_npz is not None and args.design_npz.is_file() else {}),
        },
        "rules_sha256": {
            "relational_environment_relation_rule_v0.3.json": rule_sha,
            "relational_future_family_freshness_amendment_v0.1.json": freshness_sha,
            "relational_environment_transport_execution_v0.6.json": sha256_path(args.transport_execution),
            "relational_environment_relation_producer_v0.1.json": sha256_path(args.relation_producer),
        },
        "transport_integrity": {
            "exact_species_ledgers": int(occurrence["species"]),
            "request_error_count": request_errors,
            "status_counts": occurrence.get("status_counts", {}),
            "relation_producer_workflow_run_id": int(execution["workflow_run_id"]),
            "relation_producer_head_sha": str(execution["workflow_head_sha"]),
            "corrected_transport_core_git_blobs": dict(
                transport["scientific_contract"]["corrected_transport_core_git_blobs"]
            ),
            "base_run_id": int(transport["frozen_base_component"]["workflow_run_id"]),
            "cutover_run_id": int(transport["cutover_component"]["workflow_run_id"]),
            "timeout_repair_run_id": int(transport["timeout_repair_component"]["workflow_run_id"]),
            "base_species": int(transport["final_partition"]["base_species"]),
            "cutover_success_species": int(transport["final_partition"]["cutover_success_species"]),
            "timeout_repair_species": int(transport["final_partition"]["timeout_repair_species"]),
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
