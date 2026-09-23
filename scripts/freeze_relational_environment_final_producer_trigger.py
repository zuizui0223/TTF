#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


AUDIT_SCHEMA = "ttf_relational_environment_request_error_recovery_audit_v0.3"
PROMOTION_SCHEMA = "ttf_relational_environment_final_producer_promotion_rule_v0.2"
CONTRACT_SCHEMA = "ttf_relational_environment_final_producer_contract_v0.2"
PRODUCER_SCHEMA = "ttf_relational_environment_relation_producer_v0.1"
TRIGGER_SCHEMA = "ttf_relational_environment_final_producer_trigger_v0.2"


def load(path: Path, schema: str) -> dict:
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def closed_firewall(payload: dict, *, label: str) -> None:
    firewall = payload.get("response_firewall")
    if not isinstance(firewall, dict) or not firewall:
        raise RuntimeError(f"{label} response firewall missing")
    if any(value is not False for value in firewall.values()):
        raise RuntimeError(f"{label} response firewall is open")


def build_trigger(
    *,
    audit: dict,
    promotion: dict,
    contract: dict,
    producer: dict,
    audit_artifact_id: int,
    audit_artifact_digest: str,
    frozen_on: str,
) -> dict:
    closed_firewall(audit, label="audit")
    closed_firewall(promotion, label="promotion")
    closed_firewall(contract, label="contract")
    closed_firewall(producer, label="producer")

    if audit.get("status") != "PASS_FINAL_ROUND_ZERO_REQUEST_ERROR":
        raise RuntimeError("final producer trigger requires zero-error final-round audit PASS")
    if int(audit.get("retry_round", -1)) != 4:
        raise RuntimeError("final-round audit retry round drift")
    if int(audit.get("request_error_count_after_final_round", -1)) != 0:
        raise RuntimeError("final-round audit still contains REQUEST_ERROR")
    if list(audit.get("unresolved_request_error_species", [])) != []:
        raise RuntimeError("final-round audit unresolved species must be empty")
    if audit.get("network_retrieval_performed") is not False:
        raise RuntimeError("audit-only stage unexpectedly performed network retrieval")
    if audit.get("no_fifth_round_authorized") is not True:
        raise RuntimeError("audit does not preserve no-fifth-round rule")
    if audit.get("relation_result_seen") is not False:
        raise RuntimeError("relation result was seen before producer trigger")
    if audit.get("genetic_response_used") is not False:
        raise RuntimeError("genetic response was used before producer trigger")

    if promotion.get("status") != "FROZEN_BEFORE_FINAL_ROUND_RESULT_AND_BEFORE_RELATION_RESULT":
        raise RuntimeError("promotion rule is not the frozen pre-result rule")
    if contract.get("status") != "FROZEN_BEFORE_FINAL_ROUND_RESULT_AND_BEFORE_RELATION_RESULT":
        raise RuntimeError("producer contract is not the frozen pre-result contract")
    if promotion.get("no_fifth_round_authorized") is not True:
        raise RuntimeError("promotion rule does not forbid a fifth round")
    if contract.get("merge_contract", {}).get("no_fifth_round_authorized") is not True:
        raise RuntimeError("producer contract does not forbid a fifth round")

    recovery_run = int(promotion["final_recovery_run_id"])
    if int(contract["final_recovery_run_id"]) != recovery_run:
        raise RuntimeError("promotion/contract recovery-run mismatch")
    if int(audit["recovery_workflow_run_id"]) != recovery_run:
        raise RuntimeError("audit recovery-run mismatch")
    audit_run = int(audit["audit_workflow_run_id"])
    if audit_run <= 0:
        raise RuntimeError("audit workflow-run id is invalid")
    promotion_authority = promotion.get("final_recovery_audit_authority", {})
    contract_authority = contract.get("final_recovery_audit_authority", {})
    if int(promotion_authority.get("recovery_run_id", -1)) != recovery_run:
        raise RuntimeError("promotion audit-authority recovery-run mismatch")
    if int(contract_authority.get("recovery_run_id", -1)) != recovery_run:
        raise RuntimeError("contract audit-authority recovery-run mismatch")
    obsolete = set(map(int, promotion.get("obsolete_audit_runs", [])))
    obsolete.update(map(int, contract_authority.get("obsolete_audit_runs", [])))
    if audit_run in obsolete:
        raise RuntimeError("audit workflow run is explicitly obsolete")
    if int(audit.get("species", -1)) != 556:
        raise RuntimeError("audit species count drift")
    if int(audit.get("batch_count", -1)) != 186:
        raise RuntimeError("audit batch count drift")

    artifact_name = str(promotion["final_recovery_audit_artifact_name"])
    if str(contract["final_recovery_audit_artifact_name"]) != artifact_name:
        raise RuntimeError("promotion/contract audit artifact-name mismatch")
    digest = str(audit_artifact_digest)
    if not digest.startswith("sha256:") or len(digest) != 71:
        raise RuntimeError("audit artifact digest must be sha256:<64 hex>")
    try:
        int(digest.split(":", 1)[1], 16)
    except ValueError as exc:
        raise RuntimeError("audit artifact digest is not hexadecimal") from exc

    if producer.get("status") != "FROZEN_RELATION_PRODUCER_BEFORE_RESULT":
        raise RuntimeError("canonical relation producer is not frozen pre-result")
    if int(producer["workflow_run_id"]) != int(promotion["canonical_relation_producer_run_id"]):
        raise RuntimeError("canonical producer run drift")
    if producer.get("relation_result_seen") is not False:
        raise RuntimeError("canonical producer has seen relation result")
    if producer.get("genetic_response_used") is not False:
        raise RuntimeError("canonical producer has used genetic response")

    return {
        "schema": TRIGGER_SCHEMA,
        "status": "AUTHORIZE_FINAL_PRODUCER_AFTER_ZERO_ERROR_AUDIT",
        "frozen_on": str(frozen_on),
        "recovery_run_id": recovery_run,
        "recovery_audit_run_id": audit_run,
        "recovery_audit_artifact_id": int(audit_artifact_id),
        "recovery_audit_artifact_name": artifact_name,
        "recovery_audit_artifact_digest": digest,
        "promotion_rule": "docs/supporting/relational_environment_final_producer_promotion_rule_v0.2.json",
        "producer_contract": "benchmarks/frozen/relational_environment_final_producer_contract_v0.2.json",
        "relation_result_seen": False,
        "genetic_response_used": False,
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Freeze the unique final Study-B producer trigger after zero-error round-4 audit."
    )
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--promotion-rule", type=Path, required=True)
    parser.add_argument("--producer-contract", type=Path, required=True)
    parser.add_argument("--producer-receipt", type=Path, required=True)
    parser.add_argument("--audit-artifact-id", type=int, required=True)
    parser.add_argument("--audit-artifact-digest", required=True)
    parser.add_argument("--frozen-on", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = build_trigger(
        audit=load(args.audit, AUDIT_SCHEMA),
        promotion=load(args.promotion_rule, PROMOTION_SCHEMA),
        contract=load(args.producer_contract, CONTRACT_SCHEMA),
        producer=load(args.producer_receipt, PRODUCER_SCHEMA),
        audit_artifact_id=args.audit_artifact_id,
        audit_artifact_digest=args.audit_artifact_digest,
        frozen_on=args.frozen_on,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
