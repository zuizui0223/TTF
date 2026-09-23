from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path("scripts/freeze_relational_environment_final_producer_trigger.py")


def load_module():
    spec = importlib.util.spec_from_file_location("final_trigger", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def firewall():
    return {
        "Study_B_sequence_identity_opened": False,
        "Study_B_pairwise_genetic_distances_opened": False,
        "Study_B_T_st_computed": False,
        "Study_B_beta_R_computed": False,
    }


def fixtures():
    audit = {
        "schema": "ttf_relational_environment_request_error_recovery_audit_v0.3",
        "status": "PASS_FINAL_ROUND_ZERO_REQUEST_ERROR",
        "retry_round": 4,
        "recovery_workflow_run_id": 35827478403,
        "audit_workflow_run_id": 35850000001,
        "species": 556,
        "batch_count": 186,
        "request_error_count_after_final_round": 0,
        "unresolved_request_error_species": [],
        "network_retrieval_performed": False,
        "no_fifth_round_authorized": True,
        "relation_result_seen": False,
        "genetic_response_used": False,
        "response_firewall": firewall(),
    }
    promotion = {
        "schema": "ttf_relational_environment_final_producer_promotion_rule_v0.2",
        "status": "FROZEN_BEFORE_FINAL_ROUND_RESULT_AND_BEFORE_RELATION_RESULT",
        "final_recovery_run_id": 35827478403,
        "final_recovery_audit_authority": {"recovery_run_id": 35827478403},
        "obsolete_audit_runs": [35828444286],
        "final_recovery_audit_artifact_name": "relational-environment-request-error-recovery-audit-v0.3",
        "canonical_relation_producer_run_id": 35795821824,
        "no_fifth_round_authorized": True,
        "response_firewall": firewall(),
    }
    contract = {
        "schema": "ttf_relational_environment_final_producer_contract_v0.2",
        "status": "FROZEN_BEFORE_FINAL_ROUND_RESULT_AND_BEFORE_RELATION_RESULT",
        "final_recovery_run_id": 35827478403,
        "final_recovery_audit_authority": {"recovery_run_id": 35827478403, "obsolete_audit_runs": [35828444286]},
        "final_recovery_audit_artifact_name": "relational-environment-request-error-recovery-audit-v0.3",
        "merge_contract": {"no_fifth_round_authorized": True},
        "response_firewall": firewall(),
    }
    producer = {
        "schema": "ttf_relational_environment_relation_producer_v0.1",
        "status": "FROZEN_RELATION_PRODUCER_BEFORE_RESULT",
        "workflow_run_id": 35795821824,
        "relation_result_seen": False,
        "genetic_response_used": False,
        "response_firewall": firewall(),
    }
    return audit, promotion, contract, producer


def test_final_trigger_is_deterministic_and_binds_exact_runs():
    m = load_module()
    audit, promotion, contract, producer = fixtures()
    kwargs = dict(
        audit=audit,
        promotion=promotion,
        contract=contract,
        producer=producer,
        audit_artifact_id=123456,
        audit_artifact_digest="sha256:" + "a" * 64,
        frozen_on="2026-09-23",
    )
    first = m.build_trigger(**kwargs)
    second = m.build_trigger(**kwargs)
    assert first == second
    assert first["status"] == "AUTHORIZE_FINAL_PRODUCER_AFTER_ZERO_ERROR_AUDIT"
    assert first["recovery_run_id"] == 35827478403
    assert first["recovery_audit_run_id"] == 35850000001
    assert first["recovery_audit_artifact_id"] == 123456
    assert first["recovery_audit_artifact_digest"] == "sha256:" + "a" * 64
    assert all(value is False for value in first["response_firewall"].values())


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda x: x["audit"].update(status="INCOMPLETE_TECHNICAL_EXECUTION_AFTER_FINAL_ROUND"), "zero-error"),
        (lambda x: x["audit"].update(request_error_count_after_final_round=1), "REQUEST_ERROR"),
        (lambda x: x["audit"].update(unresolved_request_error_species=["Species x"]), "unresolved"),
        (lambda x: x["audit"].update(no_fifth_round_authorized=False), "no-fifth"),
        (lambda x: x["promotion"].update(final_recovery_run_id=1), "recovery-run"),
        (lambda x: x["producer"].update(relation_result_seen=True), "relation result"),
    ],
)
def test_final_trigger_rejects_non_authorized_states(mutation, message):
    m = load_module()
    audit, promotion, contract, producer = fixtures()
    state = {
        "audit": audit,
        "promotion": promotion,
        "contract": contract,
        "producer": producer,
    }
    mutation(state)
    with pytest.raises(RuntimeError, match=message):
        m.build_trigger(
            audit=audit,
            promotion=promotion,
            contract=contract,
            producer=producer,
            audit_artifact_id=123456,
            audit_artifact_digest="sha256:" + "a" * 64,
            frozen_on="2026-09-23",
        )


def test_final_trigger_rejects_bad_artifact_digest():
    m = load_module()
    audit, promotion, contract, producer = fixtures()
    with pytest.raises(RuntimeError, match="sha256"):
        m.build_trigger(
            audit=audit,
            promotion=promotion,
            contract=contract,
            producer=producer,
            audit_artifact_id=123456,
            audit_artifact_digest="not-a-digest",
            frozen_on="2026-09-23",
        )


def test_final_trigger_rejects_obsolete_polling_audit_run():
    m = load_module()
    audit, promotion, contract, producer = fixtures()
    audit["audit_workflow_run_id"] = 35828444286
    with pytest.raises(RuntimeError, match="obsolete"):
        m.build_trigger(
            audit=audit,
            promotion=promotion,
            contract=contract,
            producer=producer,
            audit_artifact_id=123456,
            audit_artifact_digest="sha256:" + "a" * 64,
            frozen_on="2026-09-23",
        )
