import json
import subprocess
import sys
import importlib.util
from pathlib import Path


SCRIPT = Path("scripts/freeze_relational_environment_artifact_binding.py")


def load_module():
    spec = importlib.util.spec_from_file_location("env_binding", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_binding_schema_and_artifact_name_are_fixed():
    module = load_module()
    assert module.SCHEMA == "ttf_relational_environment_relation_artifact_binding_v0.3"
    assert module.ARTIFACT == "relational-environment-relation-v0.3"


def test_qualification_workflow_verifies_bound_relation_hashes():
    workflow = Path(".github/workflows/relational-environment-qualification-v03.yml").read_text()
    assert "relational_environment_relation_artifact_binding_v0.3.json" in workflow
    assert "files_sha256" in workflow
    assert "relational_environment_design_v0.3.json" in workflow
    assert "relational_environment_design_v0.3.npz" in workflow
    assert "occurrence_ledger_v0.2.json" in workflow


def test_binding_rejects_any_transport_request_error(tmp_path):
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "schema": "ttf_relational_environment_relation_design_v0.3",
        "status": "NOT_EVALUABLE_ENVIRONMENT_FRESHNESS",
        "response_firewall": {
            "study_B_sequence_identity_opened": False,
            "study_B_pairwise_genetic_distances_opened": False,
            "study_B_T_st_computed": False,
            "study_B_beta_R_computed": False,
        },
    }))
    design = tmp_path / "design.npz"
    design.write_bytes(b"placeholder")
    ledger = tmp_path / "ledger.json"
    ledger.write_text(json.dumps({
        "schema": "ttf_relational_environment_occurrence_acquisition_v0.2",
        "species": 1000,
        "status_counts": {"REQUEST_ERROR": 1, "PASS_OCCURRENCE_GEOMETRY": 999},
        "response_firewall": {
            "study_B_sequence_identity_opened": False,
            "study_B_pairwise_genetic_distances_opened": False,
            "study_B_T_st_computed": False,
            "study_B_beta_R_computed": False,
        },
    }))
    rule = tmp_path / "rule.json"
    rule.write_text("{}")
    freshness = tmp_path / "freshness.json"
    freshness.write_text("{}")
    transport = tmp_path / "transport.json"
    producer = tmp_path / "producer.json"
    transport.write_text(json.dumps({
        "schema": "ttf_relational_environment_transport_execution_v0.7",
        "status": "FROZEN_FINAL_AUTHORITATIVE_COMPOSITE_TRANSPORT_BEFORE_RELATION_RESULT",
        "relation_producer_receipt": str(producer),
        "frozen_base_component": {
            "workflow_run_id": 1,
        },
        "cutover_component": {
            "workflow_run_id": 2,
        },
        "short_timeout_repair_component": {
            "workflow_run_id": 3,
        },
        "long_timeout_repair_component": {
            "workflow_run_id": 4,
        },
        "final_partition": {
            "base_species": 165,
            "cutover_success_species": 827,
            "short_repair_species": 5,
            "long_repair_species": 3,
        },
        "scientific_contract": {
            "corrected_transport_core_git_blobs": {},
        },
        "acceptance_gate": {
            "exact_species_ledgers": 1000,
            "request_error_count_must_equal": 0,
        },
    }))
    producer.write_text(json.dumps({
        "schema": "ttf_relational_environment_relation_producer_v0.1",
        "status": "FROZEN_RELATION_PRODUCER_BEFORE_RESULT",
        "transport_execution": str(transport),
        "workflow_run_id": 123,
        "workflow_head_sha": "a" * 40,
        "relation_artifact_name": "relational-environment-relation-v0.3",
    }))
    completed = subprocess.run([
        sys.executable,
        str(SCRIPT),
        "--workflow-run-id", "123",
        "--artifact-id", "456",
        "--head-sha", "a" * 40,
        "--summary", str(summary),
        "--design-npz", str(design),
        "--occurrence-ledger", str(ledger),
        "--relation-rule", str(rule),
        "--freshness-rule", str(freshness),
        "--transport-execution", str(transport),
        "--relation-producer", str(producer),
        "--output", str(tmp_path / "binding.json"),
    ], capture_output=True, text=True)
    assert completed.returncode != 0
    assert "REQUEST_ERROR=1" in completed.stderr
