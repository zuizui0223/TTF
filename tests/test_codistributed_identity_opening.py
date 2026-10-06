import json
from pathlib import Path
import subprocess
import sys


ROOT=Path(__file__).resolve().parents[1]


def test_codistributed_empirical_opening_rule_is_frozen_pre_survivor_result():
    p=json.loads(
        (ROOT/"docs/supporting/genetic_codistributed_recurrence_empirical_opening_rule_v0.1.json").read_text()
    )
    assert p["schema"]=="ttf_genetic_codistributed_recurrence_empirical_opening_rule_v0.1"
    assert p["status"]=="FROZEN_BEFORE_SURVIVOR_FIXED_TST_REQUALIFICATION_RESULT"
    assert p["frozen_empirical_domain"]["species"]==326
    assert p["frozen_empirical_domain"]["sources"]==166
    assert p["frozen_empirical_domain"]["targets"]==160
    assert p["frozen_empirical_domain"]["directed_dyads"]==26560
    assert p["primary_inference"]["alpha"]==0.05
    assert p["primary_inference"]["asymptotic_two_way_cluster_p_value_authority"] is False
    assert all(v is False for v in p["response_firewall"].values())


def _help(path):
    completed=subprocess.run(
        [sys.executable,str(ROOT/path),"--help"],
        cwd=ROOT,capture_output=True,text=True,check=False,
    )
    assert completed.returncode==0,completed.stderr
    return completed.stdout


def test_identity_authorization_cli_imports_without_opening_response():
    out=_help("scripts/authorize_codistributed_recurrence_identity_opening.py")
    assert "--survivor-qualification" in out
    assert "--reference-dir" in out
    assert "--full-mask-ledger" in out


def test_one_shot_empirical_cli_imports_without_opening_response():
    out=_help("scripts/run_codistributed_recurrence_confirmatory_empirical.py")
    assert "--authorization" in out
    assert "--center-npz" in out
    assert "--reference-dir" in out


def test_full_confirmatory_mask_ledger_is_exactly_frozen_and_response_blind():
    import hashlib
    path=ROOT/"benchmarks/frozen/genetic_codistributed_recurrence_confirmatory_mask_full_v0.1.json"
    raw=path.read_bytes()
    assert hashlib.sha256(raw).hexdigest()=="0f3035fbc8a88a18ba274f4e97ada800b5688b8d40c26568063fbba68118e917"
    p=json.loads(raw)
    assert p["schema"]=="ttf_genetic_codistributed_recurrence_confirmatory_mask_result_v0.1"
    assert p["status"]=="PASS_TO_SURVIVOR_INFORMATION_GATE"
    assert len(p["species"])==381
    assert p["surviving_species"]==326
    assert all(row.get("nucleotide_identity_persisted") is False for row in p["species"])
    assert all(v is False for v in p["response_firewall"].values())


def test_empirical_execution_contract_binds_authoritative_resume_before_result():
    p=json.loads(
        (ROOT/"docs/supporting/genetic_codistributed_recurrence_empirical_execution_contract_v0.1.json").read_text()
    )
    assert p["schema"]=="ttf_genetic_codistributed_recurrence_empirical_execution_contract_v0.1"
    assert p["status"]=="FROZEN_BEFORE_AUTHORITATIVE_SURVIVOR_REQUALIFICATION_RESULT"
    q=p["authoritative_survivor_qualification"]
    assert q["workflow_run_id"]==37406111501
    assert q["qualification_artifact"]=="codistributed-survivor-resume-qualification-v0.1"
    assert q["center_artifact"]=="codistributed-survivor-resume-center-v0.1"
    assert p["frozen_inputs"]["full_mask_ledger_sha256"]=="0f3035fbc8a88a18ba274f4e97ada800b5688b8d40c26568063fbba68118e917"
    assert all(v is False for v in p["response_firewall"].values())


def test_empirical_workflow_is_push_only_and_blocks_reruns():
    path=ROOT/".github/workflows/codistributed-recurrence-confirmatory-empirical.yml"
    text=path.read_text()
    assert "workflow_dispatch" not in text
    assert "genetic_codistributed_recurrence_empirical_execution_v0.1.json" in text
    assert 'test "$GITHUB_EVENT_NAME" = "push"' in text
    assert 'test "$GITHUB_RUN_ATTEMPT" = "1"' in text
    assert 'run-id: 37406111501' in text
    assert "cmp \"$AUTHORIZATION\" results/pre-opening/reproduced_authorization.json" in text
