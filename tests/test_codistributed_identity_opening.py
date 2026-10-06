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
