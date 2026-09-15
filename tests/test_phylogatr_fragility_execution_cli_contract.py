from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path


EXECUTION_RULE = Path(
    "docs/supporting/genetic_phylogatr_gate_d_fragility_execution_v0.1.json"
)
EXPECTED_EXECUTION_RULE_BLOB = "a06eac22f47d90a29328a837c401a6c8b938c6ab"
SCRIPTS = [
    "scripts/plan_phylogatr_gate_d_fragility_synthetic.py",
    "scripts/run_phylogatr_gate_d_fragility_reference_shard.py",
    "scripts/aggregate_phylogatr_gate_d_fragility_references.py",
    "scripts/run_phylogatr_gate_d_fragility_observed_shard.py",
    "scripts/aggregate_phylogatr_gate_d_fragility_curve.py",
]


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def test_frozen_execution_rule_blob_matches_cli_pin() -> None:
    assert _git_blob_sha1(EXECUTION_RULE) == EXPECTED_EXECUTION_RULE_BLOB


def test_fragility_execution_clis_import_and_expose_help() -> None:
    for script in SCRIPTS:
        completed = subprocess.run(
            [sys.executable, script, "--help"],
            check=False,
            text=True,
            capture_output=True,
        )
        assert completed.returncode == 0, (script, completed.stderr)


def test_planner_rejects_same_schema_execution_rule_tampering(tmp_path: Path) -> None:
    tampered = tmp_path / "tampered_execution_rule.json"
    tampered.write_text(EXECUTION_RULE.read_text().replace("0.875", "0.876", 1))
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/plan_phylogatr_gate_d_fragility_synthetic.py",
            "--geometry-plan-receipt",
            str(tmp_path / "unused_geometry_plan.json"),
            "--phase2-manifest",
            str(tmp_path / "unused_phase2.json"),
            "--phase3-rule",
            str(tmp_path / "unused_phase3.json"),
            "--formal-qualification",
            str(tmp_path / "unused_formal.json"),
            "--execution-rule",
            str(tampered),
            "--output",
            str(tmp_path / "unused_output.json"),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "frozen fragility execution rule Git blob SHA drift" in completed.stderr


def test_production_fragility_outputs_do_not_define_secondary_gate_fields() -> None:
    production = [Path(path) for path in SCRIPTS] + [
        Path("src/ttf/phylogatr_fragility_execution.py")
    ]
    joined = "\n".join(path.read_text() for path in production)
    assert "would_pass" not in joined
    assert "phase4_identity_opening_eligible" not in joined
