from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


RULE = Path("docs/supporting/genetic_phylogatr_phase3_self_detectability_rule_v0.1.json")


def _write_shard(path: Path, cell: str, p_values: list[float]) -> None:
    payload = {
        "schema": "ttf_genetic_phylogatr_phase3_self_evaluation_shard_v0.1",
        "status": "synthetic_self_detectability_evaluation_shard",
        "geometry_fingerprint_sha256": "f" * 64,
        "cell": cell,
        "start": 0,
        "stop": len(p_values),
        "n_worlds": len(p_values),
        "alpha": 0.05,
        "rows": [
            {"replicate": index, "statistic": float(index), "p_value": float(p_value)}
            for index, p_value in enumerate(p_values)
        ],
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "qualification_claim_made": False,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _run(tmp_path: Path, null_rejections: int) -> dict:
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    _write_shard(
        input_dir / "null.json",
        "null",
        [0.01] * null_rejections + [0.5] * (500 - null_rejections),
    )
    _write_shard(input_dir / "private.json", "private_A2", [0.01] * 500)
    auth = tmp_path / "auth.json"
    auth.write_text(
        json.dumps(
            {
                "schema": "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1",
                "geometry_fingerprint_sha256": "f" * 64,
            }
        )
        + "\n"
    )
    output = tmp_path / "result.json"
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/aggregate_phylogatr_phase3_self_qualification.py",
            "--input-dir",
            str(input_dir),
            "--self-rule",
            str(RULE),
            "--phase3-authorization",
            str(auth),
            "--output",
            str(output),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0, completed.stderr
    return json.loads(output.read_text())


def test_fresh_self_qualification_passes_clean_null_and_private_A2(tmp_path: Path) -> None:
    result = _run(tmp_path, null_rejections=0)
    assert result["status"] == "PASS"
    assert result["passed"] is True
    assert result["type1_gate"]["pass"] is True
    assert result["power_gate"]["pass"] is True
    assert result["confirmatory_sequence_identity_opened"] is False


def test_fresh_self_qualification_does_not_round_borderline_type1_rescue(tmp_path: Path) -> None:
    result = _run(tmp_path, null_rejections=37)
    assert result["status"] == "SELF_DETECTABILITY_NOT_QUALIFIED"
    assert result["passed"] is False
    assert result["type1_gate"]["observed_wilson95_upper"] > 0.10
