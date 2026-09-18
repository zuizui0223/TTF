from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


SCRIPT = Path("scripts/aggregate_genetic_conditional_order_qualification.py")


def _write(path: Path, payload: dict) -> str:
    path.write_text(json.dumps(payload, sort_keys=True) + "\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rule(tmp_path: Path) -> tuple[Path, str]:
    path = tmp_path / "rule.json"
    payload = {
        "schema": "ttf_genetic_conditional_order_qualification_rule_v0.1",
        "frozen_estimator": {"alpha": 0.05},
        "synthetic_worlds": {
            "private_null_cells": [
                {"cell": "private_A2", "worlds": 100, "residual_amplitude": 2.0}
            ],
            "positive_control": {
                "cell": "same_order_A2",
                "worlds": 100,
                "residual_amplitude": 2.0,
            },
        },
        "qualification": {
            "type1_wilson95_upper_ceiling": 0.10,
            "power_wilson95_lower_floor": 0.80,
        },
    }
    return path, _write(path, payload)


def _shard(
    tmp_path: Path,
    *,
    name: str,
    rule_sha: str,
    rejects: set[int],
    missing_last: bool = False,
) -> Path:
    path = tmp_path / f"{name}.json"
    count = 99 if missing_last else 100
    rows = []
    for i in range(count):
        reject = i in rejects
        rows.append(
            {
                "absolute_replicate_index": i,
                "statistic": 0.2 if reject else -0.01,
                "p_value": 0.01 if reject else 0.50,
                "reject": reject,
            }
        )
    payload = {
        "schema": "ttf_genetic_conditional_order_qualification_shard_v0.1",
        "status": "SYNTHETIC_QUALIFICATION_SHARD_COMPLETE",
        "cell": name,
        "rule_sha256": rule_sha,
        "geometry_manifest_sha256": "a" * 64,
        "geometry_csv_sha256": "b" * 64,
        "n_eval_species": 79,
        "rows": rows,
        "outcome_firewall": {
            "development_sequence_identity_opened": False,
            "confirmatory_sequence_identity_opened": False,
            "conditional_empirical_result_opened": False,
        },
    }
    _write(path, payload)
    return path


def test_aggregator_requires_complete_cells_and_keeps_identity_closed(tmp_path: Path) -> None:
    rule, digest = _rule(tmp_path)
    private = _shard(tmp_path, name="private_A2", rule_sha=digest, rejects=set())
    positive = _shard(
        tmp_path,
        name="same_order_A2",
        rule_sha=digest,
        rejects=set(range(100)),
    )
    output = tmp_path / "result.json"
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--rule",
            str(rule),
            "--shard",
            str(private),
            "--shard",
            str(positive),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert done.returncode == 0, done.stderr
    result = json.loads(output.read_text())
    assert result["status"] == "PASS"
    assert result["passed"] is True
    assert result["type1_gate"]["pass"] is True
    assert result["power_gate"]["pass"] is True
    assert result["confirmatory_identity_opening_eligible"] is False
    assert all(value is False for value in result["outcome_firewall"].values())


def test_aggregator_rejects_missing_replicate_and_rule_hash_drift(tmp_path: Path) -> None:
    rule, digest = _rule(tmp_path)
    incomplete = _shard(
        tmp_path,
        name="private_A2",
        rule_sha=digest,
        rejects=set(),
        missing_last=True,
    )
    positive = _shard(
        tmp_path,
        name="same_order_A2",
        rule_sha=digest,
        rejects=set(range(100)),
    )
    output = tmp_path / "missing.json"
    failed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--rule",
            str(rule),
            "--shard",
            str(incomplete),
            "--shard",
            str(positive),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert failed.returncode != 0
    assert not output.exists()

    private = _shard(tmp_path, name="private_A2", rule_sha="c" * 64, rejects=set())
    output = tmp_path / "hash.json"
    failed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--rule",
            str(rule),
            "--shard",
            str(private),
            "--shard",
            str(positive),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert failed.returncode != 0
    assert "rule hash" in failed.stderr
    assert not output.exists()
