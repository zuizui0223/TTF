from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path


PLANNER = Path("scripts/plan_phylogatr_gate_d_fragility_synthetic.py")
AGGREGATOR = Path("scripts/aggregate_phylogatr_gate_d_fragility_curve.py")
PHASE3_RULE = Path("docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json")


def _load_planner_module():
    spec = importlib.util.spec_from_file_location("fragility_planner_script", PLANNER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _anchor() -> dict:
    return {
        "retention_fraction": 1.0,
        "source": "formal_phase3_qualification_receipt",
        "formal_status": "PASS",
        "A3_private_null_rejection_rate": 0.06,
        "A3_private_null_Wilson95_upper": 0.099,
        "max_private_null_Wilson95_upper": 0.099,
        "shared_A2_rejection_rate": 0.86,
        "shared_A2_Wilson95_lower": 0.82,
        "distance_of_max_private_null_upper_to_formal_0.10_ceiling": -0.001,
        "distance_of_shared_A2_lower_to_formal_0.80_floor": 0.02,
        "formal_gate_d_decision_authority": True,
        "diagnostic_decision_authority": False,
    }


def test_identical_geometry_policy_creates_zero_job_plateaus() -> None:
    planner = _load_planner_module()
    fingerprint = "f" * 64
    plan = {
        "full_geometry_fingerprint_sha256": fingerprint,
        "levels": [
            {
                "retention_fraction": retention,
                "geometry_fingerprint_sha256": fingerprint,
                "status": "SYNTHETIC_DIAGNOSTIC_RUNNABLE",
                "reference_jobs": [{"configuration": "A0", "start": 0, "stop": 1}],
                "observed_jobs": [{"shared_fraction": 0.0, "residual_amplitude": 0.5, "start": 0, "stop": 1}],
                "reference_job_count": 1,
                "observed_job_count": 1,
            }
            for retention in (0.875, 0.75, 0.625, 0.5)
        ],
    }
    out = planner._apply_identical_geometry_policy(plan)
    assert [level["status"] for level in out["levels"]] == [
        "IDENTICAL_TO_FULL_GEOMETRY_NO_RERUN",
        "IDENTICAL_TO_FULL_GEOMETRY_NO_RERUN",
        "IDENTICAL_TO_FULL_GEOMETRY_NO_RERUN",
        "IDENTICAL_TO_FULL_GEOMETRY_NO_RERUN",
    ]
    assert all(level["inherits_metrics_from_retention_fraction"] == 1.0 for level in out["levels"])
    assert all(level["reference_job_count"] == 0 for level in out["levels"])
    assert all(level["observed_job_count"] == 0 for level in out["levels"])
    assert all(level["reference_jobs"] == [] for level in out["levels"])
    assert all(level["observed_jobs"] == [] for level in out["levels"])


def test_curve_aggregator_reuses_anchor_without_synthetic_files(tmp_path: Path) -> None:
    rule = json.loads(PHASE3_RULE.read_text())
    fingerprint = "f" * 64
    levels = []
    for retention in (0.875, 0.75, 0.625, 0.5):
        levels.append(
            {
                "retention_fraction": retention,
                "status": "IDENTICAL_TO_FULL_GEOMETRY_NO_RERUN",
                "inherits_metrics_from_retention_fraction": 1.0,
                "geometry_fingerprint_sha256": fingerprint,
                "metrics": {"minimum_endpoint_disjoint_ibd_training_edges": 10},
                "reference_jobs": [],
                "observed_jobs": [],
                "reference_job_count": 0,
                "observed_job_count": 0,
                "formal_gate_d_decision_authority": False,
                "phase4_identity_opening_authority": False,
            }
        )
    plan = {
        "schema": "ttf_genetic_phylogatr_gate_d_fragility_execution_plan_v0.1",
        "status": "DIAGNOSTIC_SYNTHETIC_EXECUTION_PLAN_ONLY",
        "phase3_rule_sha256": __import__("hashlib").sha256(PHASE3_RULE.read_bytes()).hexdigest(),
        "formal_qualification_sha256": "q" * 64,
        "execution_rule_sha256": "e" * 64,
        "dataset_digest_sha256": "d" * 64,
        "full_geometry_fingerprint_sha256": fingerprint,
        "formal_anchor": _anchor(),
        "levels": levels,
        "authority_firewall": {
            "diagnostic_can_change_formal_gate_d_status": False,
            "diagnostic_can_authorize_phase4": False,
            "diagnostic_can_rescue_not_evaluable": False,
            "diagnostic_can_reclassify_not_evaluable_as_negative": False,
            "diagnostic_can_change_thresholds": False,
            "diagnostic_can_change_world_counts": False,
            "diagnostic_can_change_species_or_split": False,
            "diagnostic_can_change_bandwidth_or_graph_rule": False,
        },
    }
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    empty_input = tmp_path / "empty"
    empty_input.mkdir()
    output = tmp_path / "curve.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(AGGREGATOR),
            "--input-dir",
            str(empty_input),
            "--execution-plan",
            str(plan_path),
            "--phase3-rule",
            str(PHASE3_RULE),
            "--output",
            str(output),
        ],
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0, completed.stderr
    curve = json.loads(output.read_text())
    assert curve["status"] == "DIAGNOSTIC_FRAGILITY_CURVE_COMPLETE"
    assert [row["retention_fraction"] for row in curve["curve"]] == [1.0, 0.875, 0.75, 0.625, 0.5]
    anchor = curve["curve"][0]
    for row in curve["curve"][1:]:
        assert row["status"] == "DIAGNOSTIC_SYNTHETIC_SUMMARY"
        assert row["execution_status"] == "IDENTICAL_TO_FULL_GEOMETRY_NO_RERUN"
        assert row["source"] == "diagnostic_identical_geometry_plateau"
        assert row["inherits_metrics_from_retention_fraction"] == 1.0
        assert row["A3_private_null_Wilson95_upper"] == anchor["A3_private_null_Wilson95_upper"]
        assert row["max_private_null_Wilson95_upper"] == anchor["max_private_null_Wilson95_upper"]
        assert row["shared_A2_Wilson95_lower"] == anchor["shared_A2_Wilson95_lower"]
