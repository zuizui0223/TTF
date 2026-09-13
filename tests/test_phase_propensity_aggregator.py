from __future__ import annotations

import importlib.util
import json
from pathlib import Path


def _load_module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "aggregate_ttfc_phase_propensity.py"
    spec = importlib.util.spec_from_file_location("phase_propensity_aggregate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rule():
    return {
        "status": "frozen_before_phase_propensity_outcomes",
        "fresh_worlds": {
            "worlds_per_cell": 4,
            "aggregate_bootstrap_resamples": 99,
            "aggregate_bootstrap_seed": 1234,
            "cells": [
                {"id": "DPP_HOMOGENEOUS", "geometry_mode": "identical_uniform60", "role": "negative_control"},
                {"id": "DPP_MATCHED", "geometry_mode": "q5_matched", "role": "primary"},
                {"id": "DPP_SHIFTED", "geometry_mode": "q5_shifted", "role": "secondary"},
            ],
        },
        "terminal_rule": {"candidate_selected": None, "formal_pass_fail": None},
        "claim_firewall": {"diagnostic_only": True, "no_repair_selection": True},
    }


def _row(rep: int, *, mode: str):
    x = float(rep + 1)
    if mode == "homogeneous":
        baseline = 0.0002 * (x - 2.5)
        centered = baseline
        rms = 0.001
    else:
        baseline = 0.02 + 0.001 * x
        centered = 0.0002 * (x - 2.5)
        rms = 0.02
    return {
        "replicate": rep,
        "conditional_mean_baseline_C": baseline,
        "conditional_mean_propensity_centered_C": centered,
        "conditional_mean_baseline_minus_centered_C": baseline - centered,
        "within_geometry_phase_sd_baseline_C": 0.05,
        "within_geometry_phase_sd_propensity_centered_C": 0.05,
        "phase_fraction_positive_baseline_C": 0.6,
        "phase_fraction_positive_propensity_centered_C": 0.5,
        "mean_training_expected_residual_rms": rms,
    }


def test_phase_propensity_aggregator_executes_without_selecting_candidate(tmp_path):
    module = _load_module()
    rule = _rule()
    paths = []
    for cell in rule["fresh_worlds"]["cells"]:
        mode = "homogeneous" if cell["id"] == "DPP_HOMOGENEOUS" else "heterogeneous"
        payload = {
            "schema": "ttfc_phase_propensity_diagnostic_batch_v0.1",
            "cell_id": cell["id"],
            "rows": [_row(rep, mode=mode) for rep in range(4)],
            "candidate_selected": None,
            "formal_pass_fail": None,
        }
        path = tmp_path / f"{cell['id']}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        paths.append(path)

    result = module.aggregate(rule, paths)
    assert result["status"] == "diagnostic_complete_nonqualifying"
    assert result["candidate_selected"] is None
    assert result["formal_pass_fail"] is None
    assert result["global_gate_predeclared"] is False
    assert set(result["mechanism_tests"]) == {
        "matched_baseline_bias_replication",
        "matched_phase_propensity_contribution",
        "matched_propensity_centered_recentered",
        "primary_phase_propensity_mechanism",
        "shifted_baseline_bias_replication",
        "shifted_phase_propensity_contribution",
        "shifted_propensity_centered_recentered",
        "shifted_joint_phase_propensity_mechanism",
        "homogeneous_negative_control",
    }
    assert result["mechanism_tests"]["primary_phase_propensity_mechanism"]["supported"] is True
    assert result["mechanism_tests"]["shifted_joint_phase_propensity_mechanism"]["supported"] is True
