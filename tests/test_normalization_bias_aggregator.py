from __future__ import annotations

import importlib.util
import json
from pathlib import Path


def _load_module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "aggregate_ttfc_normalization_bias.py"
    spec = importlib.util.spec_from_file_location("normalization_bias_aggregate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rule():
    return {
        "status": "frozen_before_normalization_bias_outcomes",
        "fresh_worlds": {
            "worlds_per_cell": 4,
            "aggregate_bootstrap_resamples": 99,
            "aggregate_bootstrap_seed": 117,
            "cells": [
                {"id": "DNB_HOMOGENEOUS", "geometry_mode": "identical_uniform60", "role": "negative_control"},
                {"id": "DNB_MATCHED", "geometry_mode": "q5_matched", "role": "primary"},
                {"id": "DNB_SHIFTED", "geometry_mode": "q5_shifted", "role": "secondary"},
            ],
        },
        "terminal_rule": {"candidate_selected": None, "formal_pass_fail": None},
        "claim_firewall": {"diagnostic_only": True, "no_repair_selection": True},
    }


def _row(rep: int, cell_id: str):
    jitter = 0.001 * rep
    if cell_id == "DNB_HOMOGENEOUS":
        normalized = -0.0015 + jitter
        unscaled = -0.001 + 0.0007 * rep
    else:
        normalized = 0.03 + jitter
        unscaled = -0.0015 + 0.001 * rep
    return {
        "replicate": rep,
        "conditional_mean_normalized_C": normalized,
        "conditional_mean_unscaled_C": unscaled,
        "conditional_mean_normalized_minus_unscaled_C": normalized - unscaled,
        "within_geometry_phase_sd_normalized_C": 0.06,
        "within_geometry_phase_sd_unscaled_C": 0.04,
        "phase_fraction_positive_normalized_C": 0.65,
        "phase_fraction_positive_unscaled_C": 0.50,
    }


def test_normalization_bias_aggregator_is_nonqualifying_and_reports_frozen_tests(tmp_path):
    module = _load_module()
    rule = _rule()
    paths = []
    for cell in rule["fresh_worlds"]["cells"]:
        payload = {
            "schema": "ttfc_normalization_bias_diagnostic_batch_v0.1",
            "cell_id": cell["id"],
            "rows": [_row(rep, cell["id"]) for rep in range(4)],
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
        "matched_normalized_bias_replication",
        "matched_normalization_contribution",
        "matched_unscaled_recentered",
        "primary_normalization_bias_mechanism",
        "shifted_normalized_bias_replication",
        "shifted_normalization_contribution",
        "shifted_unscaled_recentered",
        "shifted_joint_normalization_bias_mechanism",
        "homogeneous_negative_control",
    }
    assert isinstance(result["mechanism_tests"]["primary_normalization_bias_mechanism"]["supported"], bool)
    assert isinstance(result["mechanism_tests"]["shifted_joint_normalization_bias_mechanism"]["supported"], bool)
