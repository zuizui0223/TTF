from __future__ import annotations

import importlib.util
import json
from pathlib import Path


def _load_module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "aggregate_ttfc_training_leverage.py"
    spec = importlib.util.spec_from_file_location("training_leverage_aggregate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rule():
    cells = [
        {"id": "DTL_NULL_COMPONENT", "response_mode": "null_component", "geometry_profile": "matched"},
        {"id": "DTL_PRIVATE_M", "response_mode": "null_private_mismatch", "geometry_profile": "matched"},
        {"id": "DTL_PRIVATE_C", "response_mode": "null_private_relation", "geometry_profile": "matched"},
        {"id": "DTL_POWER_M_MATCHED", "response_mode": "power_shared_mismatch", "geometry_profile": "matched"},
        {"id": "DTL_POWER_C_MATCHED", "response_mode": "power_shared_relation", "geometry_profile": "matched"},
        {"id": "DTL_POWER_M_SHIFTED", "response_mode": "power_shared_mismatch", "geometry_profile": "shifted"},
        {"id": "DTL_POWER_C_SHIFTED", "response_mode": "power_shared_relation", "geometry_profile": "shifted"},
    ]
    return {
        "status": "frozen_before_training_leverage_diagnostic_outcomes",
        "diagnostic_worlds": {
            "worlds_per_cell": 4,
            "aggregate_bootstrap_resamples": 99,
            "aggregate_bootstrap_seed": 12345,
        },
        "cells": cells,
        "terminal_rule": {"candidate_selected": None, "formal_pass_fail": None},
        "claim_firewall": {"diagnostic_only": True, "no_repair_selection": True},
    }


def _world(rep: int, shifted: bool = False):
    x = float(rep + 1)
    return {
        "raw_M": 0.01 * x,
        "raw_C": 0.02 * x - (0.01 if shifted else 0.0),
        "loo_score_sd_M": 0.002 + 0.0001 * x,
        "loo_score_sd_C": 0.004 + 0.0003 * x + (0.002 if shifted else 0.0),
        "loo_score_mean_M": 0.01 * x,
        "loo_score_mean_C": 0.02 * x,
        "loo_score_min_M": 0.01 * x - 0.003,
        "loo_score_min_C": 0.02 * x - 0.006,
        "mean_abs_loo_score_shift_M": 0.001,
        "mean_abs_loo_score_shift_C": 0.002,
        "mean_prediction_loo_sd_M": 0.001 + 0.0001 * x,
        "mean_prediction_loo_sd_C": 0.002 + 0.0002 * x,
        "mean_effective_training_system_count": 8.0 - 0.1 * x,
        "mean_dominant_training_system_share": 0.20 + 0.01 * x,
    }


def test_training_leverage_aggregator_executes_without_selecting_candidate(tmp_path):
    module = _load_module()
    rule = _rule()
    paths = []
    for cell in rule["cells"]:
        rows = []
        for rep in range(4):
            rows.append({"replicate": rep, "world": _world(rep, shifted=cell["geometry_profile"] == "shifted")})
        payload = {
            "schema": "ttfc_training_leverage_diagnostic_batch_v0.1",
            "cell_id": cell["id"],
            "rows": rows,
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
        "private_relation_field_leverage",
        "private_relation_geometry_concentration",
        "private_mismatch_field_leverage",
        "shifted_power_training_instability",
        "estimand_specificity",
        "shared_relation_negative_control",
    }
    for name in (
        "private_relation_field_leverage",
        "private_relation_geometry_concentration",
        "private_mismatch_field_leverage",
        "shifted_power_training_instability",
        "estimand_specificity",
    ):
        assert isinstance(result["mechanism_tests"][name]["supported"], bool)
