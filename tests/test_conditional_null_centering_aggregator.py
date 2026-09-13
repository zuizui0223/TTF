from __future__ import annotations

import importlib.util
import json
from pathlib import Path


def _load_module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "aggregate_ttfc_conditional_null_centering.py"
    spec = importlib.util.spec_from_file_location("conditional_null_centering_aggregate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rule():
    return {
        "status": "frozen_before_conditional_null_centering_outcomes",
        "diagnostic_worlds": {
            "geometry_worlds_per_cell": 4,
            "aggregate_bootstrap_resamples": 99,
            "aggregate_bootstrap_seed": 4401,
        },
        "cells": [
            {"id": "DCN_HOMOGENEOUS", "geometry_mode": "identical_uniform60"},
            {"id": "DCN_MATCHED", "geometry_mode": "q5_matched"},
            {"id": "DCN_SHIFTED", "geometry_mode": "q5_shifted"},
        ],
        "terminal_rule": {"candidate_selected": None, "formal_pass_fail": None},
        "claim_firewall": {"diagnostic_only": True, "no_repair_selection": True},
    }


def _row(rep: int, base: float):
    x = float(rep + 1)
    return {
        "replicate": rep,
        "conditional_mean_C": base + 0.001 * x,
        "within_geometry_phase_sd_C": 0.02 + 0.001 * x,
        "conditional_mean_mc_se_C": 0.002 + 0.0001 * x,
        "phase_fraction_positive_C": 0.5 + 0.01 * x,
        "phase_q05_C": base - 0.03,
        "phase_q50_C": base,
        "phase_q95_C": base + 0.03,
    }


def test_conditional_null_aggregator_executes_without_selection(tmp_path):
    module = _load_module()
    rule = _rule()
    bases = {
        "DCN_HOMOGENEOUS": 0.0,
        "DCN_MATCHED": 0.03,
        "DCN_SHIFTED": 0.04,
    }
    paths = []
    for cell in rule["cells"]:
        payload = {
            "schema": "ttfc_conditional_null_centering_batch_v0.1",
            "cell_id": cell["id"],
            "rows": [_row(rep, bases[cell["id"]]) for rep in range(4)],
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
    assert result["mechanism_tests"]["matched_conditional_null_bias"]["supported"] is True
    assert result["mechanism_tests"]["matched_heterogeneity_amplification"]["supported"] is True
    assert result["mechanism_tests"]["primary_geometry_conditioned_null_bias"]["supported"] is True
