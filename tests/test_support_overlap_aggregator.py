from __future__ import annotations

import importlib.util
import json
from pathlib import Path


def _load_aggregator():
    path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "aggregate_ttfc_support_overlap.py"
    )
    spec = importlib.util.spec_from_file_location("support_overlap_aggregate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rule():
    ids = [
        ("DSO_NULL_COMPONENT", "null_component", "matched"),
        ("DSO_PRIVATE_M", "null_private_mismatch", "matched"),
        ("DSO_PRIVATE_C", "null_private_relation", "matched"),
        ("DSO_POWER_M_MATCHED", "power_shared_mismatch", "matched"),
        ("DSO_POWER_C_MATCHED", "power_shared_relation", "matched"),
        ("DSO_POWER_M_SHIFTED", "power_shared_mismatch", "shifted"),
        ("DSO_POWER_C_SHIFTED", "power_shared_relation", "shifted"),
    ]
    return {
        "status": "frozen_before_support_overlap_diagnostic_outcomes",
        "diagnostic_worlds": {
            "worlds_per_cell": 8,
            "aggregate_bootstrap_resamples": 99,
            "aggregate_bootstrap_seed": 12345,
        },
        "cells": [
            {"id": cell_id, "response_mode": mode, "geometry_profile": profile}
            for cell_id, mode, profile in ids
        ],
        "terminal_rule": {"candidate_selected": None, "formal_pass_fail": None},
        "claim_firewall": {"diagnostic_only": True, "no_repair_selection": True},
    }


def _batch(cell_id: str, response_mode: str, geometry_profile: str):
    shifted = geometry_profile == "shifted"
    private_c = cell_id == "DSO_PRIVATE_C"
    private_m = cell_id == "DSO_PRIVATE_M"
    rows = []
    for rep in range(8):
        x = (rep + 1) / 9.0
        alignment_c = x
        alignment_m = 1.0 - x
        raw_c = 0.2 + 0.5 * x
        raw_m = 0.1 + 0.1 * x
        partial_c = raw_c - (0.08 + 0.02 * x if private_c else 0.01)
        partial_m = raw_m - (0.01 if private_c else 0.03 if private_m else 0.005)
        if cell_id == "DSO_POWER_M_MATCHED":
            raw_c = 0.70 + 0.01 * rep
        if cell_id == "DSO_POWER_M_SHIFTED":
            raw_c = 0.45 + 0.01 * rep
        rows.append(
            {
                "replicate": rep,
                "world": {
                    "raw_M": raw_m,
                    "partial_M": partial_m,
                    "alignment_M": alignment_m,
                    "raw_C": raw_c,
                    "partial_C": partial_c,
                    "alignment_C": alignment_c,
                    "mean_opportunity": 0.8 - (0.2 if shifted else 0.0) + 0.001 * rep,
                    "mean_prior_fraction": 0.55 if shifted else 0.25,
                    "low_support_fraction": 0.65 if shifted else 0.20,
                },
            }
        )
    return {
        "schema": "ttfc_support_overlap_diagnostic_batch_v0.1",
        "cell_id": cell_id,
        "response_mode": response_mode,
        "geometry_profile": geometry_profile,
        "candidate_selected": None,
        "formal_pass_fail": None,
        "rows": rows,
    }


def test_support_overlap_aggregator_executes_without_candidate_selection(tmp_path):
    module = _load_aggregator()
    rule = _rule()
    paths = []
    for cell in rule["cells"]:
        payload = _batch(
            cell["id"], cell["response_mode"], cell["geometry_profile"]
        )
        path = tmp_path / f"{cell['id']}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        paths.append(path)

    result = module.aggregate(rule, paths)

    assert result["status"] == "diagnostic_complete_nonqualifying"
    assert result["candidate_selected"] is None
    assert result["formal_pass_fail"] is None
    assert result["global_gate_predeclared"] is False
    assert set(result["cell_summaries"]) == {cell["id"] for cell in rule["cells"]}
    assert "private_relation_opportunity_mediation" in result["mechanism_tests"]
    assert "shifted_power_support_loss" in result["mechanism_tests"]
    assert "estimand_specificity" in result["mechanism_tests"]
    assert result["mechanism_tests"]["shifted_power_support_loss"]["supported"] is True
