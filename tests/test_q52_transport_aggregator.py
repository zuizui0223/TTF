from __future__ import annotations

import importlib.util
import json
from pathlib import Path


def _load_aggregator():
    path = Path(__file__).resolve().parents[1] / "scripts" / "aggregate_ttfm_q5_2_transport_development.py"
    spec = importlib.util.spec_from_file_location("q52_dev_aggregate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rule():
    return {
        "status": "frozen_before_q5_2_transport_development_outcomes",
        "development_inference": {"worlds_per_cell": 2, "alpha": 0.05},
        "candidate_preference_order": [
            "pooled_inverse_density",
            "self_inverse_density",
            "target_density_ratio",
        ],
        "development_cells": [
            {"id": "N", "response_mode": "null_component", "geometry_profile": "matched", "role": "type1"},
            {"id": "PM", "response_mode": "power_shared_mismatch", "geometry_profile": "matched", "role": "shared_mismatch_matched_power"},
            {"id": "PS", "response_mode": "power_shared_mismatch", "geometry_profile": "shifted", "role": "shared_mismatch_shifted_power"},
            {"id": "RC", "response_mode": "power_shared_relation", "geometry_profile": "matched", "role": "shared_relation_matched_power"},
            {"id": "RS", "response_mode": "power_shared_relation", "geometry_profile": "shifted", "role": "shared_relation_shifted_power"},
        ],
        "selection_rule": {
            "type1_max_rate": 0.08,
            "shared_mismatch_matched_min_power": 0.85,
            "shared_mismatch_shifted_min_power": 0.82,
            "shared_relation_matched_min_power": 0.95,
            "shared_relation_shifted_min_power": 0.95,
        },
        "next_if_candidate_selected": "qualify_fresh",
        "next_if_no_candidate": "stop_family",
        "claim_firewall": {"development_only": True},
    }


def _batch(cell, role, response_mode, geometry_profile, p_by_mode):
    modes = [
        "q5_1_uniform",
        "pooled_inverse_density",
        "self_inverse_density",
        "target_density_ratio",
    ]
    rows = []
    for rep in range(2):
        rows.append(
            {
                "replicate": rep,
                "modes": {
                    mode: {
                        "p_value": p_by_mode[mode],
                        "statistic": 0.25,
                        "mean_train_effective_edge_count": 12.0,
                        "min_train_effective_edge_count": 8.0,
                    }
                    for mode in modes
                },
            }
        )
    return {
        "schema": "ttfm_q5_2_geometry_transport_development_batch_v0.1",
        "cell_id": cell,
        "role": role,
        "response_mode": response_mode,
        "geometry_profile": geometry_profile,
        "rows": rows,
    }


def test_q52_aggregator_executes_and_returns_python_booleans(tmp_path):
    module = _load_aggregator()
    rule = _rule()
    paths = []
    for cell in rule["development_cells"]:
        if cell["role"] == "type1":
            p = {
                "q5_1_uniform": 0.5,
                "pooled_inverse_density": 0.5,
                "self_inverse_density": 0.5,
                "target_density_ratio": 0.5,
            }
        else:
            p = {
                "q5_1_uniform": 0.01,
                "pooled_inverse_density": 0.01,
                "self_inverse_density": 0.01,
                "target_density_ratio": 0.01,
            }
        payload = _batch(
            cell["id"],
            cell["role"],
            cell["response_mode"],
            cell["geometry_profile"],
            p,
        )
        path = tmp_path / f"{cell['id']}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        paths.append(path)

    result = module.aggregate(rule, paths)
    assert result["formal_q5_2_qualified"] is False
    assert result["q5_remains_immutable_fail"] is True
    assert result["q5_1_remains_immutable_fail"] is True
    assert result["candidate_selected"] == "pooled_inverse_density"
    assert result["status"] == "development_complete_nonqualifying"
