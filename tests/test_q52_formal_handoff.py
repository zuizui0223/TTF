from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_ttfm_q5_2_formal_batch.py"
spec = importlib.util.spec_from_file_location("run_ttfm_q5_2_formal_batch", SCRIPT)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
selected_candidate = module.selected_candidate


def _protocol() -> dict:
    return {
        "status": "frozen_before_q5_2_development_aggregate_and_formal_outcomes",
        "candidate_handoff": {
            "allowed_values": [
                "pooled_inverse_density",
                "self_inverse_density",
                "target_density_ratio",
            ]
        },
    }


def test_null_development_selection_forbids_formal_run() -> None:
    development = {
        "schema": "ttfm_q5_2_geometry_transport_development_result_v0.1",
        "candidate_selected": None,
        "candidate_checks": {},
    }
    with pytest.raises(ValueError, match="formal qualification is forbidden"):
        selected_candidate(_protocol(), development)


def test_selected_candidate_must_have_development_pass() -> None:
    development = {
        "schema": "ttfm_q5_2_geometry_transport_development_result_v0.1",
        "candidate_selected": "pooled_inverse_density",
        "candidate_checks": {
            "pooled_inverse_density": {"development_pass": False}
        },
    }
    with pytest.raises(ValueError, match="did not pass"):
        selected_candidate(_protocol(), development)


def test_selected_candidate_handoff_is_exact() -> None:
    development = {
        "schema": "ttfm_q5_2_geometry_transport_development_result_v0.1",
        "candidate_selected": "self_inverse_density",
        "candidate_checks": {
            "self_inverse_density": {"development_pass": True}
        },
    }
    assert selected_candidate(_protocol(), development) == "self_inverse_density"
