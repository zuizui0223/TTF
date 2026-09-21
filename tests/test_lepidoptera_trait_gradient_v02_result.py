from __future__ import annotations

import json
from pathlib import Path


RESULT = Path("benchmarks/frozen/lepidoptera_trait_gradient_v02_qualification_result.json")


def test_v02_trait_gradient_qualification_passes_all_frozen_gates() -> None:
    payload = json.loads(RESULT.read_text())
    assert payload["schema"] == "ttf_lepidoptera_trait_gradient_v02_qualification_result"
    assert payload["status"] == "PASS"
    assert payload["design_npz_sha256"] == "e90568668f8263e3b454ddab21a281dc6fe2e3b33bd88324b8fb1cc012082b3f"
    gates = payload["gates"]
    assert gates == {
        "all_three_pass": True,
        "geometry_trap_type1_pass": True,
        "private_type1_pass": True,
        "trait_gradient_positive_power_pass": True,
    }
    assert payload["evaluation"]["private"]["wilson95_upper"] <= 0.10
    assert payload["evaluation"]["geometry_confounded_trap"]["wilson95_upper"] <= 0.10
    assert payload["evaluation"]["trait_gradient_positive"]["wilson95_lower"] >= 0.80
    assert payload["empirical_opening_authorized"] is True
    assert all(value is False for value in payload["outcome_firewall"].values())
