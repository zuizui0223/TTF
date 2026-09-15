from __future__ import annotations

import json
from pathlib import Path

import pytest

from ttf.phylogatr_fragility_formal_anchor import validate_formal_phase3_qualification
from ttf.precision import wilson_interval


PHASE3_RULE = Path("docs/supporting/genetic_phylogatr_phase3_gate_d_rule_v0.1.json")
FINGERPRINT = "f" * 64


def _formal_receipt(tmp_path: Path) -> Path:
    counts = {
        (0.0, 0.5): 20,
        (0.0, 1.0): 22,
        (0.0, 2.0): 25,
        (0.0, 3.0): 37,
        (1.0, 2.0): 434,
    }
    cells = []
    for shared, amplitude in [(0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0)]:
        rejections = counts[(shared, amplitude)]
        interval = wilson_interval(rejections, 500)
        cells.append(
            {
                "shared_fraction": shared,
                "residual_amplitude": amplitude,
                "n_worlds": 500,
                "rejections": rejections,
                "rejection_rate": rejections / 500,
                "wilson95_low": float(interval.low),
                "wilson95_high": float(interval.high),
                "selected_pair_counts": {},
            }
        )

    a3 = next(
        cell
        for cell in cells
        if cell["shared_fraction"] == 0.0 and cell["residual_amplitude"] == 3.0
    )
    shared_a2 = next(
        cell
        for cell in cells
        if cell["shared_fraction"] == 1.0 and cell["residual_amplitude"] == 2.0
    )
    payload = {
        "schema": "ttf_genetic_phylogatr_phase3_qualification_v0.1",
        "status": "NOT_EVALUABLE",
        "geometry_fingerprint_sha256": FINGERPRINT,
        "cells": cells,
        "type1_gate": {
            "wilson95_upper_ceiling": 0.10,
            "max_observed_wilson95_upper": a3["wilson95_high"],
            "pass": False,
        },
        "power_gate": {
            "cell": {"shared_fraction": 1.0, "residual_amplitude": 2.0},
            "wilson95_lower_floor": 0.80,
            "observed_wilson95_lower": shared_a2["wilson95_low"],
            "pass": True,
        },
        "passed": False,
        "failure_interpretation": "NOT_EVALUABLE_under_this_fresh_geometry_not_absence_of_transferable_phylogeography",
        "phase4_identity_opening_eligible": False,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
    }
    path = tmp_path / "formal.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def test_formal_anchor_recomputes_narrow_decker_like_boundary(tmp_path: Path) -> None:
    path = _formal_receipt(tmp_path)
    payload = validate_formal_phase3_qualification(
        PHASE3_RULE,
        path,
        expected_geometry_fingerprint_sha256=FINGERPRINT,
    )
    assert payload["status"] == "NOT_EVALUABLE"
    assert payload["type1_gate"]["max_observed_wilson95_upper"] == pytest.approx(
        0.10033475332223055
    )
    assert payload["power_gate"]["observed_wilson95_lower"] == pytest.approx(
        0.8355052092686175
    )


def test_formal_anchor_rejects_hand_edited_wilson_value(tmp_path: Path) -> None:
    path = _formal_receipt(tmp_path)
    payload = json.loads(path.read_text())
    payload["cells"][3]["wilson95_high"] += 0.001
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with pytest.raises(RuntimeError, match="Wilson upper"):
        validate_formal_phase3_qualification(
            PHASE3_RULE,
            path,
            expected_geometry_fingerprint_sha256=FINGERPRINT,
        )


def test_formal_anchor_rejects_hand_edited_final_status(tmp_path: Path) -> None:
    path = _formal_receipt(tmp_path)
    payload = json.loads(path.read_text())
    payload["status"] = "PASS"
    payload["passed"] = True
    payload["phase4_identity_opening_eligible"] = True
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    with pytest.raises(RuntimeError, match="final passed boolean drift"):
        validate_formal_phase3_qualification(
            PHASE3_RULE,
            path,
            expected_geometry_fingerprint_sha256=FINGERPRINT,
        )
