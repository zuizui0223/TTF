from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Mapping

from .precision import wilson_interval


_PHASE3_RULE_SCHEMA = "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1"
_FORMAL_QUALIFICATION_SCHEMA = "ttf_genetic_phylogatr_phase3_qualification_v0.1"


def _load(path: Path, schema: str) -> dict:
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _close(observed: object, expected: float, label: str) -> None:
    try:
        value = float(observed)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"formal qualification lacks numeric {label}") from exc
    if not math.isclose(value, float(expected), rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(
            f"formal qualification {label} drift: {value!r} != {float(expected)!r}"
        )


def validate_formal_phase3_qualification(
    phase3_rule_path: Path,
    formal_qualification_path: Path,
    *,
    expected_geometry_fingerprint_sha256: str,
) -> dict:
    """Recompute the formal receipt's synthetic Gate-D summary from frozen counts.

    This validation is response-blind: it reads only the already-completed synthetic
    Phase-3 qualification receipt. It prevents a same-schema or hand-edited receipt
    from becoming the retention=1.0 anchor of the descriptive fragility curve.
    """
    rule = _load(phase3_rule_path, _PHASE3_RULE_SCHEMA)
    formal = _load(formal_qualification_path, _FORMAL_QUALIFICATION_SCHEMA)

    if formal.get("geometry_fingerprint_sha256") != str(
        expected_geometry_fingerprint_sha256
    ):
        raise RuntimeError("formal qualification geometry fingerprint drift")
    for key in (
        "confirmatory_sequence_identity_opened",
        "confirmatory_pairwise_genetic_distances_opened",
        "confirmatory_ttf_statistic_opened",
    ):
        if formal.get(key) is not False:
            raise RuntimeError(f"formal qualification firewall is open for {key}")

    qualification = rule["qualification"]
    expected_n = int(qualification["observed_worlds_per_cell"])
    expected_cells = [
        (float(cell[0]), float(cell[1]))
        for cell in qualification["mandatory_primary_cells"]
    ]
    cells = formal.get("cells")
    if not isinstance(cells, list) or len(cells) != len(expected_cells):
        raise RuntimeError("formal qualification mandatory-cell count drift")

    indexed: dict[tuple[float, float], Mapping[str, object]] = {}
    for cell in cells:
        if not isinstance(cell, Mapping):
            raise RuntimeError("formal qualification contains a non-object cell")
        key = (float(cell["shared_fraction"]), float(cell["residual_amplitude"]))
        if key in indexed:
            raise RuntimeError(f"duplicate formal qualification cell {key}")
        indexed[key] = cell
    if set(indexed) != set(expected_cells):
        raise RuntimeError("formal qualification mandatory-cell identity drift")

    recomputed: dict[tuple[float, float], tuple[float, float, float]] = {}
    for key in expected_cells:
        cell = indexed[key]
        n_worlds = int(cell["n_worlds"])
        rejections = int(cell["rejections"])
        if n_worlds != expected_n:
            raise RuntimeError(f"formal qualification world-count drift for cell {key}")
        if not 0 <= rejections <= n_worlds:
            raise RuntimeError(f"invalid formal qualification rejection count for cell {key}")
        rate = rejections / n_worlds
        interval = wilson_interval(rejections, n_worlds)
        _close(cell.get("rejection_rate"), rate, f"rejection_rate for cell {key}")
        _close(cell.get("wilson95_low"), interval.low, f"Wilson lower for cell {key}")
        _close(cell.get("wilson95_high"), interval.high, f"Wilson upper for cell {key}")
        recomputed[key] = (rate, float(interval.low), float(interval.high))

    type1_ceiling = float(qualification["type1_wilson95_upper_ceiling"])
    power_floor = float(qualification["shared_A2_wilson95_lower_floor"])
    private_cells = [key for key in expected_cells if key[0] == 0.0]
    max_private_upper = max(recomputed[key][2] for key in private_cells)
    type1_pass = all(recomputed[key][2] <= type1_ceiling for key in private_cells)
    shared_a2_key = (1.0, 2.0)
    if shared_a2_key not in recomputed:
        raise RuntimeError("formal qualification lacks frozen shared-A2 power cell")
    shared_a2_lower = recomputed[shared_a2_key][1]
    power_pass = shared_a2_lower >= power_floor
    passed = bool(type1_pass and power_pass)

    type1 = formal.get("type1_gate")
    power = formal.get("power_gate")
    if not isinstance(type1, Mapping) or not isinstance(power, Mapping):
        raise RuntimeError("formal qualification gate summaries missing")
    _close(type1.get("wilson95_upper_ceiling"), type1_ceiling, "type-I ceiling")
    _close(
        type1.get("max_observed_wilson95_upper"),
        max_private_upper,
        "max private-null Wilson upper",
    )
    if type1.get("pass") is not type1_pass:
        raise RuntimeError("formal qualification type-I pass boolean drift")
    _close(power.get("wilson95_lower_floor"), power_floor, "shared-A2 power floor")
    _close(
        power.get("observed_wilson95_lower"),
        shared_a2_lower,
        "shared-A2 Wilson lower",
    )
    if power.get("pass") is not power_pass:
        raise RuntimeError("formal qualification power pass boolean drift")

    expected_status = "PASS" if passed else "NOT_EVALUABLE"
    if formal.get("passed") is not passed:
        raise RuntimeError("formal qualification final passed boolean drift")
    if formal.get("status") != expected_status:
        raise RuntimeError("formal qualification final status drift")
    if formal.get("phase4_identity_opening_eligible") is not passed:
        raise RuntimeError("formal qualification Phase-4 eligibility drift")

    return formal


__all__ = ["validate_formal_phase3_qualification"]
