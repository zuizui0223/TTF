#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from ttf.genetic_geometry_io import sha256_path
from ttf.precision import wilson_interval


def _collect_level(
    input_dir: Path,
    *,
    retention: float,
    level: dict,
    rule: dict,
) -> dict:
    qualification = rule["qualification"]
    expected_cells = [
        (float(cell[0]), float(cell[1]))
        for cell in qualification["mandatory_primary_cells"]
    ]
    expected_n = int(qualification["observed_worlds_per_cell"])
    alpha = float(qualification["alpha"])
    grouped: dict[tuple[float, float], dict[int, dict]] = {cell: {} for cell in expected_cells}

    for path in sorted(input_dir.rglob("*.json")):
        shard = json.loads(path.read_text())
        if shard.get("schema") != "ttf_genetic_phylogatr_gate_d_fragility_observed_shard_v0.1":
            continue
        if float(shard["retention_fraction"]) != retention:
            continue
        if shard.get("geometry_fingerprint_sha256") != level["geometry_fingerprint_sha256"]:
            raise RuntimeError(f"fragility observed geometry fingerprint drift in {path}")
        for key in (
            "confirmatory_sequence_identity_opened",
            "confirmatory_pairwise_genetic_distances_opened",
            "confirmatory_ttf_statistic_opened",
            "formal_gate_d_decision_authority",
            "phase4_identity_opening_authority",
        ):
            if shard.get(key) is not False:
                raise RuntimeError(f"fragility observed authority/firewall drift in {path}: {key}")
        cell = (float(shard["shared_fraction"]), float(shard["residual_amplitude"]))
        if cell not in grouped:
            raise RuntimeError(f"unexpected fragility observed cell {cell}")
        for row in shard["rows"]:
            replicate = int(row["replicate"])
            if replicate in grouped[cell]:
                raise RuntimeError(f"duplicate fragility observed replicate {cell}:{replicate}")
            grouped[cell][replicate] = row

    cells: list[dict] = []
    for shared, amplitude in expected_cells:
        rows = grouped[(shared, amplitude)]
        if sorted(rows) != list(range(expected_n)):
            missing = sorted(set(range(expected_n)) - set(rows))[:10]
            raise RuntimeError(
                f"incomplete fragility observed cell retention={retention} "
                f"cell={(shared, amplitude)}; first missing={missing}"
            )
        ordered = [rows[index] for index in range(expected_n)]
        rejections = sum(float(row["p_value"]) <= alpha for row in ordered)
        interval = wilson_interval(rejections, expected_n)
        selected_counts = Counter(
            "+".join(row["selected_configurations"]) for row in ordered
        )
        cells.append(
            {
                "shared_fraction": shared,
                "residual_amplitude": amplitude,
                "n_worlds": expected_n,
                "rejections": int(rejections),
                "rejection_rate": float(rejections / expected_n),
                "wilson95_low": float(interval.low),
                "wilson95_high": float(interval.high),
                "selected_pair_counts": dict(sorted(selected_counts.items())),
            }
        )

    private_cells = [cell for cell in cells if cell["shared_fraction"] == 0.0]
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
    max_private_upper = max(float(cell["wilson95_high"]) for cell in private_cells)
    type1_ceiling = float(qualification["type1_wilson95_upper_ceiling"])
    power_floor = float(qualification["shared_A2_wilson95_lower_floor"])
    return {
        "retention_fraction": retention,
        "source": "diagnostic_thinned_geometry_synthetic_execution",
        "status": "DIAGNOSTIC_SYNTHETIC_SUMMARY",
        "geometry_fingerprint_sha256": level["geometry_fingerprint_sha256"],
        "minimum_endpoint_disjoint_ibd_training_edges": int(
            level["metrics"]["minimum_endpoint_disjoint_ibd_training_edges"]
        ),
        "A3_private_null_rejection_rate": float(a3["rejection_rate"]),
        "A3_private_null_Wilson95_upper": float(a3["wilson95_high"]),
        "max_private_null_Wilson95_upper": max_private_upper,
        "shared_A2_rejection_rate": float(shared_a2["rejection_rate"]),
        "shared_A2_Wilson95_lower": float(shared_a2["wilson95_low"]),
        "distance_of_max_private_null_upper_to_formal_0.10_ceiling": (
            max_private_upper - type1_ceiling
        ),
        "distance_of_shared_A2_lower_to_formal_0.80_floor": (
            float(shared_a2["wilson95_low"]) - power_floor
        ),
        "cells": cells,
        "formal_gate_d_decision_authority": False,
        "phase4_identity_opening_authority": False,
    }


def _plateau_row(level: dict, source_row: dict) -> dict:
    fields = (
        "A3_private_null_rejection_rate",
        "A3_private_null_Wilson95_upper",
        "max_private_null_Wilson95_upper",
        "shared_A2_rejection_rate",
        "shared_A2_Wilson95_lower",
        "distance_of_max_private_null_upper_to_formal_0.10_ceiling",
        "distance_of_shared_A2_lower_to_formal_0.80_floor",
    )
    out = {
        "retention_fraction": float(level["retention_fraction"]),
        "source": "diagnostic_identical_geometry_plateau",
        "status": str(level["status"]),
        "inherits_metrics_from_retention_fraction": float(
            level["inherits_metrics_from_retention_fraction"]
        ),
        "geometry_fingerprint_sha256": level["geometry_fingerprint_sha256"],
        "minimum_endpoint_disjoint_ibd_training_edges": int(
            level["metrics"]["minimum_endpoint_disjoint_ibd_training_edges"]
        ),
        "formal_gate_d_decision_authority": False,
        "phase4_identity_opening_authority": False,
    }
    for field in fields:
        if field not in source_row:
            raise RuntimeError(f"identical-geometry plateau source lacks {field}")
        out[field] = source_row[field]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--execution-plan", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    plan = json.loads(args.execution_plan.read_text())
    rule = json.loads(args.phase3_rule.read_text())
    if plan.get("schema") != "ttf_genetic_phylogatr_gate_d_fragility_execution_plan_v0.1":
        raise RuntimeError("fragility execution plan schema drift")
    if rule.get("schema") != "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1":
        raise RuntimeError("formal Phase-3 rule schema drift")
    if sha256_path(args.phase3_rule) != plan["phase3_rule_sha256"]:
        raise RuntimeError("formal Phase-3 rule provenance drift")
    formal_anchor = plan.get("formal_anchor")
    if not isinstance(formal_anchor, dict):
        raise RuntimeError("fragility execution plan lacks formal anchor")
    if float(formal_anchor.get("retention_fraction", -1.0)) != 1.0:
        raise RuntimeError("fragility formal anchor retention drift")
    if formal_anchor.get("source") != "formal_phase3_qualification_receipt":
        raise RuntimeError("fragility full-geometry anchor is not the formal Phase-3 receipt")

    curve: list[dict] = [dict(formal_anchor)]
    by_retention: dict[float, dict] = {1.0: curve[0]}
    identical_statuses = {
        "IDENTICAL_TO_FULL_GEOMETRY_NO_RERUN",
        "IDENTICAL_TO_HIGHER_RETENTION_GEOMETRY_NO_RERUN",
    }
    for level in plan["levels"]:
        retention = float(level["retention_fraction"])
        status = str(level["status"])
        if status in identical_statuses:
            inherited = float(level["inherits_metrics_from_retention_fraction"])
            source_row = by_retention.get(inherited)
            if source_row is None:
                raise RuntimeError(
                    f"identical-geometry plateau source retention={inherited} is unavailable"
                )
            row = _plateau_row(level, source_row)
            curve.append(row)
            by_retention[retention] = row
            continue
        if status == "STRUCTURAL_SUPPORT_BELOW_FORMAL_MINIMUM":
            row = {
                "retention_fraction": retention,
                "source": "diagnostic_thinned_geometry_support_only",
                "status": "STRUCTURAL_SUPPORT_BELOW_FORMAL_MINIMUM",
                "geometry_fingerprint_sha256": level[
                    "geometry_fingerprint_sha256"
                ],
                "minimum_endpoint_disjoint_ibd_training_edges": int(
                    level["metrics"]["minimum_endpoint_disjoint_ibd_training_edges"]
                ),
                "A3_private_null_rejection_rate": None,
                "A3_private_null_Wilson95_upper": None,
                "max_private_null_Wilson95_upper": None,
                "shared_A2_rejection_rate": None,
                "shared_A2_Wilson95_lower": None,
                "distance_of_max_private_null_upper_to_formal_0.10_ceiling": None,
                "distance_of_shared_A2_lower_to_formal_0.80_floor": None,
                "formal_gate_d_decision_authority": False,
                "phase4_identity_opening_authority": False,
            }
            curve.append(row)
            by_retention[retention] = row
            continue
        if status != "SYNTHETIC_DIAGNOSTIC_RUNNABLE":
            raise RuntimeError(f"unexpected fragility execution level status: {status}")
        row = _collect_level(
            args.input_dir,
            retention=retention,
            level=level,
            rule=rule,
        )
        curve.append(row)
        by_retention[retention] = row

    expected_retentions = [1.0] + [float(level["retention_fraction"]) for level in plan["levels"]]
    observed_retentions = [float(row["retention_fraction"]) for row in curve]
    if observed_retentions != expected_retentions:
        raise RuntimeError("fragility curve retention order/completeness drift")

    out = {
        "schema": "ttf_genetic_phylogatr_gate_d_fragility_curve_v0.1",
        "status": "DIAGNOSTIC_FRAGILITY_CURVE_COMPLETE",
        "execution_plan_sha256": sha256_path(args.execution_plan),
        "phase3_rule_sha256": plan["phase3_rule_sha256"],
        "formal_qualification_sha256": plan["formal_qualification_sha256"],
        "execution_rule_sha256": plan["execution_rule_sha256"],
        "dataset_digest_sha256": plan["dataset_digest_sha256"],
        "full_geometry_fingerprint_sha256": plan["full_geometry_fingerprint_sha256"],
        "expected_retention_fractions": expected_retentions,
        "formal_full_geometry_status": plan["formal_anchor"]["formal_status"],
        "curve": curve,
        "threshold_reference_lines": {
            "private_null_Wilson95_upper_ceiling": float(
                rule["qualification"]["type1_wilson95_upper_ceiling"]
            ),
            "shared_A2_Wilson95_lower_floor": float(
                rule["qualification"]["shared_A2_wilson95_lower_floor"]
            ),
            "read_only": True,
        },
        "authority_firewall": plan["authority_firewall"],
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "formal_gate_d_decision_made_by_this_curve": False,
        "phase4_identity_opening_authorized_by_this_curve": False,
        "diagnostic_completion_is_procedural_prerequisite_only": True,
        "claim_boundary": (
            "Only retention=1.0 imports the formal Phase-3 Gate-D decision. Thinned levels "
            "are descriptive synthetic sensitivity results and cannot rescue, alter, or reinterpret it. "
            "Identical geometry fingerprints reuse the prior summary exactly rather than adding Monte Carlo noise. "
            "Curve completion is required procedurally before Phase-4 opening but its diagnostic metrics "
            "do not enter the Phase-4 scientific PASS/FAIL decision."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": out["status"],
                "formal_full_geometry_status": out["formal_full_geometry_status"],
                "retention_levels": observed_retentions,
                "formal_gate_d_decision_made_by_this_curve": False,
                "phase4_identity_opening_authorized_by_this_curve": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
