#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell
from ttf.precision import qualify_calibration_precision
from ttf.profiled_private_null import profiled_private_pvalue

LABELS = ("A0", "A0p5", "A1", "A2", "A3", "A5", "A10", "infinite_snr")
MANDATORY = ((0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0))
MATCH = {0.5: "A0p5", 1.0: "A1", 2.0: "A2", 3.0: "A3"}


def main() -> int:
    ap = argparse.ArgumentParser(description="Aggregate frozen v0.11 fresh density-scaled qualification.")
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--source-ledger", type=Path, required=True)
    ap.add_argument("--execution-audit", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    source = json.loads(args.source_ledger.read_text())
    audit = json.loads(args.execution_audit.read_text())
    if rule.get("schema") != "ttf_v11_fresh_density_scaled_graph_rule_v0.1" or rule.get("status") != "frozen_before_v11_fresh_plant_geometry_acquisition_and_before_any_v11_synthetic_worlds":
        raise RuntimeError("v0.11 rule drift")
    if source.get("schema") != "ttf_v11_fresh_plant_geometry_source_v0.1" or source.get("fresh_geometry_freeze_pass") is not True:
        raise RuntimeError("v0.11 source is not a frozen fresh geometry pass")
    if source.get("synthetic_worlds_run_at_freeze") != 0 or source.get("empirical_trait_or_colour_fields_read") is not False:
        raise RuntimeError("v0.11 geometry was not outcome-unopened at freeze")
    if audit.get("schema") != "ttf_v11_execution_geometry_audit_v0.1" or audit.get("status") != "geometry_only_before_v11_synthetic_qualification":
        raise RuntimeError("v0.11 execution audit drift")
    if audit.get("synthetic_worlds_run") != 0 or audit.get("empirical_colour_or_pixels_opened") is not False:
        raise RuntimeError("v0.11 execution audit opened future outcomes")

    arch = rule["locked_graph_architecture"]
    exe = rule["prospective_execution"]
    if int(arch["k"]) != 15 or not np.isclose(float(arch["bandwidth_km"]), 500.0):
        raise RuntimeError("v0.11 architecture drift")
    if int(exe["seed"]) != 20260923 or int(exe["profile_strength_draws"]) != 999 or int(exe["calibration_statistic_draws"]) != 1000:
        raise RuntimeError("v0.11 execution/reference partition drift")
    if int(exe["selected_private_configurations"]) != 2:
        raise RuntimeError("v0.11 profiling count drift")
    if int(exe["reference_worlds_per_configuration"]) != 1999 or int(exe["observed_worlds_per_cell"]) != 500:
        raise RuntimeError("v0.11 world-count drift")
    if not np.isclose(float(exe["type1_wilson_upper_ceiling"]), 0.10) or not np.isclose(float(exe["power_wilson_lower_floor"]), 0.80):
        raise RuntimeError("v0.11 precision threshold drift")
    if [list(x) for x in MANDATORY] != exe["mandatory_cells"]:
        raise RuntimeError("v0.11 mandatory cell drift")

    refs: dict[str, dict] = {}
    obs: dict[tuple[float, float], dict] = {}
    fingerprint = audit["geometry"]["fingerprint_sha256"]
    expected_train = audit["split"]["train_species"]
    expected_eval = audit["split"]["eval_species"]
    for path in sorted(args.input_dir.glob("*.json")):
        d = json.loads(path.read_text())
        if d.get("schema") != "ttf_v11_fresh_density_scaled_statistics_v0.1":
            continue
        if d.get("status") != "fresh_density_scaled_synthetic_qualification_only":
            raise RuntimeError(f"unexpected v0.11 statistics status in {path}")
        if d.get("empirical_colour_outcome_opened") is not False or d.get("candidate_image_pixels_opened") is not False:
            raise RuntimeError(f"v0.11 empirical firewall failed in {path}")
        if d.get("failed_v08_or_v09_panel_reused") is not False:
            raise RuntimeError(f"failed predecessor panel was reused in {path}")
        if d.get("heldout_transfer_used_for_nuisance_selection") is not False:
            raise RuntimeError(f"heldout T leaked into nuisance profiling in {path}")
        if d["geometry"]["fingerprint_sha256"] != fingerprint:
            raise RuntimeError(f"geometry fingerprint drift in {path}")
        if d["split"]["train_species"] != expected_train or d["split"]["eval_species"] != expected_eval:
            raise RuntimeError(f"split drift in {path}")
        cfg = d["config"]
        if int(cfg["k"]) != 15 or not np.isclose(float(cfg["bandwidth"]), 500.0):
            raise RuntimeError(f"graph/bandwidth drift in {path}")
        execution = d["execution"]
        if execution.get("operator") != "exact_dense_gaussian_opportunity_corrected_segment_integrated_field":
            raise RuntimeError(f"operator drift in {path}")
        if execution.get("distance_cutoff") is not None or execution.get("kernel_approximation") is not False or execution.get("dtype_reduction") is not False or execution.get("response_dependent_pruning") is not False:
            raise RuntimeError(f"chunked exactness firewall failed in {path}")

        statistic = d.get("statistics", [])
        strength = d.get("training_strength", [])
        if d["mode"] == "reference":
            label = cfg["configuration_label"]
            if label not in LABELS or label in refs:
                raise RuntimeError(f"bad/duplicate v0.11 reference {label}")
            if len(statistic) != 1999 or len(strength) != 1999:
                raise RuntimeError(f"reference size drift for {label}")
            if not np.isclose(float(cfg["shared_fraction"]), 0.0):
                raise RuntimeError(f"reference shared fraction drift for {label}")
            refs[label] = d
        elif d["mode"] == "observed":
            key = (float(cfg["shared_fraction"]), float(cfg["amplitude"]))
            if key not in MANDATORY or key in obs:
                raise RuntimeError(f"bad/duplicate v0.11 observed cell {key}")
            if len(statistic) != 500 or len(strength) != 500:
                raise RuntimeError(f"observed size drift for {key}")
            if not np.isclose(float(cfg["noise_sd"]), 0.8):
                raise RuntimeError(f"observed noise drift for {key}")
            obs[key] = d
        else:
            raise RuntimeError(f"unknown statistics mode in {path}")

    if set(refs) != set(LABELS):
        raise RuntimeError(("reference set drift", sorted(refs)))
    if set(obs) != set(MANDATORY):
        raise RuntimeError(("observed set drift", sorted(obs)))

    reference = {
        label: (
            np.asarray(d["training_strength"], dtype=float),
            np.asarray(d["statistics"], dtype=float),
        )
        for label, d in refs.items()
    }
    cells: list[CalibrationCell] = []
    recovery: dict[str, float | None] = {}
    pair_counts: dict[str, dict[str, int]] = {}
    distance_summary: dict[str, dict[str, dict[str, float]]] = {}

    for shared, amp in MANDATORY:
        d = obs[(shared, amp)]
        t = np.asarray(d["statistics"], dtype=float)
        u = np.asarray(d["training_strength"], dtype=float)
        p = np.empty(500, dtype=float)
        pairs: list[str] = []
        distances_by_label = {label: [] for label in LABELS}
        match_label = MATCH.get(amp) if shared == 0.0 else None
        matched = 0
        for i, (ti, ui) in enumerate(zip(t, u)):
            pv, selected, _, distances = profiled_private_pvalue(
                float(ti),
                float(ui),
                reference,
                profile_draws=999,
                selected_configs=2,
                scale_floor=1e-6,
            )
            p[i] = pv
            pairs.append("+".join(selected))
            if match_label is not None and match_label in selected:
                matched += 1
            for label, value in distances.items():
                distances_by_label[label].append(float(value))

        cells.append(
            CalibrationCell(
                shared_fraction=shared,
                amplitude=amp,
                n_replicates=500,
                alpha=0.05,
                rejection_rate=float(np.mean(p <= 0.05)),
                mean_statistic=float(t.mean()),
                mean_null_statistic=float("nan"),
                median_p_value=float(np.median(p)),
            )
        )
        key = f"s{shared}-a{amp}"
        pair_counts[key] = {name: int(count) for name, count in sorted(Counter(pairs).items())}
        recovery[key] = None if match_label is None else float(matched / len(t))
        distance_summary[key] = {
            label: {
                "median": float(np.median(values)),
                "q05": float(np.quantile(values, 0.05)),
                "q95": float(np.quantile(values, 0.95)),
            }
            for label, values in distances_by_label.items()
        }

    precision = qualify_calibration_precision(
        cells,
        moderate_amplitude=2.0,
        type1_upper_ceiling=0.10,
        power_lower_floor=0.80,
    ).to_dict()
    passed = bool(precision["passed"])

    private_reference_summary = {
        label: {
            "strength_median": float(np.median(np.asarray(d["training_strength"], dtype=float))),
            "strength_mad": float(np.median(np.abs(
                np.asarray(d["training_strength"], dtype=float)
                - np.median(np.asarray(d["training_strength"], dtype=float))
            ))),
            "strength_q05": float(np.quantile(np.asarray(d["training_strength"], dtype=float), 0.05)),
            "strength_q95": float(np.quantile(np.asarray(d["training_strength"], dtype=float), 0.95)),
            "mean_T": float(np.mean(np.asarray(d["statistics"], dtype=float))),
        }
        for label, d in sorted(refs.items())
    }

    out = {
        "schema": "ttf_v11_fresh_density_scaled_qualification_v0.1",
        "status": "synthetic_qualification_on_fresh_species_disjoint_250x100_density_scaled_graph",
        "rule": str(args.rule),
        "source_ledger": str(args.source_ledger),
        "execution_audit": str(args.execution_audit),
        "geometry_fingerprint_sha256": fingerprint,
        "locked_architecture": {
            "species": 250,
            "records_per_species": 100,
            "train_species": 125,
            "eval_species": 125,
            "neighbor_fraction": 0.15,
            "k": 15,
            "bandwidth_km": 500.0,
            "master_seed": 20260923,
        },
        "cells": [cell.to_dict() for cell in cells],
        "precision_qualification": precision,
        "v11_pass": passed,
        "matching_private_label_recovery": recovery,
        "selected_pair_counts": pair_counts,
        "profile_distance_summary": distance_summary,
        "private_reference_summary": private_reference_summary,
        "empirical_colour_outcome_opened": False,
        "candidate_image_pixels_opened": False,
        "failed_v08_or_v09_panel_reused": False,
        "claim_ready": False,
        "interpretation_boundary": (
            "This is a synthetic adequacy test on one prospectively frozen fresh plant geometry with the single predeclared density-scaled k=15 architecture. "
            "A pass qualifies only that synthetic fresh-geometry architecture and does not establish an empirical flower-colour boundary. Before empirical inference, the separately frozen classifiability-only actual-geometry gate and requalification must still be applied. "
            "A failure freezes the k=15 design on this panel; it may not be tuned on these outcomes."
        ),
        "next_gate": (
            "If v11_pass is true, keep empirical colour values unopened and execute the already-frozen classifiability-only actual-geometry gate before any empirical TTF conclusion. "
            "If false, stop this design without k/bandwidth/species/record/profiling/threshold tuning."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"v11_pass": passed, "precision": precision}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
