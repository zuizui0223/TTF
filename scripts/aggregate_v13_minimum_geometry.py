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
    ap = argparse.ArgumentParser(description="Aggregate frozen v0.13 minimum-geometry synthetic requalification.")
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--source-ledger", type=Path, required=True)
    ap.add_argument("--v11-execution-audit", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    source = json.loads(args.source_ledger.read_text())
    audit = json.loads(args.v11_execution_audit.read_text())
    if rule.get("schema") != "ttf_v13_minimum_classifiability_geometry_rule_v0.1":
        raise RuntimeError("v0.13 rule drift")
    if source.get("schema") != "ttf_v13_minimum_classifiability_geometry_result_v0.1" or source.get("classifiability_gate_pass") is not True:
        raise RuntimeError("v0.13 minimum geometry is not frozen-pass")
    if source.get("species_specific_classifiability_counts_persisted") is not False:
        raise RuntimeError("v0.13 species-specific classifiability surface leaked")
    if source.get("colour_vector_values_read_by_ttf") is not False or source.get("palette_fractions_read_by_ttf") is not False or source.get("pairwise_colour_distances_computed") is not False or source.get("ttf_empirical_statistic_computed") is not False:
        raise RuntimeError("v0.13 outcome firewall failed")
    if audit.get("schema") != "ttf_v11_execution_geometry_audit_v0.1":
        raise RuntimeError("v0.11 split authority missing")

    arch = rule["actual_geometry_architecture"]
    exe = rule["synthetic_requalification"]
    if int(arch["k"]) != 2 or not np.isclose(float(arch["bandwidth_km"]), 500.0) or int(arch["records_per_species"]) != 15:
        raise RuntimeError("v0.13 architecture drift")
    if int(exe["master_seed"]) != 20260927 or int(exe["reference_worlds_per_configuration"]) != 1999 or int(exe["observed_worlds_per_cell"]) != 500:
        raise RuntimeError("v0.13 world contract drift")
    if int(exe["profile_strength_draws"]) != 999 or int(exe["calibration_statistic_draws"]) != 1000 or int(exe["selected_private_configurations"]) != 2:
        raise RuntimeError("v0.13 profiling drift")
    if not np.isclose(float(exe["type1_wilson_upper_ceiling"]), 0.10) or not np.isclose(float(exe["power_wilson_lower_floor"]), 0.80):
        raise RuntimeError("v0.13 thresholds drift")
    if [list(x) for x in MANDATORY] != exe["mandatory_cells"]:
        raise RuntimeError("v0.13 mandatory cells drift")

    expected_train = list(audit["split"]["train_species"])
    expected_eval = list(audit["split"]["eval_species"])
    refs: dict[str, dict] = {}
    obs: dict[tuple[float, float], dict] = {}
    fingerprint = None

    for path in sorted(args.input_dir.glob("*.json")):
        d = json.loads(path.read_text())
        if d.get("schema") != "ttf_v13_minimum_geometry_statistics_v0.1":
            continue
        if d.get("status") != "synthetic_requalification_on_minimum_classifiability_geometry_only":
            raise RuntimeError(f"unexpected v0.13 statistics status: {path}")
        if d.get("empirical_colour_values_read") is not False or d.get("palette_fractions_read") is not False or d.get("pairwise_colour_distances_computed") is not False or d.get("empirical_ttf_statistic_computed") is not False:
            raise RuntimeError(f"empirical firewall failed: {path}")
        if d.get("heldout_transfer_used_for_nuisance_selection") is not False:
            raise RuntimeError(f"heldout leakage: {path}")
        if d["split"]["train_species"] != expected_train or d["split"]["eval_species"] != expected_eval:
            raise RuntimeError(f"split drift: {path}")
        cfg = d["config"]
        if int(cfg["k"]) != 2 or not np.isclose(float(cfg["bandwidth"]), 500.0) or int(cfg["master_seed"]) != 20260927:
            raise RuntimeError(f"architecture drift: {path}")
        execution = d["execution"]
        if execution.get("operator") != "exact_dense_gaussian_opportunity_corrected_segment_integrated_field" or execution.get("kernel_approximation") is not False or execution.get("distance_cutoff") is not None or execution.get("dtype_reduction") is not False or execution.get("response_dependent_pruning") is not False:
            raise RuntimeError(f"exactness drift: {path}")
        geom = d["geometry"]
        if int(geom["species_count"]) != 250 or int(geom["record_count"]) != 3750 or int(geom["records_per_species"]) != 15:
            raise RuntimeError(f"geometry size drift: {path}")
        fp = geom["fingerprint_sha256"]
        if fingerprint is None:
            fingerprint = fp
        elif fingerprint != fp:
            raise RuntimeError("v0.13 geometry fingerprint drift across cells")
        statistic = d["statistics"]
        strength = d["training_strength"]
        if d["mode"] == "reference":
            label = cfg["configuration_label"]
            if label not in LABELS or label in refs or len(statistic) != 1999 or len(strength) != 1999:
                raise RuntimeError(f"bad/duplicate reference: {label}")
            refs[label] = d
        elif d["mode"] == "observed":
            key = (float(cfg["shared_fraction"]), float(cfg["amplitude"]))
            if key not in MANDATORY or key in obs or len(statistic) != 500 or len(strength) != 500:
                raise RuntimeError(f"bad/duplicate observed: {key}")
            obs[key] = d
        else:
            raise RuntimeError("unknown statistics mode")

    if set(refs) != set(LABELS) or set(obs) != set(MANDATORY):
        raise RuntimeError(("incomplete v0.13 input set", sorted(refs), sorted(obs)))

    reference = {
        label: (np.asarray(d["training_strength"], dtype=float), np.asarray(d["statistics"], dtype=float))
        for label, d in refs.items()
    }
    cells: list[CalibrationCell] = []
    recovery: dict[str, float | None] = {}
    pair_counts: dict[str, dict[str, int]] = {}

    for shared, amp in MANDATORY:
        d = obs[(shared, amp)]
        t = np.asarray(d["statistics"], dtype=float)
        u = np.asarray(d["training_strength"], dtype=float)
        p = np.empty(500, dtype=float)
        pairs: list[str] = []
        match_label = MATCH.get(amp) if shared == 0.0 else None
        matched = 0
        for i, (ti, ui) in enumerate(zip(t, u)):
            pv, selected, _, _ = profiled_private_pvalue(
                float(ti), float(ui), reference,
                profile_draws=999, selected_configs=2, scale_floor=1e-6,
            )
            p[i] = pv
            pairs.append("+".join(selected))
            if match_label is not None and match_label in selected:
                matched += 1
        cells.append(CalibrationCell(
            shared_fraction=shared,
            amplitude=amp,
            n_replicates=500,
            alpha=0.05,
            rejection_rate=float(np.mean(p <= 0.05)),
            mean_statistic=float(t.mean()),
            mean_null_statistic=float("nan"),
            median_p_value=float(np.median(p)),
        ))
        key = f"s{shared}-a{amp}"
        pair_counts[key] = {name: int(count) for name, count in sorted(Counter(pairs).items())}
        recovery[key] = None if match_label is None else float(matched / len(t))

    precision = qualify_calibration_precision(
        cells,
        moderate_amplitude=2.0,
        type1_upper_ceiling=0.10,
        power_lower_floor=0.80,
    ).to_dict()
    passed = bool(precision["passed"])
    out = {
        "schema": "ttf_v13_minimum_geometry_requalification_v0.1",
        "status": "synthetic_requalification_on_minimum_classifiability_geometry",
        "classifiability_gate_pass": True,
        "geometry_fingerprint_sha256": fingerprint,
        "locked_architecture": {
            "species": 250,
            "records_per_species": 15,
            "records": 3750,
            "train_species": 125,
            "eval_species": 125,
            "neighbor_fraction": 0.15,
            "k": 2,
            "bandwidth_km": 500.0,
            "master_seed": 20260927
        },
        "cells": [cell.to_dict() for cell in cells],
        "precision_qualification": precision,
        "v13_pass": passed,
        "matching_private_label_recovery": recovery,
        "selected_pair_counts": pair_counts,
        "empirical_colour_values_read": False,
        "palette_fractions_read": False,
        "pairwise_colour_distances_computed": False,
        "empirical_ttf_statistic_computed": False,
        "claim_ready": False,
        "interpretation_boundary": "A pass establishes only synthetic adequacy on the exact response-blind 250x15 minimum-classifiability actual geometry. It is not an empirical flower-colour result.",
        "next_gate": "If v13_pass is true, only a separately blob-pinned empirical authorization may read the frozen colour representation for these exact retained records. If false, stop empirical TTF without tuning or replacement."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"v13_pass": passed, "precision_qualification": precision}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
