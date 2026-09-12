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
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--source-ledger", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    source = json.loads(args.source_ledger.read_text())
    if rule.get("status") != "frozen_before_v09_reserve_b_geometry_is_materialized":
        raise RuntimeError("v0.9 rule drift")
    if source.get("schema") != "ttf_v09_rgfca_reserve_b_max_density_geometry_v0.1":
        raise RuntimeError("wrong v0.9 source")
    if source.get("synthetic_worlds_run_at_freeze") != 0 or source.get("empirical_trait_or_colour_fields_read") is not False:
        raise RuntimeError("v0.9 geometry was not outcome-unopened")
    if source.get("predecessor", {}).get("species_overlap") != [] or source.get("predecessor", {}).get("photo_overlap") != []:
        raise RuntimeError("v0.9 successor overlaps failed Gate-D")

    refs: dict[str, dict] = {}
    obs: dict[tuple[float, float], dict] = {}
    for path in sorted(args.input_dir.glob("*.json")):
        d = json.loads(path.read_text())
        if d.get("schema") != "ttf_v09_gate_d_statistics_v0.1":
            continue
        if d.get("empirical_colour_outcome_opened") is not False or d.get("failed_v08_gate_d_panel_reused") is not False:
            raise RuntimeError("v0.9 evidence firewall failed")
        cfg = d["config"]
        if d["mode"] == "reference":
            label = cfg["configuration_label"]
            if label in refs or len(d["statistics"]) != 1999:
                raise RuntimeError("v0.9 reference set drift")
            refs[label] = d
        else:
            key = (float(cfg["shared_fraction"]), float(cfg["amplitude"]))
            if key in obs or len(d["statistics"]) != 500:
                raise RuntimeError("v0.9 observed set drift")
            obs[key] = d

    if set(refs) != set(LABELS):
        raise RuntimeError(("reference set", sorted(refs)))
    if set(obs) != set(MANDATORY):
        raise RuntimeError(("observed set", sorted(obs)))

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
        p_values: list[float] = []
        pairs: list[str] = []
        for ti, ui in zip(t, u):
            pv, selected, _, _ = profiled_private_pvalue(
                float(ti), float(ui), reference,
                profile_draws=999, selected_configs=2, scale_floor=1e-6,
            )
            p_values.append(pv)
            pairs.append("+".join(selected))
        p = np.asarray(p_values, dtype=float)
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
        pair_counts[key] = {k: int(v) for k, v in sorted(Counter(pairs).items())}
        recovery[key] = None if shared > 0 else float(np.mean([MATCH[amp] in x.split("+") for x in pairs]))

    precision = qualify_calibration_precision(
        cells,
        moderate_amplitude=2.0,
        type1_upper_ceiling=0.10,
        power_lower_floor=0.80,
    ).to_dict()
    out = {
        "schema": "ttf_v09_gate_d_rgfca_reserve_b_v0.1",
        "status": "synthetic_deployment_adequacy_on_species_disjoint_max_density_reserve_b",
        "rule": str(args.rule),
        "source_ledger": str(args.source_ledger),
        "cells": [c.to_dict() for c in cells],
        "precision_qualification": precision,
        "gate_d_pass": bool(precision["passed"]),
        "matching_private_label_recovery": recovery,
        "selected_pair_counts": pair_counts,
        "empirical_colour_outcome_opened": False,
        "failed_v08_gate_d_panel_reused": False,
        "claim_ready": False,
        "interpretation_boundary": (
            "v0.9 Gate-D is a pre-outcome adequacy test of a new species- and photo-disjoint maximal-support design. "
            "A pass only authorizes separately frozen empirical measurement/analysis; a failure freezes this design without tuning it."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"gate_d_pass": out["gate_d_pass"], "precision": precision}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
