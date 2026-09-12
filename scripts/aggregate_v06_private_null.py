#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell
from ttf.precision import qualify_calibration_precision
from ttf.private_null_inference import envelope_upper_pvalues, upper_monte_carlo_pvalue

PANELS = ("v03", "v04")
MANDATORY = ((0.0,0.5),(0.0,1.0),(0.0,2.0),(0.0,3.0),(1.0,2.0))
ORACLE_LABEL = {0.5:"A0p5", 1.0:"A1", 2.0:"A2", 3.0:"A3"}


def _cell_from_p(
    *,
    shared: float,
    amplitude: float,
    observed: np.ndarray,
    pvalues: np.ndarray,
    null_mean: float,
) -> CalibrationCell:
    return CalibrationCell(
        shared_fraction=float(shared),
        amplitude=float(amplitude),
        n_replicates=int(len(observed)),
        alpha=0.05,
        rejection_rate=float(np.mean(pvalues <= 0.05)),
        mean_statistic=float(np.mean(observed)),
        mean_null_statistic=float(null_mean),
        median_p_value=float(np.median(pvalues)),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate v0.6 geometry-conditioned private-null development")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("status") != "frozen_before_v06_private_null_diagnostic_outcomes":
        raise RuntimeError("v0.6 private-null rule was not prospectively frozen")
    configurations = rule["private_null_family_v0.1"]["configurations"]
    expected_labels = {item["label"] for item in configurations}
    if len(expected_labels) != 8:
        raise RuntimeError("v0.6 private-null configuration grid drifted")
    if rule["private_null_family_v0.1"]["reference_worlds_per_configuration"] != 1999:
        raise RuntimeError("private reference world count drifted")
    if rule["observed_worlds_per_cell"] != 500:
        raise RuntimeError("observed world count drifted")

    paths = sorted(args.input_dir.glob("*.json"))
    reference: dict[str, dict[str, dict]] = {panel: {} for panel in PANELS}
    observed: dict[str, dict[tuple[float,float], dict]] = {panel: {} for panel in PANELS}
    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "ttf_v06_private_null_statistics_v0.1":
            continue
        if payload.get("status") != "development_only":
            raise RuntimeError(f"unexpected status in {path}")
        if payload.get("confirmatory_performance_opened") is not False or payload.get("rgfca_reserve_opened") is not False:
            raise RuntimeError(f"firewall failure in {path}")
        panel = payload["panel"]
        if panel not in PANELS:
            raise RuntimeError(f"unexpected panel {panel}")
        cfg = payload["config"]
        if payload["mode"] == "reference":
            label = cfg["configuration_label"]
            if label in reference[panel]:
                raise RuntimeError(f"duplicate reference {panel}/{label}")
            if len(payload["statistics"]) != 1999:
                raise RuntimeError(f"reference size drift for {panel}/{label}")
            reference[panel][label] = payload
        elif payload["mode"] == "observed":
            if payload.get("baseline_v05_reproduced") is not True:
                raise RuntimeError(f"observed raw baseline not reproduced in {path}")
            key = (float(cfg["shared_fraction"]), float(cfg["amplitude"]))
            if key in observed[panel]:
                raise RuntimeError(f"duplicate observed cell {panel}/{key}")
            if len(payload["statistics"]) != 500:
                raise RuntimeError(f"observed size drift for {panel}/{key}")
            observed[panel][key] = payload
        else:
            raise RuntimeError(f"unexpected mode in {path}")

    for panel in PANELS:
        if set(reference[panel]) != expected_labels:
            raise RuntimeError(f"{panel} private reference set drifted: {sorted(reference[panel])}")
        if set(observed[panel]) != set(MANDATORY):
            raise RuntimeError(f"{panel} observed cell set drifted: {sorted(observed[panel])}")

    panel_reports = {}
    envelope_all_pass = True
    for panel in PANELS:
        refs = {
            label: np.asarray(payload["statistics"], dtype=float)
            for label, payload in reference[panel].items()
        }
        ref_means = {label: float(values.mean()) for label, values in refs.items()}
        envelope_cells = []
        oracle_cells = []
        least_counts = {}
        for shared, amplitude in MANDATORY:
            obs = np.asarray(observed[panel][(shared, amplitude)]["statistics"], dtype=float)
            p_env, least = envelope_upper_pvalues(obs, refs)
            counts = Counter(least)
            least_counts[f"s{shared}-a{amplitude}"] = {
                label: int(counts.get(label, 0)) for label in sorted(refs)
            }
            # For reporting only, the envelope has no single null mean. Use the
            # maximum component mean as a conservative descriptive reference.
            envelope_cells.append(
                _cell_from_p(
                    shared=shared,
                    amplitude=amplitude,
                    observed=obs,
                    pvalues=p_env,
                    null_mean=max(ref_means.values()),
                )
            )
            oracle_label = ORACLE_LABEL[float(amplitude)]
            oracle_ref = refs[oracle_label]
            p_oracle = np.asarray(
                [upper_monte_carlo_pvalue(float(value), oracle_ref) for value in obs],
                dtype=float,
            )
            oracle_cells.append(
                _cell_from_p(
                    shared=shared,
                    amplitude=amplitude,
                    observed=obs,
                    pvalues=p_oracle,
                    null_mean=float(oracle_ref.mean()),
                )
            )

        envelope_precision = qualify_calibration_precision(
            envelope_cells,
            moderate_amplitude=2.0,
            type1_upper_ceiling=0.10,
            power_lower_floor=0.80,
        ).to_dict()
        oracle_precision = qualify_calibration_precision(
            oracle_cells,
            moderate_amplitude=2.0,
            type1_upper_ceiling=0.10,
            power_lower_floor=0.80,
        ).to_dict()
        envelope_pass = bool(envelope_precision["passed"])
        envelope_all_pass = envelope_all_pass and envelope_pass
        panel_reports[panel] = {
            "reference_summaries": {
                label: payload["summary"] for label, payload in sorted(reference[panel].items())
            },
            "envelope": {
                "cells": [cell.to_dict() for cell in envelope_cells],
                "precision": envelope_precision,
                "passed": envelope_pass,
                "least_favourable_configuration_counts": least_counts,
            },
            "oracle_same_amplitude": {
                "selectable": False,
                "cells": [cell.to_dict() for cell in oracle_cells],
                "precision": oracle_precision,
                "passed_descriptively": bool(oracle_precision["passed"]),
            },
        }

    selected = "private_envelope_raw_transfer" if envelope_all_pass else None
    out = {
        "schema": "ttf_v06_private_null_diagnostic_v0.1",
        "status": "development_only_after_v05_no_candidate",
        "rule": str(args.rule),
        "base_statistic": rule["base_statistic"],
        "panels": panel_reports,
        "selected_development_candidate": selected,
        "selection_rule_applied_mechanically": True,
        "oracle_same_amplitude_is_nonselectable": True,
        "claim_ready": False,
        "confirmatory_performance_opened": False,
        "rgfca_reserve_opened": False,
        "next_gate": (
            "If selected_development_candidate is non-null, prospectively broaden the "
            "private-null nuisance family across transition width and private-boundary "
            "offset before locking any birds/butterflies confirmation. Otherwise preserve "
            "the failure and select no v0.6 candidate from this family."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "selected_development_candidate": selected,
        "panels": {
            panel: {
                "envelope": panel_reports[panel]["envelope"]["precision"],
                "oracle": panel_reports[panel]["oracle_same_amplitude"]["precision"],
            }
            for panel in PANELS
        },
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
