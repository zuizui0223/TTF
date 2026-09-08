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

PANELS = ("v03", "v04")
MANDATORY = ((0.0,0.5),(0.0,1.0),(0.0,2.0),(0.0,3.0),(1.0,2.0))
MATCH_LABEL = {0.5:"A0p5", 1.0:"A1", 2.0:"A2", 3.0:"A3"}


def _cell(shared: float, amplitude: float, statistic: np.ndarray, p: np.ndarray) -> CalibrationCell:
    return CalibrationCell(
        shared_fraction=float(shared),
        amplitude=float(amplitude),
        n_replicates=int(len(statistic)),
        alpha=0.05,
        rejection_rate=float(np.mean(p <= 0.05)),
        mean_statistic=float(np.mean(statistic)),
        mean_null_statistic=float("nan"),
        median_p_value=float(np.median(p)),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate frozen v0.7 profiled-private development")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--v06", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("status") != "frozen_before_v07_profiled_private_null_outcomes":
        raise RuntimeError("v0.7 rule was not prospectively frozen")
    labels = set(rule["private_reference_family"]["labels"])
    if labels != {"A0","A0p5","A1","A2","A3","A5","A10","infinite_snr"}:
        raise RuntimeError("v0.7 private configuration family drifted")
    if rule["reference_partition"]["profile_strength_draws"] != 999:
        raise RuntimeError("v0.7 profile split drifted")
    if rule["reference_partition"]["calibration_statistic_draws"] != 1000:
        raise RuntimeError("v0.7 calibration split drifted")
    if rule["profiling"]["selected_configurations"] != 2:
        raise RuntimeError("v0.7 profiled config count drifted")

    reference = {panel: {} for panel in PANELS}
    observed = {panel: {} for panel in PANELS}
    for path in sorted(args.input_dir.glob("*.json")):
        d = json.loads(path.read_text())
        if d.get("schema") != "ttf_v07_profiled_statistics_v0.1":
            continue
        if d.get("status") != "development_only":
            raise RuntimeError(f"unexpected status in {path}")
        if d.get("heldout_trait_used_for_strength") is not False:
            raise RuntimeError(f"held-out trait leaked into nuisance strength in {path}")
        if d.get("confirmatory_performance_opened") is not False or d.get("rgfca_reserve_opened") is not False:
            raise RuntimeError(f"future-layer firewall failure in {path}")
        panel = d["panel"]
        cfg = d["config"]
        if d["mode"] == "reference":
            label = cfg["configuration_label"]
            if label in reference[panel]:
                raise RuntimeError(f"duplicate reference {panel}/{label}")
            if len(d["statistics"]) != 1999 or len(d["training_strength"]) != 1999:
                raise RuntimeError(f"reference size drift {panel}/{label}")
            reference[panel][label] = d
        else:
            if d.get("baseline_v05_reproduced") is not True:
                raise RuntimeError(f"raw baseline reproduction failed in {path}")
            key=(float(cfg["shared_fraction"]),float(cfg["amplitude"]))
            if key in observed[panel]:
                raise RuntimeError(f"duplicate observed cell {panel}/{key}")
            if len(d["statistics"]) != 500 or len(d["training_strength"]) != 500:
                raise RuntimeError(f"observed size drift {panel}/{key}")
            observed[panel][key]=d

    for panel in PANELS:
        if set(reference[panel]) != labels:
            raise RuntimeError(f"{panel} reference labels drifted")
        if set(observed[panel]) != set(MANDATORY):
            raise RuntimeError(f"{panel} observed cell set drifted")

    v06 = json.loads(args.v06.read_text())
    reports = {}
    all_pass = True
    for panel in PANELS:
        refs = {
            label: (
                np.asarray(d["training_strength"], dtype=float),
                np.asarray(d["statistics"], dtype=float),
            )
            for label,d in reference[panel].items()
        }
        cells=[]
        selected_counts={}
        matching_recovery={}
        distance_summary={}
        for shared,amp in MANDATORY:
            d=observed[panel][(shared,amp)]
            t=np.asarray(d["statistics"],dtype=float)
            u=np.asarray(d["training_strength"],dtype=float)
            p=np.empty(len(t),dtype=float)
            pair_counter=Counter()
            distance_acc={label:[] for label in sorted(refs)}
            match=0
            match_label=MATCH_LABEL.get(float(amp)) if shared==0.0 else None
            for i,(ti,ui) in enumerate(zip(t,u)):
                pi,selected,_,distances=profiled_private_pvalue(
                    float(ti),float(ui),refs,
                    profile_draws=999,
                    selected_configs=2,
                    scale_floor=1e-6,
                )
                p[i]=pi
                pair_counter["+".join(selected)] += 1
                if match_label is not None and match_label in selected:
                    match += 1
                for label,value in distances.items():
                    distance_acc[label].append(float(value))
            key=f"s{shared}-a{amp}"
            selected_counts[key]=dict(sorted(pair_counter.items()))
            matching_recovery[key]=(None if match_label is None else float(match/len(t)))
            distance_summary[key]={
                label:{
                    "median":float(np.median(values)),
                    "q05":float(np.quantile(values,0.05)),
                    "q95":float(np.quantile(values,0.95)),
                }
                for label,values in sorted(distance_acc.items())
            }
            cells.append(_cell(shared,amp,t,p))

        precision=qualify_calibration_precision(
            cells,
            moderate_amplitude=2.0,
            type1_upper_ceiling=0.10,
            power_lower_floor=0.80,
        ).to_dict()
        passed=bool(precision["passed"])
        all_pass=all_pass and passed

        all_u=np.concatenate([np.asarray(d["training_strength"],dtype=float) for d in reference[panel].values()])
        all_t=np.concatenate([np.asarray(d["statistics"],dtype=float) for d in reference[panel].values()])
        strength_summaries={
            label:{
                "median":float(np.median(np.asarray(d["training_strength"],dtype=float))),
                "mad":float(np.median(np.abs(np.asarray(d["training_strength"],dtype=float)-np.median(np.asarray(d["training_strength"],dtype=float))))),
                "q05":float(np.quantile(np.asarray(d["training_strength"],dtype=float),0.05)),
                "q95":float(np.quantile(np.asarray(d["training_strength"],dtype=float),0.95)),
                "mean_T":float(np.mean(np.asarray(d["statistics"],dtype=float))),
            }
            for label,d in sorted(reference[panel].items())
        }
        reports[panel]={
            "cells":[cell.to_dict() for cell in cells],
            "precision":precision,
            "passed":passed,
            "selected_pair_counts":selected_counts,
            "matching_private_label_recovery":matching_recovery,
            "profile_distance_summary":distance_summary,
            "private_reference_strength":strength_summaries,
            "pooled_private_strength_T_correlation":float(np.corrcoef(all_u,all_t)[0,1]),
            "v06_oracle_nonselectable":v06["panels"][panel]["oracle_same_amplitude"]["precision"],
        }

    selected="profiled_private_raw_transfer" if all_pass else None
    out={
        "schema":"ttf_v07_profiled_private_null_diagnostic_v0.1",
        "status":"development_only_after_v06_unrestricted_envelope_failure",
        "rule":str(args.rule),
        "panels":reports,
        "selected_development_candidate":selected,
        "selection_rule_applied_mechanically":True,
        "heldout_transfer_used_for_nuisance_selection":False,
        "claim_ready":False,
        "confirmatory_performance_opened":False,
        "rgfca_reserve_opened":False,
        "next_gate":(
            "If selected, broaden private nuisance prospectively across transition width "
            "and boundary offset before any birds/butterflies confirmation. If null, "
            "preserve failure and do not open confirmation or reserve."
        ),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"selected":selected,"panels":{p:reports[p]["precision"] for p in PANELS}},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
