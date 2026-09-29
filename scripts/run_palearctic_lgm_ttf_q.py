#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.relational_qualification import (
    calibrated_detectability_envelope,
    detectability_surface,
    dyadic_signal_support,
    residualize_primary_relation,
)


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20),b""):
            h.update(chunk)
    return h.hexdigest()


def zscore(x):
    x=np.asarray(x,dtype=float)
    sd=float(x.std())
    if not np.isfinite(sd) or sd<=np.finfo(float).eps:
        raise ValueError("constant/non-finite continuous predictor")
    return (x-float(x.mean()))/sd


def remap(values):
    x=np.asarray(values)
    levels=sorted(np.unique(x).tolist())
    lookup={v:i for i,v in enumerate(levels)}
    return np.asarray([lookup[v] for v in x],dtype=np.int64)


def evaluate_opening(info,support,envelope,rule):
    gates={}
    gates["unique_information"]=(
        info.total_unique_variance_fraction
        >= float(rule["minimum_total_unique_variance_fraction"])
    )
    gates["source_signal_breadth"]=(
        support.signal_effective_sources
        >= float(rule["minimum_signal_effective_sources"])
        and support.max_source_signal_share
        <= float(rule["maximum_single_source_signal_share"])
    )
    gates["target_signal_breadth"]=(
        support.signal_effective_targets
        >= float(rule["minimum_signal_effective_targets"])
        and support.max_target_signal_share
        <= float(rule["maximum_single_target_signal_share"])
    )

    required=[float(x) for x in rule["required_null_qualified_amplitudes"]]
    gates["null_qualification"]=all(
        bool(envelope[f"{a:g}"]["calibration_pass"])
        for a in required
    )
    a2=envelope["2"]["evaluable_grid_mde"]
    gates["A2_detectability"]=(
        a2 is not None
        and float(a2) <= float(rule["required_A2_evaluable_mde_max"])
    )
    gates["overall_pass"]=all(gates.values())
    return gates


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--design-npz",type=Path,required=True)
    ap.add_argument("--design-summary",type=Path,required=True)
    ap.add_argument("--execution",type=Path,required=True)
    ap.add_argument("--parent-rule",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--block-size",type=int,default=100)
    args=ap.parse_args()

    execution=json.loads(args.execution.read_text())
    if execution.get("schema")!="ttf_genetic_palearctic_lgm_ttf_q_execution_v0.1":
        raise RuntimeError("unexpected execution rule")
    parent=json.loads(args.parent_rule.read_text())
    if parent.get("schema")!="ttf_genetic_palearctic_lgm_subpanel_rule_v0.1":
        raise RuntimeError("unexpected parent rule")
    summary=json.loads(args.design_summary.read_text())
    if summary.get("status")!="PASS_TO_TTF_Q_RESPONSE_BLIND_CHARACTERIZATION":
        raise RuntimeError("LGM relation design is not authorized for TTF-Q characterization")
    if any(bool(v) for v in summary["response_firewall"].values()):
        raise RuntimeError("subpanel response firewall is open")

    data=np.load(args.design_npz,allow_pickle=False)
    source=remap(np.asarray(data["source_index"],dtype=np.int64))
    target=remap(np.asarray(data["target_index"],dtype=np.int64))
    r_lgm=zscore(data["R_LGM"])
    r_present=zscore(data["R_present"])
    geo=zscore(data["geographic_coverage"])
    same_order=np.asarray(data["same_order"],float)
    same_family=np.asarray(data["same_family"],float)
    ratio=zscore(np.abs(np.log(np.asarray(data["retained_occurrence_count_ratio"],float))))

    same_order=same_order-float(same_order.mean())
    same_family=same_family-float(same_family.mean())
    x=np.column_stack((r_lgm,r_present,geo,same_order,same_family,ratio))

    residual,info=residualize_primary_relation(
        source,target,r_lgm,x[:,1:]
    )
    support=dyadic_signal_support(source,target,residual)

    synth=execution["synthetic_world"]
    common=dict(
        source=source,
        target=target,
        predictors=x,
        private_amplitudes=synth["private_amplitudes"],
        primary_index=int(execution["primary_index"]),
        alpha=float(synth["alpha_one_sided"]),
        source_intercept_sd=float(synth["source_intercept_sd"]),
        target_intercept_sd=float(synth["target_intercept_sd"]),
        dyad_noise_sd=float(synth["dyad_noise_sd"]),
        fixed_coefficients=synth["fixed_coefficients"],
        source_random_slope_index=2,
        target_random_slope_index=3,
        random_slope_sd_per_amplitude=float(synth["random_slope_sd_per_amplitude"]),
        response_transform=str(synth["response_transform"]),
        seed_namespace=str(synth["seed_namespace"]),
        block_size=int(args.block_size),
    )

    null_surface=detectability_surface(
        effects=[0.0],
        worlds_per_cell=int(synth["null_worlds_per_amplitude"]),
        **common,
    )
    positive_surface=detectability_surface(
        effects=synth["positive_effects"],
        worlds_per_cell=int(synth["positive_worlds_per_cell"]),
        **common,
    )
    surface=list(null_surface)+list(positive_surface)
    envelope=calibrated_detectability_envelope(
        surface,
        target_power=float(synth["target_power"]),
        type1_wilson_upper_max=float(synth["type1_wilson95_upper_max"]),
    )

    gates=evaluate_opening(info,support,envelope,execution["opening_rule"])
    status=(
        execution["opening_rule"]["if_passed"]
        if gates["overall_pass"]
        else execution["opening_rule"]["if_failed"]
    )

    payload={
        "schema":"ttf_genetic_palearctic_lgm_ttf_q_result_v0.1",
        "status":status,
        "design_npz_sha256":sha256_path(args.design_npz),
        "design_summary_sha256":sha256_path(args.design_summary),
        "execution_sha256":sha256_path(args.execution),
        "parent_rule_sha256":sha256_path(args.parent_rule),
        "information":{
            "dyads":info.dyads,
            "source_clusters":info.source_clusters,
            "target_clusters":info.target_clusters,
            "source_target_fe_retained_variance_fraction":info.source_target_fe_retained_variance_fraction,
            "control_unique_variance_fraction_after_fe":info.control_unique_variance_fraction_after_fe,
            "total_unique_variance_fraction":info.total_unique_variance_fraction,
            "partial_sd_in_primary_sd_units":info.partial_sd_in_primary_sd_units,
            "vif_like_after_fe":info.vif_like_after_fe,
        },
        "signal_concentration":{
            "effective_dyads":support.signal_effective_dyads,
            "effective_sources":support.signal_effective_sources,
            "effective_targets":support.signal_effective_targets,
            "max_source_share":support.max_source_signal_share,
            "max_target_share":support.max_target_signal_share,
        },
        "calibrated_evaluable_envelope":envelope,
        "opening_gates":gates,
        "A3_is_diagnostic_only":True,
        "response_firewall":{
            "species_level_genetic_scores_used":False,
            "pairwise_subpanel_T_st_computed":False,
            "beta_LGM_computed":False,
        },
        "next_step":(
            "Freeze one-shot subgroup empirical opening authorization, then compute pair-specific T_st and beta_LGM exactly once."
            if gates["overall_pass"]
            else "STOP. Do not compute pair-specific subgroup genetic response."
        ),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "information":payload["information"],
        "signal_concentration":payload["signal_concentration"],
        "opening_gates":gates,
        "envelope":envelope,
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
