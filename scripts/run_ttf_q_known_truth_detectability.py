#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.relational_benchmark import known_truth_design
from ttf.relational_qualification import (
    detectability_surface,
    dyadic_signal_support,
    minimum_detectable_effects,
    residualize_primary_relation,
)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--scenario",required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    contract=json.loads(args.contract.read_text())
    if contract.get("schema")!="ttf_q_known_truth_detectability_benchmark_v0.1":
        raise RuntimeError("unexpected detectability benchmark contract")
    if contract.get("status")!="FROZEN_BEFORE_DETECTABILITY_RESULTS":
        raise RuntimeError("detectability benchmark was not frozen before results")
    if args.scenario not in contract["scenarios"]:
        raise RuntimeError("unknown frozen scenario")

    geometry=contract["geometry"]
    spec=contract["scenarios"][args.scenario]
    synth=contract["synthetic_response"]
    design=known_truth_design(
        n_source=int(geometry["n_source"]),
        n_target=int(geometry["n_target"]),
        endpoint_retained_fraction=float(spec["endpoint_retained_fraction"]),
        baseline_survival_after_fe=float(spec["baseline_survival_after_fe"]),
        geography_survival_after_baseline=float(spec["geography_survival_after_baseline"]),
        concentration=str(spec["concentration"]),
        seed=int(spec["seed"]),
    )

    predictors=np.column_stack((
        design.primary,
        design.baseline_control,
        design.geographic_control,
    ))
    residual, info=residualize_primary_relation(
        design.source,
        design.target,
        design.primary,
        predictors[:,1:],
    )
    support=dyadic_signal_support(
        design.source,
        design.target,
        residual,
    )

    surface=detectability_surface(
        design.source,
        design.target,
        predictors,
        effects=synth["effects"],
        private_amplitudes=synth["private_amplitudes"],
        primary_index=0,
        alpha=float(synth["alpha_one_sided"]),
        worlds_per_cell=int(synth["worlds_per_cell"]),
        source_intercept_sd=float(synth["source_intercept_sd"]),
        target_intercept_sd=float(synth["target_intercept_sd"]),
        dyad_noise_sd=float(synth["dyad_noise_sd"]),
        fixed_coefficients=[
            0.0,
            float(synth["fixed_baseline_control_effect"]),
            float(synth["fixed_geographic_control_effect"]),
        ],
        source_random_slope_index=1,
        target_random_slope_index=2,
        random_slope_sd_per_amplitude=float(synth["random_slope_sd_per_amplitude"]),
        response_transform=str(synth["response_transform"]),
        seed_namespace=f"{synth['seed_namespace']}|{args.scenario}",
        block_size=int(synth["block_size"]),
    )

    mde=minimum_detectable_effects(
        surface,
        target_power=float(synth["target_power"]),
    )
    reference=contract["single_cell_reference"]
    null_rows=[row for row in surface if float(row["effect"])==0.0]
    type1_pass=all(float(row["wilson95_upper"])<=0.05 for row in null_rows)
    target=[
        row for row in surface
        if float(row["effect"])==float(reference["effect"])
        and float(row["private_amplitude"])==float(reference["private_amplitude"])
    ]
    if len(target)!=1:
        raise RuntimeError("single-cell reference not found exactly once")
    power_pass=float(target[0]["wilson95_lower"])>=float(synth["target_power"])

    payload={
        "schema":"ttf_q_known_truth_detectability_scenario_v0.1",
        "status":"COMPLETE_KNOWN_TRUTH_DETECTABILITY_SCENARIO",
        "scenario":args.scenario,
        "scenario_spec":spec,
        "truth":design.truth,
        "observed_information":{
            "total_unique_variance_fraction":info.total_unique_variance_fraction,
            "control_unique_variance_fraction_after_fe":info.control_unique_variance_fraction_after_fe,
            "signal_effective_sources":support.signal_effective_sources,
            "signal_effective_targets":support.signal_effective_targets,
            "max_source_signal_share":support.max_source_signal_share,
            "max_target_signal_share":support.max_target_signal_share,
        },
        "minimum_detectable_effect_by_private_amplitude":mde,
        "single_cell_reference":{
            "effect":float(reference["effect"]),
            "private_amplitude":float(reference["private_amplitude"]),
            "type1_pass":bool(type1_pass),
            "power_pass":bool(power_pass),
            "binary_qualified":bool(type1_pass and power_pass),
            "reference_cell":target[0],
        },
        "null_cells":null_rows,
        "surface":surface,
        "boundary":contract["boundary"],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "scenario":args.scenario,
        "total_unique_variance_fraction":info.total_unique_variance_fraction,
        "signal_effective_sources":support.signal_effective_sources,
        "mde":mde,
        "binary_qualified":payload["single_cell_reference"]["binary_qualified"],
        "reference_power_lower":target[0]["wilson95_lower"],
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
