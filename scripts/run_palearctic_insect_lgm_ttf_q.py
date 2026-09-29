#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.relational_qualification import (
    calibrated_detectability_envelope,
    detectability_surface,
    dyadic_signal_support,
    relation_repeatability,
    residualize_primary_relation,
)


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    a=np.asarray(x,dtype=float)
    sd=float(a.std())
    if not np.isfinite(a).all() or sd<=np.finfo(float).eps:
        raise ValueError("non-finite or constant predictor")
    return (a-float(a.mean()))/sd


def remap(x: np.ndarray) -> np.ndarray:
    vals=sorted(map(int,np.unique(x)))
    lookup={value:i for i,value in enumerate(vals)}
    return np.asarray([lookup[int(v)] for v in x],dtype=np.int64)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",type=Path,required=True)
    ap.add_argument("--relation",type=Path,required=True)
    ap.add_argument("--replicates",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    contract=json.loads(args.contract.read_text())
    if contract.get("schema")!="ttf_palearctic_insect_lgm_refugia_v0.4":
        raise RuntimeError("unexpected contract")

    data=np.load(args.relation,allow_pickle=False)
    source=remap(np.asarray(data["source_index"],dtype=np.int64))
    target=remap(np.asarray(data["target_index"],dtype=np.int64))
    primary=zscore(np.asarray(data["R_LGM"],dtype=float))
    current=zscore(np.asarray(data["R_current"],dtype=float))
    geography=zscore(np.asarray(data["geographic_coop"],dtype=float))
    same_order=np.asarray(data["same_order"],dtype=float)
    same_order=same_order-float(same_order.mean())
    same_family=np.asarray(data["same_family"],dtype=float)
    same_family=same_family-float(same_family.mean())
    count_ratio=zscore(np.abs(np.log(np.asarray(data["presence_count_ratio"],dtype=float))))

    names=[
        "z_R_LGM",
        "z_R_current",
        "z_geographic_coop",
        "centered_same_order",
        "centered_same_family",
        "z_abs_log_GBIF_count_ratio",
    ]
    predictors=np.column_stack((
        primary,current,geography,same_order,same_family,count_ratio
    ))
    residual,nonredundancy=residualize_primary_relation(
        source,target,primary,predictors[:,1:]
    )
    support=dyadic_signal_support(source,target,residual)

    matrix=np.load(args.replicates,allow_pickle=False)
    if matrix.shape!=(100,len(primary)):
        raise RuntimeError(
            f"repeatability matrix shape drift: {matrix.shape}, expected {(100,len(primary))}"
        )
    repeatability=relation_repeatability(matrix)

    q=contract["ttf_q_qualification"]
    fixed_map=q["fixed_control_effects"]
    fixed=np.zeros(len(names),dtype=float)
    lookup={name:i for i,name in enumerate(names)}
    for name,value in fixed_map.items():
        fixed[lookup[name]]=float(value)

    common=dict(
        source=source,
        target=target,
        predictors=predictors,
        private_amplitudes=q["private_amplitudes"],
        primary_index=0,
        alpha=float(q["alpha_one_sided"]),
        source_intercept_sd=float(q["source_intercept_sd"]),
        target_intercept_sd=float(q["target_intercept_sd"]),
        dyad_noise_sd=float(q["dyad_noise_sd"]),
        fixed_coefficients=fixed,
        source_random_slope_index=lookup[q["source_random_slope_predictor"]],
        target_random_slope_index=lookup[q["target_random_slope_predictor"]],
        random_slope_sd_per_amplitude=float(q["random_slope_sd_per_amplitude"]),
        response_transform=str(q["response_transform"]),
        block_size=100,
    )
    relation_sha=sha256_path(args.relation)
    null_surface=detectability_surface(
        **common,
        effects=[0.0],
        worlds_per_cell=int(q["worlds_per_null_cell"]),
        seed_namespace=(
            "palearctic-insect-lgm-qualification-null-v0.1|" + relation_sha
        ),
    )
    positive_effects=[float(x) for x in q["effects"] if float(x)>0]
    positive_surface=detectability_surface(
        **common,
        effects=positive_effects,
        worlds_per_cell=int(q["worlds_per_positive_cell"]),
        seed_namespace=(
            "palearctic-insect-lgm-qualification-positive-v0.1|" + relation_sha
        ),
    )
    surface=sorted(
        [*null_surface,*positive_surface],
        key=lambda row:(float(row["private_amplitude"]),float(row["effect"])),
    )
    envelope=calibrated_detectability_envelope(
        surface,
        target_power=float(q["target_power"]),
        type1_wilson_upper_max=float(q["type1_wilson95_upper_ceiling"]),
    )

    domain=q["response_opening_domain"]
    calibration_all=all(
        bool(envelope[f"A{int(float(a))}"]["calibration_qualified"])
        for a in q["private_amplitudes"]
    )
    def mde(A: str):
        return envelope[A]["evaluable_grid_mde"]
    a2=mde("A2")
    a3=mde("A3")
    opening=(
        float(repeatability.icc)>=float(domain["measurement_repeatability_icc_min"])
        and calibration_all
        and a2 is not None and float(a2)<=float(domain["evaluable_grid_mde_A2_max"])
        and a3 is not None and float(a3)<=float(domain["evaluable_grid_mde_A3_max"])
        and all(np.isfinite(float(row["rejection_rate"])) for row in surface)
    )
    status=(
        "PASS_TTF_Q_AUTHORIZE_SINGLE_SUBPANEL_RESPONSE_OPENING"
        if opening
        else "NOT_EVALUABLE_PALEARCTIC_INSECT_LGM"
    )

    payload={
        "schema":"ttf_palearctic_insect_lgm_ttf_q_result_v0.1",
        "status":status,
        "inputs":{
            "contract_sha256":sha256_path(args.contract),
            "relation_sha256":sha256_path(args.relation),
            "replicates_sha256":sha256_path(args.replicates),
        },
        "predictors_in_order":names,
        "relation_repeatability":asdict(repeatability),
        "nonredundancy":asdict(nonredundancy),
        "dyadic_signal_support":asdict(support),
        "calibrated_evaluable_envelope":envelope,
        "detectability_surface":surface,
        "simulation_precision":{"null_worlds_per_cell":int(q["worlds_per_null_cell"]),"positive_worlds_per_cell":int(q["worlds_per_positive_cell"])},
        "response_opening_checks":{
            "repeatability_pass":float(repeatability.icc)>=float(domain["measurement_repeatability_icc_min"]),
            "all_A0_to_A3_calibration_qualified":calibration_all,
            "A2_mde":a2,
            "A2_mde_pass":a2 is not None and float(a2)<=float(domain["evaluable_grid_mde_A2_max"]),
            "A3_mde":a3,
            "A3_mde_pass":a3 is not None and float(a3)<=float(domain["evaluable_grid_mde_A3_max"]),
            "all_worlds_finite":True,
            "all_required":bool(opening),
        },
        "genetic_response_used":False,
        "subpanel_response_authorized":bool(opening),
        "next_step":(
            "freeze authorization receipt, then compute exactly one subpanel T_st/beta_LGM"
            if opening
            else "STOP_WITHOUT_SUBPANEL_GENETIC_RESPONSE"
        ),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status,
        "repeatability_icc":repeatability.icc,
        "A2_mde":a2,
        "A3_mde":a3,
        "subpanel_response_authorized":bool(opening),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
