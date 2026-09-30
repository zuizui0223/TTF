#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


INPUT = Path("manuscript/generated/ttf_q_scalar_results_v0.1.json")
C2_RECEIPT = Path("benchmarks/frozen/ttf_q_c2_result_receipt_v0.1.json")
STOP_CASE = Path("benchmarks/frozen/ttf_q_palearctic_prospective_stop_case_v0.1.json")
OUTPUT = Path("manuscript/figures/ttf_q_v0.1")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    data=json.loads(INPUT.read_text())
    c2=json.loads(C2_RECEIPT.read_text())
    stop=json.loads(STOP_CASE.read_text())
    OUTPUT.mkdir(parents=True,exist_ok=True)
    mpl.rcParams["svg.fonttype"]="none"
    mpl.rcParams["svg.hashsalt"]="ttf-q-v0.1"

    # Figure 1: response-blind information ladder on the shared development geometry.
    stages=["Endpoint\nidentity","Baseline\ncontrols","Geography"]
    current=[
        data["current_climate"]["development_endpoint_survival"],
        data["current_climate"]["development_unique_before_geography"],
        data["current_climate"]["development_unique_after_geography"],
    ]
    historical=[
        data["historical_climate"]["development_endpoint_survival"],
        data["historical_climate"]["development_unique_before_geography"],
        data["historical_climate"]["development_unique_after_geography"],
    ]
    fig,ax=plt.subplots(figsize=(7.4,4.7))
    ax.plot(stages,current,marker="o",label="Current climate")
    ax.plot(stages,historical,marker="o",label="Historical climate")
    ax.set_ylim(0,0.75)
    ax.set_ylabel("Fraction of total focal-relation variance retained")
    ax.set_title("Response-blind relational information survival")
    ax.legend()
    fig.tight_layout()
    p1=OUTPUT/"figure1_information_ladder.svg"
    fig.savefig(p1,format="svg",metadata={"Date":None})
    plt.close(fig)

    # Figure 2: known-truth calibrated detectability.
    A=[0,1,2,3]
    fig,ax=plt.subplots(figsize=(7.4,4.7))
    reference=data["known_truth_detectability"]["reference"]["evaluable_mde_A0_A1_A2_A3"]
    endpoint=data["known_truth_detectability"]["endpoint_loss"]["evaluable_mde_A0_A1_A2_A3"]
    redundant=data["known_truth_detectability"]["control_redundancy"]["evaluable_mde_A0_A1_A2_A3"]
    ax.plot(A,reference,marker="o",label="Reference")
    ax.plot([x-0.035 for x in A],endpoint,marker="s",label="Endpoint loss")
    ax.plot([x+0.035 for x in A],redundant,marker="^",linestyle="--",label="Control redundancy")
    ax.text(
        0.98,0.05,
        "Source-concentrated design: null qualification FAIL at all A",
        transform=ax.transAxes,ha="right",va="bottom",fontsize=9,
    )
    ax.set_xticks(A)
    ax.set_xlabel("Private-heterogeneity amplitude A")
    ax.set_ylabel("Calibration-qualified grid MDE")
    ax.set_title("Known-truth designs share a binary label but not an envelope")
    ax.legend()
    fig.tight_layout()
    p2=OUTPUT/"figure2_known_truth_mde.svg"
    fig.savefig(p2,format="svg",metadata={"Date":None})
    plt.close(fig)

    # Figure 3: empirical calibrated detectability.
    fig,ax=plt.subplots(figsize=(7.4,4.7))
    ax.plot(
        A,
        data["current_climate"]["development_evaluable_mde_A0_A1_A2_A3"],
        marker="o",
        label="Current climate",
    )
    ax.plot(
        A,
        data["historical_climate"]["development_evaluable_mde_A0_A1_A2_A3"],
        marker="o",
        label="Historical climate | current",
    )
    ax.set_xticks(A)
    ax.set_xlabel("Private-heterogeneity amplitude A")
    ax.set_ylabel("Calibration-qualified grid MDE")
    ax.set_title("Empirical response-blind climate relation envelopes")
    ax.legend()
    fig.tight_layout()
    p3=OUTPUT/"figure3_empirical_mde.svg"
    fig.savefig(p3,format="svg",metadata={"Date":None})
    plt.close(fig)

    # Figure 4: concentration is distinct from information amount.
    ref=data["known_truth_detectability"]["reference"]
    con=data["known_truth_detectability"]["source_concentration"]
    labels=["Reference","Source-concentrated"]
    effective_sources=[ref["signal_effective_sources"],con["signal_effective_sources"]]
    unique=[ref["unique_fraction"],con["unique_fraction"]]

    fig,ax=plt.subplots(figsize=(7.4,4.7))
    x=[0,1]
    ax.bar(x,effective_sources)
    ax.set_xticks(x,labels)
    ax.set_ylabel("Inverse-Herfindahl effective source count")
    ax.set_title("Equal unique information can have very different endpoint concentration")
    qualification=[
        ref["calibration_pass_all"],
        con["calibration_pass_all"],
    ]
    for i,(n,u,passed) in enumerate(zip(effective_sources,unique,qualification)):
        state="PASS" if passed else "FAIL"
        ax.text(
            i,n+0.5,
            f"unique fraction = {u:.2f}\nnull qualification: {state}",
            ha="center",va="bottom",
        )
    ax.set_ylim(0,max(effective_sources)*1.22)
    fig.tight_layout()
    p4=OUTPUT/"figure4_source_concentration.svg"
    fig.savefig(p4,format="svg",metadata={"Date":None})
    plt.close(fig)

    # Figure 5: a prospectively frozen application can stop before outcome opening.
    broad=[
        c2["development"]["C2_total_unique_variance_fraction"],
        c2["development"]["max_source_signal_share_after_geography"],
        c2["development"]["max_target_signal_share_after_geography"],
    ]
    pal=[
        stop["authoritative_deterministic_ttf_q"]["unique_information"]["observed_total_unique_variance_fraction"],
        stop["authoritative_deterministic_ttf_q"]["source_signal_breadth"]["maximum_single_source_share"],
        stop["authoritative_deterministic_ttf_q"]["target_signal_breadth"]["maximum_single_target_share"],
    ]
    labels=["Unique relation\nfraction","Max source\nsignal share","Max target\nsignal share"]
    x=[0,1,2]
    width=0.34
    fig,ax=plt.subplots(figsize=(7.4,4.9))
    ax.bar([v-width/2 for v in x],broad,width=width,label="Broad climate geometry")
    ax.bar([v+width/2 for v in x],pal,width=width,label="Palearctic prospective case")
    ax.axhline(0.15,linestyle="--",linewidth=1.2)
    ax.text(
        2.48,0.154,
        "Palearctic predeclared 0.15 boundary",
        ha="right",va="bottom",fontsize=8.5,
    )
    ax.set_xticks(x,labels)
    ax.set_ylim(0,0.34)
    ax.set_ylabel("Response-blind relation diagnostic")
    ax.set_title("Prospective Palearctic stop relative to broad climate geometry")
    ax.legend()
    ax.text(
        0.02,0.98,
        "Palearctic rule: unique fraction ≥ 0.15; endpoint shares ≤ 0.15",
        transform=ax.transAxes,ha="left",va="top",fontsize=8.5,
    )
    fig.tight_layout()
    p5=OUTPUT/"figure5_prospective_stop.svg"
    fig.savefig(p5,format="svg",metadata={"Date":None})
    plt.close(fig)

    manifest={
        "schema":"ttf_q_figure_manifest_v0.1",
        "status":"DESCRIPTIVE_RENDERING_FROM_FROZEN_SCALAR_HANDOFF",
        "input":str(INPUT),
        "input_sha256":sha256(INPUT),
        "additional_frozen_inputs":{
            str(C2_RECEIPT):sha256(C2_RECEIPT),
            str(STOP_CASE):sha256(STOP_CASE),
        },
        "new_inference_performed":False,
        "figures":{
            p.name:sha256(p)
            for p in (p1,p2,p3,p4,p5)
        },
    }
    manifest_path=OUTPUT/"figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":manifest["status"],
        "figures":len(manifest["figures"]),
        "output":str(OUTPUT),
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
