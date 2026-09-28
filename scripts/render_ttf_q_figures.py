#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


INPUT = Path("manuscript/generated/ttf_q_scalar_results_v0.1.json")
OUTPUT = Path("manuscript/figures/ttf_q_v0.1")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    data=json.loads(INPUT.read_text())
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
    ax.plot(stages,historical,marker="o",label="Historical climate | current")
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
    for key,label in [
        ("reference","Reference"),
        ("endpoint_loss","Endpoint loss"),
        ("control_redundancy","Control redundancy"),
    ]:
        y=data["known_truth_detectability"][key]["evaluable_mde_A0_A1_A2_A3"]
        ax.plot(A,y,marker="o",label=label)
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
    for i,(n,u) in enumerate(zip(effective_sources,unique)):
        ax.text(i,n+0.5,f"unique fraction = {u:.2f}",ha="center",va="bottom")
    ax.set_ylim(0,max(effective_sources)*1.22)
    fig.tight_layout()
    p4=OUTPUT/"figure4_source_concentration.svg"
    fig.savefig(p4,format="svg",metadata={"Date":None})
    plt.close(fig)

    manifest={
        "schema":"ttf_q_figure_manifest_v0.1",
        "status":"DESCRIPTIVE_RENDERING_FROM_FROZEN_SCALAR_HANDOFF",
        "input":str(INPUT),
        "input_sha256":sha256(INPUT),
        "new_inference_performed":False,
        "figures":{
            p.name:sha256(p)
            for p in (p1,p2,p3,p4)
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
