#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    rows={}
    for path in sorted(args.input_dir.glob("*.json")):
        payload=json.loads(path.read_text())
        if payload.get("schema") in {
            "ttf_q_known_truth_detectability_scenario_v0.1",
            "ttf_q_known_truth_detectability_scenario_v0.2",
        }:
            rows[payload["scenario"]]=payload
    required={"reference","endpoint_loss","control_redundancy","source_concentration"}
    if set(rows)!=required:
        raise RuntimeError(f"detectability scenario set drift: {sorted(rows)}")

    scenario_schemas={p["schema"] for p in rows.values()}
    if len(scenario_schemas)!=1:
        raise RuntimeError(f"mixed detectability scenario schema versions: {sorted(scenario_schemas)}")
    scenario_schema=next(iter(scenario_schemas))
    aggregate_schema=(
        "ttf_q_known_truth_detectability_benchmark_result_v0.2"
        if scenario_schema.endswith("_v0.2")
        else "ttf_q_known_truth_detectability_benchmark_result_v0.1"
    )

    pairs=[]
    for a,b in itertools.combinations(sorted(rows),2):
        pa,pb=rows[a],rows[b]
        same_binary=(
            pa["single_cell_reference"]["binary_qualified"]
            == pb["single_cell_reference"]["binary_qualified"]
        )
        calibrated_a=pa.get("calibrated_evaluable_envelope")
        calibrated_b=pb.get("calibrated_evaluable_envelope")
        if calibrated_a is not None and calibrated_b is not None:
            envelope_a={
                amplitude: cell.get("evaluable_grid_mde")
                for amplitude,cell in calibrated_a.items()
            }
            envelope_b={
                amplitude: cell.get("evaluable_grid_mde")
                for amplitude,cell in calibrated_b.items()
            }
        else:
            envelope_a=pa["minimum_detectable_effect_by_private_amplitude"]
            envelope_b=pb["minimum_detectable_effect_by_private_amplitude"]
        different_envelope=(envelope_a!=envelope_b)
        pairs.append({
            "a":a,
            "b":b,
            "same_binary_classification":same_binary,
            "different_evaluable_envelope":different_envelope,
            "demonstrates_binary_information_loss":bool(same_binary and different_envelope),
        })

    demonstrations=[row for row in pairs if row["demonstrates_binary_information_loss"]]
    payload={
        "schema":aggregate_schema,
        "scenario_schema":scenario_schema,
        "status":(
            "COMPLETE_CALIBRATED_KNOWN_TRUTH_DETECTABILITY_BENCHMARK"
            if aggregate_schema.endswith("_v0.2")
            else "COMPLETE_KNOWN_TRUTH_DETECTABILITY_BENCHMARK"
        ),
        "scenarios":{
            name:{
                "truth":p["truth"],
                "observed_information":p["observed_information"],
                "minimum_detectable_effect_by_private_amplitude":p["minimum_detectable_effect_by_private_amplitude"],
                "calibrated_evaluable_envelope":p.get("calibrated_evaluable_envelope"),
                "single_cell_reference":p["single_cell_reference"],
                "null_rejection_rate_by_private_amplitude":{
                    f"{float(row['private_amplitude']):g}":float(row["rejection_rate"])
                    for row in p["null_cells"]
                },
            }
            for name,p in rows.items()
        },
        "pairwise_comparisons":pairs,
        "binary_information_loss_demonstrations":len(demonstrations),
        "required_comparison_pass":bool(demonstrations),
        "boundary":next(iter(rows.values()))["boundary"],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":payload["status"],
        "binary_information_loss_demonstrations":len(demonstrations),
        "required_comparison_pass":payload["required_comparison_pass"],
        "scenario_summary":{
            name:{
                "binary":x["single_cell_reference"]["binary_qualified"],
                "raw_mde":x["minimum_detectable_effect_by_private_amplitude"],
                "calibrated_envelope":x.get("calibrated_evaluable_envelope"),
            }
            for name,x in rows.items()
        },
    },sort_keys=True))
    if not demonstrations:
        raise RuntimeError("pre-frozen scenarios did not demonstrate binary information loss")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
