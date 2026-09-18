#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import seed_for
from ttf.phylogatr_compact_execution import (
    prepare_phylogatr_compact_cached_transfer,
    prepare_phylogatr_compact_ttf_design,
    score_phylogatr_compact_world_batch,
)
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_fragility_execution import load_fragility_execution_context


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--execution-plan", type=Path, required=True)
    ap.add_argument("--geometry-dir", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--retention", type=float, required=True)
    ap.add_argument("--config-label", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    context = load_fragility_execution_context(
        args.execution_plan,
        args.geometry_dir,
        args.phase2_manifest,
        args.phase3_rule,
        args.retention,
    )
    rule = context.phase3_rule
    configs = rule["synthetic_world"]["private_reference_configurations"]
    label = str(args.config_label)
    if label not in configs:
        raise RuntimeError(f"unknown diagnostic private reference configuration: {label}")
    start, stop = int(args.start), int(args.stop)
    requested = {"configuration": label, "start": start, "stop": stop}
    if requested not in context.level["reference_jobs"]:
        raise RuntimeError("requested fragility reference shard is not in frozen execution plan")

    core = rule["core_method"]
    design = prepare_phylogatr_compact_ttf_design(
        context.geometries,
        train_species=context.train_species,
        eval_species=context.eval_species,
        bandwidth=float(core["bandwidth_km"]),
        prior_strength=float(core["prior_strength"]),
        segment_points=int(core["segment_points"]),
        min_training_edges=int(
            rule["geometry_contract"]["minimum_endpoint_disjoint_ibd_training_edges"]
        ),
        strength_neighbours=int(core["strength_neighbours"]),
    )
    execution = rule["execution_defaults"]
    cached = prepare_phylogatr_compact_cached_transfer(
        design,
        edge_chunk_size=int(execution["edge_chunk_size"]),
        train_chunk_size=int(execution["train_chunk_size"]),
    )
    amplitude, noise_sd = map(float, configs[label])
    world_rule = rule["synthetic_world"]
    batch_width = int(execution["batch_width"])
    rows: list[dict[str, float | int]] = []
    for batch_start in range(start, stop, batch_width):
        batch_stop = min(batch_start + batch_width, stop)
        worlds = [
            simulate_genetic_distance_world(
                context.geometries,
                shared_fraction=0.0,
                residual_amplitude=amplitude,
                ibd_strength=float(world_rule["ibd_strength_primary"]),
                noise_sd=noise_sd,
                transition_width=float(world_rule["transition_width"]),
                noise_dimensions=int(world_rule["noise_dimensions"]),
                seed=seed_for(
                    context.master_seed,
                    "genetic_fragility_private_reference",
                    context.retention_label,
                    label,
                    replicate,
                ),
            )
            for replicate in range(batch_start, batch_stop)
        ]
        scored = score_phylogatr_compact_world_batch(design, worlds, cached)
        for column, replicate in enumerate(range(batch_start, batch_stop)):
            if int(scored.n_eval_species[column]) != len(design.eval_species):
                raise RuntimeError("non-finite held-out species count in fragility reference shard")
            rows.append(
                {
                    "replicate": int(replicate),
                    "training_strength": float(scored.training_strengths[column]),
                    "statistic": float(scored.statistics[column]),
                }
            )

    out = {
        "schema": "ttf_genetic_phylogatr_gate_d_fragility_reference_shard_v0.1",
        "status": "DIAGNOSTIC_SYNTHETIC_PRIVATE_REFERENCE_SHARD",
        "retention_fraction": context.retention_fraction,
        "retention_label": context.retention_label,
        "geometry_fingerprint_sha256": context.level["geometry_fingerprint_sha256"],
        "configuration": label,
        "residual_amplitude": amplitude,
        "noise_sd": noise_sd,
        "start": start,
        "stop": stop,
        "n_worlds": stop - start,
        "rows": rows,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "formal_gate_d_decision_authority": False,
        "phase4_identity_opening_authority": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in out.items() if key != "rows"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
