#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import seed_for
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_compact_execution import (
    prepare_phylogatr_compact_cached_transfer,
    prepare_phylogatr_compact_ttf_design,
    score_phylogatr_compact_world_batch,
)
from ttf.phylogatr_phase3 import load_phylogatr_phase3_context


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--authorization", type=Path, required=True)
    ap.add_argument("--config-label", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    context = load_phylogatr_phase3_context(
        args.geometry,
        args.phase2_manifest,
        args.phase3_rule,
        args.authorization,
    )
    rule = context.rule
    auth = context.authorization
    configs = rule["synthetic_world"]["private_reference_configurations"]
    label = str(args.config_label)
    if label not in configs:
        raise RuntimeError(f"unknown private reference configuration: {label}")
    total = int(rule["qualification"]["reference_worlds_per_configuration"])
    start, stop = int(args.start), int(args.stop)
    if not 0 <= start < stop <= total:
        raise RuntimeError("reference shard range is outside frozen replicate range")

    batch_width = int(auth["execution"]["batch_width"])
    if batch_width < 1:
        raise RuntimeError("invalid authorized batch width")
    core = rule["core_method"]
    design = prepare_phylogatr_compact_ttf_design(
        context.geometries,
        train_species=context.train_species,
        eval_species=context.eval_species,
        bandwidth=float(core["bandwidth_km"]),
        prior_strength=float(core["prior_strength"]),
        segment_points=int(core["segment_points"]),
        min_training_edges=int(rule["geometry_contract"]["minimum_endpoint_disjoint_ibd_training_edges"]),
        strength_neighbours=int(core["strength_neighbours"]),
    )
    cached = prepare_phylogatr_compact_cached_transfer(
        design,
        edge_chunk_size=int(auth["execution"]["edge_chunk_size"]),
        train_chunk_size=int(auth["execution"]["train_chunk_size"]),
    )

    amplitude, noise_sd = map(float, configs[label])
    world_rule = rule["synthetic_world"]
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
                    "genetic_private_reference",
                    label,
                    replicate,
                ),
            )
            for replicate in range(batch_start, batch_stop)
        ]
        scored = score_phylogatr_compact_world_batch(design, worlds, cached)
        for column, replicate in enumerate(range(batch_start, batch_stop)):
            if int(scored.n_eval_species[column]) != len(design.eval_species):
                raise RuntimeError("non-finite held-out species count in fresh reference shard")
            rows.append(
                {
                    "replicate": int(replicate),
                    "training_strength": float(scored.training_strengths[column]),
                    "statistic": float(scored.statistics[column]),
                }
            )

    if [row["replicate"] for row in rows] != list(range(start, stop)):
        raise RuntimeError("fresh reference shard replicate ordering drift")
    out = {
        "schema": "ttf_genetic_phylogatr_phase3_reference_shard_v0.1",
        "status": "synthetic_private_reference_shard",
        "geometry_fingerprint_sha256": auth["geometry_fingerprint_sha256"],
        "configuration": label,
        "residual_amplitude": amplitude,
        "noise_sd": noise_sd,
        "start": start,
        "stop": stop,
        "n_worlds": stop - start,
        "rows": rows,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "qualification_claim_made": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in out.items() if key != "rows"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
