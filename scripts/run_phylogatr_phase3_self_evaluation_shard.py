#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.calibration import seed_for
from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_self_detectability import (
    prepare_genetic_self_detectability,
    score_genetic_self_world_batch,
)
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_phase3 import load_phylogatr_phase3_context
from ttf.private_null_inference import upper_monte_carlo_pvalue


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--phase3-authorization", type=Path, required=True)
    ap.add_argument("--self-rule", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--cell", choices=("null", "private_A2"), required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    context = load_phylogatr_phase3_context(
        args.geometry,
        args.phase2_manifest,
        args.phase3_rule,
        args.phase3_authorization,
    )
    rule = json.loads(args.self_rule.read_text())
    references = json.loads(args.references.read_text())
    if rule.get("schema") != "ttf_genetic_phylogatr_phase3_self_detectability_rule_v0.1":
        raise RuntimeError("fresh self rule schema drift")
    if references.get("schema") != "ttf_genetic_phylogatr_phase3_self_references_v0.1":
        raise RuntimeError("complete fresh self reference required")
    if references.get("geometry_fingerprint_sha256") != context.authorization["geometry_fingerprint_sha256"]:
        raise RuntimeError("fresh self-reference geometry fingerprint drift")
    if any(value is not False for value in rule["outcome_firewall"].values()):
        raise RuntimeError("fresh self-detectability outcome firewall is open")

    syn = rule["synthetic_worlds"]
    total = int(
        syn["null_evaluation_worlds"]
        if args.cell == "null"
        else syn["private_A2_evaluation_worlds"]
    )
    start, stop = int(args.start), int(args.stop)
    if not 0 <= start < stop <= total:
        raise RuntimeError("fresh self-evaluation range outside frozen universe")
    reference = np.asarray(references["statistics"], dtype=float)
    if len(reference) != int(syn["independent_null_reference_worlds"]):
        raise RuntimeError("fresh self-reference length drift")

    core = context.rule["core_method"]
    design = prepare_genetic_ttf_design(
        context.geometries,
        train_species=context.train_species,
        eval_species=context.eval_species,
        bandwidth=float(core["bandwidth_km"]),
        prior_strength=float(core["prior_strength"]),
        segment_points=int(core["segment_points"]),
        min_training_edges=int(context.rule["geometry_contract"]["minimum_endpoint_disjoint_ibd_training_edges"]),
        strength_neighbours=int(core["strength_neighbours"]),
    )
    geometry_rule = rule["geometry"]
    self_design = prepare_genetic_self_detectability(
        design,
        bandwidth=float(geometry_rule["bandwidth_km"]),
        prior_strength=float(geometry_rule["prior_strength"]),
        prior_mean=float(geometry_rule["prior_mean"]),
        segment_points=int(geometry_rule["segment_points"]),
    )

    world_spec = (
        syn["null_evaluation_world"]
        if args.cell == "null"
        else syn["private_A2_evaluation_world"]
    )
    seed_tag = "genetic_self_null_eval" if args.cell == "null" else "genetic_self_private_A2_eval"
    batch_width = int(rule["execution"]["world_batch_width"])
    rows: list[dict[str, float | int]] = []
    for batch_start in range(start, stop, batch_width):
        batch_stop = min(batch_start + batch_width, stop)
        worlds = [
            simulate_genetic_distance_world(
                context.geometries,
                shared_fraction=float(world_spec["shared_fraction"]),
                residual_amplitude=float(world_spec["residual_amplitude"]),
                ibd_strength=float(syn["ibd_strength"]),
                noise_sd=float(world_spec["noise_sd"]),
                transition_width=float(syn["transition_width"]),
                noise_dimensions=int(syn["noise_dimensions"]),
                seed=seed_for(context.master_seed, seed_tag, replicate),
            )
            for replicate in range(batch_start, batch_stop)
        ]
        scored = score_genetic_self_world_batch(design, self_design, worlds)
        for column, replicate in enumerate(range(batch_start, batch_stop)):
            if int(scored.n_species[column]) != len(design.eval_species):
                raise RuntimeError("non-finite species count in fresh self evaluation")
            statistic = float(scored.statistics[column])
            rows.append(
                {
                    "replicate": int(replicate),
                    "statistic": statistic,
                    "p_value": float(upper_monte_carlo_pvalue(statistic, reference)),
                }
            )
    if [row["replicate"] for row in rows] != list(range(start, stop)):
        raise RuntimeError("fresh self-evaluation ordering drift")

    out = {
        "schema": "ttf_genetic_phylogatr_phase3_self_evaluation_shard_v0.1",
        "status": "synthetic_self_detectability_evaluation_shard",
        "geometry_fingerprint_sha256": context.authorization["geometry_fingerprint_sha256"],
        "cell": args.cell,
        "start": start,
        "stop": stop,
        "n_worlds": stop - start,
        "alpha": float(rule["inference"]["alpha"]),
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
