#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import seed_for
from ttf.genetic_self_detectability import prepare_genetic_self_detectability
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_compact_execution import prepare_phylogatr_compact_ttf_design
from ttf.phylogatr_compact_self_detectability import (
    score_phylogatr_compact_self_world_batch,
)
from ttf.phylogatr_phase3 import load_phylogatr_phase3_context


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--phase3-authorization", type=Path, required=True)
    ap.add_argument("--self-rule", type=Path, required=True)
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
    self_rule = json.loads(args.self_rule.read_text())
    if self_rule.get("schema") != "ttf_genetic_phylogatr_phase3_self_detectability_rule_v0.1":
        raise RuntimeError("fresh self-detectability rule schema drift")
    if any(value is not False for value in self_rule["outcome_firewall"].values()):
        raise RuntimeError("fresh self-detectability outcome firewall is open")

    syn = self_rule["synthetic_worlds"]
    total = int(syn["independent_null_reference_worlds"])
    start, stop = int(args.start), int(args.stop)
    if not 0 <= start < stop <= total:
        raise RuntimeError("fresh self-reference range outside frozen universe")

    core = context.rule["core_method"]
    design = prepare_phylogatr_compact_ttf_design(
        context.geometries,
        train_species=context.train_species,
        eval_species=context.eval_species,
        bandwidth=float(core["bandwidth_km"]),
        prior_strength=float(core["prior_strength"]),
        segment_points=int(core["segment_points"]),
        min_training_edges=int(context.rule["geometry_contract"]["minimum_endpoint_disjoint_ibd_training_edges"]),
        strength_neighbours=int(core["strength_neighbours"]),
    )
    geometry_rule = self_rule["geometry"]
    self_design = prepare_genetic_self_detectability(
        design,
        bandwidth=float(geometry_rule["bandwidth_km"]),
        prior_strength=float(geometry_rule["prior_strength"]),
        prior_mean=float(geometry_rule["prior_mean"]),
        segment_points=int(geometry_rule["segment_points"]),
    )

    batch_width = int(self_rule["execution"]["world_batch_width"])
    spec = syn["reference_world"]
    rows: list[dict[str, float | int]] = []
    for batch_start in range(start, stop, batch_width):
        batch_stop = min(batch_start + batch_width, stop)
        worlds = [
            simulate_genetic_distance_world(
                context.geometries,
                shared_fraction=float(spec["shared_fraction"]),
                residual_amplitude=float(spec["residual_amplitude"]),
                ibd_strength=float(syn["ibd_strength"]),
                noise_sd=float(spec["noise_sd"]),
                transition_width=float(syn["transition_width"]),
                noise_dimensions=int(syn["noise_dimensions"]),
                seed=seed_for(context.master_seed, "genetic_self_reference", replicate),
            )
            for replicate in range(batch_start, batch_stop)
        ]
        scored = score_phylogatr_compact_self_world_batch(design, self_design, worlds)
        for column, replicate in enumerate(range(batch_start, batch_stop)):
            if int(scored.n_species[column]) != len(design.eval_species):
                raise RuntimeError("non-finite species count in fresh self reference")
            rows.append(
                {"replicate": int(replicate), "statistic": float(scored.statistics[column])}
            )

    if [row["replicate"] for row in rows] != list(range(start, stop)):
        raise RuntimeError("fresh self-reference ordering drift")
    out = {
        "schema": "ttf_genetic_phylogatr_phase3_self_reference_shard_v0.1",
        "status": "independent_null_reference_shard",
        "geometry_fingerprint_sha256": context.authorization["geometry_fingerprint_sha256"],
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
