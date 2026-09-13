#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.calibration import seed_for
from ttf.genetic_batch_execution import prepare_genetic_cached_transfer, score_genetic_world_batch
from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_phase3 import load_phylogatr_phase3_context
from ttf.profiled_private_null import profiled_private_pvalue


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--authorization", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--shared-fraction", type=float, required=True)
    ap.add_argument("--residual-amplitude", type=float, required=True)
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
    refs_payload = json.loads(args.references.read_text())
    if refs_payload.get("schema") != "ttf_genetic_phylogatr_phase3_references_v0.1":
        raise RuntimeError("ordered fresh phase-3 reference payload required")
    if refs_payload.get("geometry_fingerprint_sha256") != auth["geometry_fingerprint_sha256"]:
        raise RuntimeError("fresh reference payload geometry fingerprint drift")
    if refs_payload.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("fresh reference payload identity firewall drift")

    shared = float(args.shared_fraction)
    amplitude = float(args.residual_amplitude)
    qualification = rule["qualification"]
    allowed_cells = {
        (float(cell[0]), float(cell[1]))
        for cell in qualification["mandatory_primary_cells"]
    }
    if (shared, amplitude) not in allowed_cells:
        raise RuntimeError("fresh observed cell is not one of the frozen mandatory cells")
    total = int(qualification["observed_worlds_per_cell"])
    start, stop = int(args.start), int(args.stop)
    if not 0 <= start < stop <= total:
        raise RuntimeError("fresh observed shard range is outside frozen replicate range")

    references = {
        label: (
            np.asarray(data["training_strength"], dtype=float),
            np.asarray(data["statistic"], dtype=float),
        )
        for label, data in refs_payload["references"].items()
    }
    expected_labels = set(rule["synthetic_world"]["private_reference_configurations"])
    if set(references) != expected_labels:
        raise RuntimeError("fresh reference configuration set drift")

    core = rule["core_method"]
    design = prepare_genetic_ttf_design(
        context.geometries,
        train_species=context.train_species,
        eval_species=context.eval_species,
        bandwidth=float(core["bandwidth_km"]),
        prior_strength=float(core["prior_strength"]),
        segment_points=int(core["segment_points"]),
        min_training_edges=int(rule["geometry_contract"]["minimum_endpoint_disjoint_ibd_training_edges"]),
        strength_neighbours=int(core["strength_neighbours"]),
    )
    cached = prepare_genetic_cached_transfer(
        design,
        edge_chunk_size=int(auth["execution"]["edge_chunk_size"]),
        train_chunk_size=int(auth["execution"]["train_chunk_size"]),
    )

    batch_width = int(auth["execution"]["batch_width"])
    world_rule = rule["synthetic_world"]
    rows: list[dict] = []
    for batch_start in range(start, stop, batch_width):
        batch_stop = min(batch_start + batch_width, stop)
        worlds = [
            simulate_genetic_distance_world(
                context.geometries,
                shared_fraction=shared,
                residual_amplitude=amplitude,
                ibd_strength=float(world_rule["ibd_strength_primary"]),
                noise_sd=float(world_rule["locality_noise_sd"]),
                transition_width=float(world_rule["transition_width"]),
                noise_dimensions=int(world_rule["noise_dimensions"]),
                seed=seed_for(
                    context.master_seed,
                    "genetic_observed",
                    shared,
                    amplitude,
                    replicate,
                ),
            )
            for replicate in range(batch_start, batch_stop)
        ]
        scored = score_genetic_world_batch(design, worlds, cached)
        for column, replicate in enumerate(range(batch_start, batch_stop)):
            if int(scored.n_eval_species[column]) != len(design.eval_species):
                raise RuntimeError("non-finite held-out species count in fresh observed shard")
            p_value, selected, component, distances = profiled_private_pvalue(
                float(scored.statistics[column]),
                float(scored.training_strengths[column]),
                references,
                profile_draws=int(qualification["profile_strength_draws"]),
                selected_configs=int(core["profiled_private_selected_configurations"]),
            )
            rows.append(
                {
                    "replicate": int(replicate),
                    "p_value": float(p_value),
                    "statistic": float(scored.statistics[column]),
                    "training_strength": float(scored.training_strengths[column]),
                    "selected_configurations": list(selected),
                    "component_p_values": {key: float(value) for key, value in component.items()},
                    "profile_distances": {key: float(value) for key, value in distances.items()},
                }
            )

    if [row["replicate"] for row in rows] != list(range(start, stop)):
        raise RuntimeError("fresh observed shard replicate ordering drift")
    out = {
        "schema": "ttf_genetic_phylogatr_phase3_observed_shard_v0.1",
        "status": "synthetic_mandatory_observed_cell_shard",
        "geometry_fingerprint_sha256": auth["geometry_fingerprint_sha256"],
        "shared_fraction": shared,
        "residual_amplitude": amplitude,
        "start": start,
        "stop": stop,
        "n_worlds": stop - start,
        "alpha": float(qualification["alpha"]),
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
