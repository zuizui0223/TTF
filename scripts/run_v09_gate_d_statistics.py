#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.batch import score_prepared_batch
from ttf.calibration import seed_for
from ttf.core import SpeciesSample, edge_turnover, split_species
from ttf.geometry import geometry_fingerprint, simulate_fixed_geometry_boundary_world
from ttf.geometry_control import length_orthogonalized_turnover
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.private_strength import edge_midpoint_neighbor_indices, training_private_strength_from_indices
from ttf.transfer import prepare_transfer


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate v0.9 reserve-B Gate-D joint strength/transfer statistics.")
    ap.add_argument("--mode", choices=["observed", "reference"], required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--source-ledger", type=Path, required=True)
    ap.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    ap.add_argument("--min-records", type=int, default=100)
    ap.add_argument("--shared-fraction", type=float, required=True)
    ap.add_argument("--amplitude", type=float, required=True)
    ap.add_argument("--noise-sd", type=float, required=True)
    ap.add_argument("--configuration-label", default=None)
    ap.add_argument("--replicates", type=int, required=True)
    ap.add_argument("--world-batch-size", type=int, default=10)
    ap.add_argument("--strength-neighbours", type=int, default=4)
    ap.add_argument("--transition-width", type=float, default=0.2)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    source = json.loads(args.source_ledger.read_text())
    if source.get("schema") != "ttf_v09_rgfca_reserve_b_max_density_geometry_v0.1":
        raise RuntimeError("wrong v0.9 reserve-B source")
    if source.get("empirical_trait_or_colour_fields_read") is not False:
        raise RuntimeError("empirical colour already read")
    if source.get("candidate_performance_evaluated_at_freeze") is not False or source.get("synthetic_worlds_run_at_freeze") != 0:
        raise RuntimeError("v0.9 geometry opened before Gate-D")
    audit = source.get("source_audit", {})
    if audit.get("candidate_pixels_opened_by_audit") is not False or audit.get("colour_fields_parsed") is not False or audit.get("measurement_authorized") is not False:
        raise RuntimeError("v0.9 outcome firewall failed")
    predecessor = source.get("predecessor", {})
    if predecessor.get("species_overlap") != [] or predecessor.get("photo_overlap") != []:
        raise RuntimeError("v0.9 design is not disjoint from failed Gate-D")

    exe = source["prospective_execution"]
    k = int(exe["k"])
    bandwidth = float(exe["bandwidth_km"])
    master_seed = int(exe["seed"])
    eval_fraction = float(exe["eval_fraction"])
    if k != 3 or abs(bandwidth - 500.0) > 1e-12 or abs(eval_fraction - 0.5) > 1e-12 or master_seed != 20260920:
        raise RuntimeError("v0.9 execution contract drifted")

    if args.mode == "reference":
        if float(args.shared_fraction) != 0.0 or not args.configuration_label:
            raise RuntimeError("private reference must be shared=0 with a label")
    elif args.configuration_label is not None:
        raise RuntimeError("observed mode cannot carry configuration label")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=args.min_records,
    )
    if excluded or len(geometries) != 250 or sum(len(g.coordinates) for g in geometries) != 25000:
        raise RuntimeError("v0.9 reserve-B geometry size drift")
    if any(len(g.coordinates) != 100 for g in geometries):
        raise RuntimeError("v0.9 requires exactly 100 records per species")

    labels = [g.species for g in geometries]
    split_seed = seed_for(master_seed, "fixed_geometry_split")
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=split_seed)
    if len(train) != 125 or len(evaluation) != 125:
        raise RuntimeError("v0.9 split drift")

    gmap = {g.species: g for g in geometries}
    templates = [
        SpeciesSample(species=n, coordinates=gmap[n].coordinates, trait=np.arange(len(gmap[n].coordinates), dtype=float))
        for n in train + evaluation
    ]
    graphs = fixed_graphs(templates, k=k)
    template_edges = edges_on_fixed_graphs(templates, graphs)
    train_length = {n: template_edges[n].length for n in train}
    strength_index = {
        n: edge_midpoint_neighbor_indices(template_edges[n].midpoint, k=args.strength_neighbours)
        for n in train
    }
    prepared = prepare_transfer(
        [template_edges[n] for n in train],
        [template_edges[n] for n in evaluation],
        bandwidth=bandwidth,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )

    stats: list[float] = []
    strengths: list[float] = []
    for start in range(0, args.replicates, args.world_batch_size):
        stop = min(start + args.world_batch_size, args.replicates)
        train_cols = {n: [] for n in train}
        eval_cols = {n: [] for n in evaluation}
        batch_u: list[float] = []
        for replicate in range(start, stop):
            if args.mode == "reference":
                seed = seed_for(master_seed, "v09_gate_d_private_reference", args.configuration_label, replicate)
            else:
                seed = seed_for(master_seed, "v09_gate_d_observed", float(args.shared_fraction), float(args.amplitude), replicate)
            world = simulate_fixed_geometry_boundary_world(
                geometries,
                shared_fraction=float(args.shared_fraction),
                amplitude=float(args.amplitude),
                noise_sd=float(args.noise_sd),
                transition_width=float(args.transition_width),
                seed=seed,
            )
            smap = {s.species: s for s in world.samples}
            turnover = {n: edge_turnover(smap[n], graphs[n]) for n in train + evaluation}
            tr = {n: length_orthogonalized_turnover(turnover[n], train_length[n]) for n in train}
            u, _ = training_private_strength_from_indices(tr, strength_index, train)
            batch_u.append(float(u))
            for n in train:
                train_cols[n].append(tr[n])
            for n in evaluation:
                eval_cols[n].append(turnover[n])
        scored = score_prepared_batch(
            prepared,
            {n: np.column_stack(train_cols[n]) for n in train},
            {n: np.column_stack(eval_cols[n]) for n in evaluation},
        )
        stats.extend(map(float, scored.statistics))
        strengths.extend(batch_u)

    t = np.asarray(stats, dtype=float)
    u = np.asarray(strengths, dtype=float)
    if t.shape != (args.replicates,) or u.shape != t.shape or not np.isfinite(t).all() or not np.isfinite(u).all():
        raise RuntimeError("incomplete v0.9 Gate-D statistics")

    payload = {
        "schema": "ttf_v09_gate_d_statistics_v0.1",
        "status": "reserve_b_max_density_synthetic_deployment_adequacy_only",
        "mode": args.mode,
        "geometry": {
            "input": str(args.input),
            "source_ledger": str(args.source_ledger),
            "fingerprint_sha256": geometry_fingerprint(geometries),
            "species_count": 250,
            "record_count": 25000,
            "records_per_species": 100,
        },
        "split": {"split_seed": split_seed, "train_species": list(train), "eval_species": list(evaluation)},
        "config": {
            "shared_fraction": float(args.shared_fraction),
            "amplitude": float(args.amplitude),
            "noise_sd": float(args.noise_sd),
            "configuration_label": args.configuration_label,
            "replicates": args.replicates,
            "k": k,
            "bandwidth": bandwidth,
            "master_seed": master_seed,
            "strength_neighbours": args.strength_neighbours,
            "transition_width": args.transition_width,
        },
        "statistics": stats,
        "training_strength": strengths,
        "summary": {
            "statistic_mean": float(t.mean()),
            "statistic_sd": float(t.std(ddof=1)),
            "strength_mean": float(u.mean()),
            "strength_sd": float(u.std(ddof=1)),
            "strength_median": float(np.median(u)),
            "strength_mad": float(np.median(np.abs(u - np.median(u)))),
        },
        "empirical_colour_outcome_opened": False,
        "failed_v08_gate_d_panel_reused": False,
        "claim_ready": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"mode": args.mode, "shared": args.shared_fraction, "amplitude": args.amplitude, "label": args.configuration_label, "summary": payload["summary"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
