#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import seed_for
from ttf.chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from ttf.core import SpeciesSample, edge_turnover
from ttf.geometry import geometry_fingerprint, simulate_fixed_geometry_boundary_world
from ttf.geometry_control import length_orthogonalized_turnover
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.private_strength import edge_midpoint_neighbor_indices, training_private_strength_from_indices

REFERENCE_CONFIG = {
    "A0": (0.0, 0.8),
    "A0p5": (0.5, 0.8),
    "A1": (1.0, 0.8),
    "A2": (2.0, 0.8),
    "A3": (3.0, 0.8),
    "A5": (5.0, 0.8),
    "A10": (10.0, 0.8),
    "infinite_snr": (1.0, 0.0),
}
MANDATORY = {(0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0)}


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate frozen v0.13 minimum-geometry synthetic requalification statistics.")
    ap.add_argument("--mode", choices=["observed", "reference"], required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--source-ledger", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--v11-execution-audit", type=Path, required=True)
    ap.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    ap.add_argument("--shared-fraction", type=float, default=None)
    ap.add_argument("--amplitude", type=float, default=None)
    ap.add_argument("--configuration-label", choices=tuple(REFERENCE_CONFIG), default=None)
    ap.add_argument("--edge-chunk-size", type=int, default=32)
    ap.add_argument("--train-chunk-size", type=int, default=4096)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    source = json.loads(args.source_ledger.read_text())
    rule = json.loads(args.rule.read_text())
    v11_audit = json.loads(args.v11_execution_audit.read_text())
    if source.get("schema") != "ttf_v13_minimum_classifiability_geometry_result_v0.1" or source.get("classifiability_gate_pass") is not True:
        raise RuntimeError("v0.13 minimum geometry source is not a frozen pass")
    if source.get("species_specific_classifiability_counts_persisted") is not False:
        raise RuntimeError("v0.13 classifiability surface leaked species-specific counts")
    if source.get("colour_vector_values_read_by_ttf") is not False or source.get("palette_fractions_read_by_ttf") is not False or source.get("pairwise_colour_distances_computed") is not False or source.get("ttf_empirical_statistic_computed") is not False:
        raise RuntimeError("v0.13 empirical outcome firewall failed")
    if source.get("new_image_acquisition_performed") is not False or source.get("species_replacement_performed") is not False:
        raise RuntimeError("v0.13 inherited denominator was modified")
    if rule.get("schema") != "ttf_v13_minimum_classifiability_geometry_rule_v0.1" or rule.get("status") != "frozen_after_v12_classifiability_fail_before_any_empirical_colour_values_are_opened":
        raise RuntimeError("v0.13 rule drift")
    if v11_audit.get("schema") != "ttf_v11_execution_geometry_audit_v0.1":
        raise RuntimeError("v0.11 split authority missing")

    arch = rule["actual_geometry_architecture"]
    exe = rule["synthetic_requalification"]
    k = int(arch["k"])
    bandwidth = float(arch["bandwidth_km"])
    master_seed = int(exe["master_seed"])
    if k != 2 or not np.isclose(bandwidth, 500.0) or master_seed != 20260927:
        raise RuntimeError("v0.13 architecture drift")
    if int(arch["records_per_species"]) != 15 or int(arch["train_species"]) != 125 or int(arch["eval_species"]) != 125:
        raise RuntimeError("v0.13 retained design drift")
    if int(exe["profile_strength_draws"]) != 999 or int(exe["calibration_statistic_draws"]) != 1000 or int(exe["selected_private_configurations"]) != 2:
        raise RuntimeError("v0.13 profiling rule drift")

    if args.mode == "reference":
        if args.configuration_label is None or args.shared_fraction is not None or args.amplitude is not None:
            raise RuntimeError("reference mode requires only a frozen configuration label")
        amplitude, noise_sd = REFERENCE_CONFIG[args.configuration_label]
        shared_fraction = 0.0
        replicates = int(exe["reference_worlds_per_configuration"])
        seed_namespace = str(exe["reference_seed_namespace"])
        if replicates != 1999:
            raise RuntimeError("v0.13 reference count drift")
    else:
        if args.configuration_label is not None or args.shared_fraction is None or args.amplitude is None:
            raise RuntimeError("observed mode requires shared fraction/amplitude only")
        shared_fraction = float(args.shared_fraction)
        amplitude = float(args.amplitude)
        if (shared_fraction, amplitude) not in MANDATORY:
            raise RuntimeError("observed cell outside frozen v0.13 mandatory set")
        noise_sd = 0.8
        replicates = int(exe["observed_worlds_per_cell"])
        seed_namespace = str(exe["observed_seed_namespace"])
        if replicates != 500:
            raise RuntimeError("v0.13 observed count drift")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=15,
    )
    if excluded or len(geometries) != 250 or sum(len(g.coordinates) for g in geometries) != 3750:
        raise RuntimeError("v0.13 actual geometry size drift")
    if any(len(g.coordinates) != 15 for g in geometries):
        raise RuntimeError("v0.13 requires exactly 15 records/species")

    labels = [g.species for g in geometries]
    expected_train = list(v11_audit["split"]["train_species"])
    expected_eval = list(v11_audit["split"]["eval_species"])
    if len(expected_train) != 125 or len(expected_eval) != 125 or set(expected_train + expected_eval) != set(labels):
        raise RuntimeError("v0.13 species set does not exactly support the frozen v0.11 split")
    train = expected_train
    evaluation = expected_eval

    gmap = {g.species: g for g in geometries}
    templates = [
        SpeciesSample(species=name, coordinates=gmap[name].coordinates, trait=np.arange(15, dtype=float))
        for name in train + evaluation
    ]
    graphs = fixed_graphs(templates, k=k)
    template_edges = edges_on_fixed_graphs(templates, graphs)
    train_length = {name: template_edges[name].length for name in train}
    strength_neighbours = 4
    strength_index = {
        name: edge_midpoint_neighbor_indices(template_edges[name].midpoint, k=strength_neighbours)
        for name in train
    }
    prepared = prepare_chunked_transfer(
        [template_edges[name] for name in train],
        [template_edges[name] for name in evaluation],
        bandwidth=bandwidth,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )

    train_values = {
        name: np.empty((template_edges[name].n_edges, replicates), dtype=np.float64)
        for name in train
    }
    eval_values = {
        name: np.empty((template_edges[name].n_edges, replicates), dtype=np.float64)
        for name in evaluation
    }
    strengths = np.empty(replicates, dtype=np.float64)
    transition_width = 0.2

    for replicate in range(replicates):
        if args.mode == "reference":
            use_seed = seed_for(master_seed, seed_namespace, args.configuration_label, replicate)
        else:
            use_seed = seed_for(master_seed, seed_namespace, shared_fraction, amplitude, replicate)
        world = simulate_fixed_geometry_boundary_world(
            geometries,
            shared_fraction=shared_fraction,
            amplitude=amplitude,
            noise_sd=noise_sd,
            transition_width=transition_width,
            seed=use_seed,
        )
        smap = {sample.species: sample for sample in world.samples}
        train_response: dict[str, np.ndarray] = {}
        for name in train:
            turnover = edge_turnover(smap[name], graphs[name])
            response = length_orthogonalized_turnover(turnover, train_length[name])
            train_values[name][:, replicate] = response
            train_response[name] = response
        for name in evaluation:
            eval_values[name][:, replicate] = edge_turnover(smap[name], graphs[name])
        strength, _ = training_private_strength_from_indices(train_response, strength_index, train)
        strengths[replicate] = float(strength)
        if (replicate + 1) % 100 == 0 or replicate + 1 == replicates:
            print(json.dumps({"stage":"response_worlds","mode":args.mode,"completed":replicate+1,"replicates":replicates}), flush=True)

    scored = score_chunked_batch(
        prepared,
        train_values,
        eval_values,
        edge_chunk_size=int(args.edge_chunk_size),
        train_chunk_size=int(args.train_chunk_size),
    )
    statistic = np.asarray(scored.statistics, dtype=float)
    if statistic.shape != (replicates,) or strengths.shape != statistic.shape or not np.isfinite(statistic).all() or not np.isfinite(strengths).all():
        raise RuntimeError("v0.13 statistic generation incomplete/nonfinite")

    payload = {
        "schema": "ttf_v13_minimum_geometry_statistics_v0.1",
        "status": "synthetic_requalification_on_minimum_classifiability_geometry_only",
        "mode": args.mode,
        "geometry": {
            "input": str(args.input),
            "source_ledger": str(args.source_ledger),
            "fingerprint_sha256": geometry_fingerprint(geometries),
            "species_count": 250,
            "record_count": 3750,
            "records_per_species": 15
        },
        "split": {"authority":"results/v11_execution_geometry_audit_v0.1.json","train_species":train,"eval_species":evaluation},
        "config": {
            "shared_fraction": shared_fraction,
            "amplitude": amplitude,
            "noise_sd": noise_sd,
            "configuration_label": args.configuration_label,
            "replicates": replicates,
            "k": k,
            "bandwidth": bandwidth,
            "master_seed": master_seed,
            "seed_namespace": seed_namespace,
            "strength_neighbours": strength_neighbours,
            "transition_width": transition_width
        },
        "execution": {
            "operator": "exact_dense_gaussian_opportunity_corrected_segment_integrated_field",
            "storage": "chunked_exact_only",
            "edge_chunk_size": int(args.edge_chunk_size),
            "train_chunk_size": int(args.train_chunk_size),
            "kernel_approximation": False,
            "distance_cutoff": None,
            "dtype_reduction": False,
            "response_dependent_pruning": False
        },
        "statistics": statistic.tolist(),
        "training_strength": strengths.tolist(),
        "empirical_colour_values_read": False,
        "palette_fractions_read": False,
        "pairwise_colour_distances_computed": False,
        "empirical_ttf_statistic_computed": False,
        "heldout_transfer_used_for_nuisance_selection": False,
        "claim_ready": False
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"stage":"synthetic_cell_complete","mode":args.mode,"configuration_label":args.configuration_label,"shared_fraction":shared_fraction,"amplitude":amplitude}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
