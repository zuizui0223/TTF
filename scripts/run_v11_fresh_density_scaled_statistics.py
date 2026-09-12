#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import seed_for
from ttf.chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from ttf.core import SpeciesSample, edge_turnover, split_species
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
    ap = argparse.ArgumentParser(description="Generate frozen v0.11 fresh density-scaled qualification statistics.")
    ap.add_argument("--mode", choices=["observed", "reference"], required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--source-ledger", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--execution-audit", type=Path, required=True)
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
    audit = json.loads(args.execution_audit.read_text())
    if source.get("schema") != "ttf_v11_fresh_plant_geometry_source_v0.1" or source.get("fresh_geometry_freeze_pass") is not True:
        raise RuntimeError("v0.11 fresh geometry source is not frozen-pass")
    if source.get("candidate_performance_evaluated_at_freeze") is not False or source.get("synthetic_worlds_run_at_freeze") != 0:
        raise RuntimeError("v0.11 geometry was performance-opened before qualification")
    if source.get("empirical_trait_or_colour_fields_read") is not False or source.get("candidate_image_pixels_opened") is not False:
        raise RuntimeError("v0.11 empirical outcome firewall failed")
    if source.get("species_replacement_performed") is not False or source.get("additional_pages_after_result_performed") is not False:
        raise RuntimeError("v0.11 fresh geometry acquisition was post-hoc modified")
    if rule.get("schema") != "ttf_v11_fresh_density_scaled_graph_rule_v0.1" or rule.get("status") != "frozen_before_v11_fresh_plant_geometry_acquisition_and_before_any_v11_synthetic_worlds":
        raise RuntimeError("v0.11 rule drift")
    if audit.get("schema") != "ttf_v11_execution_geometry_audit_v0.1" or audit.get("status") != "geometry_only_before_v11_synthetic_qualification":
        raise RuntimeError("v0.11 execution audit missing or wrong")
    if audit.get("synthetic_worlds_run") != 0 or audit.get("empirical_colour_or_pixels_opened") is not False or audit.get("candidate_selectable") is not False:
        raise RuntimeError("v0.11 execution audit firewall failed")

    arch = rule["locked_graph_architecture"]
    exe = rule["prospective_execution"]
    k = int(arch["k"])
    bandwidth = float(arch["bandwidth_km"])
    master_seed = int(exe["seed"])
    eval_fraction = float(exe["eval_fraction"])
    strength_neighbours = int(exe["strength_neighbours"])
    transition_width = float(exe["transition_width"])
    if k != 15 or not np.isclose(bandwidth, 500.0) or master_seed != 20260923 or not np.isclose(eval_fraction, 0.5):
        raise RuntimeError("v0.11 locked architecture drift")
    if strength_neighbours != 4 or not np.isclose(transition_width, 0.2):
        raise RuntimeError("v0.11 nuisance/simulator constants drift")
    if int(exe["profile_strength_draws"]) != 999 or int(exe["calibration_statistic_draws"]) != 1000:
        raise RuntimeError("v0.11 reference partition drift")
    if int(exe["selected_private_configurations"]) != 2:
        raise RuntimeError("v0.11 profiled configuration count drift")
    if args.edge_chunk_size != int(audit["chunked_exact_execution"]["edge_chunk_size"]):
        raise RuntimeError("edge chunk size differs from outcome-free audit")
    if args.train_chunk_size != int(audit["chunked_exact_execution"]["train_chunk_size"]):
        raise RuntimeError("train chunk size differs from outcome-free audit")

    if args.mode == "reference":
        if args.configuration_label is None:
            raise RuntimeError("reference mode requires a configuration label")
        if args.shared_fraction is not None or args.amplitude is not None:
            raise RuntimeError("reference mode derives amplitude/noise only from the frozen label")
        amplitude, noise_sd = REFERENCE_CONFIG[args.configuration_label]
        shared_fraction = 0.0
        replicates = int(exe["reference_worlds_per_configuration"])
        if replicates != 1999:
            raise RuntimeError("v0.11 reference replicate count drift")
    else:
        if args.configuration_label is not None:
            raise RuntimeError("observed mode cannot carry a reference label")
        if args.shared_fraction is None or args.amplitude is None:
            raise RuntimeError("observed mode requires shared fraction and amplitude")
        shared_fraction = float(args.shared_fraction)
        amplitude = float(args.amplitude)
        if (shared_fraction, amplitude) not in MANDATORY:
            raise RuntimeError("observed cell is outside the frozen v0.11 mandatory set")
        noise_sd = float(exe["noise_sd"])
        if not np.isclose(noise_sd, 0.8):
            raise RuntimeError("v0.11 observed noise drift")
        replicates = int(exe["observed_worlds_per_cell"])
        if replicates != 500:
            raise RuntimeError("v0.11 observed replicate count drift")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=100,
    )
    if excluded or len(geometries) != 250 or sum(len(g.coordinates) for g in geometries) != 25000:
        raise RuntimeError("v0.11 fresh geometry size drift")
    if any(len(g.coordinates) != 100 for g in geometries):
        raise RuntimeError("v0.11 requires exactly 100 records per species")
    fingerprint = geometry_fingerprint(geometries)
    if fingerprint != audit["geometry"]["fingerprint_sha256"]:
        raise RuntimeError("v0.11 geometry fingerprint differs from execution audit")

    labels = [g.species for g in geometries]
    split_seed = seed_for(master_seed, "fixed_geometry_split")
    train, evaluation = split_species(labels, eval_fraction=eval_fraction, seed=split_seed)
    if list(train) != audit["split"]["train_species"] or list(evaluation) != audit["split"]["eval_species"]:
        raise RuntimeError("v0.11 split differs from outcome-free execution audit")
    if len(train) != 125 or len(evaluation) != 125:
        raise RuntimeError("v0.11 split size drift")

    gmap = {g.species: g for g in geometries}
    templates = [
        SpeciesSample(
            species=name,
            coordinates=gmap[name].coordinates,
            trait=np.arange(len(gmap[name].coordinates), dtype=float),
        )
        for name in train + evaluation
    ]
    graphs = fixed_graphs(templates, k=k)
    template_edges = edges_on_fixed_graphs(templates, graphs)
    if sum(template_edges[n].n_edges for n in train) != int(audit["graph"]["train_edge_count"]):
        raise RuntimeError("v0.11 train edge count drift")
    if sum(template_edges[n].n_edges for n in evaluation) != int(audit["graph"]["eval_edge_count"]):
        raise RuntimeError("v0.11 eval edge count drift")

    train_length = {name: template_edges[name].length for name in train}
    strength_index = {
        name: edge_midpoint_neighbor_indices(
            template_edges[name].midpoint,
            k=strength_neighbours,
        )
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

    # All response worlds are materialized once so the expensive geometry-only
    # Gaussian operator is traversed only once for this cell/configuration.
    train_values = {
        name: np.empty((template_edges[name].n_edges, replicates), dtype=np.float64)
        for name in train
    }
    eval_values = {
        name: np.empty((template_edges[name].n_edges, replicates), dtype=np.float64)
        for name in evaluation
    }
    strengths = np.empty(replicates, dtype=np.float64)

    for replicate in range(replicates):
        if args.mode == "reference":
            use_seed = seed_for(
                master_seed,
                "v11_fresh_density_scaled_private_reference",
                args.configuration_label,
                replicate,
            )
        else:
            use_seed = seed_for(
                master_seed,
                "v11_fresh_density_scaled_observed",
                shared_fraction,
                amplitude,
                replicate,
            )
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
        strength, _ = training_private_strength_from_indices(
            train_response,
            strength_index,
            train,
        )
        strengths[replicate] = float(strength)
        if (replicate + 1) % 100 == 0 or replicate + 1 == replicates:
            print(json.dumps({
                "stage": "response_worlds",
                "mode": args.mode,
                "label": args.configuration_label,
                "shared": shared_fraction,
                "amplitude": amplitude,
                "completed": replicate + 1,
                "replicates": replicates,
            }), flush=True)

    scored = score_chunked_batch(
        prepared,
        train_values,
        eval_values,
        edge_chunk_size=args.edge_chunk_size,
        train_chunk_size=args.train_chunk_size,
    )
    statistic = np.asarray(scored.statistics, dtype=float)
    if statistic.shape != (replicates,) or strengths.shape != statistic.shape:
        raise RuntimeError("v0.11 statistic generation incomplete")
    if not np.isfinite(statistic).all() or not np.isfinite(strengths).all():
        raise RuntimeError("v0.11 statistics contain non-finite values")

    payload = {
        "schema": "ttf_v11_fresh_density_scaled_statistics_v0.1",
        "status": "fresh_density_scaled_synthetic_qualification_only",
        "mode": args.mode,
        "geometry": {
            "input": str(args.input),
            "source_ledger": str(args.source_ledger),
            "execution_audit": str(args.execution_audit),
            "fingerprint_sha256": fingerprint,
            "species_count": 250,
            "record_count": 25000,
            "records_per_species": 100,
        },
        "split": {
            "master_seed": master_seed,
            "split_seed": int(split_seed),
            "train_species": list(train),
            "eval_species": list(evaluation),
        },
        "config": {
            "shared_fraction": shared_fraction,
            "amplitude": amplitude,
            "noise_sd": noise_sd,
            "configuration_label": args.configuration_label,
            "replicates": replicates,
            "k": k,
            "neighbor_fraction": float(arch["neighbor_fraction"]),
            "bandwidth": bandwidth,
            "strength_neighbours": strength_neighbours,
            "transition_width": transition_width,
            "training_response": "v03_length_rank_orthogonalized_turnover",
            "evaluation_score": "raw_spearman",
            "private_strength": "equal_species_mean_training_edge_midpoint_neighbor_coherence",
            "profile_strength_draws": int(exe["profile_strength_draws"]),
            "calibration_statistic_draws": int(exe["calibration_statistic_draws"]),
            "selected_private_configurations": int(exe["selected_private_configurations"]),
        },
        "execution": {
            "operator": "exact_dense_gaussian_opportunity_corrected_segment_integrated_field",
            "storage": "memory_safe_chunked_execution_only",
            "edge_chunk_size": int(args.edge_chunk_size),
            "train_chunk_size": int(args.train_chunk_size),
            "distance_cutoff": None,
            "kernel_approximation": False,
            "dtype_reduction": False,
            "response_dependent_pruning": False,
        },
        "statistics": statistic.tolist(),
        "training_strength": strengths.tolist(),
        "summary": {
            "statistic_mean": float(statistic.mean()),
            "statistic_sd": float(statistic.std(ddof=1)),
            "strength_mean": float(strengths.mean()),
            "strength_sd": float(strengths.std(ddof=1)),
            "strength_median": float(np.median(strengths)),
            "strength_mad": float(np.median(np.abs(strengths - np.median(strengths)))),
        },
        "empirical_colour_outcome_opened": False,
        "candidate_image_pixels_opened": False,
        "failed_v08_or_v09_panel_reused": False,
        "heldout_transfer_used_for_nuisance_selection": False,
        "claim_ready": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "mode": args.mode,
        "label": args.configuration_label,
        "shared": shared_fraction,
        "amplitude": amplitude,
        "replicates": replicates,
        "summary": payload["summary"],
    }, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
