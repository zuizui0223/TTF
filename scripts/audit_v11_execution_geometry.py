#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import seed_for
from ttf.core import SpeciesSample, split_species
from ttf.geometry import geometry_fingerprint
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs


def summary(values: list[int]) -> dict[str, float | int]:
    x = np.asarray(values, dtype=float)
    return {
        "min": int(np.min(x)),
        "q10": float(np.quantile(x, 0.10)),
        "median": float(np.median(x)),
        "q90": float(np.quantile(x, 0.90)),
        "max": int(np.max(x)),
        "mean": float(np.mean(x)),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Audit outcome-free v0.11 k=15 execution geometry.")
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--source-ledger", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    ap.add_argument("--edge-chunk-size", type=int, default=32)
    ap.add_argument("--train-chunk-size", type=int, default=4096)
    ap.add_argument("--reference-worlds", type=int, default=1999)
    ap.add_argument("--observed-worlds", type=int, default=500)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    source = json.loads(args.source_ledger.read_text())
    rule = json.loads(args.rule.read_text())
    if source.get("schema") != "ttf_v11_fresh_plant_geometry_source_v0.1":
        raise RuntimeError("wrong v0.11 source ledger")
    if source.get("fresh_geometry_freeze_pass") is not True:
        raise RuntimeError("v0.11 geometry did not freeze successfully")
    if source.get("candidate_performance_evaluated_at_freeze") is not False:
        raise RuntimeError("candidate performance already opened")
    if source.get("empirical_trait_or_colour_fields_read") is not False:
        raise RuntimeError("empirical colour already opened")
    if source.get("synthetic_worlds_run_at_freeze") != 0:
        raise RuntimeError("synthetic worlds already ran at freeze")
    if rule.get("status") != "frozen_before_v11_fresh_plant_geometry_acquisition_and_before_any_v11_synthetic_worlds":
        raise RuntimeError("v0.11 rule status drift")

    arch = rule["locked_graph_architecture"]
    exe = rule["prospective_execution"]
    k = int(arch["k"])
    if k != 15 or int(arch["records_per_species"]) != 100:
        raise RuntimeError("v0.11 graph architecture drift")
    if not np.isclose(float(arch["bandwidth_km"]), 500.0):
        raise RuntimeError("v0.11 bandwidth drift")
    master_seed = int(exe["seed"])
    if master_seed != 20260923 or not np.isclose(float(exe["eval_fraction"]), 0.5):
        raise RuntimeError("v0.11 execution seed/split drift")
    if args.reference_worlds != int(exe["reference_worlds_per_configuration"]):
        raise RuntimeError("reference world count drift")
    if args.observed_worlds != int(exe["observed_worlds_per_cell"]):
        raise RuntimeError("observed world count drift")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=100,
    )
    if excluded or len(geometries) != 250:
        raise RuntimeError("v0.11 geometry species count drift")
    if sum(len(g.coordinates) for g in geometries) != 25000:
        raise RuntimeError("v0.11 geometry record count drift")
    if any(len(g.coordinates) != 100 for g in geometries):
        raise RuntimeError("v0.11 requires exactly 100 records per species")

    labels = [g.species for g in geometries]
    split_seed = seed_for(master_seed, "fixed_geometry_split")
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=split_seed)
    if len(train) != 125 or len(evaluation) != 125 or set(train) & set(evaluation):
        raise RuntimeError("v0.11 split drift")

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
    edges = edges_on_fixed_graphs(templates, graphs)
    train_counts = [int(edges[name].n_edges) for name in train]
    eval_counts = [int(edges[name].n_edges) for name in evaluation]
    n_train_edges = int(sum(train_counts))
    n_eval_edges = int(sum(eval_counts))
    segment_points = 5

    dense_projection_elements = int(n_train_edges * n_eval_edges)
    dense_projection_bytes = int(dense_projection_elements * 8)
    segment_kernel_pairs = int(dense_projection_elements * segment_points)
    edge_chunk = int(args.edge_chunk_size)
    train_chunk = int(args.train_chunk_size)
    if edge_chunk < 1 or train_chunk < 1:
        raise RuntimeError("invalid chunk sizes")
    max_eval_edges = int(max(eval_counts))

    # Conservative peak arrays for the exact two-pass chunked implementation.
    kernel_bytes = int(edge_chunk * segment_points * train_chunk * 8)
    projection_chunk_bytes = int(edge_chunk * train_chunk * 8)
    reference_train_value_chunk_bytes = int(train_chunk * args.reference_worlds * 8)
    reference_predicted_species_bytes = int(max_eval_edges * args.reference_worlds * 8)
    observed_train_value_chunk_bytes = int(train_chunk * args.observed_worlds * 8)
    observed_predicted_species_bytes = int(max_eval_edges * args.observed_worlds * 8)

    payload = {
        "schema": "ttf_v11_execution_geometry_audit_v0.1",
        "status": "geometry_only_before_v11_synthetic_qualification",
        "candidate_selectable": False,
        "synthetic_worlds_run": 0,
        "empirical_colour_or_pixels_opened": False,
        "geometry": {
            "input": str(args.input),
            "source_ledger": str(args.source_ledger),
            "fingerprint_sha256": geometry_fingerprint(geometries),
            "species": 250,
            "records": 25000,
            "records_per_species": 100,
        },
        "split": {
            "master_seed": master_seed,
            "split_seed": int(split_seed),
            "train_species": list(train),
            "eval_species": list(evaluation),
            "train_species_count": 125,
            "eval_species_count": 125,
        },
        "graph": {
            "k": k,
            "train_edge_count": n_train_edges,
            "eval_edge_count": n_eval_edges,
            "train_edges_per_species": summary(train_counts),
            "eval_edges_per_species": summary(eval_counts),
            "segment_points": segment_points,
        },
        "dense_execution": {
            "stored_projection_elements": dense_projection_elements,
            "stored_projection_bytes_float64": dense_projection_bytes,
            "stored_projection_gib_float64": float(dense_projection_bytes / (1024 ** 3)),
            "segment_kernel_pairs_per_geometry_pass": segment_kernel_pairs,
            "reference_projection_multiply_scalar_products_per_configuration": int(dense_projection_elements * args.reference_worlds),
            "observed_projection_multiply_scalar_products_per_cell": int(dense_projection_elements * args.observed_worlds),
        },
        "chunked_exact_execution": {
            "edge_chunk_size": edge_chunk,
            "train_chunk_size": train_chunk,
            "distance_cutoff": None,
            "kernel_approximation": False,
            "dtype_reduction": False,
            "response_dependent_pruning": False,
            "kernel_chunk_bytes_float64": kernel_bytes,
            "projection_chunk_bytes_float64": projection_chunk_bytes,
            "reference_train_value_chunk_bytes_float64": reference_train_value_chunk_bytes,
            "reference_max_predicted_species_bytes_float64": reference_predicted_species_bytes,
            "observed_train_value_chunk_bytes_float64": observed_train_value_chunk_bytes,
            "observed_max_predicted_species_bytes_float64": observed_predicted_species_bytes,
            "largest_listed_transient_gib": float(max(
                kernel_bytes,
                projection_chunk_bytes,
                reference_train_value_chunk_bytes,
                reference_predicted_species_bytes,
                observed_train_value_chunk_bytes,
                observed_predicted_species_bytes,
            ) / (1024 ** 3)),
            "mathematical_operator": "identical Gaussian opportunity-corrected segment-integrated field; execution order only",
        },
        "claim_boundary": "This audit measures only the frozen v0.11 graph/operator size and deterministic split. It contains no synthetic response performance and cannot qualify the method.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "train_edges": n_train_edges,
        "eval_edges": n_eval_edges,
        "dense_gib": payload["dense_execution"]["stored_projection_gib_float64"],
        "chunk_transient_gib": payload["chunked_exact_execution"]["largest_listed_transient_gib"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
