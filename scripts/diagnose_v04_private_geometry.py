#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import seed_for
from ttf.core import SpeciesSample, edge_turnover, spearman_rho
from ttf.geometry import simulate_fixed_geometry_boundary_world
from ttf.inference import centered_species_bootstrap_mean_test
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.private_geometry_control import (
    midpoint_radial_distance,
    partial_spearman_controls,
    pooled_affine_frame,
    random_private_transition_propensity,
    rank_residualized_turnover,
)
from ttf.transfer import prepare_transfer


CELLS = ((0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0))
VARIANTS = (
    "v03_raw",
    "eval_geometry_partial",
    "train_geometry_residualized",
    "train_and_eval_geometry",
)


def mean_finite(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    return float(x.mean()) if len(x) else float("nan")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Development-only diagnosis of v0.3 private-boundary geometry leakage."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--source-ledger", type=Path, required=True)
    parser.add_argument("--v03-result", type=Path, required=True)
    parser.add_argument("--reserve-ledger", type=Path, required=True)
    parser.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    parser.add_argument("--min-records", type=int, default=80)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--bandwidth", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--resamples", type=int, default=1999)
    parser.add_argument("--propensity-directions", type=int, default=2048)
    parser.add_argument("--transition-width", type=float, default=0.2)
    parser.add_argument("--noise-sd", type=float, default=0.8)
    parser.add_argument("--segment-points", type=int, default=5)
    parser.add_argument("--prior-strength", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    v03 = json.loads(args.v03_result.read_text())
    source = json.loads(args.source_ledger.read_text())
    reserve = json.loads(args.reserve_ledger.read_text())
    if v03.get("v03_fresh_external_pass") is not False:
        raise RuntimeError("v0.4 diagnosis requires the immutable v0.3 fresh failure")
    if v03["precision_qualification"]["type1_pass"] is not False:
        raise RuntimeError("v0.3 failure is not the expected private-boundary type-I failure")
    if reserve.get("synthetic_worlds_run_at_freeze") != 0:
        raise RuntimeError("reserve holdout has been opened; stop")
    if reserve.get("candidate_performance_evaluated_at_freeze") is not False:
        raise RuntimeError("reserve candidate performance was unexpectedly evaluated")
    if source.get("confirmatory_holdout") is not True:
        raise RuntimeError("fresh animal source ledger drifted")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=args.min_records,
    )
    if excluded:
        raise RuntimeError(f"fresh geometry exclusions drifted: {excluded}")
    geometry_map = {item.species: item for item in geometries}
    train = tuple(v03["split"]["train_species"])
    evaluation = tuple(v03["split"]["eval_species"])
    used = train + evaluation
    if set(used) != set(geometry_map) or set(train) & set(evaluation):
        raise RuntimeError("v0.3 split no longer matches fresh geometry")

    templates = [
        SpeciesSample(
            species=name,
            coordinates=geometry_map[name].coordinates,
            trait=np.arange(len(geometry_map[name].coordinates), dtype=float),
        )
        for name in used
    ]
    graphs = fixed_graphs(templates, k=args.k)
    template_edges = edges_on_fixed_graphs(templates, graphs)
    prepared = prepare_transfer(
        [template_edges[name] for name in train],
        [template_edges[name] for name in evaluation],
        bandwidth=args.bandwidth,
        prior_strength=args.prior_strength,
        prior_mean=0.0,
        segment_points=args.segment_points,
    )

    center, scale, active = pooled_affine_frame(geometries)
    length = {name: template_edges[name].length for name in used}
    propensity = {}
    radial = {}
    for name in used:
        g = geometry_map[name]
        propensity[name] = random_private_transition_propensity(
            g.coordinates,
            graphs[name],
            center=center,
            scale=scale,
            active=active,
            transition_width=args.transition_width,
            n_directions=args.propensity_directions,
            seed=seed_for(args.seed, "v04_private_propensity_directions"),
        )
        radial[name] = midpoint_radial_distance(
            g.coordinates,
            graphs[name],
            center=center,
            scale=scale,
        )

    geometry_rows = []
    for name in evaluation:
        projection = prepared.eval_projection[name]
        leverage_l2 = np.sqrt(np.sum(projection * projection, axis=1))
        opportunity = prepared.eval_opportunity[name]
        geometry_rows.append(
            {
                "species": name,
                "rho_propensity_length": spearman_rho(propensity[name], length[name]),
                "rho_propensity_radial": spearman_rho(propensity[name], radial[name]),
                "rho_propensity_opportunity": spearman_rho(propensity[name], opportunity),
                "rho_propensity_projection_l2": spearman_rho(propensity[name], leverage_l2),
                "rho_radial_opportunity": spearman_rho(radial[name], opportunity),
            }
        )
    geometry_summary = {
        key: mean_finite([row[key] for row in geometry_rows])
        for key in geometry_rows[0]
        if key != "species"
    }

    frozen_v03_cells = {
        (float(cell["shared_fraction"]), float(cell["amplitude"])): cell
        for cell in v03["cells"]
    }
    cell_results = []
    for shared, amplitude in CELLS:
        variant_p = defaultdict(list)
        variant_stat = defaultdict(list)
        target_prop_rho = []
        target_length_rho = []
        prediction_prop_rho = []
        prediction_opportunity_rho = []

        for replicate in range(args.replicates):
            world = simulate_fixed_geometry_boundary_world(
                geometries,
                shared_fraction=shared,
                amplitude=amplitude,
                noise_sd=args.noise_sd,
                transition_width=args.transition_width,
                seed=seed_for(
                    args.seed,
                    "fixed_geometry_world",
                    shared,
                    amplitude,
                    replicate,
                ),
            )
            sample_map = {sample.species: sample for sample in world.samples}
            turnover = {
                name: edge_turnover(sample_map[name], graphs[name]) for name in used
            }

            v03_train = {
                name: rank_residualized_turnover(turnover[name], [length[name]])
                for name in train
            }
            geometry_train = {
                name: rank_residualized_turnover(
                    turnover[name], [length[name], propensity[name]]
                )
                for name in train
            }
            n_train_edges = max(sl.stop for sl in prepared.train_slices.values())
            v03_values = np.empty(n_train_edges, dtype=float)
            geometry_values = np.empty(n_train_edges, dtype=float)
            for name in train:
                sl = prepared.train_slices[name]
                v03_values[sl] = v03_train[name]
                geometry_values[sl] = geometry_train[name]

            scores = {variant: [] for variant in VARIANTS}
            for name in evaluation:
                projection = prepared.eval_projection[name]
                prior = prepared.eval_prior_offset[name]
                pred_v03 = projection @ v03_values + prior
                pred_geometry = projection @ geometry_values + prior
                target = turnover[name]
                controls = [length[name], propensity[name]]

                scores["v03_raw"].append(spearman_rho(pred_v03, target))
                scores["eval_geometry_partial"].append(
                    partial_spearman_controls(pred_v03, target, controls)
                )
                scores["train_geometry_residualized"].append(
                    spearman_rho(pred_geometry, target)
                )
                scores["train_and_eval_geometry"].append(
                    partial_spearman_controls(pred_geometry, target, controls)
                )
                target_prop_rho.append(spearman_rho(target, propensity[name]))
                target_length_rho.append(spearman_rho(target, length[name]))
                prediction_prop_rho.append(spearman_rho(pred_v03, propensity[name]))
                prediction_opportunity_rho.append(
                    spearman_rho(pred_v03, prepared.eval_opportunity[name])
                )

            null_seed = seed_for(
                args.seed,
                "fixed_geometry_bootstrap",
                shared,
                amplitude,
                replicate,
            )
            for variant in VARIANTS:
                arr = np.asarray(scores[variant], dtype=float)
                bootstrap = centered_species_bootstrap_mean_test(
                    arr,
                    n_bootstrap=args.resamples,
                    seed=null_seed,
                )
                variant_stat[variant].append(float(arr.mean()))
                variant_p[variant].append(float(bootstrap.p_value))

        variants = {}
        for variant in VARIANTS:
            p = np.asarray(variant_p[variant], dtype=float)
            stats = np.asarray(variant_stat[variant], dtype=float)
            variants[variant] = {
                "rejection_rate": float(np.mean(p <= 0.05)),
                "mean_statistic": float(stats.mean()),
                "median_p": float(np.median(p)),
            }
        expected = frozen_v03_cells[(shared, amplitude)]
        if not np.isclose(
            variants["v03_raw"]["rejection_rate"],
            float(expected["rejection_rate"]),
            atol=0.0,
            rtol=0.0,
        ):
            raise RuntimeError(
                f"v0.4 diagnostic failed to reproduce frozen v0.3 cell {(shared, amplitude)}"
            )
        if not np.isclose(
            variants["v03_raw"]["mean_statistic"],
            float(expected["mean_statistic"]),
            atol=1e-12,
            rtol=0.0,
        ):
            raise RuntimeError(
                f"v0.4 diagnostic statistic drift for {(shared, amplitude)}"
            )

        cell_results.append(
            {
                "shared_fraction": shared,
                "amplitude": amplitude,
                "variants": variants,
                "mechanism": {
                    "mean_target_private_propensity_rho": mean_finite(target_prop_rho),
                    "mean_target_length_rho": mean_finite(target_length_rho),
                    "mean_v03_prediction_private_propensity_rho": mean_finite(prediction_prop_rho),
                    "mean_v03_prediction_opportunity_rho": mean_finite(prediction_opportunity_rho),
                },
            }
        )

    payload = {
        "schema": "ttf_v04_private_geometry_diagnostic_v0.1",
        "status": "development_only_after_v03_fresh_failure",
        "claim_ready": False,
        "candidate_selected": False,
        "reserve_opened": False,
        "v03_failure": {
            "path": str(args.v03_result),
            "max_private_rejection": v03["precision_qualification"]["max_zero_shared_estimate"],
            "max_private_upper95": v03["precision_qualification"]["max_zero_shared_upper95"],
            "shared_power": v03["precision_qualification"]["full_shared_moderate_estimate"],
            "shared_power_lower95": v03["precision_qualification"]["full_shared_moderate_lower95"],
        },
        "diagnostic_nuisance": {
            "name": "expected_random_private_median_hyperplane_tanh_edge_contrast",
            "outcome_free": True,
            "directions": args.propensity_directions,
            "transition_width": args.transition_width,
            "frame": "exact_pooled_affine_standardization_used_by_gate_i_simulator",
            "amplitude_used": False,
        },
        "geometry_only_eval_summary": geometry_summary,
        "geometry_only_eval_species": geometry_rows,
        "design": {
            "replicates": args.replicates,
            "resamples": args.resamples,
            "k": args.k,
            "bandwidth": args.bandwidth,
            "seed": args.seed,
            "paired_worlds_with_frozen_v03": True,
            "variants": list(VARIANTS),
            "selection_rule": None,
        },
        "cells": cell_results,
        "interpretation_boundary": (
            "This reuses the failed fresh animal panel only as development data to diagnose why "
            "strong species-private transitions produce positive transfer. It cannot qualify v0.4, "
            "and the untouched RGFCA reserve remains closed."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "geometry_summary": geometry_summary,
        "cells": [
            {
                "shared": row["shared_fraction"],
                "amplitude": row["amplitude"],
                "v03": row["variants"]["v03_raw"]["rejection_rate"],
                "eval_partial": row["variants"]["eval_geometry_partial"]["rejection_rate"],
                "train_geometry": row["variants"]["train_geometry_residualized"]["rejection_rate"],
                "both": row["variants"]["train_and_eval_geometry"]["rejection_rate"],
            }
            for row in cell_results
        ],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
