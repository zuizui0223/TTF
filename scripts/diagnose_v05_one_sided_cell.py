#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import CalibrationCell, seed_for
from ttf.core import SpeciesSample, edge_turnover
from ttf.geometry import simulate_fixed_geometry_boundary_world
from ttf.geometry_control import length_orthogonalized_turnover
from ttf.inference import centered_species_bootstrap_mean_test
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.one_sided_geometry_control import residual_rank_correlation
from ttf.private_geometry_control import (
    pooled_affine_frame,
    random_private_transition_propensity,
)
from ttf.transfer import prepare_transfer

VARIANTS = (
    "prediction_propensity_only",
    "prediction_length_and_propensity",
    "target_propensity_only",
    "two_sided_propensity_only",
    "v04_two_sided_length_and_propensity",
)


def reference_v04_cell(panel: str, result: dict, shared: float, amplitude: float) -> dict:
    if panel == "v04":
        cells = result["cells"]
    elif panel == "v03":
        cells = result["formal_seed_cells"]
    else:
        raise ValueError(panel)
    for cell in cells:
        if np.isclose(float(cell["shared_fraction"]), shared) and np.isclose(float(cell["amplitude"]), amplitude):
            return cell
    raise RuntimeError(f"reference cell missing for {panel} {(shared, amplitude)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one paired v0.5 side-specific geometry diagnostic cell.")
    parser.add_argument("--panel", choices=["v03", "v04"], required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--split-result", type=Path, required=True)
    parser.add_argument("--baseline-result", type=Path, required=True)
    parser.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    parser.add_argument("--min-records", type=int, default=80)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--resamples", type=int, default=1999)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--bandwidth", type=float, required=True)
    parser.add_argument("--world-seed", type=int, required=True)
    parser.add_argument("--propensity-seed", type=int, default=20260908)
    parser.add_argument("--propensity-directions", type=int, default=2048)
    parser.add_argument("--transition-width", type=float, default=0.2)
    parser.add_argument("--noise-sd", type=float, default=0.8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    split_result = json.loads(args.split_result.read_text())
    baseline_result = json.loads(args.baseline_result.read_text())
    if args.panel == "v03":
        if split_result.get("v03_fresh_external_pass") is not False:
            raise RuntimeError("v03 panel is not a failed development panel")
        if baseline_result.get("candidate_lock_valid_on_development_rule") is not True:
            raise RuntimeError("formal v04 seed baseline missing for v03 panel")
    else:
        if split_result.get("v04_fresh_external_pass") is not False:
            raise RuntimeError("v04 panel is not a failed development panel")
        if split_result["precision_qualification"]["power_pass"] is not False:
            raise RuntimeError("v04 failure is not the expected power failure")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=args.min_records,
    )
    if excluded:
        raise RuntimeError(f"unexpected excluded taxa: {excluded}")
    gmap = {g.species: g for g in geometries}
    train = tuple(split_result["split"]["train_species"])
    evaluation = tuple(split_result["split"]["eval_species"])
    if set(train) & set(evaluation) or set(train + evaluation) != set(gmap):
        raise RuntimeError("frozen split does not match geometry")

    templates = [
        SpeciesSample(
            species=name,
            coordinates=gmap[name].coordinates,
            trait=np.arange(len(gmap[name].coordinates), dtype=float),
        )
        for name in train + evaluation
    ]
    graphs = fixed_graphs(templates, k=args.k)
    template_edges = edges_on_fixed_graphs(templates, graphs)
    train_length = {name: template_edges[name].length for name in train}
    eval_length = {name: template_edges[name].length for name in evaluation}
    center, scale, active = pooled_affine_frame(geometries)
    eval_propensity = {
        name: random_private_transition_propensity(
            gmap[name].coordinates,
            graphs[name],
            center=center,
            scale=scale,
            active=active,
            transition_width=args.transition_width,
            n_directions=args.propensity_directions,
            seed=args.propensity_seed,
        )
        for name in evaluation
    }
    prepared = prepare_transfer(
        [template_edges[name] for name in train],
        [template_edges[name] for name in evaluation],
        bandwidth=args.bandwidth,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )
    n_train_edges = max(sl.stop for sl in prepared.train_slices.values())

    p_values = {variant: [] for variant in VARIANTS}
    statistics = {variant: [] for variant in VARIANTS}
    null_means = {variant: [] for variant in VARIANTS}

    for replicate in range(args.replicates):
        world = simulate_fixed_geometry_boundary_world(
            geometries,
            shared_fraction=args.shared_fraction,
            amplitude=args.amplitude,
            noise_sd=args.noise_sd,
            transition_width=args.transition_width,
            seed=seed_for(
                args.world_seed,
                "fixed_geometry_world",
                args.shared_fraction,
                args.amplitude,
                replicate,
            ),
        )
        sample_map = {sample.species: sample for sample in world.samples}
        turnover = {
            name: edge_turnover(sample_map[name], graphs[name])
            for name in train + evaluation
        }
        flat = np.empty(n_train_edges, dtype=float)
        for name in train:
            sl = prepared.train_slices[name]
            flat[sl] = length_orthogonalized_turnover(turnover[name], train_length[name])

        scores = {variant: [] for variant in VARIANTS}
        for name in evaluation:
            prediction = prepared.eval_projection[name] @ flat + prepared.eval_prior_offset[name]
            target = turnover[name]
            length = eval_length[name]
            propensity = eval_propensity[name]
            scores["prediction_propensity_only"].append(
                residual_rank_correlation(
                    prediction, target,
                    prediction_nuisances=[propensity],
                    target_nuisances=[],
                )
            )
            scores["prediction_length_and_propensity"].append(
                residual_rank_correlation(
                    prediction, target,
                    prediction_nuisances=[length, propensity],
                    target_nuisances=[],
                )
            )
            scores["target_propensity_only"].append(
                residual_rank_correlation(
                    prediction, target,
                    prediction_nuisances=[],
                    target_nuisances=[propensity],
                )
            )
            scores["two_sided_propensity_only"].append(
                residual_rank_correlation(
                    prediction, target,
                    prediction_nuisances=[propensity],
                    target_nuisances=[propensity],
                )
            )
            scores["v04_two_sided_length_and_propensity"].append(
                residual_rank_correlation(
                    prediction, target,
                    prediction_nuisances=[length, propensity],
                    target_nuisances=[length, propensity],
                )
            )

        null_seed = seed_for(
            args.world_seed,
            "fixed_geometry_bootstrap",
            args.shared_fraction,
            args.amplitude,
            replicate,
        )
        for variant in VARIANTS:
            arr = np.asarray(scores[variant], dtype=float)
            if len(arr) != len(evaluation) or not np.isfinite(arr).all():
                raise RuntimeError(f"non-finite side-specific scores for {variant}")
            boot = centered_species_bootstrap_mean_test(
                arr,
                n_bootstrap=args.resamples,
                seed=null_seed,
            )
            statistics[variant].append(float(arr.mean()))
            p_values[variant].append(float(boot.p_value))
            null_means[variant].append(float(boot.null_means.mean()))

    variants = {}
    for variant in VARIANTS:
        p = np.asarray(p_values[variant], dtype=float)
        variants[variant] = CalibrationCell(
            shared_fraction=float(args.shared_fraction),
            amplitude=float(args.amplitude),
            n_replicates=args.replicates,
            alpha=0.05,
            rejection_rate=float(np.mean(p <= 0.05)),
            mean_statistic=float(np.mean(statistics[variant])),
            mean_null_statistic=float(np.mean(null_means[variant])),
            median_p_value=float(np.median(p)),
        ).to_dict()

    reference = reference_v04_cell(
        args.panel,
        baseline_result,
        float(args.shared_fraction),
        float(args.amplitude),
    )
    baseline = variants["v04_two_sided_length_and_propensity"]
    if not np.isclose(float(baseline["rejection_rate"]), float(reference["rejection_rate"]), atol=0.0, rtol=0.0):
        raise RuntimeError(
            f"v05 diagnostic did not reproduce frozen/formal v04 rejection: {baseline['rejection_rate']} != {reference['rejection_rate']}"
        )
    if not np.isclose(float(baseline["mean_statistic"]), float(reference["mean_statistic"]), atol=1e-12, rtol=0.0):
        raise RuntimeError("v05 diagnostic did not reproduce frozen/formal v04 mean statistic")

    payload = {
        "schema": "ttf_v05_one_sided_diagnostic_cell_v0.1",
        "status": "development_only",
        "panel": args.panel,
        "geometry": str(args.input),
        "split_result": str(args.split_result),
        "baseline_result": str(args.baseline_result),
        "shared_fraction": float(args.shared_fraction),
        "amplitude": float(args.amplitude),
        "config": {
            "replicates": args.replicates,
            "resamples": args.resamples,
            "k": args.k,
            "bandwidth": args.bandwidth,
            "world_seed": args.world_seed,
            "propensity_seed": args.propensity_seed,
            "propensity_directions": args.propensity_directions,
            "transition_width": args.transition_width,
            "noise_sd": args.noise_sd,
        },
        "variants": variants,
        "baseline_reproduced_exactly": True,
        "claim_ready": False,
        "rgfca_reserve_opened": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "panel": args.panel,
        "shared": args.shared_fraction,
        "amplitude": args.amplitude,
        "rejections": {k:v["rejection_rate"] for k,v in variants.items()},
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
