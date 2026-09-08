#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.batch import score_prepared_batch
from ttf.calibration import seed_for
from ttf.core import SpeciesSample, edge_turnover
from ttf.geometry import geometry_fingerprint, simulate_fixed_geometry_boundary_world
from ttf.geometry_control import length_orthogonalized_turnover
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.transfer import prepare_transfer


def _baseline_cell(data: dict, panel: str, shared: float, amplitude: float) -> dict:
    cells = data["nonselectable_v03_raw_benchmark"][panel]["cells"]
    for cell in cells:
        if np.isclose(float(cell["shared_fraction"]), shared) and np.isclose(
            float(cell["amplitude"]), amplitude
        ):
            return cell
    raise RuntimeError(f"v0.5 raw baseline missing for {panel} {(shared, amplitude)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate raw v0.3 transfer statistics for v0.6 private-null development."
    )
    parser.add_argument("--mode", choices=["observed", "reference"], required=True)
    parser.add_argument("--panel", choices=["v03", "v04"], required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--split-result", type=Path, required=True)
    parser.add_argument("--baseline-v05", type=Path, required=True)
    parser.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    parser.add_argument("--min-records", type=int, default=80)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--noise-sd", type=float, required=True)
    parser.add_argument("--configuration-label", default=None)
    parser.add_argument("--replicates", type=int, required=True)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--bandwidth", type=float, required=True)
    parser.add_argument("--world-seed", type=int, required=True)
    parser.add_argument("--world-batch-size", type=int, default=25)
    parser.add_argument("--transition-width", type=float, default=0.2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.mode == "reference":
        if not np.isclose(args.shared_fraction, 0.0):
            raise RuntimeError("v0.6 private reference must have shared_fraction=0")
        if not args.configuration_label:
            raise RuntimeError("reference mode requires --configuration-label")
    elif args.configuration_label is not None:
        raise RuntimeError("observed mode must not carry a private configuration label")

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

    split_result = json.loads(args.split_result.read_text())
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
    prepared = prepare_transfer(
        [template_edges[name] for name in train],
        [template_edges[name] for name in evaluation],
        bandwidth=args.bandwidth,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )

    statistics: list[float] = []
    for batch_start in range(0, args.replicates, args.world_batch_size):
        batch_stop = min(batch_start + args.world_batch_size, args.replicates)
        train_columns = {name: [] for name in train}
        eval_columns = {name: [] for name in evaluation}
        for replicate in range(batch_start, batch_stop):
            if args.mode == "observed":
                world_seed = seed_for(
                    args.world_seed,
                    "fixed_geometry_world",
                    float(args.shared_fraction),
                    float(args.amplitude),
                    replicate,
                )
            else:
                world_seed = seed_for(
                    args.world_seed,
                    "v06_private_reference",
                    args.configuration_label,
                    replicate,
                )
            world = simulate_fixed_geometry_boundary_world(
                geometries,
                shared_fraction=float(args.shared_fraction),
                amplitude=float(args.amplitude),
                noise_sd=float(args.noise_sd),
                transition_width=float(args.transition_width),
                seed=world_seed,
            )
            sample_map = {sample.species: sample for sample in world.samples}
            turnover = {
                name: edge_turnover(sample_map[name], graphs[name])
                for name in train + evaluation
            }
            for name in train:
                train_columns[name].append(
                    length_orthogonalized_turnover(turnover[name], train_length[name])
                )
            for name in evaluation:
                eval_columns[name].append(turnover[name])

        scored = score_prepared_batch(
            prepared,
            {name: np.column_stack(train_columns[name]) for name in train},
            {name: np.column_stack(eval_columns[name]) for name in evaluation},
        )
        statistics.extend(map(float, scored.statistics))

    values = np.asarray(statistics, dtype=float)
    if len(values) != args.replicates or not np.isfinite(values).all():
        raise RuntimeError("v0.6 statistic generation was incomplete")

    baseline_reproduced = None
    if args.mode == "observed":
        baseline_data = json.loads(args.baseline_v05.read_text())
        baseline = _baseline_cell(
            baseline_data,
            args.panel,
            float(args.shared_fraction),
            float(args.amplitude),
        )
        baseline_reproduced = bool(
            np.isclose(values.mean(), float(baseline["mean_statistic"]), atol=1e-12, rtol=0.0)
        )
        if not baseline_reproduced:
            raise RuntimeError(
                f"raw statistic drift: {values.mean()} != {baseline['mean_statistic']}"
            )

    payload = {
        "schema": "ttf_v06_private_null_statistics_v0.1",
        "status": "development_only",
        "mode": args.mode,
        "panel": args.panel,
        "geometry": {
            "input": str(args.input),
            "fingerprint_sha256": geometry_fingerprint(geometries),
            "species_count": len(geometries),
            "record_count": int(sum(len(g.coordinates) for g in geometries)),
        },
        "split": {
            "train_species": list(train),
            "eval_species": list(evaluation),
        },
        "config": {
            "shared_fraction": float(args.shared_fraction),
            "amplitude": float(args.amplitude),
            "noise_sd": float(args.noise_sd),
            "configuration_label": args.configuration_label,
            "replicates": args.replicates,
            "k": args.k,
            "bandwidth": args.bandwidth,
            "world_seed": args.world_seed,
            "world_batch_size": args.world_batch_size,
            "transition_width": args.transition_width,
            "training_response": "v03_length_rank_orthogonalized_turnover",
            "evaluation_score": "raw_spearman",
            "statistic": "equal_species_mean_heldout_spearman",
        },
        "statistics": statistics,
        "summary": {
            "mean": float(values.mean()),
            "sd": float(values.std(ddof=1)),
            "median": float(np.median(values)),
            "q05": float(np.quantile(values, 0.05)),
            "q50": float(np.quantile(values, 0.50)),
            "q95": float(np.quantile(values, 0.95)),
            "min": float(values.min()),
            "max": float(values.max()),
        },
        "baseline_v05_reproduced": baseline_reproduced,
        "candidate_selectable_by_this_file": False,
        "confirmatory_performance_opened": False,
        "rgfca_reserve_opened": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "mode": args.mode,
        "panel": args.panel,
        "label": args.configuration_label,
        "shared": args.shared_fraction,
        "amplitude": args.amplitude,
        "noise_sd": args.noise_sd,
        "summary": payload["summary"],
        "baseline_reproduced": baseline_reproduced,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
