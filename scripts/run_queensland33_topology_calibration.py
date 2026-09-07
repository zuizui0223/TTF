#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from ttf.calibration import CalibrationCell, seed_for
from ttf.core import SpeciesSample, split_species
from ttf.inference import heldout_species_bootstrap_test
from ttf.semisynthetic import geometry_bandwidth, normalize_geometry

BASINS = ("GOL", "LOG", "BRI", "PIN", "MCY", "NOO", "TCB", "MRY")


def load_layout(path: Path, layout_index: int) -> dict[str, np.ndarray]:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "ttf_queensland33_feasible_layouts_v0.1":
        raise RuntimeError("unexpected feasible-layout schema")
    layout = payload["layouts"][int(layout_index)]
    geometry: dict[str, np.ndarray] = {}
    for row in layout["rows"]:
        if not row["eligible"]:
            continue
        supported = [basin for basin in BASINS if int(row["basin_counts"][basin]) >= 3]
        if len(supported) < 4:
            raise RuntimeError(f"eligibility drift for {row['code']}")
        x = np.asarray([BASINS.index(basin) for basin in supported], dtype=float)
        geometry[str(row["code"])] = np.column_stack([x, np.zeros_like(x)])
    if len(geometry) != 21:
        raise RuntimeError(f"expected 21 eligible species, got {len(geometry)}")
    return geometry


def crossing_species(geometry: dict[str, np.ndarray], cut: float) -> tuple[str, ...]:
    return tuple(sorted(
        species for species, coords in geometry.items()
        if float(coords[:, 0].min()) < cut < float(coords[:, 0].max())
    ))


def simulate_ordered_world(
    geometry: dict[str, np.ndarray],
    *,
    shared_fraction: float,
    amplitude: float,
    noise_sd: float,
    transition_width: float,
    min_shared_crossing_species: int,
    k: int,
    seed: int,
) -> tuple[tuple[SpeciesSample, ...], float, float | None, tuple[str, ...]]:
    if not 0.0 <= shared_fraction <= 1.0 or amplitude < 0 or noise_sd < 0 or transition_width <= 0:
        raise ValueError("invalid synthetic-world parameters")
    normalized = normalize_geometry(geometry)
    coords = normalized.coordinates
    names = tuple(sorted(coords))
    rng = np.random.default_rng(int(seed))

    # Published basin order is the only topology supplied.  Boundaries live
    # between adjacent basin positions; Mary-Brisbane is never privileged.
    raw_cuts = np.arange(0.5, len(BASINS) - 0.5, 1.0)
    normalized_cuts = (raw_cuts - float(normalized.center[0])) / float(normalized.scale)
    cut_crossings = [crossing_species({s: np.column_stack([
        (geometry[s][:, 0] - normalized.center[0]) / normalized.scale,
        np.zeros(len(geometry[s]))
    ]) for s in geometry}, float(cut)) for cut in normalized_cuts]

    shared_cut = None
    shared_crossing: tuple[str, ...] = ()
    n_shared = int(round(shared_fraction * len(names)))
    shared_names: set[str] = set()
    if n_shared > 0:
        eligible_cut_ids = [i for i, crossing in enumerate(cut_crossings) if len(crossing) >= min_shared_crossing_species]
        if not eligible_cut_ids:
            raise RuntimeError("no topology cut satisfies the prespecified shared-crossing requirement")
        cut_id = int(rng.choice(eligible_cut_ids))
        shared_cut = float(normalized_cuts[cut_id])
        shared_crossing = cut_crossings[cut_id]
        crossing_order = list(shared_crossing)
        rng.shuffle(crossing_order)
        chosen = crossing_order[:min(n_shared, len(crossing_order))]
        if len(chosen) < n_shared:
            remaining = [name for name in names if name not in chosen]
            rng.shuffle(remaining)
            chosen.extend(remaining[:n_shared-len(chosen)])
        shared_names = set(chosen)

    samples = []
    for name in names:
        x = coords[name][:, 0]
        if name in shared_names:
            assert shared_cut is not None
            cut = shared_cut
        else:
            # Private worlds retain strong within-species transitions but choose
            # their break independently among global inter-basin cuts that the
            # species actually spans.  This is deliberately adversarial under
            # uneven sampling support.
            valid = [float(cut) for cut in normalized_cuts if float(x.min()) < cut < float(x.max())]
            if not valid:
                raise RuntimeError(f"no private cut available for {name}")
            cut = float(rng.choice(valid))
        latent = float(amplitude) * np.tanh((x - cut) / float(transition_width))
        trait = latent + rng.normal(0.0, float(noise_sd), size=len(x))
        samples.append(SpeciesSample(species=name, coordinates=coords[name], trait=trait))

    return tuple(samples), geometry_bandwidth(coords, k=k), shared_cut, shared_crossing


def main() -> int:
    parser = argparse.ArgumentParser(description="Calibrate TTF on basin-layout uncertainty for the Queensland 33-species panel.")
    parser.add_argument("--layouts", type=Path, required=True)
    parser.add_argument("--layout-index", type=int, required=True)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--bootstrap", type=int, default=999)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--noise-sd", type=float, default=0.8)
    parser.add_argument("--transition-width", type=float, default=0.20)
    parser.add_argument("--min-shared-crossing-species", type=int, default=14)
    parser.add_argument("--split-seed", type=int, default=20260907)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geometry = load_layout(args.layouts, args.layout_index)
    labels = tuple(sorted(geometry))
    train, evaluation = split_species(labels, eval_fraction=11.0 / len(labels), seed=args.split_seed)
    if len(train) != 10 or len(evaluation) != 11:
        raise RuntimeError("Queensland topology benchmark requires frozen 10/11 split")

    p_values = []
    statistics = []
    null_means = []
    bandwidths = []
    crossing_counts = []
    for replicate in range(args.replicates):
        samples, bandwidth, _shared_cut, crossing = simulate_ordered_world(
            geometry,
            shared_fraction=float(args.shared_fraction),
            amplitude=float(args.amplitude),
            noise_sd=float(args.noise_sd),
            transition_width=float(args.transition_width),
            min_shared_crossing_species=int(args.min_shared_crossing_species),
            k=int(args.k),
            seed=seed_for(args.seed, "qld33-topology-world", args.layout_index, args.shared_fraction, args.amplitude, replicate),
        )
        result = heldout_species_bootstrap_test(
            samples,
            train_species=train,
            eval_species=evaluation,
            k=int(args.k),
            bandwidth=float(bandwidth),
            n_bootstrap=int(args.bootstrap),
            seed=seed_for(args.seed, "qld33-topology-bootstrap", args.layout_index, args.shared_fraction, args.amplitude, replicate),
        )
        p_values.append(result.p_value)
        statistics.append(result.observed.statistic)
        null_means.append(result.null_mean)
        bandwidths.append(bandwidth)
        crossing_counts.append(len(crossing))

    p = np.asarray(p_values, dtype=float)
    cell = CalibrationCell(
        shared_fraction=float(args.shared_fraction),
        amplitude=float(args.amplitude),
        n_replicates=int(args.replicates),
        alpha=float(args.alpha),
        rejection_rate=float(np.mean(p <= args.alpha)),
        mean_statistic=float(np.mean(statistics)),
        mean_null_statistic=float(np.mean(null_means)),
        median_p_value=float(np.median(p)),
    )
    payload = {
        "schema": "ttf_queensland33_topology_cell_v0.1",
        "layout_index": int(args.layout_index),
        "cell": asdict(cell),
        "design": {
            "n_eligible_species": len(labels),
            "train_species": list(train),
            "eval_species": list(evaluation),
            "split": "10 train / 11 evaluation species",
            "split_seed": int(args.split_seed),
            "geometry": "ordered_basin_topology_v0.1: x=published south-to-north basin column index, y=0",
            "named_Mary_Brisbane_boundary_used": False,
            "genetic_outcomes_used": False,
            "layout_selected_by_outcome": False,
            "k": int(args.k),
            "bandwidth_rule": "pooled median species-local kNN edge length after normalization",
            "bandwidth": float(np.mean(bandwidths)),
            "noise_sd": float(args.noise_sd),
            "transition_width": float(args.transition_width),
            "min_shared_crossing_species": int(args.min_shared_crossing_species),
            "bootstrap_resamples": int(args.bootstrap),
            "master_seed": int(args.seed),
        },
        "diagnostics": {
            "mean_crossing_species": float(np.mean(crossing_counts)) if crossing_counts else 0.0,
            "min_crossing_species": int(min(crossing_counts)) if crossing_counts else 0,
            "max_crossing_species": int(max(crossing_counts)) if crossing_counts else 0,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"layout": args.layout_index, "cell": payload["cell"], "diagnostics": payload["diagnostics"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
