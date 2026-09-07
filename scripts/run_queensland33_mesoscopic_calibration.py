#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from run_queensland33_topology_calibration import (
    BASINS,
    load_layout as load_legacy_layout,
    simulate_ordered_world,
)
from ttf.calibration import CalibrationCell, seed_for
from ttf.core import SpeciesSample, split_species
from ttf.mesoscopic import mesoscopic_bootstrap_test
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.semisynthetic import normalize_geometry


def load_mesoscopic_layout(path: Path, layout_index: int) -> tuple[dict[str, np.ndarray], dict]:
    """Read either the generator artifact schema or the compact frozen snapshot."""

    payload = json.loads(path.read_text())
    if payload.get("schema") != "ttf_queensland33_topology_layouts_frozen_v0.2":
        return load_legacy_layout(path, layout_index)
    if tuple(payload.get("basin_column_order", ())) != BASINS:
        raise RuntimeError("frozen basin-column order drift")
    matches = [layout for layout in payload["layouts"] if int(layout["i"]) == int(layout_index)]
    if len(matches) != 1:
        raise RuntimeError(f"expected one frozen layout {layout_index}, got {len(matches)}")
    layout = matches[0]
    geometry: dict[str, np.ndarray] = {}
    for code, counts in layout["c"].items():
        if len(counts) != len(BASINS):
            raise RuntimeError(f"basin-count length drift for {code}")
        supported_ids = [i for i, count in enumerate(counts) if int(count) >= 3]
        if len(supported_ids) < 4:
            raise RuntimeError(f"eligibility drift for {code}")
        x = np.asarray(supported_ids, dtype=float)
        geometry[str(code)] = np.column_stack([x, np.zeros_like(x)])
    expected = set(payload.get("eligibility", {}).get("eligible_codes", ()))
    if expected and set(geometry) != expected:
        raise RuntimeError("frozen eligible-species set drift")
    if len(geometry) != 21:
        raise RuntimeError(f"expected 21 eligible species, got {len(geometry)}")
    return geometry, {
        "layout_schema": payload["schema"],
        "layout_mode": payload.get("layout_mode", "ordered_nonzero_counts"),
        "solver": {},
        "source_schema": payload.get("source_schema"),
        "lineage": payload.get("lineage", {}),
    }


def candidate_cuts(geometry: dict[str, np.ndarray]) -> np.ndarray:
    normalized = normalize_geometry(geometry)
    raw = np.arange(0.5, len(BASINS) - 0.5, 1.0)
    return (raw - float(normalized.center[0])) / float(normalized.scale)


def geometry_graphs(geometry: dict[str, np.ndarray], *, k: int) -> dict[str, np.ndarray]:
    normalized = normalize_geometry(geometry)
    samples = [
        SpeciesSample(
            species=species,
            coordinates=normalized.coordinates[species],
            trait=np.zeros(len(normalized.coordinates[species]), dtype=float),
        )
        for species in sorted(normalized.coordinates)
    ]
    return fixed_graphs(samples, k=int(k))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Calibrate the mesoscopic TTF estimator on the frozen Queensland "
            "8-basin topology envelope without opening genetic outcomes."
        )
    )
    parser.add_argument("--layouts", type=Path, required=True)
    parser.add_argument("--layout-index", type=int, required=True)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--bootstrap", type=int, default=999)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--k", type=int, default=2)
    parser.add_argument("--noise-sd", type=float, default=0.8)
    parser.add_argument("--transition-width", type=float, default=0.20)
    parser.add_argument("--min-shared-crossing-species", type=int, default=14)
    parser.add_argument("--prior-strength", type=float, default=1.0)
    parser.add_argument("--split-seed", type=int, default=20260907)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geometry, layout_meta = load_mesoscopic_layout(args.layouts, args.layout_index)
    labels = tuple(sorted(geometry))
    train, evaluation = split_species(
        labels,
        eval_fraction=11.0 / len(labels),
        seed=int(args.split_seed),
    )
    if len(train) != 10 or len(evaluation) != 11:
        raise RuntimeError("Queensland mesoscopic benchmark requires frozen 10/11 split")

    cuts = candidate_cuts(geometry)
    if len(cuts) != 7:
        raise RuntimeError("Queensland mesoscopic benchmark requires seven inter-basin cuts")
    graphs = geometry_graphs(geometry, k=int(args.k))

    p_values: list[float] = []
    statistics: list[float] = []
    field_entropy: list[float] = []
    max_field_probability: list[float] = []
    finite_eval: list[int] = []
    crossing_counts: list[int] = []

    for replicate in range(int(args.replicates)):
        world_seed = seed_for(
            int(args.seed),
            "qld33-mesoscopic-world",
            int(args.layout_index),
            float(args.shared_fraction),
            float(args.amplitude),
            replicate,
        )
        samples, _bandwidth, _shared_cut, crossing = simulate_ordered_world(
            geometry,
            shared_fraction=float(args.shared_fraction),
            amplitude=float(args.amplitude),
            noise_sd=float(args.noise_sd),
            transition_width=float(args.transition_width),
            min_shared_crossing_species=int(args.min_shared_crossing_species),
            k=int(args.k),
            seed=world_seed,
        )
        sample_map = {sample.species: sample for sample in samples}
        edge_map = edges_on_fixed_graphs(samples, graphs)
        result = mesoscopic_bootstrap_test(
            [edge_map[s] for s in train],
            [edge_map[s] for s in evaluation],
            cuts,
            prior_strength=float(args.prior_strength),
            n_bootstrap=int(args.bootstrap),
            seed=seed_for(
                int(args.seed),
                "qld33-mesoscopic-bootstrap",
                int(args.layout_index),
                float(args.shared_fraction),
                float(args.amplitude),
                replicate,
            ),
        )
        # Fail loudly if simulation and graph maps ever drift in labels.
        if set(sample_map) != set(edge_map):
            raise RuntimeError("sample/edge label drift")
        p_values.append(float(result.p_value))
        statistics.append(float(result.observed.statistic))
        probability = np.asarray(result.observed.field.probability, dtype=float)
        field_entropy.append(float(-np.sum(probability * np.log(probability))))
        max_field_probability.append(float(np.max(probability)))
        finite_eval.append(int(result.observed.n_eval_species))
        crossing_counts.append(len(crossing))

    p = np.asarray(p_values, dtype=float)
    cell = CalibrationCell(
        shared_fraction=float(args.shared_fraction),
        amplitude=float(args.amplitude),
        n_replicates=int(args.replicates),
        alpha=float(args.alpha),
        rejection_rate=float(np.mean(p <= float(args.alpha))),
        mean_statistic=float(np.mean(statistics)),
        mean_null_statistic=0.0,
        median_p_value=float(np.median(p)),
    )
    payload = {
        "schema": "ttf_queensland33_mesoscopic_cell_v0.3",
        "layout_index": int(args.layout_index),
        "cell": asdict(cell),
        "design": {
            "estimator": "mesoscopic boundary-location log-score gain over geometry-conditioned private-cut baseline",
            "n_candidate_cuts": int(len(cuts)),
            "n_eligible_species": int(len(labels)),
            "train_species": list(train),
            "eval_species": list(evaluation),
            "split": "10 train / 11 evaluation species",
            "split_seed": int(args.split_seed),
            "layout_schema": layout_meta["layout_schema"],
            "layout_mode": layout_meta["layout_mode"],
            "exact_frozen_layout_envelope": True,
            "layout_selected_by_outcome": False,
            "named_Mary_Brisbane_boundary_used": False,
            "genetic_outcomes_used": False,
            "k": int(args.k),
            "noise_sd": float(args.noise_sd),
            "transition_width": float(args.transition_width),
            "min_shared_crossing_species": int(args.min_shared_crossing_species),
            "prior_strength": float(args.prior_strength),
            "bootstrap_resamples": int(args.bootstrap),
            "master_seed": int(args.seed),
        },
        "diagnostics": {
            "mean_field_entropy": float(np.mean(field_entropy)),
            "mean_max_field_probability": float(np.mean(max_field_probability)),
            "mean_finite_eval_species": float(np.mean(finite_eval)),
            "min_finite_eval_species": int(min(finite_eval)),
            "mean_shared_crossing_species": float(np.mean(crossing_counts)) if crossing_counts else 0.0,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "layout": int(args.layout_index),
                "shared_fraction": float(args.shared_fraction),
                "amplitude": float(args.amplitude),
                "rejection_rate": cell.rejection_rate,
                "mean_statistic": cell.mean_statistic,
                "mean_max_field_probability": payload["diagnostics"]["mean_max_field_probability"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
