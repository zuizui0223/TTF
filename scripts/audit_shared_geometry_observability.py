#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from run_geometry_calibration import columns, load_geometry
from ttf import knn_edges
from ttf.geometry import geometry_fingerprint
from ttf.shared_geometry_observability import shared_transition_observability


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Audit outcome-free observability of the Gate-I shared hyperplane on "
            "one fixed sampling geometry. This is interpretive only and cannot "
            "select candidates or discard geometries."
        )
    )
    parser.add_argument("--panel", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    parser.add_argument("--min-records", type=int, default=80)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--directions", type=int, default=2048)
    parser.add_argument("--seed", type=int, default=20260912)
    parser.add_argument("--transition-width", type=float, default=0.2)
    parser.add_argument("--nontrivial-edge-contrast", type=float, default=0.25)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=args.min_records,
    )
    if excluded:
        raise RuntimeError(f"observability audit unexpectedly excluded taxa: {excluded}")
    graphs = {
        g.species: knn_edges(g.coordinates, k=args.k)
        for g in geometries
    }
    audit = shared_transition_observability(
        geometries,
        graphs,
        transition_width=args.transition_width,
        n_directions=args.directions,
        seed=args.seed,
        nontrivial_edge_contrast=args.nontrivial_edge_contrast,
    )

    payload = {
        "schema": "ttf_shared_geometry_observability_v0.1",
        "panel": args.panel,
        "role": "geometry_only_interpretation_diagnostic_not_candidate_selection",
        "candidate_selectable": False,
        "qualification_gate": False,
        "geometry": {
            "input": str(args.input),
            "fingerprint_sha256": geometry_fingerprint(geometries),
            "species_count": len(geometries),
            "record_count": int(sum(len(g.coordinates) for g in geometries)),
            "records_per_species": {
                g.species: len(g.coordinates) for g in geometries
            },
            "coordinate_columns": args.coordinate_columns,
        },
        "config": {
            "k": args.k,
            "directions": args.directions,
            "seed": args.seed,
            "transition_width": args.transition_width,
            "nontrivial_edge_contrast": args.nontrivial_edge_contrast,
            "amplitude_used": False,
            "trait_or_candidate_outcomes_used": False,
            "shared_hyperplane_offset": "pooled_median_projection_exactly_as_gate_i_simulator",
        },
        "summary": audit.summary(),
        "species_mean_contrast": dict(sorted(audit.species_mean_contrast.items())),
        "species_median_contrast": dict(sorted(audit.species_median_contrast.items())),
        "species_nontrivial_fraction": dict(sorted(audit.species_nontrivial_fraction.items())),
        "interpretation_boundary": (
            "This measures how strongly the simulator's shared transition is visible on "
            "the fixed sampling graph before noise and before any TTF estimator is fit. "
            "It cannot qualify a method, tune a candidate, or justify dropping a hard panel."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"panel": args.panel, **payload["summary"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
