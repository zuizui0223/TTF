#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import qualify_calibration, seed_for
from ttf.core import SpeciesSample, spearman_rho
from ttf.geometry_batch_v04 import run_geometry_calibration_v04_locked_batched
from ttf.nulls import fixed_graphs
from ttf.precision import qualify_calibration_precision
from ttf.private_geometry_control import (
    pooled_affine_frame,
    random_private_transition_propensity,
)


CELLS = ((0.0, 0.5), (0.0, 1.0), (0.0, 2.0), (0.0, 3.0), (1.0, 2.0))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit the explicit v0.4 propensity seed lock on development data."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--v03-result", type=Path, required=True)
    parser.add_argument("--diagnostic", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--reserve-ledger", type=Path, required=True)
    parser.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    parser.add_argument("--min-records", type=int, default=80)
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--bandwidth", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--resamples", type=int, default=1999)
    parser.add_argument("--propensity-directions", type=int, default=2048)
    parser.add_argument("--formal-propensity-seed", type=int, default=20260908)
    parser.add_argument("--transition-width", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    v03 = json.loads(args.v03_result.read_text())
    diagnostic = json.loads(args.diagnostic.read_text())
    selection = json.loads(args.selection.read_text())
    reserve = json.loads(args.reserve_ledger.read_text())
    if v03.get("v03_fresh_external_pass") is not False:
        raise RuntimeError("seed audit requires failed v0.3 development panel")
    if selection.get("selected_candidate") != "eval_geometry_partial":
        raise RuntimeError("unexpected v0.4 candidate selection")
    if selection["locked_definition"]["private_propensity_seed"] != args.formal_propensity_seed:
        raise RuntimeError("formal seed disagrees with frozen candidate definition")
    if reserve.get("synthetic_worlds_run_at_freeze") != 0:
        raise RuntimeError("RGFCA reserve was opened")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=args.min_records,
    )
    if excluded:
        raise RuntimeError(f"unexpected geometry exclusions: {excluded}")
    train = tuple(v03["split"]["train_species"])
    evaluation = tuple(v03["split"]["eval_species"])
    gmap = {g.species: g for g in geometries}

    cells = run_geometry_calibration_v04_locked_batched(
        geometries,
        shared_fractions=(0.0, 1.0),
        amplitudes=(0.5, 1.0, 2.0, 3.0),
        n_replicates=args.replicates,
        n_bootstrap=args.resamples,
        train_species=train,
        eval_species=evaluation,
        k=args.k,
        bandwidth=args.bandwidth,
        prior_strength=0.25,
        segment_points=5,
        noise_sd=0.8,
        transition_width=args.transition_width,
        alpha=0.05,
        seed=args.seed,
        world_batch_size=25,
        propensity_directions=args.propensity_directions,
        propensity_seed=args.formal_propensity_seed,
        propensity_transition_width=args.transition_width,
    )
    cells = tuple(
        cell for cell in cells
        if (float(cell.shared_fraction), float(cell.amplitude)) in set(CELLS)
    )
    if len(cells) != len(CELLS):
        raise RuntimeError("seed audit mandatory cell set drifted")
    point = qualify_calibration(
        cells,
        moderate_amplitude=2.0,
        type1_ceiling=0.10,
        power_floor=0.80,
    )
    precision = qualify_calibration_precision(
        cells,
        moderate_amplitude=2.0,
        type1_upper_ceiling=0.10,
        power_lower_floor=0.80,
    )

    diag_cells = {
        (float(row["shared_fraction"]), float(row["amplitude"])):
            row["variants"]["eval_geometry_partial"]
        for row in diagnostic["cells"]
    }
    shifts = []
    for cell in cells:
        key = (float(cell.shared_fraction), float(cell.amplitude))
        old = diag_cells[key]
        shifts.append({
            "shared_fraction": key[0],
            "amplitude": key[1],
            "diagnostic_derived_seed_rejection": float(old["rejection_rate"]),
            "formal_seed_rejection": float(cell.rejection_rate),
            "absolute_rejection_shift": abs(float(cell.rejection_rate) - float(old["rejection_rate"])),
            "diagnostic_derived_seed_mean_statistic": float(old["mean_statistic"]),
            "formal_seed_mean_statistic": float(cell.mean_statistic),
            "absolute_mean_statistic_shift": abs(float(cell.mean_statistic) - float(old["mean_statistic"])),
        })

    templates = [
        SpeciesSample(
            name,
            gmap[name].coordinates,
            np.arange(len(gmap[name].coordinates), dtype=float),
        )
        for name in train + evaluation
    ]
    graphs = fixed_graphs(templates, k=args.k)
    center, scale, active = pooled_affine_frame(geometries)
    derived_seed = seed_for(args.seed, "v04_private_propensity_directions")
    propensity_rho = []
    for name in evaluation:
        formal = random_private_transition_propensity(
            gmap[name].coordinates,
            graphs[name],
            center=center,
            scale=scale,
            active=active,
            transition_width=args.transition_width,
            n_directions=args.propensity_directions,
            seed=args.formal_propensity_seed,
        )
        derived = random_private_transition_propensity(
            gmap[name].coordinates,
            graphs[name],
            center=center,
            scale=scale,
            active=active,
            transition_width=args.transition_width,
            n_directions=args.propensity_directions,
            seed=derived_seed,
        )
        propensity_rho.append(spearman_rho(formal, derived))

    payload = {
        "schema": "ttf_v04_propensity_seed_lock_audit_v0.1",
        "status": "development_only_seed_definition_audit",
        "claim_ready": False,
        "formal_candidate_seed": args.formal_propensity_seed,
        "diagnostic_derived_seed": derived_seed,
        "directions": args.propensity_directions,
        "propensity_rank_agreement": {
            "mean_spearman": float(np.mean(propensity_rho)),
            "min_spearman": float(np.min(propensity_rho)),
            "max_spearman": float(np.max(propensity_rho)),
        },
        "cell_shifts": shifts,
        "formal_seed_cells": [cell.to_dict() for cell in cells],
        "point_qualification": point.to_dict(),
        "precision_qualification": precision.to_dict(),
        "candidate_lock_valid_on_development_rule": bool(point.passed),
        "reserve_opened": False,
        "interpretation_boundary": (
            "This audit resolves a seed-serialization mismatch on already failed development data. "
            "It does not qualify v0.4. A confirmatory panel remains required."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "candidate_lock_valid": payload["candidate_lock_valid_on_development_rule"],
        "propensity_rank_agreement": payload["propensity_rank_agreement"],
        "point": payload["point_qualification"],
        "precision": payload["precision_qualification"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
