#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from run_geometry_calibration import columns, load_geometry
from ttf.calibration import seed_for
from ttf.core import split_species
from ttf.geometry import geometry_fingerprint
from ttf.geometry_batch_v03 import (
    V03_LOCKED_CANDIDATE,
    run_geometry_calibration_v03_locked_batched,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run one locked v0.3 RGFCA-reserve intended-domain qualification cell."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--source-ledger", type=Path, required=True)
    parser.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    parser.add_argument("--min-records", type=int, default=20)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--resamples", type=int, default=1999)
    parser.add_argument("--eval-fraction", type=float, default=0.5)
    parser.add_argument("--k", type=int, default=3)
    parser.add_argument("--bandwidth", type=float, default=500.0)
    parser.add_argument("--world-batch-size", type=int, default=25)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.source_ledger.read_text())
    if source.get("schema") != "ttf_v03_rgfca_reserve_geometry_v0.1":
        raise RuntimeError("source ledger is not v0.3 RGFCA reserve geometry")
    if source.get("validation_role") != "prospectively_frozen_fresh_intended_domain_geometry_for_v03":
        raise RuntimeError("reserve geometry role drifted")
    if source.get("empirical_trait_or_colour_fields_read") is not False:
        raise RuntimeError("reserve geometry read empirical trait/colour fields")
    if int(source.get("synthetic_worlds_run_at_freeze", -1)) != 0:
        raise RuntimeError("reserve geometry opened synthetic outcomes before qualification")
    if source.get("candidate_performance_evaluated_at_freeze") is not False:
        raise RuntimeError("reserve candidate performance was evaluated at freeze")
    audit = source.get("source_audit", {})
    if audit.get("discovery_taxon_overlap") != 0:
        raise RuntimeError("reserve source overlaps discovery taxa")
    if audit.get("candidate_pixels_opened_by_audit") is not False:
        raise RuntimeError("reserve source pixels were opened")
    if audit.get("colour_fields_parsed") is not False:
        raise RuntimeError("reserve source parsed colour fields")
    if audit.get("measurement_authorized") is not False:
        raise RuntimeError("reserve measurement was authorized at geometry freeze")

    prospective = source.get("prospective_execution", {})
    if int(prospective.get("k")) != args.k:
        raise RuntimeError("reserve k differs from frozen geometry ledger")
    if abs(float(prospective.get("bandwidth_km")) - args.bandwidth) > 1e-12:
        raise RuntimeError("reserve bandwidth differs from frozen geometry ledger")
    if abs(float(prospective.get("eval_fraction")) - args.eval_fraction) > 1e-12:
        raise RuntimeError("reserve eval fraction differs from frozen geometry ledger")
    if int(prospective.get("seed")) != args.seed:
        raise RuntimeError("reserve seed differs from frozen geometry ledger")

    development = json.loads(Path("results/v03_length_control_ablation.json").read_text())
    if development.get("selected_candidate") != V03_LOCKED_CANDIDATE:
        raise RuntimeError("locked candidate differs from frozen development selection")

    external = json.loads(Path("results/v03_fresh_external_precision_v0.1.json").read_text())
    if external.get("candidate") != V03_LOCKED_CANDIDATE:
        raise RuntimeError("fresh external candidate differs from reserve candidate")
    if external.get("v03_fresh_external_pass") is not True:
        raise RuntimeError("reserve intended-domain qualification is forbidden unless fresh external passes")
    if external.get("precision_qualification", {}).get("passed") is not True:
        raise RuntimeError("fresh external Wilson qualification did not pass")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=args.min_records,
    )
    if excluded:
        raise RuntimeError(f"reserve qualification unexpectedly excluded taxa: {excluded}")
    if len(geometries) != 250:
        raise RuntimeError(f"reserve qualification requires exactly 250 taxa, found {len(geometries)}")
    if sum(len(item.coordinates) for item in geometries) != 5000:
        raise RuntimeError("reserve qualification requires exactly 5,000 records")

    labels = [item.species for item in geometries]
    split_seed = seed_for(args.seed, "fixed_geometry_split")
    train, evaluation = split_species(labels, eval_fraction=args.eval_fraction, seed=split_seed)
    if len(train) != 125 or len(evaluation) != 125:
        raise RuntimeError("reserve qualification requires fixed 125/125 split")

    cells = run_geometry_calibration_v03_locked_batched(
        geometries,
        shared_fractions=[args.shared_fraction],
        amplitudes=[args.amplitude],
        n_replicates=args.replicates,
        n_bootstrap=args.resamples,
        train_species=train,
        eval_species=evaluation,
        k=args.k,
        bandwidth=args.bandwidth,
        prior_strength=0.25,
        segment_points=5,
        noise_sd=0.8,
        transition_width=0.2,
        alpha=0.05,
        seed=args.seed,
        world_batch_size=args.world_batch_size,
    )
    if len(cells) != 1:
        raise RuntimeError("reserve v0.3 cell runner produced unexpected grid")

    counts = {item.species: len(item.coordinates) for item in geometries}
    payload = {
        "schema": "ttf_v03_reserve_geometry_cell_v0.1",
        "candidate": V03_LOCKED_CANDIDATE,
        "qualification_role": "fresh_rgfca_reserve_intended_domain_geometry",
        "geometry": {
            "input": str(args.input),
            "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
            "source_ledger": str(args.source_ledger),
            "source_ledger_sha256": hashlib.sha256(args.source_ledger.read_bytes()).hexdigest(),
            "fingerprint_sha256": geometry_fingerprint(geometries),
            "species_count": len(geometries),
            "record_count": int(sum(counts.values())),
            "records_per_species": counts,
            "coordinate_columns": args.coordinate_columns,
            "min_records": args.min_records,
        },
        "split": {
            "split_seed": split_seed,
            "train_species": list(train),
            "eval_species": list(evaluation),
        },
        "config": {
            "shared_fraction": float(args.shared_fraction),
            "amplitude": float(args.amplitude),
            "replicates": args.replicates,
            "resamples": args.resamples,
            "eval_fraction": args.eval_fraction,
            "k": args.k,
            "bandwidth": args.bandwidth,
            "prior_strength": 0.25,
            "prior_mean": 0.0,
            "segment_points": 5,
            "noise_sd": 0.8,
            "transition_width": 0.2,
            "alpha": 0.05,
            "seed": args.seed,
            "world_batch_size": args.world_batch_size,
            "inference": "heldout_species_bootstrap",
            "execution": "exact_dense_projection_batched_worlds",
            "training_response": "within_species_turnover_rank_orthogonalized_against_edge_length_rank",
            "evaluation_score": "unchanged_raw_spearman",
        },
        "cell": cells[0].to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["cell"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
