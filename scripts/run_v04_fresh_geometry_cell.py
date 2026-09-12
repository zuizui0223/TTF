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
from ttf.geometry_batch_v04 import (
    V04_LOCKED_CANDIDATE,
    V04_PROPENSITY_DIRECTIONS,
    V04_PROPENSITY_SEED,
    V04_PROPENSITY_TRANSITION_WIDTH,
    run_geometry_calibration_v04_locked_batched,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run one locked v0.4 fresh external geometry qualification cell."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--source-ledger", type=Path, required=True)
    parser.add_argument("--coordinate-columns", type=columns, default=["x_km", "y_km", "z_km"])
    parser.add_argument("--min-records", type=int, default=80)
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--resamples", type=int, default=1999)
    parser.add_argument("--eval-fraction", type=float, default=0.5)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--bandwidth", type=float, required=True)
    parser.add_argument("--world-batch-size", type=int, default=25)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.source_ledger.read_text())
    if source.get("validation_role") != "prospectively_frozen_second_external_geometry_for_v04_or_successor":
        raise RuntimeError("source ledger is not the prospective second external holdout")
    if source.get("confirmatory_holdout") is not True:
        raise RuntimeError("v0.4 fresh geometry is not marked confirmatory")
    if source.get("candidate_performance_evaluated_at_freeze") is not False:
        raise RuntimeError("v0.4 fresh outcome firewall was not preserved")
    if int(source.get("synthetic_worlds_run_at_freeze", -1)) != 0:
        raise RuntimeError("v0.4 fresh geometry had synthetic outcomes before candidate lock")
    if source.get("prior_gate_i_taxa_overlap") != []:
        raise RuntimeError("v0.4 fresh geometry overlaps prior Gate-I taxa")
    prospective = source.get("prospective_execution", {})
    if int(prospective.get("k")) != int(args.k):
        raise RuntimeError("v0.4 fresh k differs from pre-outcome ledger")
    if abs(float(prospective.get("bandwidth_km")) - float(args.bandwidth)) > 1e-12:
        raise RuntimeError("v0.4 fresh bandwidth differs from pre-outcome ledger")
    if abs(float(prospective.get("eval_fraction")) - float(args.eval_fraction)) > 1e-12:
        raise RuntimeError("v0.4 fresh eval fraction differs from pre-outcome ledger")
    if int(prospective.get("seed")) != int(args.seed):
        raise RuntimeError("v0.4 fresh world/split seed differs from pre-outcome ledger")

    selection = json.loads(Path("results/v04_candidate_selection_v0.1.json").read_text())
    if selection.get("selected_candidate") != V04_LOCKED_CANDIDATE:
        raise RuntimeError("code candidate differs from frozen v0.4 selection")
    locked = selection.get("locked_definition", {})
    if int(locked.get("private_propensity_seed")) != V04_PROPENSITY_SEED:
        raise RuntimeError("v0.4 propensity seed drift")
    if int(locked.get("private_propensity_directions")) != V04_PROPENSITY_DIRECTIONS:
        raise RuntimeError("v0.4 propensity direction-count drift")
    if abs(float(locked.get("private_propensity_transition_width")) - V04_PROPENSITY_TRANSITION_WIDTH) > 1e-12:
        raise RuntimeError("v0.4 propensity width drift")

    seed_audit = json.loads(Path("results/v04_propensity_seed_lock_audit_v0.1.json").read_text())
    if seed_audit.get("candidate_lock_valid_on_development_rule") is not True:
        raise RuntimeError("formal v0.4 seed lock did not survive development audit")

    geometries, excluded = load_geometry(
        args.input,
        species_column="species",
        coordinate_columns=args.coordinate_columns,
        block_column=None,
        min_records=args.min_records,
    )
    if excluded:
        raise RuntimeError(f"v0.4 qualification unexpectedly excluded taxa: {excluded}")
    if len(geometries) != 48:
        raise RuntimeError(f"v0.4 qualification requires exactly 48 taxa, found {len(geometries)}")

    labels = [item.species for item in geometries]
    split_seed = seed_for(args.seed, "fixed_geometry_split")
    train, evaluation = split_species(labels, eval_fraction=args.eval_fraction, seed=split_seed)
    if len(train) != 24 or len(evaluation) != 24:
        raise RuntimeError("v0.4 qualification requires fixed 24/24 species split")

    cells = run_geometry_calibration_v04_locked_batched(
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
        propensity_directions=V04_PROPENSITY_DIRECTIONS,
        propensity_seed=V04_PROPENSITY_SEED,
        propensity_transition_width=V04_PROPENSITY_TRANSITION_WIDTH,
    )
    if len(cells) != 1:
        raise RuntimeError("v0.4 fresh cell runner produced unexpected grid")

    record_counts = {item.species: len(item.coordinates) for item in geometries}
    payload = {
        "schema": "ttf_v04_fresh_geometry_cell_v0.1",
        "candidate": V04_LOCKED_CANDIDATE,
        "qualification_role": "second_fresh_external_confirmatory_geometry",
        "geometry": {
            "input": str(args.input),
            "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
            "source_ledger": str(args.source_ledger),
            "source_ledger_sha256": hashlib.sha256(args.source_ledger.read_bytes()).hexdigest(),
            "fingerprint_sha256": geometry_fingerprint(geometries),
            "species_count": len(geometries),
            "record_count": int(sum(record_counts.values())),
            "records_per_species": record_counts,
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
            "training_response": "v03_length_rank_orthogonalized_turnover",
            "evaluation_score": "partial_spearman_control_edge_length_and_random_private_geometry_propensity",
            "propensity_directions": V04_PROPENSITY_DIRECTIONS,
            "propensity_seed": V04_PROPENSITY_SEED,
            "propensity_transition_width": V04_PROPENSITY_TRANSITION_WIDTH,
        },
        "cell": cells[0].to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["cell"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
