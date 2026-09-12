#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EARTH_RADIUS_KM = 6371.0088
DEFAULT_SOURCE_COMMIT = "e6ad7aa76261906b6f757ffc324ea6cd062a0917"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ecef_km(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lat = np.radians(np.asarray(latitude, dtype=float))
    lon = np.radians(np.asarray(longitude, dtype=float))
    cos_lat = np.cos(lat)
    return EARTH_RADIUS_KM * np.column_stack(
        (
            cos_lat * np.cos(lon),
            cos_lat * np.sin(lon),
            np.sin(lat),
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Freeze the intended TTF Gate-I-B geometry from the pre-inference "
            "RGFCA measurement commit. No empirical colour columns are read."
        )
    )
    parser.add_argument("--fcp-root", type=Path, required=True)
    parser.add_argument("--expected-fcp-commit", default=DEFAULT_SOURCE_COMMIT)
    parser.add_argument("--outer-index", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    args = parser.parse_args()

    root = args.fcp_root.resolve()
    source_commit = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        text=True,
    ).strip()
    if source_commit != args.expected_fcp_commit:
        raise RuntimeError(
            f"FCP source commit drift: {source_commit} != {args.expected_fcp_commit}"
        )

    measurement_path = root / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
    execution_path = root / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
    measured_path = root / "data/derived/global_monte_carlo_measured_photos_v1.csv"
    measurement = json.loads(measurement_path.read_text(encoding="utf-8"))
    execution = json.loads(execution_path.read_text(encoding="utf-8"))

    if measurement.get("status") != "complete_global_location_blind_measurement_and_join":
        raise RuntimeError("RGFCA source is not the frozen completed measurement state")
    if measurement.get("g1_g3_inference_run") is not False:
        raise RuntimeError("RGFCA source commit is not pre-inference")
    if measurement.get("external_overlay_opened") is not False:
        raise RuntimeError("RGFCA source commit has an opened external overlay")
    if execution.get("status") != "frozen_before_capacity_outcome_before_candidate_pixels_and_before_global_colour_outcome":
        raise RuntimeError("RGFCA execution contract is not the pre-outcome freeze")

    expected_measured_sha = str(measurement["lineage"]["measured_table_sha256"])
    observed_measured_sha = sha256_file(measured_path)
    if observed_measured_sha != expected_measured_sha:
        raise RuntimeError("RGFCA measured-table checksum drift")

    # Outcome firewall: these are the only measured-table columns ever loaded.
    allowed_columns = [
        "photo_id",
        "species",
        "latitude",
        "longitude",
        "global_classifiable",
    ]
    frame = pd.read_csv(measured_path, usecols=allowed_columns)
    classifiable = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    pool = frame.loc[
        classifiable,
        ["photo_id", "species", "latitude", "longitude"],
    ].copy()
    pool["species"] = pool["species"].astype(str)
    pool["photo_id"] = pd.to_numeric(pool["photo_id"], errors="raise").astype(np.int64)
    if pool["photo_id"].duplicated().any():
        raise RuntimeError("classifiable RGFCA pool has duplicate photo IDs")

    gate = execution["input_gate"]
    minimum_pool = int(gate["minimum_classifiable_photos_per_species"])
    counts = pool.groupby("species", observed=True).size()
    eligible = tuple(sorted(counts[counts >= minimum_pool].index.astype(str)))
    expected_species = int(measurement["postmeasurement_gate"]["evaluable_species"])
    if len(eligible) != expected_species:
        raise RuntimeError(
            f"eligible species drift: {len(eligible)} != {expected_species}"
        )
    pool = pool.loc[pool["species"].isin(set(eligible))].copy()
    if len(pool) != 21424:
        raise RuntimeError(f"eligible RGFCA row count drift: {len(pool)} != 21424")

    sys.path.insert(0, str(root))
    from fcp_pipeline.global_repeated_atlas import (  # noqa: PLC0415
        build_repeated_atlas_schedule,
        schedule_audit,
    )

    outer = execution["outer_schedule"]
    schedule = build_repeated_atlas_schedule(
        pool["photo_id"].to_numpy(np.int64),
        pool["species"].to_numpy(object),
        n_outer=int(outer["observed_resamples"]),
        species_per_outer=int(outer["species_per_resample"]),
        photos_per_species=int(outer["photos_per_species"]),
        minimum_pool_photos_per_species=minimum_pool,
        species_seed=int(outer["species_seed"]),
        photo_master_seed=int(outer["photo_master_seed"]),
    )
    outer_index = int(args.outer_index)
    if outer_index < 0 or outer_index >= schedule.n_outer:
        raise ValueError("outer_index outside frozen RGFCA schedule")

    selected_species = [str(x) for x in schedule.outer_species[outer_index]]
    selected_ids = schedule.outer_photo_ids[outer_index].reshape(-1).astype(np.int64)
    if len(selected_species) != int(outer["species_per_resample"]):
        raise RuntimeError("RGFCA outer species count drift")
    if len(selected_ids) != len(selected_species) * int(outer["photos_per_species"]):
        raise RuntimeError("RGFCA outer photo count drift")
    if len(np.unique(selected_ids)) != len(selected_ids):
        raise RuntimeError("RGFCA outer realization contains duplicate photo IDs")

    selected = pool.loc[pool["photo_id"].isin(set(map(int, selected_ids)))].copy()
    if len(selected) != len(selected_ids):
        raise RuntimeError("not every frozen RGFCA outer photo ID was recovered")
    expected_set = set(selected_species)
    if set(selected["species"].astype(str)) != expected_set:
        raise RuntimeError("selected RGFCA species differ from frozen outer schedule")
    selected_counts = selected.groupby("species", observed=True).size()
    photos_per_species = int(outer["photos_per_species"])
    if not (selected_counts == photos_per_species).all():
        raise RuntimeError("frozen RGFCA outer geometry is not species-equal")

    xyz = ecef_km(
        selected["latitude"].to_numpy(float),
        selected["longitude"].to_numpy(float),
    )
    selected["x_km"] = xyz[:, 0]
    selected["y_km"] = xyz[:, 1]
    selected["z_km"] = xyz[:, 2]
    selected["rgfca_outer_index"] = outer_index
    selected = selected.sort_values(["species", "photo_id"], kind="mergesort")
    columns = [
        "species",
        "photo_id",
        "latitude",
        "longitude",
        "x_km",
        "y_km",
        "z_km",
        "rgfca_outer_index",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, columns=columns, index=False)

    ledger = {
        "schema": "ttf_gate_i_b_rgfca_geometry_source_v0.1",
        "gate": "I-B intended empirical TTF sampling geometry",
        "source_repository": "zuizui0223/fcp",
        "source_commit": source_commit,
        "source_state": "post_measurement_pre_g1_g3_inference",
        "source_measurement_status": measurement["status"],
        "source_g1_g3_inference_run": measurement["g1_g3_inference_run"],
        "source_external_overlay_opened": measurement["external_overlay_opened"],
        "source_measured_table_sha256": observed_measured_sha,
        "source_measurement_manifest_sha256": sha256_file(measurement_path),
        "source_execution_contract_sha256": sha256_file(execution_path),
        "measured_columns_read": allowed_columns,
        "empirical_colour_columns_read": False,
        "eligible_pool": {
            "species": len(eligible),
            "rows": int(len(pool)),
            "minimum_classifiable_photos_per_species": minimum_pool,
        },
        "frozen_outer": {
            "outer_index": outer_index,
            "species": len(selected_species),
            "records": int(len(selected)),
            "photos_per_species": photos_per_species,
            "species_labels": sorted(selected_species),
            "schedule_audit": schedule_audit(schedule),
        },
        "coordinates": {
            "input": "latitude/longitude",
            "ttf_execution": "Earth-centred spherical x/y/z in km",
            "earth_radius_km": EARTH_RADIUS_KM,
        },
        "output": {
            "path": str(args.output),
            "sha256": sha256_file(args.output),
        },
        "claim_boundary": (
            "This file freezes sampling geometry only. No empirical flower-colour "
            "measurement is read into TTF Gate-I-B semi-synthetic calibration."
        ),
    }
    args.ledger.parent.mkdir(parents=True, exist_ok=True)
    args.ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(json.dumps(ledger["frozen_outer"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
