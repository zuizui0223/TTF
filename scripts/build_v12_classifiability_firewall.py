#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def mid(photo_id: object, salt: str) -> str:
    payload = f"{salt}\x1fphoto\x1f{photo_id}".encode("utf-8")
    return "FCPR-" + hashlib.sha256(payload).hexdigest().upper()[:24]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--measurement-contract", type=Path, required=True)
    ap.add_argument("--execution-contract", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    measurement = json.loads(args.measurement_contract.read_text())
    execution = json.loads(args.execution_contract.read_text())
    if rule.get("schema") != "ttf_v12_classifiability_actual_geometry_rule_v0.1":
        raise RuntimeError("wrong v0.12 rule")
    if rule.get("status") != "frozen_after_v11_synthetic_pass_before_any_v11_empirical_colour_values_are_opened":
        raise RuntimeError("v0.12 rule is not frozen at the required boundary")
    if measurement.get("protocol") != "random-photo-first-measurement-v1" or measurement.get("status") != "frozen_before_any_fresh_candidate_image_pixel":
        raise RuntimeError("unexpected frozen FCP measurement contract")
    if execution.get("protocol") != "random-photo-first-measurement-execution-v1" or execution.get("status") != "prospectively_frozen_before_any_fresh_candidate_image_pixel":
        raise RuntimeError("unexpected frozen FCP execution contract")
    if tuple(execution["blinding"]["worker_allowed_fields"]) != ("measurement_id", "image_filename", "photo_license"):
        raise RuntimeError("FCP blind worker interface drifted")
    salt = str(execution["blinding"]["salt"])

    paths = sorted(args.input_dir.glob("v11_shard_*_photos.csv"))
    if len(paths) != 8:
        raise RuntimeError(f"expected 8 frozen v11 photo shards, found {len(paths)}")
    frames = [pd.read_csv(path, dtype={"photo_id": str, "observation_id": str}) for path in paths]
    photos = pd.concat(frames, ignore_index=True)
    required = {"species", "inat_taxon_id", "observation_id", "photo_id", "photo_url_large", "photo_license", "latitude", "longitude"}
    if not required.issubset(photos.columns):
        raise RuntimeError(f"v11 photo artifacts missing fields: {sorted(required-set(photos.columns))}")
    if len(photos) != 25000 or photos["photo_id"].astype(str).nunique() != 25000:
        raise RuntimeError("v11 artifact denominator is not exactly 25,000 unique photos")
    counts = photos.groupby("species", observed=True).size()
    if len(counts) != 250 or not (counts == 100).all():
        raise RuntimeError("v11 artifact is not exactly 250 species x 100 photos")

    geometry = pd.read_csv(args.geometry, dtype={"source_photo_id": str, "source_observation_id": str})
    greq = {"species", "inat_taxon_id", "source_photo_id", "source_observation_id", "latitude", "longitude", "x_km", "y_km", "z_km"}
    if not greq.issubset(geometry.columns):
        raise RuntimeError("v11 frozen geometry schema drifted")
    if len(geometry) != 25000 or geometry["source_photo_id"].astype(str).nunique() != 25000:
        raise RuntimeError("v11 geometry denominator drifted")
    if set(photos["photo_id"].astype(str)) != set(geometry["source_photo_id"].astype(str)):
        raise RuntimeError("artifact photo IDs do not exactly match frozen geometry")

    photos = photos.copy()
    photos["measurement_id"] = [mid(x, salt) for x in photos["photo_id"].astype(str)]
    if photos["measurement_id"].nunique() != len(photos):
        raise RuntimeError("measurement IDs are not unique")
    photos["image_filename"] = photos["measurement_id"] + ".jpg"
    worker = photos[["measurement_id", "image_filename", "photo_license"]].copy()
    acquisition = photos[["measurement_id", "image_filename", "photo_url_large", "photo_license"]].copy()

    meta = geometry.merge(
        photos[["measurement_id", "photo_id", "species", "inat_taxon_id"]],
        left_on=["source_photo_id", "species", "inat_taxon_id"],
        right_on=["photo_id", "species", "inat_taxon_id"],
        how="inner",
        validate="one_to_one",
    ).drop(columns=["photo_id"])
    if len(meta) != 25000 or meta["measurement_id"].nunique() != 25000:
        raise RuntimeError("metadata join key lost rows")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    worker_path = args.output_dir / "worker_manifest.csv"
    acquisition_path = args.output_dir / "acquisition_key.csv"
    metadata_path = args.output_dir / "metadata_join_key.csv"
    worker.to_csv(worker_path, index=False, lineterminator="\n")
    acquisition.to_csv(acquisition_path, index=False, lineterminator="\n")
    meta.to_csv(metadata_path, index=False, lineterminator="\n")
    manifest = {
        "schema": "ttf_v12_classifiability_firewall_v0.1",
        "status": "actual_geometry_measurement_firewall_frozen_before_pixels",
        "candidate_rows": 25000,
        "species": 250,
        "records_per_species": 100,
        "measurement_ids": 25000,
        "worker_fields": list(worker.columns),
        "worker_contains_species": False,
        "worker_contains_coordinates": False,
        "worker_contains_source_url": False,
        "metadata_join_opened_to_measurement": False,
        "colour_values_opened": False,
        "pixels_opened_at_firewall_freeze": False,
        "replacement_after_failure": False,
        "v11_geometry_sha256": sha256_file(args.geometry),
        "worker_sha256": sha256_file(worker_path),
        "acquisition_sha256": sha256_file(acquisition_path),
        "metadata_join_sha256": sha256_file(metadata_path),
        "measurement_contract_protocol": measurement["protocol"],
        "execution_contract_protocol": execution["protocol"],
        "measurement_evaluable_terminal_status": "classified_four_state_morph",
        "claim_boundary": "This firewall opens no colour value. It only prepares the exact 25,000 frozen v0.11 photos for the independently frozen location-blind FCP measurement pipeline and preserves the coordinate/species join key outside measurement workers.",
    }
    (args.output_dir / "firewall_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"firewall": "pass", "rows": 25000, "species": 250}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
