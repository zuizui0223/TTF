#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

EARTH_RADIUS_KM = 6371.0088


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def xyz_km(latitude: float, longitude: float) -> tuple[float, float, float]:
    lat = math.radians(latitude)
    lon = math.radians(longitude)
    cos_lat = math.cos(lat)
    return (
        EARTH_RADIUS_KM * cos_lat * math.cos(lon),
        EARTH_RADIUS_KM * cos_lat * math.sin(lon),
        EARTH_RADIUS_KM * math.sin(lat),
    )


def resolve_field(fieldnames: list[str], candidates: tuple[str, ...], role: str) -> str:
    folded = {name.casefold(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.casefold() in folded:
            return folded[candidate.casefold()]
    raise RuntimeError(f"cannot resolve {role} column from {fieldnames}")


def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Freeze the v0.9 RGFCA reserve-B max-density geometry from the exact "
            "species complement of the failed v0.8 Gate-D reserve. No colour or pixels are read."
        )
    )
    p.add_argument("--fcp-root", type=Path, required=True)
    p.add_argument("--source-commit", required=True)
    p.add_argument("--rule", type=Path, required=True)
    p.add_argument("--failed-reserve-ledger", type=Path, required=True)
    p.add_argument("--failed-gate-d-result", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--ledger", type=Path, required=True)
    args = p.parse_args()

    rule = json.loads(args.rule.read_text())
    failed_reserve = json.loads(args.failed_reserve_ledger.read_text())
    failed_gate = json.loads(args.failed_gate_d_result.read_text())
    if rule.get("status") != "frozen_before_v09_reserve_b_geometry_is_materialized":
        raise RuntimeError("v0.9 rule is not prospectively frozen")
    if failed_gate.get("gate_d_pass") is not False:
        raise RuntimeError("v0.9 successor is only defined after the frozen Gate-D failure")
    if failed_gate.get("empirical_colour_outcome_opened") is not False:
        raise RuntimeError("failed Gate-D opened empirical colour")
    if failed_reserve.get("schema") != "ttf_v03_rgfca_reserve_geometry_v0.1":
        raise RuntimeError("unexpected failed reserve ledger")
    if failed_reserve.get("empirical_trait_or_colour_fields_read") is not False:
        raise RuntimeError("failed reserve already read empirical outcomes")

    candidate_path = args.fcp_root / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
    reserve_audit_path = args.fcp_root / "docs/supporting/rgfca_reserve_metadata_audit_v1.json"
    geometry_audit_path = args.fcp_root / "data/frozen/rgfca_reserve_geometry_audit_v1.csv"
    reserve_contract_path = args.fcp_root / "docs/supporting/rgfca_reserve_replication_contract_v1.json"

    audit = json.loads(reserve_audit_path.read_text())
    contract = json.loads(reserve_contract_path.read_text())
    if audit.get("status") != "metadata_only_reserve_admission_passed":
        raise RuntimeError("reserve metadata audit is not admitted")
    if audit.get("rows") != 50000 or audit.get("species") != 500 or audit.get("photos_per_species") != 100:
        raise RuntimeError("reserve metadata dimensions drifted")
    if audit.get("discovery_taxon_overlap") != 0:
        raise RuntimeError("reserve cohort overlaps discovery taxa")
    if audit.get("candidate_pixels_opened_by_audit") is not False or audit.get("colour_fields_parsed") is not False or audit.get("measurement_authorized") is not False:
        raise RuntimeError("reserve outcome firewall already opened")
    if contract.get("status") != "frozen_before_reserve_pixels":
        raise RuntimeError("reserve contract is not pre-pixel")
    expected_sha = audit["source_sha256"]["data/frozen/global_monte_carlo_candidate_photos_v1.csv"]
    if sha256_file(candidate_path) != expected_sha:
        raise RuntimeError("candidate metadata hash drifted")

    with geometry_audit_path.open(newline="", encoding="utf-8") as handle:
        all_species = {str(row["species"]).strip() for row in csv.DictReader(handle)}
    if len(all_species) != 500:
        raise RuntimeError(f"expected 500 reserve species, found {len(all_species)}")

    used_species = set(map(str, failed_reserve.get("selection", {}).get("selected_species", [])))
    if len(used_species) != 250 or not used_species <= all_species:
        raise RuntimeError("failed reserve species set drifted")
    complement = tuple(sorted(all_species - used_species))
    if len(complement) != 250:
        raise RuntimeError(f"expected 250-species complement, found {len(complement)}")
    complement_set = set(complement)

    per_species: dict[str, list[dict[str, object]]] = {name: [] for name in complement}
    with candidate_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or ())
        species_col = resolve_field(fields, ("species", "scientific_name"), "species")
        photo_col = resolve_field(fields, ("photo_id", "inat_photo_id"), "photo id")
        lat_col = resolve_field(fields, ("latitude", "decimalLatitude", "lat"), "latitude")
        lon_col = resolve_field(fields, ("longitude", "decimalLongitude", "lon", "lng"), "longitude")
        for row_number, row in enumerate(reader, start=2):
            species = str(row[species_col]).strip()
            if species not in complement_set:
                continue
            photo_id = str(row[photo_col]).strip()
            if not photo_id:
                raise RuntimeError(f"empty photo id at row {row_number}")
            latitude = float(row[lat_col])
            longitude = float(row[lon_col])
            if not (-90 <= latitude <= 90 and -180 <= longitude <= 180):
                raise RuntimeError(f"invalid coordinate at row {row_number}")
            per_species[species].append(
                {"photo_id": photo_id, "latitude": latitude, "longitude": longitude}
            )

    rows_out: list[dict[str, object]] = []
    for species in complement:
        rows = per_species[species]
        if len(rows) != 100:
            raise RuntimeError(f"{species} has {len(rows)} rows, expected 100")
        if len({str(row["photo_id"]) for row in rows}) != 100:
            raise RuntimeError(f"duplicate photo ids for {species}")
        for row in rows:
            x, y, z = xyz_km(float(row["latitude"]), float(row["longitude"]))
            rows_out.append(
                {
                    "species": species,
                    "x_km": x,
                    "y_km": y,
                    "z_km": z,
                    "source_photo_id": row["photo_id"],
                    "latitude": row["latitude"],
                    "longitude": row["longitude"],
                }
            )

    rows_out.sort(key=lambda row: (str(row["species"]), str(row["source_photo_id"])))
    if len(rows_out) != 25000:
        raise RuntimeError(f"expected 25,000 records, found {len(rows_out)}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("species", "x_km", "y_km", "z_km", "source_photo_id", "latitude", "longitude"),
        )
        writer.writeheader()
        writer.writerows(rows_out)

    failed_photo_ids = set()
    with Path(failed_reserve["output"]["path"]).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            failed_photo_ids.add(str(row["source_photo_id"]))
    new_photo_ids = {str(row["source_photo_id"]) for row in rows_out}
    if failed_photo_ids & new_photo_ids:
        raise RuntimeError("v0.9 reserve-B overlaps failed Gate-D photo ids")

    payload = {
        "schema": "ttf_v09_rgfca_reserve_b_max_density_geometry_v0.1",
        "validation_role": "outcome_unopened_species_disjoint_max_density_successor_design_after_v08_gate_d_failure",
        "source_repository": "zuizui0223/fcp",
        "source_commit": args.source_commit,
        "source_files": {
            str(candidate_path.relative_to(args.fcp_root)): sha256_file(candidate_path),
            str(reserve_audit_path.relative_to(args.fcp_root)): sha256_file(reserve_audit_path),
            str(geometry_audit_path.relative_to(args.fcp_root)): sha256_file(geometry_audit_path),
            str(reserve_contract_path.relative_to(args.fcp_root)): sha256_file(reserve_contract_path),
        },
        "predecessor": {
            "failed_reserve_ledger": str(args.failed_reserve_ledger),
            "failed_reserve_ledger_sha256": sha256_file(args.failed_reserve_ledger),
            "failed_gate_d_result": str(args.failed_gate_d_result),
            "failed_gate_d_result_sha256": sha256_file(args.failed_gate_d_result),
            "species_overlap": sorted(used_species & complement_set),
            "photo_overlap": [],
        },
        "source_audit": {
            "rows": audit["rows"],
            "species": audit["species"],
            "photos_per_species": audit["photos_per_species"],
            "discovery_taxon_overlap": audit["discovery_taxon_overlap"],
            "candidate_pixels_opened_by_audit": audit["candidate_pixels_opened_by_audit"],
            "colour_fields_parsed": audit["colour_fields_parsed"],
            "measurement_authorized": audit["measurement_authorized"],
        },
        "selection": {
            "species_rule": "exact 250-species complement of failed v0.8 Gate-D reserve within the frozen 500-species reserve cohort",
            "photo_rule": "all 100 frozen metadata rows per complement species; no density search",
            "species_count": 250,
            "records_per_species": 100,
            "selected_species": list(complement),
        },
        "output": {
            "path": str(args.output),
            "sha256": sha256_file(args.output),
            "records": len(rows_out),
            "species": len(complement),
        },
        "empirical_trait_or_colour_fields_read": False,
        "synthetic_worlds_run_at_freeze": 0,
        "candidate_performance_evaluated_at_freeze": False,
        "prospective_execution": {
            "k": 3,
            "bandwidth_km": 500.0,
            "eval_fraction": 0.5,
            "seed": 20260920,
            "train_species": 125,
            "eval_species": 125,
        },
        "claim_boundary": (
            "This is a new species- and photo-disjoint max-density deployment-design layer created after the frozen v0.8 Gate-D failure. "
            "It does not rescue or retune the failed 250x20 panel. No pixels, empirical colour outcomes, or synthetic performance are opened at freeze."
        ),
    }
    args.ledger.parent.mkdir(parents=True, exist_ok=True)
    args.ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"species": 250, "records": 25000, "species_overlap": 0, "photo_overlap": 0}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
