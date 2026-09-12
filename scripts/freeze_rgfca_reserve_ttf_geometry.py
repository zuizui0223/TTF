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


def priority(seed: int, *parts: object) -> str:
    return hashlib.sha256(
        "|".join(map(str, (int(seed),) + parts)).encode("utf-8")
    ).hexdigest()


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
    by_fold = {name.casefold(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.casefold() in by_fold:
            return by_fold[candidate.casefold()]
    raise RuntimeError(f"cannot resolve {role} column from {fieldnames}")


def main() -> int:
    p = argparse.ArgumentParser(
        description=(
            "Freeze a TTF sampling geometry from the outcome-unopened RGFCA reserve "
            "metadata cohort. No colour/pixel fields are read."
        )
    )
    p.add_argument("--fcp-root", type=Path, required=True)
    p.add_argument("--source-commit", required=True)
    p.add_argument("--species-count", type=int, default=250)
    p.add_argument("--records-per-species", type=int, default=20)
    p.add_argument("--seed", type=int, default=20260908)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--ledger", type=Path, required=True)
    args = p.parse_args()

    candidate_path = args.fcp_root / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
    reserve_audit_path = args.fcp_root / "docs/supporting/rgfca_reserve_metadata_audit_v1.json"
    geometry_audit_path = args.fcp_root / "data/frozen/rgfca_reserve_geometry_audit_v1.csv"
    reserve_contract_path = args.fcp_root / "docs/supporting/rgfca_reserve_replication_contract_v1.json"

    audit = json.loads(reserve_audit_path.read_text())
    contract = json.loads(reserve_contract_path.read_text())
    if audit.get("status") != "metadata_only_reserve_admission_passed":
        raise RuntimeError("reserve metadata audit is not admitted")
    if audit.get("rows") != 50000 or audit.get("species") != 500:
        raise RuntimeError("reserve metadata audit dimensions drifted")
    if audit.get("discovery_taxon_overlap") != 0:
        raise RuntimeError("reserve cohort overlaps discovery taxa")
    if audit.get("candidate_pixels_opened_by_audit") is not False:
        raise RuntimeError("reserve audit opened candidate pixels")
    if audit.get("colour_fields_parsed") is not False:
        raise RuntimeError("reserve audit parsed colour fields")
    if audit.get("measurement_authorized") is not False:
        raise RuntimeError("reserve measurement was already authorized at geometry freeze source")
    if contract.get("status") != "frozen_before_reserve_pixels":
        raise RuntimeError("reserve replication contract is not pre-pixel")
    expected_candidate_sha = audit["source_sha256"]["data/frozen/global_monte_carlo_candidate_photos_v1.csv"]
    actual_candidate_sha = sha256_file(candidate_path)
    if actual_candidate_sha != expected_candidate_sha:
        raise RuntimeError("candidate metadata hash differs from reserve audit")

    with geometry_audit_path.open(newline="", encoding="utf-8") as handle:
        reserve_rows = list(csv.DictReader(handle))
    reserve_species = {str(row["species"]).strip() for row in reserve_rows}
    if len(reserve_species) != 500:
        raise RuntimeError(f"reserve geometry audit has {len(reserve_species)} species, expected 500")

    if not 20 <= args.species_count <= 500:
        raise ValueError("species-count must lie in [20, 500]")
    if not 8 <= args.records_per_species <= 100:
        raise ValueError("records-per-species must lie in [8, 100]")

    selected_species = tuple(
        sorted(
            reserve_species,
            key=lambda name: (priority(args.seed, "species", name), name),
        )[: args.species_count]
    )
    selected_set = set(selected_species)

    per_species: dict[str, list[dict[str, object]]] = {name: [] for name in selected_species}
    with candidate_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or ())
        species_col = resolve_field(fieldnames, ("species", "scientific_name"), "species")
        photo_col = resolve_field(fieldnames, ("photo_id", "inat_photo_id"), "photo id")
        lat_col = resolve_field(fieldnames, ("latitude", "decimalLatitude", "lat"), "latitude")
        lon_col = resolve_field(fieldnames, ("longitude", "decimalLongitude", "lon", "lng"), "longitude")
        for row_number, row in enumerate(reader, start=2):
            species = str(row[species_col]).strip()
            if species not in selected_set:
                continue
            photo_id = str(row[photo_col]).strip()
            if not photo_id:
                raise RuntimeError(f"empty photo id at source row {row_number}")
            try:
                latitude = float(row[lat_col])
                longitude = float(row[lon_col])
            except (TypeError, ValueError) as exc:
                raise RuntimeError(f"invalid coordinate at source row {row_number}") from exc
            if not (-90 <= latitude <= 90 and -180 <= longitude <= 180):
                raise RuntimeError(f"out-of-range coordinate at source row {row_number}")
            per_species[species].append(
                {
                    "photo_id": photo_id,
                    "latitude": latitude,
                    "longitude": longitude,
                }
            )

    output_rows: list[dict[str, object]] = []
    source_counts: dict[str, int] = {}
    for species in selected_species:
        rows = per_species[species]
        source_counts[species] = len(rows)
        if len(rows) != 100:
            raise RuntimeError(
                f"reserve source has {len(rows)} rows for {species}, expected exactly 100"
            )
        rows.sort(
            key=lambda row: (
                priority(args.seed, "photo", species, row["photo_id"]),
                str(row["photo_id"]),
            )
        )
        chosen = rows[: args.records_per_species]
        if len({str(row["photo_id"]) for row in chosen}) != len(chosen):
            raise RuntimeError(f"duplicate selected photo ids for {species}")
        for row in chosen:
            x, y, z = xyz_km(float(row["latitude"]), float(row["longitude"]))
            output_rows.append(
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

    output_rows.sort(key=lambda row: (str(row["species"]), str(row["source_photo_id"])))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "species", "x_km", "y_km", "z_km", "source_photo_id", "latitude", "longitude"
            ),
        )
        writer.writeheader()
        writer.writerows(output_rows)

    payload = {
        "schema": "ttf_v03_rgfca_reserve_geometry_v0.1",
        "validation_role": "prospectively_frozen_fresh_intended_domain_geometry_for_v03",
        "source_repository": "zuizui0223/fcp",
        "source_commit": args.source_commit,
        "source_files": {
            str(candidate_path.relative_to(args.fcp_root)): sha256_file(candidate_path),
            str(reserve_audit_path.relative_to(args.fcp_root)): sha256_file(reserve_audit_path),
            str(geometry_audit_path.relative_to(args.fcp_root)): sha256_file(geometry_audit_path),
            str(reserve_contract_path.relative_to(args.fcp_root)): sha256_file(reserve_contract_path),
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
            "seed": args.seed,
            "species_rule": "lowest SHA256(seed|species|species_name), ties by species name",
            "photo_rule": "within selected species, lowest SHA256(seed|photo|species|photo_id), ties by photo_id",
            "species_count": args.species_count,
            "records_per_species": args.records_per_species,
            "selected_species": list(selected_species),
        },
        "output": {
            "path": str(args.output),
            "sha256": sha256_file(args.output),
            "records": len(output_rows),
            "species": len(selected_species),
        },
        "empirical_trait_or_colour_fields_read": False,
        "synthetic_worlds_run_at_freeze": 0,
        "candidate_performance_evaluated_at_freeze": False,
        "prospective_execution": {
            "k": 3,
            "bandwidth_km": 500.0,
            "eval_fraction": 0.5,
            "seed": 20260908,
            "train_species": 125,
            "eval_species": 125,
        },
        "claim_boundary": (
            "This freezes a species- and photo-disjoint RGFCA reserve sampling geometry only. "
            "No reserve pixels, colour outcomes, synthetic TTF worlds or v0.3 performance are opened."
        ),
    }
    args.ledger.parent.mkdir(parents=True, exist_ok=True)
    args.ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"species": len(selected_species), "records": len(output_rows)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
