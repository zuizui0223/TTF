#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import pandas as pd

EARTH_RADIUS_KM = 6371.0088


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob(root: Path, rel: str) -> str:
    return subprocess.check_output(["git", "-C", str(root), "rev-parse", f"HEAD:{rel}"], text=True).strip()


def rank_key(taxon_id: int, seed: int) -> tuple[bytes, int]:
    return hashlib.sha256(f"{int(seed)}\x1f{int(taxon_id)}".encode("utf-8")).digest(), int(taxon_id)


def taxon_digest(values: list[int]) -> str:
    text = "\n".join(str(int(x)) for x in sorted(set(values))) + "\n"
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def xyz(latitude: float, longitude: float) -> tuple[float, float, float]:
    lat = math.radians(float(latitude))
    lon = math.radians(float(longitude))
    c = math.cos(lat)
    return (
        EARTH_RADIUS_KM * c * math.cos(lon),
        EARTH_RADIUS_KM * c * math.sin(lon),
        EARTH_RADIUS_KM * math.sin(lat),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--fcp-root", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--output-geometry", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_v11_fresh_density_scaled_graph_rule_v0.1":
        raise RuntimeError("wrong v0.11 rule")
    if rule.get("status") != "frozen_before_v11_fresh_plant_geometry_acquisition_and_before_any_v11_synthetic_worlds":
        raise RuntimeError("v0.11 rule is not pre-geometry")
    source = rule["fresh_geometry_source"]
    fcp_root = args.fcp_root.resolve()
    commit = subprocess.check_output(["git", "-C", str(fcp_root), "rev-parse", "HEAD"], text=True).strip()
    if commit != source["commit"]:
        raise RuntimeError("FCP source commit drift")
    for path_key, blob_key in (
        ("capacity_manifest", "capacity_manifest_git_blob_sha"),
        ("capacity_selected_species", "capacity_selected_species_git_blob_sha"),
        ("candidate_species_budget", "candidate_species_budget_git_blob_sha"),
        ("original_candidate_authorization", "original_candidate_authorization_git_blob_sha"),
        ("candidate_acquisition_contract", "candidate_acquisition_contract_git_blob_sha"),
    ):
        if git_blob(fcp_root, source[path_key]) != source[blob_key]:
            raise RuntimeError(f"FCP blob drift: {path_key}")

    capacity = pd.read_csv(fcp_root / source["capacity_selected_species"])
    ranked = sorted(capacity["inat_taxon_id"].astype(int).tolist(), key=lambda x: rank_key(x, int(source["species_hash_seed"])))
    previous = tuple(ranked[:1000])
    fresh_taxa = tuple(ranked[1000:1250])
    expected_fresh = set(fresh_taxa)
    if len(expected_fresh) != 250 or set(previous) & expected_fresh:
        raise RuntimeError("fresh rank slice drift")

    manifests = []
    audits = []
    photos = []
    for path in sorted(args.input_dir.glob("v11_shard_*.json")):
        d = json.loads(path.read_text())
        if d.get("schema") == "ttf_v11_fresh_plant_acquisition_shard_v0.1":
            manifests.append(d)
    for path in sorted(args.input_dir.glob("v11_shard_*_species.csv")):
        audits.append(pd.read_csv(path))
    for path in sorted(args.input_dir.glob("v11_shard_*_photos.csv")):
        photos.append(pd.read_csv(path))
    if len(manifests) != 8 or len(audits) != 8 or len(photos) != 8:
        raise RuntimeError(f"expected 8 complete v0.11 shards, got manifests={len(manifests)} audits={len(audits)} photos={len(photos)}")
    if sorted(int(x["shard_index"]) for x in manifests) != list(range(8)):
        raise RuntimeError("v0.11 shard index set drift")
    if any(x.get("candidate_image_pixels_opened") is not False or x.get("flower_colour_used") is not False for x in manifests):
        raise RuntimeError("v0.11 acquisition opened forbidden outcomes")
    if any(int(x.get("synthetic_worlds_run", -1)) != 0 or x.get("candidate_performance_evaluated") is not False for x in manifests):
        raise RuntimeError("v0.11 acquisition evaluated synthetic performance")

    audit = pd.concat(audits, ignore_index=True).sort_values("v11_row_index", kind="mergesort").reset_index(drop=True)
    photo = pd.concat(photos, ignore_index=True) if photos else pd.DataFrame()
    if len(audit) != 250 or audit["inat_taxon_id"].astype(int).nunique() != 250 or audit["v11_row_index"].astype(int).nunique() != 250:
        raise RuntimeError("v0.11 species audit denominator drift")
    audit_taxa = set(audit["inat_taxon_id"].astype(int).tolist())
    if audit_taxa != expected_fresh:
        raise RuntimeError("v0.11 acquired species differ from rank 1001-1250")
    if sorted(audit["hash_rank_one_based"].astype(int).tolist()) != list(range(1001, 1251)):
        raise RuntimeError("v0.11 rank positions drift")

    request_attempts = sum(int(x["request_attempts"]) for x in manifests)
    request_errors = sum(int(x["request_errors"]) for x in manifests)
    request_error_fraction = float(request_errors / request_attempts) if request_attempts else 1.0
    full_target = audit["full_target"].astype(str).str.casefold().eq("true")
    full_target_species = int(full_target.sum())
    all_full = full_target_species == 250
    request_gate = request_error_fraction <= 0.05

    photo_ok = False
    photo_reason = ""
    if all_full:
        required_photo_cols = {"species", "inat_taxon_id", "observation_id", "photo_id", "latitude", "longitude"}
        if not required_photo_cols.issubset(photo.columns):
            photo_reason = "missing_required_photo_columns"
        elif len(photo) != 25000:
            photo_reason = f"wrong_total_photo_rows:{len(photo)}"
        elif photo["photo_id"].astype(str).duplicated().any():
            photo_reason = "duplicate_photo_ids"
        elif photo["observation_id"].astype(str).duplicated().any():
            photo_reason = "duplicate_observation_ids"
        elif set(photo["inat_taxon_id"].astype(int).tolist()) != expected_fresh:
            photo_reason = "photo_taxon_set_drift"
        else:
            counts = photo.groupby(photo["inat_taxon_id"].astype(int), observed=True).size()
            if len(counts) != 250 or not (counts == 100).all():
                photo_reason = "not_exactly_100_rows_per_species"
            else:
                photo_ok = True
    else:
        photo_reason = "at_least_one_species_failed_full_target"

    freeze_pass = bool(all_full and request_gate and photo_ok)
    geometry_sha = None
    geometry_records = 0
    if freeze_pass:
        geometry_rows = []
        for row in photo.itertuples(index=False):
            x, y, z = xyz(float(row.latitude), float(row.longitude))
            geometry_rows.append({
                "species": str(row.species),
                "inat_taxon_id": int(row.inat_taxon_id),
                "x_km": x,
                "y_km": y,
                "z_km": z,
                "source_observation_id": str(row.observation_id),
                "source_photo_id": str(row.photo_id),
                "latitude": float(row.latitude),
                "longitude": float(row.longitude),
            })
        geometry = pd.DataFrame(geometry_rows).sort_values(["inat_taxon_id", "source_photo_id"], kind="mergesort").reset_index(drop=True)
        args.output_geometry.parent.mkdir(parents=True, exist_ok=True)
        geometry.to_csv(args.output_geometry, index=False, lineterminator="\n")
        geometry_sha = sha256_file(args.output_geometry)
        geometry_records = int(len(geometry))
    elif args.output_geometry.exists():
        args.output_geometry.unlink()

    failing = audit.loc[~full_target, ["hash_rank_one_based", "species", "inat_taxon_id", "after_observer_cap", "candidate_page_errors"]].to_dict(orient="records")
    ledger = {
        "schema": "ttf_v11_fresh_plant_geometry_source_v0.1",
        "status": "fresh_geometry_frozen_without_outcomes" if freeze_pass else "not_evaluable_fresh_geometry_due_frozen_metadata_acquisition_gate",
        "fresh_geometry_freeze_pass": freeze_pass,
        "rule": str(args.rule),
        "source_repository": source["repository"],
        "source_commit": commit,
        "source_blobs": {
            source["capacity_manifest"]: git_blob(fcp_root, source["capacity_manifest"]),
            source["capacity_selected_species"]: git_blob(fcp_root, source["capacity_selected_species"]),
            source["candidate_species_budget"]: git_blob(fcp_root, source["candidate_species_budget"]),
            source["original_candidate_authorization"]: git_blob(fcp_root, source["original_candidate_authorization"]),
            source["candidate_acquisition_contract"]: git_blob(fcp_root, source["candidate_acquisition_contract"]),
        },
        "selection": {
            "capacity_species": 4730,
            "raw_photo_target": 100,
            "species_hash_seed": int(source["species_hash_seed"]),
            "rank_range_one_based": [1001, 1250],
            "previous_rank_range_one_based": [1, 1000],
            "previous_1000_taxon_digest": taxon_digest(list(previous)),
            "v11_250_taxon_digest": taxon_digest(list(fresh_taxa)),
            "taxon_overlap_with_previous_1000": sorted(set(previous) & expected_fresh),
            "selected_species": 250,
        },
        "acquisition": {
            "shards": 8,
            "request_attempts": request_attempts,
            "request_errors": request_errors,
            "request_error_fraction": request_error_fraction,
            "request_error_fraction_ceiling": 0.05,
            "request_gate_pass": request_gate,
            "full_target_species": full_target_species,
            "required_full_target_species": 250,
            "all_species_full_target_pass": all_full,
            "photo_structure_pass": photo_ok,
            "photo_structure_reason": photo_reason,
            "failing_species": failing,
        },
        "output": {
            "geometry_path": str(args.output_geometry) if freeze_pass else None,
            "geometry_sha256": geometry_sha,
            "species": 250 if freeze_pass else 0,
            "records": geometry_records,
            "records_per_species": 100 if freeze_pass else None,
        },
        "prospective_execution": rule["prospective_execution"],
        "locked_graph_architecture": rule["locked_graph_architecture"],
        "candidate_image_pixels_opened": False,
        "empirical_trait_or_colour_fields_read": False,
        "synthetic_worlds_run_at_freeze": 0,
        "candidate_performance_evaluated_at_freeze": False,
        "species_replacement_performed": False,
        "additional_pages_after_result_performed": False,
        "claim_boundary": "This ledger freezes only a fresh plant sampling geometry acquisition outcome. It contains no image pixels, flower-colour values, synthetic TTF performance or empirical TTF conclusion. A failed freeze may not be rescued by species substitution, additional pages, target relaxation or a favourable rerun.",
    }
    args.output_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.output_ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"fresh_geometry_freeze_pass": freeze_pass, "full_target_species": full_target_species, "request_error_fraction": request_error_fraction, "records": geometry_records}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
