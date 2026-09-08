#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from fcp_pipeline.global_candidate_acquisition import deterministic_candidate_pages, stable_candidate_query
from fcp_pipeline.global_capacity_handoff import taxon_digest
from fcp_pipeline.random_photo_h9_pool import geographic_maximin, observer_cap, parse_h9_observation
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

PARSER_COLUMNS = [
    "species", "inat_taxon_id", "observation_id", "photo_id", "photo_url_large",
    "photo_license", "attribution", "latitude", "longitude", "positional_accuracy_m",
    "observed_on", "observer_id", "observer",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob(root: Path, rel: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", f"HEAD:{rel}"], text=True
    ).strip()


def rank_key(taxon_id: int, seed: int) -> tuple[bytes, int]:
    payload = f"{int(seed)}\x1f{int(taxon_id)}".encode("utf-8")
    return hashlib.sha256(payload).digest(), int(taxon_id)


def maximum_span_km(frame: pd.DataFrame) -> float:
    if len(frame) < 2:
        return 0.0
    lat = np.deg2rad(frame["latitude"].to_numpy(dtype=float))
    lon = np.deg2rad(frame["longitude"].to_numpy(dtype=float))
    c = np.cos(lat)
    xyz = np.column_stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return float(np.max(np.arccos(dot)) * 6371.0088)


def load_exclusions(fcp_root: Path, parent: dict[str, object]) -> tuple[set[int], set[int], dict[str, str]]:
    observation_ids: set[int] = set()
    photo_ids: set[int] = set()
    hashes: dict[str, str] = {}
    paths = parent.get("prior_experiment_exclusion_sources") or []
    if not isinstance(paths, list) or not paths:
        raise RuntimeError("FCP parent contract lacks prior experiment exclusions")
    for rel in paths:
        path = fcp_root / str(rel)
        if not path.exists():
            raise RuntimeError(f"missing exclusion source: {rel}")
        frame = pd.read_csv(path, usecols=lambda c: c in {"observation_id", "photo_id"})
        if not {"observation_id", "photo_id"}.issubset(frame.columns):
            raise RuntimeError(f"bad exclusion source: {rel}")
        observation_ids.update(frame["observation_id"].dropna().astype(int).tolist())
        photo_ids.update(frame["photo_id"].dropna().astype(int).tolist())
        hashes[str(rel)] = sha256_file(path)
    return observation_ids, photo_ids, hashes


def results(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    raw = payload.get("results") or []
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
        raise RuntimeError("iNaturalist results is not a sequence")
    return [row for row in raw if isinstance(row, Mapping)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fcp-root", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, required=True)
    ap.add_argument("--shard-count", type=int, default=8)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    fcp_root = args.fcp_root.resolve()
    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_v11_fresh_density_scaled_graph_rule_v0.1":
        raise RuntimeError("wrong v0.11 rule")
    if rule.get("status") != "frozen_before_v11_fresh_plant_geometry_acquisition_and_before_any_v11_synthetic_worlds":
        raise RuntimeError("v0.11 rule is not pre-geometry")

    source = rule["fresh_geometry_source"]
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
            raise RuntimeError(f"FCP source blob drift: {path_key}")

    capacity_manifest = json.loads((fcp_root / source["capacity_manifest"]).read_text())
    if capacity_manifest.get("status") != "complete_metadata_only_capacity_scan_target_selected_after_transport_recovery_v3":
        raise RuntimeError("capacity source is not the frozen v3 recovery")
    if capacity_manifest.get("candidate_image_pixels_opened") is not False or capacity_manifest.get("flower_colour_used") is not False:
        raise RuntimeError("capacity source opened forbidden outcomes")
    if int(capacity_manifest.get("selected_species") or 0) != 4730 or int(capacity_manifest.get("selected_raw_photo_target") or 0) != 100:
        raise RuntimeError("capacity dimensions drifted")

    selected_path = fcp_root / source["capacity_selected_species"]
    if sha256_file(selected_path) != source["capacity_selected_species_sha256"]:
        raise RuntimeError("capacity selected-species SHA256 drift")
    frame = pd.read_csv(selected_path)
    required = {"species", "inat_taxon_id", "selected_raw_photo_target"}
    if not required.issubset(frame.columns):
        raise RuntimeError("capacity species table lacks required columns")
    if len(frame) != 4730 or frame["inat_taxon_id"].nunique() != 4730:
        raise RuntimeError("capacity species denominator drift")
    if not (frame["selected_raw_photo_target"].astype(int) == 100).all():
        raise RuntimeError("capacity raw-photo target drift")

    seed = int(source["species_hash_seed"])
    ranked = sorted(frame["inat_taxon_id"].astype(int).tolist(), key=lambda x: rank_key(x, seed))
    previous = tuple(ranked[:1000])
    fresh_taxa = tuple(ranked[1000:1250])
    if len(previous) != 1000 or len(fresh_taxa) != 250 or set(previous) & set(fresh_taxa):
        raise RuntimeError("v0.11 rank slice drift")
    original_auth = json.loads((fcp_root / source["original_candidate_authorization"]).read_text())
    if taxon_digest(pd.Series(previous, dtype="int64")) != original_auth["candidate_taxon_id_sha256"]:
        raise RuntimeError("rank 1-1000 no longer reproduces original candidate authorization")

    fresh = frame.loc[frame["inat_taxon_id"].astype(int).isin(set(fresh_taxa))].copy()
    fresh["inat_taxon_id"] = fresh["inat_taxon_id"].astype(int)
    fresh["species"] = fresh["species"].astype(str)
    fresh = fresh.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)
    if len(fresh) != 250:
        raise RuntimeError("v0.11 selected species count drift")
    fresh["v11_row_index"] = np.arange(250, dtype=int)

    shard_count = int(args.shard_count)
    shard_index = int(args.shard_index)
    if shard_count != 8 or shard_index < 0 or shard_index >= shard_count:
        raise RuntimeError("v0.11 shard contract drift")
    shard = fresh.loc[(fresh["v11_row_index"] % shard_count) == shard_index].copy().reset_index(drop=True)
    if shard.empty:
        raise RuntimeError("v0.11 shard unexpectedly empty")

    acquisition_contract = json.loads((fcp_root / source["candidate_acquisition_contract"]).read_text())
    if acquisition_contract.get("status") != "frozen_after_v2_discovery_success_before_capacity_outcome_and_before_candidate_pixels":
        raise RuntimeError("FCP acquisition contract drift")
    q = acquisition_contract["query"]
    selection = acquisition_contract["selection"]
    m = rule["metadata_acquisition"]
    checks = {
        "candidate_page_seed": int(q["candidate_page_seed"]),
        "candidate_pages_per_species": int(q["candidate_pages_per_species"]),
        "per_page": int(q["per_page"]),
        "maximum_api_page": int(q["maximum_api_page"]),
        "maximum_positional_accuracy_m": int(q["maximum_positional_accuracy_m"]),
        "flowering_term_id": int(q["flowering_term_id"]),
        "flowering_term_value_id": int(q["flowering_term_value_id"]),
        "request_retries": int(q["request_retries"]),
        "observer_cap_per_species": int(selection["observer_cap_per_species"]),
    }
    for key, value in checks.items():
        if value != int(m[key]):
            raise RuntimeError(f"v0.11 acquisition rule drift: {key}")
    if sorted(map(str.casefold, q["allowed_photo_licenses"])) != sorted(map(str.casefold, m["allowed_photo_licenses"])):
        raise RuntimeError("photo license rule drift")
    if int(q["request_retries"]) != 0 or bool(q["extra_page_after_underperformance"]):
        raise RuntimeError("FCP query would permit result-dependent rescue")

    parent = json.loads((fcp_root / "docs/supporting/global_monte_carlo_capacity_scan_contract_v1.json").read_text())
    excluded_obs, excluded_photo, exclusion_hashes = load_exclusions(fcp_root, parent)
    allowed = frozenset(str(x).casefold() for x in q["allowed_photo_licenses"])
    client = InaturalistObservationClient(
        request_interval_seconds=float(q["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=0,
        user_agent=f"ttf-v11-fresh-plant-geometry-s{shard_index}/1.0 (github.com/zuizui0223/TTF)",
    )

    request_attempts = 0
    request_errors = 0
    audits: list[dict[str, object]] = []
    photo_rows: list[dict[str, object]] = []

    def request(taxon_id: int, page: int) -> tuple[Mapping[str, object] | None, str]:
        nonlocal request_attempts, request_errors
        request_attempts += 1
        params = stable_candidate_query(
            taxon_id,
            page=page,
            per_page=int(q["per_page"]),
            maximum_positional_accuracy_m=int(q["maximum_positional_accuracy_m"]),
            flowering_term_id=int(q["flowering_term_id"]),
            flowering_term_value_id=int(q["flowering_term_value_id"]),
            allowed_photo_licenses=tuple(q["allowed_photo_licenses"]),
        )
        try:
            return client.observations(params), ""
        except Exception as exc:
            request_errors += 1
            return None, f"{type(exc).__name__}:{str(exc)[:180]}"

    for row in shard.itertuples(index=False):
        taxon_id = int(row.inat_taxon_id)
        probe, probe_error = request(taxon_id, 1)
        total_results = 0
        pages: tuple[int, ...] = ()
        payload_by_page: dict[int, Mapping[str, object]] = {}
        page_errors: dict[int, str] = {}
        if probe is not None:
            try:
                total_results = int(probe.get("total_results") or len(results(probe)))
                pages = deterministic_candidate_pages(
                    taxon_id,
                    total_results,
                    per_page=int(q["per_page"]),
                    maximum_api_page=int(q["maximum_api_page"]),
                    pages_per_species=int(q["candidate_pages_per_species"]),
                    seed=int(q["candidate_page_seed"]),
                )
                if 1 in pages:
                    payload_by_page[1] = probe
            except Exception as exc:
                probe_error = f"{type(exc).__name__}:{str(exc)[:180]}"
        if probe_error:
            page_errors[1] = probe_error
        if probe is not None:
            for page in pages:
                if page == 1:
                    continue
                payload, error = request(taxon_id, page)
                if payload is not None:
                    payload_by_page[int(page)] = payload
                if error:
                    page_errors[int(page)] = error

        raw_rows: list[Mapping[str, object]] = []
        for page in pages:
            payload = payload_by_page.get(int(page))
            if payload is not None:
                raw_rows.extend(results(payload))

        parsed_rows: list[dict[str, object]] = []
        locally_eligible = 0
        prior_excluded = 0
        wrong_taxon = 0
        for observation in raw_rows:
            taxon = observation.get("taxon") or {}
            if isinstance(taxon, Mapping) and int(taxon.get("id") or -1) != taxon_id:
                wrong_taxon += 1
                continue
            parsed = parse_h9_observation(
                observation,
                expected_taxon_id=taxon_id,
                maximum_positional_accuracy_m=float(q["maximum_positional_accuracy_m"]),
                allowed_photo_licenses=allowed,
            )
            if parsed is None:
                continue
            locally_eligible += 1
            if int(parsed["observation_id"]) in excluded_obs or int(parsed["photo_id"]) in excluded_photo:
                prior_excluded += 1
                continue
            parsed_rows.append(parsed)

        metadata = pd.DataFrame(parsed_rows, columns=PARSER_COLUMNS)
        if len(metadata):
            metadata = metadata.drop_duplicates("observation_id", keep="first")
            metadata = metadata.drop_duplicates("photo_id", keep="first").reset_index(drop=True)
        deduplicated = int(len(metadata))
        capped = observer_cap(metadata, 2) if len(metadata) else metadata
        capped_n = int(len(capped))
        full_target = capped_n >= 100
        chosen = geographic_maximin(capped, 100) if full_target else capped.iloc[0:0].copy()
        if full_target and len(chosen) != 100:
            raise RuntimeError("geographic maximin did not return exactly 100 rows")
        if full_target:
            chosen = chosen.copy()
            chosen["v11_row_index"] = int(row.v11_row_index)
            chosen["shard_index"] = shard_index
            photo_rows.extend(chosen.to_dict(orient="records"))

        audits.append({
            "v11_row_index": int(row.v11_row_index),
            "shard_index": shard_index,
            "species": str(row.species),
            "inat_taxon_id": taxon_id,
            "hash_rank_one_based": int(ranked.index(taxon_id) + 1),
            "selected_raw_photo_target": 100,
            "total_results_probe": int(total_results),
            "candidate_pages": json.dumps(list(pages), separators=(",", ":")),
            "candidate_page_errors": json.dumps(page_errors, sort_keys=True, separators=(",", ":")),
            "raw_results_from_candidate_pages": int(len(raw_rows)),
            "locally_eligible": int(locally_eligible),
            "prior_colour_experiment_excluded": int(prior_excluded),
            "wrong_taxon": int(wrong_taxon),
            "deduplicated_metadata_rows": deduplicated,
            "after_observer_cap": capped_n,
            "maximum_span_km_after_observer_cap": maximum_span_km(capped),
            "full_target": bool(full_target),
            "selected_candidate_rows": int(len(chosen)),
        })

    audit = pd.DataFrame(audits).sort_values("v11_row_index", kind="mergesort").reset_index(drop=True)
    photos = pd.DataFrame(photo_rows)
    candidate_columns = PARSER_COLUMNS + ["row_hash", "v11_row_index", "shard_index"]
    for column in candidate_columns:
        if column not in photos.columns:
            photos[column] = pd.Series(dtype="object")
    photos = photos[candidate_columns]
    if len(photos):
        photos = photos.sort_values(["v11_row_index", "row_hash"], kind="mergesort").reset_index(drop=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = args.output_dir / f"v11_shard_{shard_index:02d}_species.csv"
    photos_path = args.output_dir / f"v11_shard_{shard_index:02d}_photos.csv"
    manifest_path = args.output_dir / f"v11_shard_{shard_index:02d}.json"
    audit.to_csv(audit_path, index=False, lineterminator="\n")
    photos.to_csv(photos_path, index=False, lineterminator="\n")
    manifest = {
        "schema": "ttf_v11_fresh_plant_acquisition_shard_v0.1",
        "status": "complete_metadata_only_v11_fresh_plant_acquisition_shard",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "shard_species": int(len(shard)),
        "rank_range_one_based": [1001, 1250],
        "species_hash_seed": seed,
        "previous_1000_taxon_digest": taxon_digest(pd.Series(previous, dtype="int64")),
        "v11_250_taxon_digest": taxon_digest(pd.Series(fresh_taxa, dtype="int64")),
        "request_attempts": int(request_attempts),
        "request_errors": int(request_errors),
        "full_target_species": int(audit["full_target"].sum()),
        "candidate_rows": int(len(photos)),
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "synthetic_worlds_run": 0,
        "candidate_performance_evaluated": False,
        "lineage": {
            "fcp_commit": commit,
            "capacity_manifest_blob_sha": git_blob(fcp_root, source["capacity_manifest"]),
            "capacity_selected_species_blob_sha": git_blob(fcp_root, source["capacity_selected_species"]),
            "candidate_species_budget_blob_sha": git_blob(fcp_root, source["candidate_species_budget"]),
            "candidate_acquisition_contract_blob_sha": git_blob(fcp_root, source["candidate_acquisition_contract"]),
            "prior_colour_exclusion_sha256": exclusion_hashes,
            "species_audit_sha256": sha256_file(audit_path),
            "candidate_photos_sha256": sha256_file(photos_path),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"shard": shard_index, "species": len(shard), "full_target": int(audit["full_target"].sum()), "request_errors": request_errors, "rows": len(photos)}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
