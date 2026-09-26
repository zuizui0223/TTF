#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from ttf.checkpointed_gbif_occurrence import (
    deterministic_page_offsets,
    missing_page_offsets,
    request_timeout_seconds,
    species_state_key,
)


GBIF = "https://api.gbif.org/v1"
USER_AGENT = "ttf-butterfly-resource-envelope/0.1 (https://github.com/zuizui0223/TTF)"
RETRYABLE_HTTP = {429, 500, 502, 503, 504}


def atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def get_json(
    path: str,
    params: dict[str, object],
    *,
    deadline: float,
    attempts: int = 3,
) -> dict:
    url = f"{GBIF}/{path}?{urlencode(params)}"
    last_detail = ""
    for attempt in range(attempts):
        timeout = request_timeout_seconds(deadline, per_request_cap=30.0)
        try:
            request = Request(
                url,
                headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
            )
            with urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            body = ""
            try:
                body = exc.read(512).decode("utf-8", errors="replace").replace("\n", " ")
            except Exception:
                pass
            last_detail = f"HTTP {exc.code}: {body[:300]}"
            if exc.code not in RETRYABLE_HTTP:
                raise RuntimeError(f"GBIF request failed: {url} ({last_detail})") from exc
            retry_after = exc.headers.get("Retry-After")
            try:
                delay = float(retry_after) if retry_after else float(2**attempt)
            except (TypeError, ValueError):
                delay = float(2**attempt)
        except TimeoutError:
            raise
        except Exception as exc:
            last_detail = f"{type(exc).__name__}: {exc}"
            delay = float(2**attempt)

        remaining = deadline - time.monotonic()
        if remaining <= 1.0:
            raise TimeoutError("species total deadline exhausted")
        time.sleep(min(10.0, delay, max(0.0, remaining - 1.0)))

    raise RuntimeError(f"GBIF request failed after {attempts} attempts: {url} ({last_detail})")


def load_pilot(path: Path) -> tuple[str, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "ttf_butterfly_resource_envelope_pilot_v0.1":
        raise RuntimeError("unexpected pilot schema")
    if payload.get("status") != "EXPLORATORY_RESPONSE_BLIND_PILOT_SELECTED":
        raise RuntimeError("pilot is not in response-blind selected state")
    names = tuple(map(str, payload.get("pilot_species", [])))
    if len(names) != 10 or len(set(names)) != 10:
        raise RuntimeError("expected exact ten-species exploratory pilot")
    return names


def metadata_for_species(
    species: str,
    state_dir: Path,
    *,
    deadline: float,
    maximum_pages: int,
) -> dict:
    meta_path = state_dir / "metadata.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("species") != species:
            raise RuntimeError("checkpoint species identity drift")
        return meta

    match = get_json(
        "species/match",
        {"name": species, "strict": "true"},
        deadline=deadline,
    )
    usage_key = int(match.get("usageKey") or 0)
    canonical = str(match.get("canonicalName") or "").strip()
    rank = str(match.get("rank") or "").strip().upper()
    match_type = str(match.get("matchType") or "").strip().upper()
    accepted = (
        usage_key > 0
        and rank == "SPECIES"
        and match_type in {"EXACT", "FUZZY"}
        and canonical.lower() == species.lower()
    )
    if not accepted:
        meta = {
            "species": species,
            "status": "REJECTED_GBIF_TAXON_MATCH",
            "usage_key": usage_key,
            "canonical_name": canonical,
            "rank": rank,
            "match_type": match_type,
            "page_offsets": [],
        }
        atomic_json(meta_path, meta)
        return meta

    count = get_json(
        "occurrence/search",
        {
            "taxonKey": usage_key,
            "hasCoordinate": "true",
            "hasGeospatialIssue": "false",
            "occurrenceStatus": "PRESENT",
            "year": "2010,2026",
            "limit": 1,
        },
        deadline=deadline,
    )
    total = int(count.get("count") or 0)
    offsets = deterministic_page_offsets(
        total,
        page_size=300,
        maximum_pages=maximum_pages,
    )
    meta = {
        "species": species,
        "status": "MATCHED",
        "usage_key": usage_key,
        "canonical_name": canonical,
        "rank": rank,
        "match_type": match_type,
        "total_coordinate_records_2010_2026": total,
        "page_size": 300,
        "page_offsets": list(offsets),
    }
    atomic_json(meta_path, meta)
    return meta


def fetch_species_pages(
    species: str,
    state_root: Path,
    *,
    species_seconds: float,
    maximum_pages: int,
) -> dict:
    state_dir = state_root / species_state_key(species)
    state_dir.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + float(species_seconds)

    try:
        meta = metadata_for_species(
            species,
            state_dir,
            deadline=deadline,
            maximum_pages=maximum_pages,
        )
    except Exception as exc:
        return {
            "species": species,
            "status": "PARTIAL_TRANSPORT",
            "stage": "metadata",
            "error": f"{type(exc).__name__}: {exc}",
            "completed_offsets": [],
            "missing_offsets": [],
        }

    if meta.get("status") == "REJECTED_GBIF_TAXON_MATCH":
        return {
            "species": species,
            "status": "REJECTED_GBIF_TAXON_MATCH",
            "completed_offsets": [],
            "missing_offsets": [],
        }

    offsets = tuple(map(int, meta["page_offsets"]))
    pages_dir = state_dir / "pages"
    pages_dir.mkdir(exist_ok=True)
    completed = []
    for offset in offsets:
        page_path = pages_dir / f"offset_{offset:06d}.json"
        if page_path.exists():
            completed.append(offset)
    pending = list(missing_page_offsets(offsets, completed))

    error = None
    for offset in pending:
        try:
            payload = get_json(
                "occurrence/search",
                {
                    "taxonKey": int(meta["usage_key"]),
                    "hasCoordinate": "true",
                    "hasGeospatialIssue": "false",
                    "occurrenceStatus": "PRESENT",
                    "year": "2010,2026",
                    "limit": 300,
                    "offset": offset,
                },
                deadline=deadline,
            )
            records = []
            for item in payload.get("results") or []:
                key = int(item.get("key") or 0)
                lat = item.get("decimalLatitude")
                lon = item.get("decimalLongitude")
                if key <= 0 or lat is None or lon is None:
                    continue
                records.append(
                    {
                        "key": key,
                        "decimalLatitude": float(lat),
                        "decimalLongitude": float(lon),
                        "year": item.get("year"),
                        "basisOfRecord": item.get("basisOfRecord"),
                        "datasetKey": item.get("datasetKey"),
                    }
                )
            atomic_json(
                pages_dir / f"offset_{offset:06d}.json",
                {
                    "species": species,
                    "usage_key": int(meta["usage_key"]),
                    "offset": offset,
                    "records": records,
                },
            )
            completed.append(offset)
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
            break

    missing = list(missing_page_offsets(offsets, completed))
    return {
        "species": species,
        "status": "COMPLETE" if not missing else "PARTIAL_TRANSPORT",
        "usage_key": int(meta["usage_key"]),
        "planned_offsets": list(offsets),
        "completed_offsets": sorted(set(completed)),
        "missing_offsets": missing,
        "error": error,
    }


def combine_records(state_root: Path, species: tuple[str, ...]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for name in species:
        state_dir = state_root / species_state_key(name)
        for page in sorted((state_dir / "pages").glob("offset_*.json")):
            payload = json.loads(page.read_text(encoding="utf-8"))
            for record in payload.get("records", []):
                rows.append(
                    {
                        "species": name,
                        "gbif_key": int(record["key"]),
                        "latitude": float(record["decimalLatitude"]),
                        "longitude": float(record["decimalLongitude"]),
                        "year": record.get("year"),
                        "basis_of_record": record.get("basisOfRecord"),
                        "dataset_key": record.get("datasetKey"),
                    }
                )
    unique: dict[tuple[str, int], dict[str, object]] = {}
    for row in rows:
        unique[(str(row["species"]), int(row["gbif_key"]))] = row
    return [unique[key] for key in sorted(unique)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot-json", type=Path, required=True)
    ap.add_argument("--state-dir", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--species-seconds", type=float, default=600.0)
    ap.add_argument("--maximum-pages", type=int, default=12)
    args = ap.parse_args()

    species = load_pilot(args.pilot_json)
    args.state_dir.mkdir(parents=True, exist_ok=True)
    ledgers = [
        fetch_species_pages(
            name,
            args.state_dir,
            species_seconds=args.species_seconds,
            maximum_pages=args.maximum_pages,
        )
        for name in species
    ]

    rows = combine_records(args.state_dir, species)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species",
        "gbif_key",
        "latitude",
        "longitude",
        "year",
        "basis_of_record",
        "dataset_key",
    ]
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "schema": "ttf_butterfly_resource_envelope_occurrence_pilot_v0.1",
        "status": (
            "COMPLETE"
            if all(row["status"] != "PARTIAL_TRANSPORT" for row in ledgers)
            else "RESUMABLE_PARTIAL_TRANSPORT"
        ),
        "species": ledgers,
        "occurrence_rows": len(rows),
        "transport": {
            "page_checkpointing": True,
            "resume_only_unfinished_pages": True,
            "species_total_deadline_seconds": float(args.species_seconds),
            "request_timeout_cap_seconds": 30.0,
            "request_attempts": 3,
            "maximum_pages_per_species": int(args.maximum_pages),
            "page_size": 300,
        },
        "scientific_scope": {
            "exploratory": True,
            "genetic_response_used": False,
            "phylogatr_coordinates_used_as_range_data": False,
        },
    }
    args.output_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.output_ledger.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": payload["status"], "rows": len(rows)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
