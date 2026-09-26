#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import time
import subprocess
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
    request_seconds: float = 30.0,
    attempts: int = 3,
) -> dict:
    """Fetch one GBIF JSON response with a true wall-clock request cap.

    urllib/socket timeouts do not bound the total duration of response.read().
    curl --max-time does, so every page request is externally bounded while the
    deterministic query itself remains unchanged.
    """
    url = f"{GBIF}/{path}?{urlencode(params)}"
    last_detail = ""
    for attempt in range(attempts):
        timeout = request_timeout_seconds(
            deadline,
            per_request_cap=float(request_seconds),
        )
        cmd = [
            "curl",
            "--location",
            "--silent",
            "--show-error",
            "--max-time",
            f"{timeout:.3f}",
            "--header",
            f"User-Agent: {USER_AGENT}",
            "--header",
            "Accept: application/json",
            "--write-out",
            "\\n%{http_code}",
            url,
        ]
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout + 2.0,
                check=False,
            )
        except subprocess.TimeoutExpired:
            proc = None
            last_detail = f"hard subprocess timeout after {timeout:.1f}s"
            delay = float(2**attempt)
        else:
            stdout = proc.stdout or ""
            if "\n" in stdout:
                body, code_text = stdout.rsplit("\n", 1)
            else:
                body, code_text = stdout, "000"
            try:
                status_code = int(code_text.strip())
            except ValueError:
                status_code = 0

            if proc.returncode == 0 and status_code == 200:
                try:
                    return json.loads(body)
                except json.JSONDecodeError as exc:
                    last_detail = f"JSONDecodeError: {exc}"
                    delay = float(2**attempt)
                else:
                    delay = float(2**attempt)
            else:
                stderr = (proc.stderr or "").strip().replace("\n", " ")
                last_detail = (
                    f"curl_rc={proc.returncode} HTTP {status_code}: "
                    f"{stderr[:180]} {body[:180]}"
                ).strip()
                if status_code and status_code not in RETRYABLE_HTTP:
                    raise RuntimeError(
                        f"GBIF request failed: {url} ({last_detail})"
                    )
                delay = float(2**attempt)

        remaining = deadline - time.monotonic()
        if remaining <= 1.0:
            raise TimeoutError("species total deadline exhausted")
        time.sleep(min(10.0, delay, max(0.0, remaining - 1.0)))

    raise RuntimeError(
        f"GBIF request failed after {attempts} attempts: {url} ({last_detail})"
    )

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


RESOLVER_VERSION = "gbif-exact-accepted-species-v0.2"


def _accepted_match(species: str, payload: dict) -> tuple[bool, int, str, str, str]:
    usage_key = int(payload.get("usageKey") or payload.get("key") or 0)
    canonical = str(payload.get("canonicalName") or "").strip()
    rank = str(payload.get("rank") or "").strip().upper()
    match_type = str(payload.get("matchType") or "").strip().upper()
    accepted = (
        usage_key > 0
        and rank == "SPECIES"
        and canonical.lower() == species.lower()
        and (
            not match_type
            or match_type in {"EXACT", "FUZZY"}
        )
    )
    return accepted, usage_key, canonical, rank, match_type


def resolve_species_metadata(
    species: str,
    *,
    deadline: float,
    request_seconds: float,
) -> dict:
    strict = get_json(
        "species/match",
        {"name": species, "strict": "true"},
        deadline=deadline,
        request_seconds=request_seconds,
    )
    accepted, usage_key, canonical, rank, match_type = _accepted_match(
        species, strict
    )
    if accepted:
        return {
            "species": species,
            "status": "MATCHED",
            "usage_key": usage_key,
            "canonical_name": canonical,
            "rank": rank,
            "match_type": match_type,
            "resolver_version": RESOLVER_VERSION,
            "resolver_route": "species_match_strict",
        }

    search = get_json(
        "species/search",
        {"q": species, "rank": "SPECIES", "limit": 20},
        deadline=deadline,
        request_seconds=request_seconds,
    )
    candidates = []
    for item in search.get("results") or []:
        item_rank = str(item.get("rank") or "").strip().upper()
        item_canonical = str(item.get("canonicalName") or "").strip()
        item_status = str(
            item.get("taxonomicStatus") or item.get("status") or ""
        ).strip().upper()
        item_key = int(item.get("key") or 0)
        if (
            item_key > 0
            and item_rank == "SPECIES"
            and item_canonical.lower() == species.lower()
            and item_status == "ACCEPTED"
        ):
            candidates.append(
                (item_key, item_canonical, item_rank, item_status)
            )
    unique = {}
    for candidate in candidates:
        unique[candidate[0]] = candidate
    if len(unique) == 1:
        item_key, item_canonical, item_rank, item_status = next(
            iter(unique.values())
        )
        return {
            "species": species,
            "status": "MATCHED",
            "usage_key": item_key,
            "canonical_name": item_canonical,
            "rank": item_rank,
            "match_type": "EXACT_SEARCH_FALLBACK",
            "taxonomic_status": item_status,
            "resolver_version": RESOLVER_VERSION,
            "resolver_route": "species_search_exact_accepted_unique",
        }

    return {
        "species": species,
        "status": "REJECTED_GBIF_TAXON_MATCH",
        "usage_key": usage_key,
        "canonical_name": canonical,
        "rank": rank,
        "match_type": match_type,
        "resolver_version": RESOLVER_VERSION,
        "resolver_route": "strict_then_exact_accepted_search",
        "exact_accepted_search_candidates": len(unique),
        "page_offsets": [],
    }


def metadata_for_species(
    species: str,
    state_dir: Path,
    *,
    deadline: float,
    maximum_pages: int,
    request_seconds: float,
) -> dict:
    meta_path = state_dir / "metadata.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("species") != species:
            raise RuntimeError("checkpoint species identity drift")
        # Keep every prior successful exact strict match. Re-evaluate only
        # previously rejected metadata under the new general resolver.
        if meta.get("status") != "REJECTED_GBIF_TAXON_MATCH":
            return meta
        if meta.get("resolver_version") == RESOLVER_VERSION:
            return meta

    meta = resolve_species_metadata(
        species,
        deadline=deadline,
        request_seconds=request_seconds,
    )
    if meta.get("status") == "REJECTED_GBIF_TAXON_MATCH":
        atomic_json(meta_path, meta)
        return meta

    usage_key = int(meta["usage_key"])
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
        request_seconds=request_seconds,
    )
    total = int(count.get("count") or 0)
    offsets = deterministic_page_offsets(
        total,
        page_size=300,
        maximum_pages=maximum_pages,
    )
    meta.update(
        {
            "total_coordinate_records_2010_2026": total,
            "page_size": 300,
            "page_offsets": list(offsets),
        }
    )
    atomic_json(meta_path, meta)
    return meta

def fetch_species_pages(
    species: str,
    state_root: Path,
    *,
    species_seconds: float,
    maximum_pages: int,
    request_seconds: float,
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
            request_seconds=request_seconds,
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
            "usage_key": int(meta.get("usage_key") or 0),
            "canonical_name": str(meta.get("canonical_name") or ""),
            "rank": str(meta.get("rank") or ""),
            "match_type": str(meta.get("match_type") or ""),
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
                request_seconds=request_seconds,
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
        "canonical_name": str(meta.get("canonical_name") or ""),
        "rank": str(meta.get("rank") or ""),
        "match_type": str(meta.get("match_type") or ""),
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
    ap.add_argument("--only-species", type=str, default=None)
    ap.add_argument("--state-dir", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--species-seconds", type=float, default=600.0)
    ap.add_argument("--maximum-pages", type=int, default=12)
    ap.add_argument("--request-seconds", type=float, default=30.0)
    args = ap.parse_args()

    species = load_pilot(args.pilot_json)
    if args.only_species is not None:
        selected = str(args.only_species).strip()
        if selected not in species:
            raise RuntimeError(f"requested species is not in frozen pilot: {selected}")
        species = (selected,)
    args.state_dir.mkdir(parents=True, exist_ok=True)
    ledgers = [
        fetch_species_pages(
            name,
            args.state_dir,
            species_seconds=args.species_seconds,
            maximum_pages=args.maximum_pages,
            request_seconds=args.request_seconds,
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
            "request_timeout_cap_seconds": float(args.request_seconds),
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
