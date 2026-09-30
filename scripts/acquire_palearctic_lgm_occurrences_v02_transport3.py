#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

try:
    import scripts.acquire_palearctic_lgm_occurrences_v02_transport2 as transport2
except ModuleNotFoundError:
    import acquire_palearctic_lgm_occurrences_v02_transport2 as transport2

from ttf.relational_external import shard_items


GBIF = "https://api.gbif.org/v1"
USER_AGENT = "ttf-palearctic-lgm-transport3/0.1 (https://github.com/zuizui0223/TTF)"
REQUEST_TIMEOUT_SECONDS = 25
REQUEST_RETRIES = 8
SPECIES_WORKERS = 2


def get_json_bounded(path: str, params: dict[str, object], retries: int = REQUEST_RETRIES) -> dict:
    url = f"{GBIF}/{path}?{urlencode(params)}"
    error: Exception | None = None
    detail = ""
    for attempt in range(int(retries)):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                payload = json.loads(response.read().decode("utf-8"))
            time.sleep(0.05)
            return payload
        except HTTPError as exc:
            error = exc
            try:
                body = exc.read(512).decode("utf-8", errors="replace").replace("\n", " ")
            except Exception:
                body = ""
            detail = f"HTTP {exc.code}: {body[:300]}"
            if exc.code not in {429, 500, 502, 503, 504}:
                break
            retry_after = exc.headers.get("Retry-After")
            try:
                wait = float(retry_after) if retry_after else float(2**attempt)
            except (TypeError, ValueError):
                wait = float(2**attempt)
            time.sleep(min(30.0, max(1.0, wait)))
        except Exception as exc:
            error = exc
            detail = f"{type(exc).__name__}: {exc}"
            time.sleep(min(30.0, float(2**attempt)))
    suffix = f" ({detail})" if detail else ""
    raise RuntimeError(f"GBIF request failed after bounded retries: {url}{suffix}") from error


# Patch only transport mechanics. The page set, filters, realm test, deterministic
# priority, thinning, and occurrence cap/minimum remain those in transport2.
transport2.get_json = get_json_bounded


def safe_fetch(species: str, realm):
    try:
        rows, meta = transport2.fetch_species_palearctic_parallel(species, realm)
    except Exception as exc:
        rows = []
        meta = {
            "species": species,
            "status": "REQUEST_ERROR",
            "error": f"{type(exc).__name__}: {exc}",
        }
    meta["transport3"] = {
        "page_workers_per_species": int(transport2.PAGE_WORKERS),
        "species_workers": SPECIES_WORKERS,
        "request_timeout_seconds": REQUEST_TIMEOUT_SECONDS,
        "request_retries": REQUEST_RETRIES,
        "scientific_page_offsets_changed": False,
    }
    return rows, meta


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", type=Path, required=True)
    ap.add_argument("--realm-geojson", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--shard-index", type=int, required=True)
    ap.add_argument("--shards", type=int, required=True)
    args = ap.parse_args()

    census = json.loads(args.census.read_text())
    if census.get("schema") != "ttf_genetic_palearctic_lgm_nonlepidoptera_census_v0.2":
        raise RuntimeError("unexpected non-Lepidoptera Palearctic-LGM census")
    if census.get("status") != "PASS_TO_FRESH_V02_GBIF_CENSUS":
        raise RuntimeError("v0.2 census did not authorize external occurrence acquisition")
    if any(bool(v) for v in census["response_firewall"].values()):
        raise RuntimeError("v0.2 genetic response firewall is open")

    names = [str(row["species"]) for row in census["retained_species"]]
    if len(names) != 23 or len(set(names)) != 23:
        raise RuntimeError("expected exact frozen 23-species non-Lepidoptera panel")
    shard = list(shard_items(names, shard_index=args.shard_index, shards=args.shards))
    realm = transport2.load_palearctic_geometry(args.realm_geojson)

    with ThreadPoolExecutor(max_workers=min(SPECIES_WORKERS, len(shard))) as pool:
        fetched = list(pool.map(lambda sp: safe_fetch(sp, realm), shard))

    output: list[dict[str, object]] = []
    ledger: list[dict[str, object]] = []
    for species, (rows, meta) in zip(shard, fetched):
        if str(meta.get("species")) != species:
            raise RuntimeError("transport3 species order drift")
        for row in rows:
            row["role"] = "both"
        output.extend(rows)
        meta["role"] = "both"
        ledger.append(meta)
        print(json.dumps({
            "species": species,
            "status": meta["status"],
            "retained": meta.get("retained_palearctic_after_exact_dedup_thinning_cap", 0),
        }, sort_keys=True), flush=True)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "species", "role", "source_key", "latitude", "longitude",
        "priority_rank", "priority_sha256",
    ]
    with args.output_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)

    payload = {
        "schema": "ttf_genetic_palearctic_lgm_occurrence_shard_v0.2",
        "shard_index": int(args.shard_index),
        "shards": int(args.shards),
        "species_requested": len(shard),
        "retained_occurrence_rows": len(output),
        "species": ledger,
        "transport": {
            "mode": "bounded_parallel_frozen_offsets",
            "page_workers_per_species": int(transport2.PAGE_WORKERS),
            "species_workers": SPECIES_WORKERS,
            "request_timeout_seconds": REQUEST_TIMEOUT_SECONDS,
            "request_retries": REQUEST_RETRIES,
            "scientific_page_offsets_changed": False,
        },
        "response_firewall": {
            "species_level_genetic_scores_used": False,
            "pairwise_v02_T_st_computed": False,
            "beta_LGM_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
