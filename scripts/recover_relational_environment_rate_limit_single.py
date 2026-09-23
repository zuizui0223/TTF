#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import acquire_relational_environment_occurrences as acquisition

GBIF = acquisition.GBIF
USER_AGENT = "ttf-relational-environment-rate-limit-singleton/0.1 (https://github.com/zuizui0223/TTF)"
FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def slow_get_json(path: str, params: dict[str, object], retries: int = 14) -> dict:
    url = f"{GBIF}/{path}?{urlencode(params)}"
    error: Exception | None = None
    detail = ""
    for attempt in range(retries):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urlopen(request, timeout=120) as response:
                payload = json.loads(response.read().decode("utf-8"))
            time.sleep(3.0)
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
                wait = float(retry_after) if retry_after else max(60.0, float(2 ** attempt))
            except (TypeError, ValueError):
                wait = max(60.0, float(2 ** attempt))
            time.sleep(min(600.0, wait))
        except Exception as exc:
            error = exc
            detail = f"{type(exc).__name__}: {exc}"
            time.sleep(min(180.0, max(60.0, float(2 ** attempt))))
    suffix = f" ({detail})" if detail else ""
    raise RuntimeError(f"GBIF request failed: {url}{suffix}") from error


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trigger", type=Path, required=True)
    ap.add_argument("--species-index", type=int, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    trigger = json.loads(args.trigger.read_text())
    if trigger.get("schema") != "ttf_relational_environment_rate_limit_singleton_trigger_v0.1":
        raise RuntimeError("unexpected singleton recovery trigger")
    if trigger.get("status") != "FROZEN_SINGLETON_HTTP_RECOVERY_BEFORE_RELATION_RESULT":
        raise RuntimeError("singleton recovery is not frozen")
    names = list(map(str, trigger["species"]))
    if not 0 <= args.species_index < len(names):
        raise RuntimeError("species index outside frozen singleton recovery")
    species = names[args.species_index]
    if trigger.get("relation_result_seen") is not False or trigger.get("genetic_response_used") is not False:
        raise RuntimeError("singleton recovery was frozen after result opening")

    acquisition.get_json = slow_get_json
    time.sleep(float(trigger["initial_cooldown_seconds"]))
    try:
        retained, ledger = acquisition.fetch_species(species)
    except Exception as exc:
        retained, ledger = [], {
            "species": species,
            "status": "REQUEST_ERROR",
            "error": f"{type(exc).__name__}: {exc}",
        }

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(retained)

    payload = {
        "schema": "ttf_relational_environment_rate_limit_singleton_result_v0.1",
        "species_index": args.species_index,
        "species": [ledger],
        "request_error_count": int(ledger["status"] == "REQUEST_ERROR"),
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"species": species, "status": ledger["status"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
