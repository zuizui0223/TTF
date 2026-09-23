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

import scripts.acquire_relational_environment_occurrences as acquisition

GBIF = acquisition.GBIF
USER_AGENT = "ttf-relational-environment-rate-limit-recovery/0.1 (https://github.com/zuizui0223/TTF)"
FIELDS = ["species", "source_key", "latitude", "longitude", "priority_rank", "priority_sha256"]


def slow_get_json(path: str, params: dict[str, object], retries: int = 12) -> dict:
    url = f"{GBIF}/{path}?{urlencode(params)}"
    error: Exception | None = None
    detail = ""
    for attempt in range(retries):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urlopen(request, timeout=120) as response:
                payload = json.loads(response.read().decode("utf-8"))
            time.sleep(2.0)
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
                wait = float(retry_after) if retry_after else max(30.0, float(2 ** attempt))
            except (TypeError, ValueError):
                wait = max(30.0, float(2 ** attempt))
            time.sleep(min(300.0, wait))
        except Exception as exc:
            error = exc
            detail = f"{type(exc).__name__}: {exc}"
            time.sleep(min(120.0, max(30.0, float(2 ** attempt))))
    suffix = f" ({detail})" if detail else ""
    raise RuntimeError(f"GBIF request failed: {url}{suffix}") from error


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_relational_environment_rate_limit_recovery_v0.1":
        raise RuntimeError("unexpected rate-limit recovery schema")
    if rule.get("status") != "FROZEN_HTTP_ONLY_RECOVERY_BEFORE_RELATION_RESULT":
        raise RuntimeError("rate-limit recovery is not frozen")
    if rule.get("relation_result_seen") is not False or rule.get("genetic_response_used") is not False:
        raise RuntimeError("rate-limit recovery was frozen after outcome opening")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("rate-limit recovery response firewall is open")

    names = list(map(str, rule["source_request_error_species"]))
    if len(names) != 4 or len(set(names)) != 4:
        raise RuntimeError("expected exact four frozen HTTP-429 species")

    acquisition.get_json = slow_get_json
    time.sleep(float(rule["transport_only_changes"]["initial_cooldown_seconds"]))

    retained_rows = []
    ledgers = []
    for species in names:
        try:
            retained, ledger = acquisition.fetch_species(species)
        except Exception as exc:
            retained, ledger = [], {
                "species": species,
                "status": "REQUEST_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            }
        retained_rows.extend(retained)
        ledgers.append(ledger)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(retained_rows)

    payload = {
        "schema": "ttf_relational_environment_rate_limit_recovery_result_v0.1",
        "species_requested": names,
        "request_error_count": sum(row["status"] == "REQUEST_ERROR" for row in ledgers),
        "species": ledgers,
        "response_firewall": {
            "Study_B_sequence_identity_opened": False,
            "Study_B_pairwise_genetic_distances_opened": False,
            "Study_B_T_st_computed": False,
            "Study_B_beta_R_computed": False,
        },
    }
    args.output_ledger.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "request_error_count": payload["request_error_count"],
        "statuses": {row["species"]: row["status"] for row in ledgers},
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
