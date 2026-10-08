#!/usr/bin/env python3
"""Descriptive historical host temporal-plausibility flag census only.

This script does NOT change frozen taxon membership, predictors, LGM times,
qualification rules or genetic response availability.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

FROZEN_CANDIDATES = "c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b"
# The previously frozen five-name candidate concern list was response-blind.
TEMPORAL_REVIEW_NAMES = frozenset({
    "Brassica oleracea", "Malus domestica", "Prunus domestica",
    "Prunus persica", "Zea mays",
})


def audit_frozen_hosts(path: Path) -> dict:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != FROZEN_CANDIDATES:
        raise RuntimeError("not the exact original 642-species candidate table")
    with path.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 642 or len({r["species"] for r in rows}) != 642:
        raise RuntimeError("candidate universe drift")

    flagged = []
    for row in rows:
        hosts = {name.strip() for name in row["accepted_host_names"].split(";") if name.strip()}
        concern = sorted(hosts & TEMPORAL_REVIEW_NAMES)
        if concern:
            flagged.append({"species": row["species"], "concern_host_names": concern})
    return {
        "schema": "ttf_historical_host_memory_temporal_plausibility_flags_v0.1",
        "status": "RESPONSE_BLIND_DESCRIPTIVE_TEMPORAL_REVIEW_ONLY",
        "source_candidate_sha256": digest,
        "source_panel_species": 642,
        "previously_frozen_concern_names": sorted(TEMPORAL_REVIEW_NAMES),
        "flagged_insects": len(flagged),
        "unflagged_insects": 642 - len(flagged),
        "flagged_species": sorted(flagged, key=lambda r: r["species"]),
        "interpretation": "Present HOSTS/WCVP interaction labels do not prove the same consumer-host interaction or domestic taxon existed at 21 ka; this is not an automated exclusion.",
        "frozen_panel_modified": False,
        "candidate_species_removed": 0,
        "backfill": False,
        "genetic_response_opened": False,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    payload = audit_frozen_hosts(args.candidate)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({"flagged_insects": payload["flagged_insects"],
                      "unflagged_insects": payload["unflagged_insects"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
