#!/usr/bin/env python3
"""Response-blind GBIF taxon/count pre-census for Relational TTF.

This stage never reads nucleotide identity, genetic distance, T_st, or any
empirical ecology-genetics association. It only asks whether the frozen fresh
candidate species have enough recent georeferenced occurrence availability to
justify a later environmental-niche acquisition.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

try:
    from scripts.freeze_gbif_animal_geometry import accepted_taxon_match, get_json
except ModuleNotFoundError:
    # Support both direct CLI execution and importlib-based test loading without
    # requiring scripts/ to be an importable package.
    import importlib.util

    _sibling = Path(__file__).with_name("freeze_gbif_animal_geometry.py")
    _spec = importlib.util.spec_from_file_location("_ttf_freeze_gbif_animal_geometry", _sibling)
    if _spec is None or _spec.loader is None:
        raise RuntimeError(f"cannot load GBIF helper module from {_sibling}")
    _module = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_module)
    accepted_taxon_match = _module.accepted_taxon_match
    get_json = _module.get_json

SCHEMA = "ttf_relational_gbif_precensus_v0.1"


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def occurrence_count(usage_key: int, start_year: int, end_year: int) -> int:
    payload = get_json(
        "occurrence/search",
        {
            "taxon_key": int(usage_key),
            "has_coordinate": "true",
            "has_geospatial_issue": "false",
            "occurrence_status": "present",
            "year": f"{int(start_year)},{int(end_year)}",
            "limit": 1,
        },
    )
    return int(payload.get("count") or 0)


def inspect_species(
    species: str,
    *,
    start_year: int,
    end_year: int,
    minimum_count: int,
    allow_canonical_fuzzy: bool,
) -> dict[str, object]:
    match = get_json("species/match", {"name": species, "strict": "true"})
    usage_key = int(match.get("usageKey") or 0)
    accepted, acceptance_rule = accepted_taxon_match(
        species,
        match,
        allow_canonical_fuzzy=allow_canonical_fuzzy,
    )
    row: dict[str, object] = {
        "species": species,
        "gbif_usage_key": usage_key,
        "gbif_match_type": str(match.get("matchType") or ""),
        "gbif_canonical_name": str(match.get("canonicalName") or ""),
        "gbif_rank": str(match.get("rank") or ""),
        "gbif_match_acceptance_rule": acceptance_rule,
        "match_accepted": bool(accepted and usage_key > 0),
        "occurrence_count_2010_endyear": 0,
        "count_ge_minimum": False,
        "status": "",
    }
    if not accepted or usage_key <= 0:
        row["status"] = "REJECTED_GBIF_TAXON_MATCH"
        return row

    total = occurrence_count(usage_key, start_year, end_year)
    row["occurrence_count_2010_endyear"] = int(total)
    row["count_ge_minimum"] = bool(total >= int(minimum_count))
    row["status"] = "COUNT_GE_MINIMUM" if row["count_ge_minimum"] else "COUNT_BELOW_MINIMUM"
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    ap.add_argument("--start-year", type=int, default=2010)
    ap.add_argument("--end-year", type=int, default=2026)
    ap.add_argument("--minimum-count", type=int, default=30)
    ap.add_argument("--start-index", type=int, default=0)
    ap.add_argument("--max-species", type=int, default=0, help="0 means all remaining candidates")
    ap.add_argument("--pause-seconds", type=float, default=0.10)
    ap.add_argument("--allow-canonical-fuzzy", action="store_true")
    args = ap.parse_args()

    if args.start_year > args.end_year:
        raise ValueError("start-year must be <= end-year")
    if args.minimum_count < 1:
        raise ValueError("minimum-count must be >= 1")
    if args.start_index < 0 or args.max_species < 0 or args.pause_seconds < 0:
        raise ValueError("start/max/pause arguments must be non-negative")

    rows = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    if not rows or "species" not in rows[0]:
        raise ValueError("candidate CSV must contain a species column")
    species = [str(row["species"]).strip() for row in rows]
    if any(not s for s in species) or len(species) != len(set(species)):
        raise ValueError("candidate species must be unique and non-empty")

    stop = len(species) if args.max_species == 0 else min(
        len(species), args.start_index + args.max_species
    )
    selected = species[args.start_index:stop]
    results: list[dict[str, object]] = []
    for i, name in enumerate(selected):
        try:
            result = inspect_species(
                name,
                start_year=args.start_year,
                end_year=args.end_year,
                minimum_count=args.minimum_count,
                allow_canonical_fuzzy=args.allow_canonical_fuzzy,
            )
        except Exception as exc:  # network/provider boundary; preserve row and continue
            result = {
                "species": name,
                "gbif_usage_key": 0,
                "gbif_match_type": "",
                "gbif_canonical_name": "",
                "gbif_rank": "",
                "gbif_match_acceptance_rule": "",
                "match_accepted": False,
                "occurrence_count_2010_endyear": 0,
                "count_ge_minimum": False,
                "status": "REQUEST_ERROR",
                "error": f"{type(exc).__name__}: {exc}",
            }
        results.append(result)
        if args.pause_seconds and i + 1 < len(selected):
            time.sleep(args.pause_seconds)

    fields = [
        "species",
        "gbif_usage_key",
        "gbif_match_type",
        "gbif_canonical_name",
        "gbif_rank",
        "gbif_match_acceptance_rule",
        "match_accepted",
        "occurrence_count_2010_endyear",
        "count_ge_minimum",
        "status",
        "error",
    ]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)

    counts: dict[str, int] = {}
    for row in results:
        key = str(row["status"])
        counts[key] = counts.get(key, 0) + 1
    ledger = {
        "schema": SCHEMA,
        "status": "RESPONSE_BLIND_GBIF_COUNT_PRECENSUS_COMPLETE",
        "candidate_csv_sha256": sha256_path(args.candidates),
        "output_csv_sha256": sha256_path(args.output_csv),
        "provider": "GBIF occurrence/search + species/match",
        "filters": {
            "start_year": args.start_year,
            "end_year": args.end_year,
            "has_coordinate": True,
            "has_geospatial_issue": False,
            "occurrence_status": "present",
            "minimum_count_candidate": args.minimum_count,
        },
        "slice": {
            "start_index_zero_based": args.start_index,
            "species_requested": len(selected),
            "candidate_total": len(species),
        },
        "status_counts": dict(sorted(counts.items())),
        "important_boundary": (
            "COUNT_GE_MINIMUM is only a provider-availability triage. It is not the "
            "final >=30 post-duplicate/post-thinning niche criterion and cannot authorize "
            "Schoener's D or any genetic response opening."
        ),
        "response_firewall": {
            "sequence_identity_opened": False,
            "pairwise_genetic_distances_opened": False,
            "T_st_computed": False,
            "environmental_overlap_computed": False,
            "ecology_genetics_association_computed": False,
        },
    }
    args.output_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.output_ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status_counts": ledger["status_counts"], "species": len(selected)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
