#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--expected-species", type=int, default=1000)
    args = ap.parse_args()

    shards = []
    for path in sorted(args.input_dir.glob("*.json")):
        payload = json.loads(path.read_text())
        if payload.get("schema") == "ttf_relational_gbif_occurrence_feasibility_shard_v0.1":
            shards.append(payload)
    if not shards:
        raise RuntimeError("no relational GBIF shard payloads found")
    digests = {x["source_census_digest_sha256"] for x in shards}
    if len(digests) != 1:
        raise RuntimeError("source census digest drift across shards")
    results = [row for shard in shards for row in shard["results"]]
    names = [row["species"] for row in results]
    if len(names) != len(set(names)):
        raise RuntimeError("duplicate species across shards")
    if len(results) != args.expected_species:
        raise RuntimeError(f"expected {args.expected_species} species, got {len(results)}")
    counts = Counter(str(row["status"]) for row in results)
    usable = counts.get("usable_lower_bound", 0)
    payload = {
        "schema": "ttf_relational_gbif_occurrence_feasibility_v0.1",
        "status": "PASS_OCCURRENCE_AVAILABILITY" if usable >= 500 else "NOT_YET_PASS_OCCURRENCE_AVAILABILITY",
        "source_census_digest_sha256": next(iter(digests)),
        "species": len(results),
        "usable_ge_30_thinned": usable,
        "status_counts": dict(sorted(counts.items())),
        "pass_rule": ">=500 of 1000 fresh species reach >=30 occurrence points after provisional 10-km thinning in 2010-2026",
        "interpretation": "Occurrence-availability gate only. This does not freeze the environmental variables, PCA, niche metric, or final niche sampling representation.",
        "results": sorted(results, key=lambda x: x["species"]),
        "firewall": {
            "nucleotide_identity_opened": False,
            "pairwise_genetic_distances_opened": False,
            "T_st_computed": False,
            "ecology_genetics_association_computed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: payload[k] for k in ("status", "species", "usable_ge_30_thinned", "status_counts")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
