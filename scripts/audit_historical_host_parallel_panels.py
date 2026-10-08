#!/usr/bin/env python3
"""Response-blind audit of two independently frozen historical-host panels.

This script never accesses genetic sequence characters, genetic distances,
post-IBD responses or historical-host effect estimates.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


MEMORY_SHA256 = "c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b"
CONNECTIVITY_SHA256 = "a402483497bc7bdf9b7427bfc859ce9dcc9183eebf5964d4591f352847130484"


def read_pinned_panel(path: Path, expected_sha256: str, expected_species: int) -> dict[str, frozenset[str]]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected_sha256:
        raise RuntimeError(f"frozen CSV SHA256 mismatch: {path}; got {digest}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != expected_species:
        raise RuntimeError(f"frozen species count mismatch in {path}")
    associations = {}
    for row in rows:
        name = row["species"].strip()
        raw = [x.strip() for x in row["accepted_host_ids"].split(";")]
        ids = [x[:-2] if x.endswith(".0") and x[:-2].isdigit() else x for x in raw]
        if (not name or not ids or not all(x.isdigit() for x in ids)
                or len(ids) != len(set(ids)) or name in associations):
            raise RuntimeError(f"invalid or duplicate frozen species-host identity in {path}: {name}")
        associations[name] = frozenset(ids)
    return associations


def compute_overlap(memory: dict[str, frozenset[str]], connectivity: dict[str, frozenset[str]]) -> dict:
    shared = sorted(memory.keys() & connectivity.keys())
    disagreements = sorted(name for name in shared if memory[name] != connectivity[name])
    return {
        "memory_species": len(memory),
        "connectivity_species": len(connectivity),
        "union_species": len(memory.keys() | connectivity.keys()),
        "overlap_species": len(shared),
        "overlap_species_names": shared,
        "overlap_accepted_host_ids_identical_species": len(shared) - len(disagreements),
        "overlap_accepted_host_id_disagreements": disagreements,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--memory", type=Path, required=True)
    parser.add_argument("--connectivity", type=Path, required=True)
    parser.add_argument("--frozen-receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    memory = read_pinned_panel(args.memory, MEMORY_SHA256, 642)
    connectivity = read_pinned_panel(args.connectivity, CONNECTIVITY_SHA256, 135)
    result = compute_overlap(memory, connectivity)
    frozen = json.loads(args.frozen_receipt.read_text(encoding="utf-8"))

    if (result["overlap_species"] != frozen["overlap"]["species"]
            or result["overlap_species_names"] != frozen["overlap"]["shared_species"]
            or result["overlap_accepted_host_id_disagreements"]):
        raise RuntimeError("frozen cross-program overlap / host-ID agreement drift")
    if (result["union_species"] != 704 or result["overlap_accepted_host_ids_identical_species"] != 73):
        raise RuntimeError("expected response-blind 642-by-135 audit identities not reproduced")

    payload = {
        "schema": "ttf_historical_host_parallel_panel_exact_source_replay_v0.1",
        "status": "PASS_RESPONSE_BLIND_FROZEN_SOURCE_AND_HOST_ID_REPLAY",
        "memory_sha256": MEMORY_SHA256,
        "connectivity_sha256": CONNECTIVITY_SHA256,
        "comparison": result,
        "interpretation": (
            "The 73 shared insect identities have the same recorded HOSTS/WCVP accepted "
            "plant IDs in both panels; this is database consistency, not fossil validation "
            "and not independent genetic replication."
        ),
        "genetic_response_opened": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in {"overlap_species_names", "overlap_accepted_host_id_disagreements"}},
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
