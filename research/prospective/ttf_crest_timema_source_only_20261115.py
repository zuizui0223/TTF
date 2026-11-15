"""Response-blind Timema source identity and graph preflight.

Consumes *only* declared metadata columns. Never decodes mating counts,
host choices, Fst, or other biological responses. Byte hashes must be
established independently from a pinned public archive version.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
from pathlib import Path

PAIR = "Pair"
HOST_META = ("Pair", "Species1", "Species2")


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _meta_rows(raw: bytes, required: tuple[str, ...], label: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter=",")
    cols = reader.fieldnames or ()
    if len(cols) != len(set(cols)) or not set(required).issubset(cols):
        raise ValueError(f"{label}: duplicated or missing metadata columns")
    out = []
    for row in reader:
        vals = {k: (row.get(k) or "").strip() for k in required}
        if not all(vals.values()):
            raise ValueError(f"{label}: missing required metadata key")
        out.append(vals)
    return out


def _unique_pairs(rows: list[dict[str, str]], label: str) -> set[str]:
    ids = [r[PAIR] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{label}: duplicate Pair IDs")
    return set(ids)


def species_components(host_rows: list[dict[str, str]]) -> list[list[str]]:
    adjacency: dict[str, set[str]] = collections.defaultdict(set)
    for row in host_rows:
        a, b = row["Species1"], row["Species2"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    remaining = set(adjacency)
    components = []
    while remaining:
        root = min(remaining)
        seen, pending = {root}, [root]
        while pending:
            for node in sorted(adjacency[pending.pop()]):
                if node not in seen:
                    seen.add(node)
                    pending.append(node)
        remaining -= seen
        components.append(sorted(seen))
    return sorted(components, key=lambda c: (-len(c), c))


def audit_bytes(host_bytes: bytes, count_bytes: bytes,
                expected_host_sha256: str, expected_count_sha256: str) -> dict:
    for name, raw, declared in (
        ("HostRIDat.csv", host_bytes, expected_host_sha256),
        ("CountData.csv", count_bytes, expected_count_sha256),
    ):
        if not isinstance(declared, str) or len(declared) != 64 or sha256_bytes(raw) != declared:
            raise ValueError(f"{name}: source SHA256 mismatch or not pinned")
    host = _meta_rows(host_bytes, HOST_META, "HostRIDat")
    mating = _meta_rows(count_bytes, (PAIR,), "CountData")
    hids = _unique_pairs(host, "HostRIDat")
    mids = _unique_pairs(mating, "CountData")
    if not mids <= hids:
        raise ValueError("CountData references Pair absent from HostRIDat")
    comp = species_components(host)
    return {
        "status": "SOURCE_ONLY_NO_BIOLOGICAL_RESPONSES_READ",
        "original_host_pairs": len(host),
        "original_mating_table_pair_keys": len(mating),
        "pair_id_overlap": len(hids & mids),
        "host_only_pair_ids": sorted(hids - mids),
        "species_labels": sorted({v for r in host for k, v in r.items() if k != PAIR}),
        "species_graph_components": comp,
        "species_graph_component_sizes": [len(c) for c in comp],
        "source_sha256": {"HostRIDat.csv": expected_host_sha256, "CountData.csv": expected_count_sha256},
        "mating_data_outcome_completeness": "NOT_CHECKED",
        "species_disjoint_transfer_admissibility": "NOT_QUALIFIED",
        "note": "Row-key overlap is not outcome completeness or independent dyad replication.",
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--host", required=True, type=Path)
    p.add_argument("--mating", required=True, type=Path)
    p.add_argument("--host-sha256", required=True)
    p.add_argument("--mating-sha256", required=True)
    a = p.parse_args()
    print(json.dumps(
        audit_bytes(a.host.read_bytes(), a.mating.read_bytes(), a.host_sha256, a.mating_sha256),
        ensure_ascii=False, indent=2, sort_keys=True
    ))


if __name__ == "__main__":
    main()
