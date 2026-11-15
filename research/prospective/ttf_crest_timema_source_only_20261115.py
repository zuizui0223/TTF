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
IMPUTED_GBS_PROVENANCE = "0.15897 + nuc * 19.41020"


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
        if label == "HostRIDat" and "GBS_ref" in cols:
            vals["GBS_ref"] = (row.get("GBS_ref") or "").strip()
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


def _species_disjoint_split_envelope(host_rows: list[dict[str, str]], mating_pair_ids: set[str]) -> dict:
    """Source-only split feasibility; cross-boundary dyads are discarded.

    Reports pair-key availability, NOT usable outcomes or calibrated power.
    """
    names = sorted({r[k] for r in host_rows for k in ("Species1", "Species2")})
    n = len(names)
    if n < 2 or n > 18:
        return {"status": "NOT_ENUMERABLE", "n_species": n}
    first = names[0]
    best, eligible = None, 0
    for bits in range(1 << (n - 1)):
        training = {first} | {names[i + 1] for i in range(n - 1) if bits & (1 << i)}
        testing = set(names) - training
        if not testing:
            continue
        tr = [r["Pair"] for r in host_rows if r["Pair"] in mating_pair_ids and
              r["Species1"] in training and r["Species2"] in training]
        te = [r["Pair"] for r in host_rows if r["Pair"] in mating_pair_ids and
              r["Species1"] in testing and r["Species2"] in testing]
        if tr and te:
            eligible += 1
        score = (min(len(tr), len(te)), len(tr) + len(te), -abs(len(tr) - len(te)))
        current = {"train_species": sorted(training), "test_species": sorted(testing),
                   "train_mating_pair_keys": len(tr), "test_mating_pair_keys": len(te),
                   "cross_species_pair_keys_discarded": len(mating_pair_ids) - len(tr) - len(te),
                   "score": score}
        if best is None or score > best["score"]:
            best = current
    if best is None:
        return {"status": "NOT_EVALUABLE", "n_species": n}
    best.pop("score")
    return {"status": "GEOMETRY_ONLY_NO_OUTCOME_CLAIM", "n_species": n,
            "n_species_disjoint_nonempty_splits": eligible, "best_balance_capacity": best}


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
    heldout_envelope = _species_disjoint_split_envelope(host, mids)
    has_gbs_provenance = bool(host and "GBS_ref" in host[0])
    imputed_gbs_count = sum(r.get("GBS_ref") == IMPUTED_GBS_PROVENANCE for r in host) if has_gbs_provenance else None
    missing_gbs_ref_count = sum(not r.get("GBS_ref") for r in host) if has_gbs_provenance else None
    return {
        "status": "SOURCE_ONLY_NO_BIOLOGICAL_RESPONSES_READ",
        "original_host_pairs": len(host),
        "original_mating_table_pair_keys": len(mating),
        "pair_id_overlap": len(hids & mids),
        "host_only_pair_ids": sorted(hids - mids),
        "species_labels": sorted({r[k] for r in host for k in ("Species1", "Species2")}),
        "gbs_provenance_column_present": has_gbs_provenance,
        "gbs_proxy_imputed_pair_count": imputed_gbs_count,
        "gbs_provenance_blank_pair_count": missing_gbs_ref_count,
        "gbs_provenance_rule_source": "zgompert/TimemaRI GeneticsSummary.R",
        "species_graph_components": comp,
        "species_disjoint_source_geometry": heldout_envelope,
        "species_graph_component_sizes": [len(c) for c in comp],
        "source_sha256": {"HostRIDat.csv": expected_host_sha256, "CountData.csv": expected_count_sha256},
        "mating_data_outcome_completeness": "NOT_CHECKED",
        "species_disjoint_transfer_admissibility": "NOT_QUALIFIED",
        "observed_genetic_predictor_admissibility": "NOT_QUALIFIED_EXCLUDE_PROXY_GBS",
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
