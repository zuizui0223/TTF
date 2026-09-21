#!/usr/bin/env python3
"""Freeze the fresh geometry-eligible candidate universe for Relational TTF.

The scan is deliberately response-blind. It reads taxonomy, marker labels,
FASTA headers, occurrence coordinates, and graph geometry only. It never reads
sequence characters, genetic distances, T_st, or any ecological-response
association.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from ttf.phylogatr_confirmatory import (
    _candidate_from_row,
    collapse_whitespace,
    is_coi_family_locus,
    read_genes_rows,
)

TAG = "relational-ttf-v0.1"


def species_hash(source_sha256: str, species: str) -> str:
    return hashlib.sha256(
        f"{TAG}|{source_sha256}|{species}".encode("utf-8")
    ).hexdigest()


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def species_digest(names: list[str]) -> str:
    return hashlib.sha256(("\n".join(names) + "\n").encode("utf-8")).hexdigest()


def read_exclusion_species(path: Path) -> set[str]:
    """Read species from a Phase-1 manifest, generic JSON list, CSV, or text."""
    if path.suffix.lower() == ".csv":
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if rows and "species" not in rows[0]:
            raise ValueError(f"{path} lacks species column")
        return {
            collapse_whitespace(row.get("species", ""))
            for row in rows
            if collapse_whitespace(row.get("species", ""))
        }
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, list):
            return {collapse_whitespace(x) for x in payload if collapse_whitespace(x)}
        if isinstance(payload.get("selected_panels"), dict):
            return {
                collapse_whitespace(x)
                for x in payload["selected_panels"]
                if collapse_whitespace(x)
            }
        for key in ("species", "selected_species", "excluded_species"):
            if isinstance(payload.get(key), list):
                return {
                    collapse_whitespace(x)
                    for x in payload[key]
                    if collapse_whitespace(x)
                }
        raise ValueError(f"cannot identify species list in {path}")
    return {
        collapse_whitespace(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if collapse_whitespace(line)
    }


def metadata_rows_by_species(
    root: Path,
    *,
    aliases: tuple[str, ...],
    exclusions: set[str],
    min_localities: int,
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_genes_rows(root / "genes.txt"):
        species = collapse_whitespace(row.get("species", ""))
        if (
            row.get("kingdom", "") != "Animalia"
            or row.get("class", "") == "Aves"
            or row.get("order", "") == "Chiroptera"
            or not species
            or species in exclusions
        ):
            continue
        if not is_coi_family_locus(str(row.get("gene", "")), species, aliases):
            continue
        try:
            aligned = int(float(row.get("num_seqs_aligned", "0") or 0))
        except ValueError:
            aligned = 0
        # A necessary response-blind prefilter only. Exact eligibility is still
        # determined from header-linked unique localities below.
        if aligned < int(min_localities):
            continue
        grouped[species].append(row)
    return grouped


def choose_candidate_for_species(
    root: Path,
    rows: list[dict[str, str]],
    *,
    aliases: tuple[str, ...],
    min_localities: int,
    min_endpoint_training_edges: int,
    neighbor_fraction: float,
):
    eligible = []
    statuses: list[str] = []
    for row in rows:
        candidate, status = _candidate_from_row(
            root,
            row,
            aliases=aliases,
            excluded_species=set(),
            min_localities=int(min_localities),
            min_endpoint_training_edges=int(min_endpoint_training_edges),
            neighbor_fraction=float(neighbor_fraction),
        )
        statuses.append(status)
        if candidate is not None:
            eligible.append(candidate)
    if not eligible:
        status = Counter(statuses).most_common(1)[0][0] if statuses else "no_metadata_rows"
        return None, status
    eligible.sort(key=lambda x: (-x.n_localities, -x.n_headers, x.raw_gene))
    return eligible[0], "eligible"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--protocol", type=Path, required=True)
    ap.add_argument("--execution-rule", type=Path, required=True)
    ap.add_argument("--source-sha256", required=True)
    ap.add_argument("--exclude", type=Path, action="append", default=[])
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-receipt", type=Path, required=True)
    args = ap.parse_args()

    if args.n < 1:
        raise ValueError("--n must be >= 1")
    if len(args.source_sha256) != 64:
        raise ValueError("--source-sha256 must be a 64-character SHA-256")

    root = args.root.resolve()
    if not (root / "genes.txt").is_file() or not (root / "cite.txt").is_file():
        raise FileNotFoundError("root must contain genes.txt and cite.txt")

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    execution = json.loads(args.execution_rule.read_text(encoding="utf-8"))
    aliases = tuple(protocol["marker_contract"]["normalized_aliases"])
    p1 = execution["phase1"]

    exclusions: set[str] = set()
    exclusion_receipts: list[dict[str, object]] = []
    for path in args.exclude:
        names = read_exclusion_species(path)
        exclusions |= names
        exclusion_receipts.append(
            {"path": str(path), "species": len(names), "sha256": sha256_path(path)}
        )

    grouped = metadata_rows_by_species(
        root,
        aliases=aliases,
        exclusions=exclusions,
        min_localities=int(p1["minimum_unique_localities"]),
    )
    ranked = sorted(
        grouped,
        key=lambda name: (species_hash(args.source_sha256, name), name),
    )

    selected = []
    status_counts: Counter[str] = Counter()
    checked = 0
    for species in ranked:
        checked += 1
        candidate, status = choose_candidate_for_species(
            root,
            grouped[species],
            aliases=aliases,
            min_localities=int(p1["minimum_unique_localities"]),
            min_endpoint_training_edges=int(
                p1["minimum_endpoint_disjoint_ibd_training_edges"]
            ),
            neighbor_fraction=float(p1["neighbor_fraction"]),
        )
        status_counts[status] += 1
        if candidate is None:
            continue
        selected.append((candidate, checked))
        if len(selected) == args.n:
            break

    if len(selected) < args.n:
        raise RuntimeError(
            f"only {len(selected)} geometry-eligible fresh species found; requested {args.n}"
        )

    fields = [
        "index",
        "species",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "raw_gene",
        "raw_dir",
        "n_localities",
        "n_headers",
        "graph_k",
        "edges",
        "min_endpoint_disjoint_training_edges",
        "rank_after_support",
        "rank_before_support",
    ]
    rows_out = []
    for index, (candidate, pre_rank) in enumerate(selected):
        rows_out.append(
            {
                "index": index,
                "species": candidate.species,
                "phylum": candidate.phylum,
                "class": candidate.class_name,
                "order": candidate.order,
                "family": candidate.family,
                "genus": candidate.genus,
                "raw_gene": candidate.raw_gene,
                "raw_dir": candidate.raw_dir,
                "n_localities": candidate.n_localities,
                "n_headers": candidate.n_headers,
                "graph_k": candidate.geometry.graph_k,
                "edges": candidate.geometry.n_edges,
                "min_endpoint_disjoint_training_edges": (
                    candidate.geometry.min_endpoint_disjoint_training_edges
                ),
                "rank_after_support": index + 1,
                "rank_before_support": pre_rank,
            }
        )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows_out)

    names = [str(row["species"]) for row in rows_out]
    order_counts = Counter(str(row["order"]) for row in rows_out)
    class_counts = Counter(str(row["class"]) for row in rows_out)
    family_counts = Counter(str(row["family"]) for row in rows_out)
    localities = sorted(int(row["n_localities"]) for row in rows_out)
    edges = sorted(int(row["edges"]) for row in rows_out)
    min_disjoint = sorted(
        int(row["min_endpoint_disjoint_training_edges"]) for row in rows_out
    )

    def quantile_nearest(values: list[int], q: float) -> float:
        if not values:
            return float("nan")
        pos = round(q * (len(values) - 1))
        return float(values[int(pos)])

    receipt = {
        "schema": "ttf_relational_fresh_candidate_census_v0.1",
        "status": "FRESH_GEOMETRY_ELIGIBLE_RESPONSE_BLIND",
        "source_archive_sha256": args.source_sha256,
        "ranking_rule": (
            "SHA256('relational-ttf-v0.1|<source_sha256>|<species>'), species tie-break; "
            "scan in rank order and retain only frozen Phase-1 geometry-eligible species"
        ),
        "exclusions": exclusion_receipts,
        "unique_excluded_species": len(exclusions),
        "metadata_ranked_candidates_after_exclusions": len(ranked),
        "ranked_candidates_checked": checked,
        "geometry_status_counts": dict(sorted(status_counts.items())),
        "selected_species_count": len(names),
        "selected_species_digest_sha256": species_digest(names),
        "selected_overlap_with_exclusions": len(set(names) & exclusions),
        "candidate_csv_sha256": sha256_path(args.output_csv),
        "taxonomic_summary": {
            "classes": len(class_counts),
            "orders": len(order_counts),
            "families": len(family_counts),
            "insecta": class_counts.get("Insecta", 0),
            "top_orders": dict(order_counts.most_common(12)),
            "largest_order_fraction": (
                max(order_counts.values()) / len(names) if names else None
            ),
        },
        "geometry_summary": {
            "localities_median": quantile_nearest(localities, 0.5),
            "localities_q10": quantile_nearest(localities, 0.1),
            "localities_q90": quantile_nearest(localities, 0.9),
            "edges_median": quantile_nearest(edges, 0.5),
            "min_endpoint_disjoint_median": quantile_nearest(min_disjoint, 0.5),
        },
        "response_firewall": {
            "sequence_characters_used": False,
            "pairwise_genetic_distances_opened": False,
            "T_st_computed": False,
            "ecology_genetics_association_computed": False,
        },
    }
    args.output_receipt.parent.mkdir(parents=True, exist_ok=True)
    args.output_receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "selected": len(names),
                "checked": checked,
                "species_digest": receipt["selected_species_digest_sha256"],
                "candidate_csv_sha256": receipt["candidate_csv_sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
