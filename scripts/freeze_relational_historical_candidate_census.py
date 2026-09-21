#!/usr/bin/env python3
"""Freeze the Study-C geometry-eligible candidate universe.

This script is response-blind. It verifies and unions the frozen S1, S2, and
Study-B design universes, excludes them before ranking, and then reads only
taxonomy, marker labels, FASTA headers, occurrence coordinates, and graph
geometry. It never reads nucleotide character identity or any empirical T_st.
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

TAG = "relational-history-candidate-v0.1"


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def species_digest(names) -> str:
    payload = "\n".join(map(str, names)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def species_hash(source_sha256: str, species: str) -> str:
    return hashlib.sha256(
        f"{TAG}|{source_sha256}|{species}".encode("utf-8")
    ).hexdigest()


def load_b_candidates(path: Path) -> list[str]:
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    if not rows or "species" not in rows[0]:
        raise RuntimeError("Study-B candidate CSV lacks species column")
    names = [collapse_whitespace(row["species"]) for row in rows]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise RuntimeError("Study-B candidate species are blank or duplicated")
    return names


def verified_prior_exclusions(
    history_rule: dict,
    *,
    s1_design_path: Path,
    s2_census_path: Path,
    b_candidates_path: Path,
) -> tuple[set[str], dict]:
    contract = history_rule["independent_species_domain"]["prior_universe_reproduction"]

    s1_rule = contract["S1_butterfly_trait"]
    if sha256_path(s1_design_path) != str(s1_rule["expected_design_sha256"]):
        raise RuntimeError("S1 design SHA drift")
    s1 = json.loads(s1_design_path.read_text(encoding="utf-8"))
    s1_names = list(map(str, s1["species_order"]))
    if len(s1_names) != int(s1_rule["expected_species"]):
        raise RuntimeError("S1 species-count drift")

    s2_rule = contract["S2_host_resource_geography"]
    s2 = json.loads(s2_census_path.read_text(encoding="utf-8"))
    if "species" not in s2:
        raise RuntimeError("S2 full census must contain the selected species rows")
    s2_names = [str(row["species"]) for row in s2["species"]]
    if len(s2_names) != int(s2_rule["expected_species"]):
        raise RuntimeError("S2 species-count drift")
    if species_digest(s2_names) != str(s2_rule["expected_selected_species_sha256"]):
        raise RuntimeError("S2 selected-species digest drift")

    b_rule = contract["B_environment_candidate_universe"]
    if sha256_path(b_candidates_path) != str(b_rule["expected_sha256"]):
        raise RuntimeError("Study-B candidate CSV SHA drift")
    b_names = load_b_candidates(b_candidates_path)
    if len(b_names) != int(b_rule["expected_species"]):
        raise RuntimeError("Study-B candidate species-count drift")

    union = set(s1_names) | set(s2_names) | set(b_names)
    return union, {
        "S1_species": len(s1_names),
        "S2_species": len(s2_names),
        "B_species": len(b_names),
        "union_species": len(union),
        "S1_digest": species_digest(s1_names),
        "S2_digest": species_digest(s2_names),
        "B_digest": species_digest(b_names),
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
    statuses = []
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
        return None, (Counter(statuses).most_common(1)[0][0] if statuses else "no_metadata_rows")
    eligible.sort(key=lambda x: (-x.n_localities, -x.n_headers, x.raw_gene))
    return eligible[0], "eligible"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--source-sha256", required=True)
    ap.add_argument("--protocol", type=Path, required=True)
    ap.add_argument("--history-rule", type=Path, required=True)
    ap.add_argument("--s1-design-json", type=Path, required=True)
    ap.add_argument("--s2-full-census-json", type=Path, required=True)
    ap.add_argument("--b-candidates", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-receipt", type=Path, required=True)
    args = ap.parse_args()

    if len(args.source_sha256) != 64:
        raise ValueError("--source-sha256 must be a 64-character SHA-256")
    root = args.root.resolve()
    if not (root / "genes.txt").is_file() or not (root / "cite.txt").is_file():
        raise FileNotFoundError("root must contain genes.txt and cite.txt")

    protocol = json.loads(args.protocol.read_text(encoding="utf-8"))
    history = json.loads(args.history_rule.read_text(encoding="utf-8"))
    if history.get("schema") != "ttf_relational_historical_climate_exposure_rule_v0.1":
        raise RuntimeError("unexpected Study-C history rule")
    if args.source_sha256 != history["independent_species_domain"]["source_archive_sha256"]:
        raise RuntimeError("source archive SHA drift")

    exclusions, exclusion_receipt = verified_prior_exclusions(
        history,
        s1_design_path=args.s1_design_json,
        s2_census_path=args.s2_full_census_json,
        b_candidates_path=args.b_candidates,
    )

    aliases = tuple(protocol["marker_contract"]["normalized_aliases"])
    geometry = history["independent_species_domain"]["geometry_eligibility"]
    grouped = metadata_rows_by_species(
        root,
        aliases=aliases,
        exclusions=exclusions,
        min_localities=int(geometry["minimum_unique_localities"]),
    )
    ranked = sorted(
        grouped,
        key=lambda name: (species_hash(args.source_sha256, name), name),
    )

    cap = int(geometry["candidate_cap"])
    minimum = int(history["independent_species_domain"]["minimum_candidate_species"])
    selected = []
    status_counts: Counter[str] = Counter()
    checked = 0
    for species in ranked:
        checked += 1
        candidate, status = choose_candidate_for_species(
            root,
            grouped[species],
            aliases=aliases,
            min_localities=int(geometry["minimum_unique_localities"]),
            min_endpoint_training_edges=int(
                geometry["minimum_endpoint_disjoint_ibd_training_edges"]
            ),
            neighbor_fraction=float(geometry["neighbor_fraction"]),
        )
        status_counts[status] += 1
        if candidate is not None:
            selected.append((candidate, checked))
            if len(selected) == cap:
                break

    if len(selected) < minimum:
        receipt = {
            "schema": "ttf_relational_historical_candidate_census_v0.1",
            "status": "NOT_EVALUABLE_HISTORICAL_CANDIDATE_GEOMETRY",
            "source_archive_sha256": args.source_sha256,
            "verified_prior_exclusions": exclusion_receipt,
            "unique_excluded_species": len(exclusions),
            "metadata_ranked_candidates_after_exclusions": len(ranked),
            "ranked_candidates_checked": checked,
            "geometry_status_counts": dict(sorted(status_counts.items())),
            "selected_species_count": len(selected),
            "minimum_required_species": minimum,
            "response_firewall": {
                "sequence_characters_used": False,
                "pairwise_genetic_distances_opened": False,
                "T_st_computed": False,
                "beta_hist_computed": False,
            },
        }
        args.output_receipt.parent.mkdir(parents=True, exist_ok=True)
        args.output_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        print(json.dumps(receipt, sort_keys=True))
        return 0

    fields = [
        "index", "species", "phylum", "class", "order", "family", "genus",
        "raw_gene", "raw_dir", "n_localities", "n_headers", "graph_k", "edges",
        "min_endpoint_disjoint_training_edges", "rank_after_support",
        "rank_before_support",
    ]
    rows_out = []
    for index, (candidate, pre_rank) in enumerate(selected):
        rows_out.append({
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
        })

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows_out)

    names = [str(row["species"]) for row in rows_out]
    receipt = {
        "schema": "ttf_relational_historical_candidate_census_v0.1",
        "status": "FRESH_HISTORICAL_GEOMETRY_ELIGIBLE_RESPONSE_BLIND",
        "source_archive_sha256": args.source_sha256,
        "ranking_rule": (
            "SHA256('relational-history-candidate-v0.1|<source_sha256>|<species>'), "
            "species tie-break; exclusions applied before ranking"
        ),
        "verified_prior_exclusions": exclusion_receipt,
        "unique_excluded_species": len(exclusions),
        "metadata_ranked_candidates_after_exclusions": len(ranked),
        "ranked_candidates_checked": checked,
        "geometry_status_counts": dict(sorted(status_counts.items())),
        "selected_species_count": len(names),
        "selected_species_digest_sha256": species_digest(names),
        "selected_overlap_with_exclusions": len(set(names) & exclusions),
        "candidate_csv_sha256": sha256_path(args.output_csv),
        "response_firewall": {
            "sequence_characters_used": False,
            "pairwise_genetic_distances_opened": False,
            "T_st_computed": False,
            "beta_hist_computed": False,
        },
    }
    args.output_receipt.parent.mkdir(parents=True, exist_ok=True)
    args.output_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": receipt["status"],
        "selected": len(names),
        "checked": checked,
        "excluded": len(exclusions),
        "species_digest": receipt["selected_species_digest_sha256"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
