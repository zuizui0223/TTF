#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.core import SpeciesEdges
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.phylogatr_confirmatory import (
    Phase1PanelCandidate,
    apply_species_cap,
    choose_one_panel_per_species,
    deterministic_species_split,
    phase1_dataset_digest,
    scan_phylogatr_phase1,
    sha256_path,
    sha256_text,
)


def species_edges(panel: Phase1PanelCandidate) -> SpeciesEdges:
    coordinates = np.asarray(panel.geometry.coordinates, dtype=float)
    nodes = np.asarray(panel.geometry.edge_nodes, dtype=np.int64)
    start = coordinates[nodes[:, 0]]
    end = coordinates[nodes[:, 1]]
    return SpeciesEdges(
        species=panel.species,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=np.zeros(len(nodes), dtype=float),
    )


def nearest_training_support(
    train_edges: dict[str, SpeciesEdges],
    eval_edges: dict[str, SpeciesEdges],
    *,
    bandwidth_km: float,
    chunk_size: int = 256,
) -> dict:
    if not train_edges or not eval_edges:
        return {
            "bandwidth_km": float(bandwidth_km),
            "eval_species": len(eval_edges),
            "status": "not_computed_empty_split",
        }
    train_midpoint = np.vstack([train_edges[name].midpoint for name in sorted(train_edges)])
    per_species: dict[str, dict] = {}
    for name in sorted(eval_edges):
        midpoint = eval_edges[name].midpoint
        nearest = np.empty(len(midpoint), dtype=float)
        for start in range(0, len(midpoint), int(chunk_size)):
            stop = min(start + int(chunk_size), len(midpoint))
            delta = midpoint[start:stop, None, :] - train_midpoint[None, :, :]
            nearest[start:stop] = np.sqrt(np.min(np.sum(delta * delta, axis=2), axis=1))
        per_species[name] = {
            "edges": int(len(nearest)),
            "median_nearest_training_edge_km": float(np.median(nearest)),
            "p90_nearest_training_edge_km": float(np.percentile(nearest, 90)),
            "fraction_within_bandwidth": float(np.mean(nearest <= float(bandwidth_km))),
        }
    medians = np.asarray(
        [row["median_nearest_training_edge_km"] for row in per_species.values()], dtype=float
    )
    fractions = np.asarray(
        [row["fraction_within_bandwidth"] for row in per_species.values()], dtype=float
    )
    return {
        "bandwidth_km": float(bandwidth_km),
        "eval_species": len(per_species),
        "species_median_nearest_training_edge_km_median": float(np.median(medians)),
        "species_median_nearest_training_edge_km_p90": float(np.percentile(medians, 90)),
        "species_fraction_within_bandwidth_median": float(np.median(fractions)),
        "species_fraction_within_bandwidth_p10": float(np.percentile(fractions, 10)),
        "species_with_zero_edges_within_bandwidth": int(np.count_nonzero(fractions == 0.0)),
        "per_eval_species": per_species,
    }


def _load_json(path: Path, schema: str) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _assert_firewall(payload: dict, key: str) -> None:
    firewall = payload[key]
    if not isinstance(firewall, dict) or not firewall:
        raise RuntimeError(f"missing frozen outcome firewall in {key}")
    if any(value is not False for value in firewall.values()):
        raise RuntimeError(f"outcome firewall is open in {key}")


def _species_list_digest(species: list[str]) -> str:
    return hashlib.sha256(("\n".join(species) + "\n").encode("utf-8")).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Freeze response-blind phase-1 geometry from a fresh phylogatR download."
    )
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--protocol", type=Path, required=True)
    ap.add_argument("--parser-rule", type=Path, required=True)
    ap.add_argument("--digest-rule", type=Path, required=True)
    ap.add_argument("--execution-rule", type=Path, required=True)
    ap.add_argument("--decker-exclusion", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-manifest", type=Path, required=True)
    args = ap.parse_args()

    protocol = _load_json(args.protocol, "ttf_genetic_phylogatr_confirmatory_protocol_v0.1")
    parser_rule = _load_json(args.parser_rule, "ttf_genetic_phylogatr_phase1_parser_rule_v0.1")
    digest_rule = _load_json(args.digest_rule, "ttf_genetic_phylogatr_phase1_digest_rule_v0.1")
    execution = _load_json(args.execution_rule, "ttf_genetic_phylogatr_phase1_execution_rule_v0.1")
    exclusion = _load_json(args.decker_exclusion, "ttf_genetic_decker_species_exclusion_v0.1")
    _assert_firewall(protocol, "outcome_firewall_at_freeze")
    _assert_firewall(parser_rule, "outcome_firewall")
    _assert_firewall(digest_rule, "firewall")
    _assert_firewall(execution, "outcome_firewall")
    _assert_firewall(exclusion, "outcome_firewall")

    root = args.root.resolve()
    genes_path = root / "genes.txt"
    cite_path = root / "cite.txt"
    if not genes_path.is_file() or not cite_path.is_file():
        raise FileNotFoundError("root must be the phylogatr-results directory containing genes.txt and cite.txt")

    exclusion_species = sorted(str(name) for name in exclusion["species"])
    if len(exclusion_species) != int(exclusion["species_count"]):
        raise RuntimeError("Decker exclusion species count drift")
    if _species_list_digest(exclusion_species) != exclusion["species_list_sha256"]:
        raise RuntimeError("Decker exclusion species-list digest drift")

    p1 = execution["phase1"]
    aliases = tuple(protocol["marker_contract"]["normalized_aliases"])
    scan = scan_phylogatr_phase1(
        root,
        aliases=aliases,
        excluded_species=exclusion_species,
        min_localities=int(p1["minimum_unique_localities"]),
        min_endpoint_training_edges=int(p1["minimum_endpoint_disjoint_ibd_training_edges"]),
        neighbor_fraction=float(p1["neighbor_fraction"]),
    )
    one_panel = choose_one_panel_per_species(scan.candidates)
    dataset_digest, digest_payload = phase1_dataset_digest(root, one_panel)
    selected = apply_species_cap(
        one_panel,
        dataset_digest=dataset_digest,
        maximum_species=int(p1["maximum_species"]),
    )
    train, evaluation = deterministic_species_split(selected, dataset_digest=dataset_digest)
    if set(train) & set(evaluation):
        raise RuntimeError("fresh confirmatory split overlap")
    if set(train) | set(evaluation) != {panel.species for panel in selected}:
        raise RuntimeError("fresh confirmatory split is not a partition")

    geometry_by_species = {
        panel.species: SpeciesGeometry(
            species=panel.species,
            coordinates=np.asarray(panel.geometry.coordinates, dtype=float),
        )
        for panel in selected
    }
    edges_by_species = {panel.species: species_edges(panel) for panel in selected}
    selected_by_species = {panel.species: panel for panel in selected}

    fieldnames = [
        "species",
        "phylum",
        "class",
        "order",
        "family",
        "raw_gene",
        "locality_index",
        "latitude",
        "longitude",
        "x_km",
        "y_km",
        "z_km",
        "graph_k",
    ]
    rows: list[dict] = []
    for species in sorted(selected_by_species):
        panel = selected_by_species[species]
        for index, ((lat, lon), (x, y, z)) in enumerate(
            zip(panel.canonical_latlon, panel.geometry.coordinates)
        ):
            rows.append(
                {
                    "species": species,
                    "phylum": panel.phylum,
                    "class": panel.class_name,
                    "order": panel.order,
                    "family": panel.family,
                    "raw_gene": panel.raw_gene,
                    "locality_index": int(index),
                    "latitude": float(lat),
                    "longitude": float(lon),
                    "x_km": float(x),
                    "y_km": float(y),
                    "z_km": float(z),
                    "graph_k": int(panel.geometry.graph_k),
                }
            )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    bandwidth = float(p1["support_census_bandwidth_km"])
    support = nearest_training_support(
        {name: edges_by_species[name] for name in train},
        {name: edges_by_species[name] for name in evaluation},
        bandwidth_km=bandwidth,
    )
    final_species_count = int(len(selected))
    minimum_panel = int(p1["minimum_panel_size"])
    panel_pass = final_species_count >= minimum_panel

    panels = {}
    for species in sorted(selected_by_species):
        panel = selected_by_species[species]
        panels[species] = {
            "taxonomy": {
                "kingdom": panel.kingdom,
                "phylum": panel.phylum,
                "class": panel.class_name,
                "order": panel.order,
                "family": panel.family,
                "genus": panel.genus,
            },
            "raw_gene": panel.raw_gene,
            "raw_dir": panel.raw_dir,
            "aligned_fasta_relative": str(panel.fasta_path.relative_to(root)),
            "occurrence_relative": str(panel.occurrence_path.relative_to(root)),
            "aligned_header_sha256": panel.fasta_header_sha256,
            "occurrence_sha256": panel.occurrence_sha256,
            "aligned_header_count": panel.n_headers,
            "matched_header_count": int(panel.matched_header_count),
            "unique_localities": panel.n_localities,
            "graph_k": int(panel.geometry.graph_k),
            "edge_count": int(panel.geometry.n_edges),
            "min_endpoint_disjoint_training_edges": int(
                panel.geometry.min_endpoint_disjoint_training_edges
            ),
            "cap_key_sha256": sha256_text(f"{dataset_digest}|{species}"),
            "split_key_sha256": sha256_text(f"{dataset_digest}|split|{species}"),
        }

    manifest = {
        "schema": "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1",
        "status": (
            "FROZEN_RESPONSE_BLIND_PHASE1_GEOMETRY"
            if panel_pass
            else "NOT_EVALUABLE_PHASE1_PANEL_TOO_SMALL"
        ),
        "dataset_digest_sha256": dataset_digest,
        "digest_payload_sha256": hashlib.sha256(
            json.dumps(
                digest_payload,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest(),
        "provenance": {
            "cite_sha256": sha256_path(cite_path),
            "genes_sha256": sha256_path(genes_path),
            "protocol_sha256": sha256_path(args.protocol),
            "parser_rule_sha256": sha256_path(args.parser_rule),
            "digest_rule_sha256": sha256_path(args.digest_rule),
            "execution_rule_sha256": sha256_path(args.execution_rule),
            "decker_exclusion_manifest_sha256": sha256_path(args.decker_exclusion),
            "decker_exclusion_species_list_sha256": exclusion["species_list_sha256"],
        },
        "response_blind": {
            "sequence_characters_used": False,
            "sequence_characters_hashed": False,
            "pairwise_genetic_distances_opened": False,
            "ttf_statistic_opened": False,
            "decker_empirical_genetic_outcomes_opened": False,
            "allowed_fasta_information": "header strings only",
        },
        "rules": {
            "minimum_unique_localities": int(p1["minimum_unique_localities"]),
            "neighbor_fraction": float(p1["neighbor_fraction"]),
            "minimum_endpoint_disjoint_ibd_training_edges": int(
                p1["minimum_endpoint_disjoint_ibd_training_edges"]
            ),
            "maximum_species": int(p1["maximum_species"]),
            "minimum_panel_size": minimum_panel,
            "support_census_bandwidth_km": bandwidth,
            "support_census_is_descriptive_only": True,
        },
        "census": {
            "genes_rows": int(scan.genes_rows),
            "row_status_counts": scan.status_counts,
            "eligible_panel_candidates": int(len(scan.candidates)),
            "eligible_species_before_cap": int(len(one_panel)),
            "panel_cap_applied": bool(len(one_panel) > int(p1["maximum_species"])),
            "final_species": final_species_count,
            "minimum_panel_pass": bool(panel_pass),
        },
        "geometry_csv_sha256": sha256_path(args.output_csv),
        "geometry_fingerprint_sha256": geometry_fingerprint(
            [geometry_by_species[name] for name in sorted(geometry_by_species)]
        ) if geometry_by_species else None,
        "localities_total": int(sum(panel.n_localities for panel in selected)),
        "edges_total": int(sum(panel.geometry.n_edges for panel in selected)),
        "split": {
            "train_species": list(train),
            "eval_species": list(evaluation),
            "train_count": len(train),
            "eval_count": len(evaluation),
        },
        "support_at_500km_descriptive_only": support,
        "selected_panels": panels,
        "next_gate": (
            "phase_2 canonical-valid/noncanonical character-mask admissibility may be opened"
            if panel_pass
            else "STOP; do not open character masks or nucleotide identities"
        ),
        "qualification_claim_made": False,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
    }
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "dataset_digest_sha256": dataset_digest,
                "geometry_fingerprint_sha256": manifest["geometry_fingerprint_sha256"],
                "species": final_species_count,
                "train": len(train),
                "eval": len(evaluation),
                "localities": manifest["localities_total"],
                "edges": manifest["edges_total"],
                "support": {
                    key: value for key, value in support.items() if key != "per_eval_species"
                },
                "sequence_characters_used": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
