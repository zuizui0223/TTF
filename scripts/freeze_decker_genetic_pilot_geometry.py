#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from census_decker_genetic_geometry import (
    fasta_ids,
    latlon_to_ecef_km,
    occurrence_coordinates,
)
from ttf.calibration import seed_for
from ttf.core import SpeciesEdges
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.geometry import SpeciesGeometry, geometry_fingerprint


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def choose_one_panel_per_species(census: dict) -> list[dict]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in census["panels_ledger"]:
        if row.get("precalibration_eligible"):
            grouped[str(row["species"])].append(row)
    selected: list[dict] = []
    for species, rows in grouped.items():
        rows = sorted(
            rows,
            key=lambda row: (
                -int(row["unique_localities"]),
                -int(row["matched_records"]),
                str(row["locus"]),
            ),
        )
        selected.append(rows[0])
    return sorted(selected, key=lambda row: str(row["species"]))


def species_edges(species: str, coordinates: np.ndarray, nodes: np.ndarray) -> SpeciesEdges:
    start = coordinates[nodes[:, 0]]
    end = coordinates[nodes[:, 1]]
    return SpeciesEdges(
        species=species,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=np.zeros(len(nodes), dtype=float),
    )


def class_stratified_split(selected: list[dict], master_seed: int) -> tuple[tuple[str, ...], tuple[str, ...]]:
    train: list[str] = []
    evaluation: list[str] = []
    for class_name in sorted({str(row["class"]) for row in selected}):
        labels = sorted(str(row["species"]) for row in selected if row["class"] == class_name)
        rng = np.random.default_rng(seed_for(master_seed, "genetic_decker_split", class_name))
        permutation = rng.permutation(len(labels))
        n_eval = int(round(0.5 * len(labels)))
        evaluation.extend(labels[int(i)] for i in permutation[:n_eval])
        train.extend(labels[int(i)] for i in permutation[n_eval:])
    return tuple(sorted(train)), tuple(sorted(evaluation))


def nearest_training_support(
    train_edges: dict[str, SpeciesEdges],
    eval_edges: dict[str, SpeciesEdges],
    *,
    bandwidth_km: float,
    chunk_size: int = 256,
) -> dict:
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--rule", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path, required=True)
    args = parser.parse_args()

    census = json.loads(args.census.read_text())
    rule = json.loads(args.rule.read_text())
    if census.get("schema") != "ttf_genetic_decker_geometry_census_v0.1":
        raise RuntimeError("unexpected census schema")
    if rule.get("schema") != "ttf_genetic_decker_pilot_rule_v0.1":
        raise RuntimeError("unexpected pilot rule schema")
    firewall = rule["outcome_firewall"]
    if any(bool(value) for value in firewall.values()):
        raise RuntimeError("outcome firewall is not closed")

    selected = choose_one_panel_per_species(census)
    expected = int(rule["one_panel_per_species"]["expected_species"])
    if len(selected) != expected:
        raise RuntimeError(f"selected species drift: {len(selected)} != {expected}")

    geometry_by_species = {}
    metadata_by_species: dict[str, dict] = {}
    edge_by_species: dict[str, SpeciesEdges] = {}
    csv_rows: list[dict] = []
    fraction = float(rule["geometry_filter"]["neighbor_fraction"])

    for panel in selected:
        species = str(panel["species"])
        fasta = args.root / str(panel["fasta"]).split("external/decker/", 1)[-1]
        occurrence = args.root / str(panel["occurrence"]).split("external/decker/", 1)[-1]
        ids = set(fasta_ids(fasta))
        latlon, _ = occurrence_coordinates(occurrence, ids)
        unique_latlon = np.unique(latlon, axis=0)
        ecef = latlon_to_ecef_km(unique_latlon)
        geometry = prepare_density_scaled_genetic_geometry(
            ecef,
            neighbor_fraction=fraction,
        )
        if geometry.n_localities != int(panel["unique_localities"]):
            raise RuntimeError(f"locality-count drift for {species}")
        if geometry.graph_k != int(panel["graph_k"]):
            raise RuntimeError(f"graph-k drift for {species}")
        if geometry.n_edges != int(panel["edge_count"]):
            raise RuntimeError(f"edge-count drift for {species}")

        geometry_by_species[species] = SpeciesGeometry(species=species, coordinates=ecef)
        edge_by_species[species] = species_edges(species, ecef, geometry.edge_nodes)
        metadata_by_species[species] = {
            "class": str(panel["class"]),
            "locus": str(panel["locus"]),
            "fasta_sequences": int(panel["fasta_sequences"]),
            "matched_records": int(panel["matched_records"]),
            "unique_localities": int(panel["unique_localities"]),
            "graph_k": int(panel["graph_k"]),
            "edge_count": int(panel["edge_count"]),
            "min_endpoint_disjoint_training_edges": int(
                panel["min_endpoint_disjoint_training_edges"]
            ),
            "source_fasta": str(panel["fasta"]),
            "source_occurrence": str(panel["occurrence"]),
        }
        for index, ((lat, lon), (x, y, z)) in enumerate(zip(unique_latlon, ecef)):
            csv_rows.append(
                {
                    "species": species,
                    "class": str(panel["class"]),
                    "locus": str(panel["locus"]),
                    "locality_index": int(index),
                    "latitude": float(lat),
                    "longitude": float(lon),
                    "x_km": float(x),
                    "y_km": float(y),
                    "z_km": float(z),
                    "graph_k": int(panel["graph_k"]),
                }
            )

    train, evaluation = class_stratified_split(selected, int(rule["split"]["master_seed"]))
    if len(train) != int(rule["split"]["expected_train_species"]):
        raise RuntimeError("train split size drift")
    if len(evaluation) != int(rule["split"]["expected_eval_species"]):
        raise RuntimeError("evaluation split size drift")
    if set(train) & set(evaluation) or set(train) | set(evaluation) != set(geometry_by_species):
        raise RuntimeError("species split is not a partition")

    fieldnames = [
        "species", "class", "locus", "locality_index", "latitude", "longitude",
        "x_km", "y_km", "z_km", "graph_k",
    ]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)

    bandwidth = float(rule["core_v011_inheritance"]["bandwidth_km"])
    support = nearest_training_support(
        {name: edge_by_species[name] for name in train},
        {name: edge_by_species[name] for name in evaluation},
        bandwidth_km=bandwidth,
    )
    selected_edges = [metadata_by_species[name]["edge_count"] for name in sorted(metadata_by_species)]
    selected_localities = [
        metadata_by_species[name]["unique_localities"] for name in sorted(metadata_by_species)
    ]
    selected_edge_medians = [
        float(np.median(edge_by_species[name].length)) for name in sorted(edge_by_species)
    ]
    manifest = {
        "schema": "ttf_genetic_decker_pilot_geometry_v0.1",
        "status": "frozen_response_blind_development_geometry",
        "source_repository": rule["source"]["repository"],
        "source_commit": rule["source"]["commit"],
        "rule_sha256": sha256_path(args.rule),
        "census_sha256": sha256_path(args.census),
        "geometry_csv_sha256": sha256_path(args.output_csv),
        "geometry_fingerprint_sha256": geometry_fingerprint(
            [geometry_by_species[name] for name in sorted(geometry_by_species)]
        ),
        "genetic_outcomes_opened": False,
        "monmonier_results_opened": False,
        "species_count": len(geometry_by_species),
        "records_localities": int(sum(selected_localities)),
        "edge_count": int(sum(selected_edges)),
        "class_counts": {
            class_name: int(sum(meta["class"] == class_name for meta in metadata_by_species.values()))
            for class_name in sorted({meta["class"] for meta in metadata_by_species.values()})
        },
        "locality_count": {
            "min": int(min(selected_localities)),
            "median": float(np.median(selected_localities)),
            "p90": float(np.percentile(selected_localities, 90)),
            "max": int(max(selected_localities)),
        },
        "edge_count_per_species": {
            "min": int(min(selected_edges)),
            "median": float(np.median(selected_edges)),
            "p90": float(np.percentile(selected_edges, 90)),
            "max": int(max(selected_edges)),
        },
        "species_median_edge_length_km": {
            "median": float(np.median(selected_edge_medians)),
            "p10": float(np.percentile(selected_edge_medians, 10)),
            "p90": float(np.percentile(selected_edge_medians, 90)),
        },
        "split": {
            "train_species": list(train),
            "eval_species": list(evaluation),
            "train_count": len(train),
            "eval_count": len(evaluation),
        },
        "support_at_inherited_bandwidth": support,
        "selected_panels": metadata_by_species,
    }
    args.output_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "geometry_fingerprint_sha256": manifest["geometry_fingerprint_sha256"],
                "geometry_csv_sha256": manifest["geometry_csv_sha256"],
                "species_count": manifest["species_count"],
                "records_localities": manifest["records_localities"],
                "edge_count": manifest["edge_count"],
                "support": {
                    key: value
                    for key, value in support.items()
                    if key != "per_eval_species"
                },
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
