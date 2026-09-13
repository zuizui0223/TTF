#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry

SOURCE_COMMIT = "11534ba70fbe3705e667651030edef5b849fa949"
FASTA_SUFFIX = "_missingfiltered_finalfiltered_final.fasta"
EARTH_RADIUS_KM = 6371.0088


def fasta_ids(path: Path) -> tuple[str, ...]:
    out: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                label = line[1:].strip().split()[0]
                if label:
                    out.append(label)
    return tuple(out)


def occurrence_coordinates(path: Path, allowed_ids: set[str]) -> tuple[np.ndarray, int]:
    rows: list[tuple[float, float]] = []
    total = 0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"phylogatr_id", "latitude", "longitude"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"occurrence schema missing required columns: {path}")
        for row in reader:
            total += 1
            if row["phylogatr_id"] not in allowed_ids:
                continue
            try:
                lat = float(row["latitude"])
                lon = float(row["longitude"])
            except (TypeError, ValueError):
                continue
            if np.isfinite(lat) and np.isfinite(lon):
                rows.append((lat, lon))
    if not rows:
        return np.empty((0, 2), dtype=float), total
    return np.asarray(rows, dtype=float), total


def latlon_to_ecef_km(latlon: np.ndarray) -> np.ndarray:
    x = np.asarray(latlon, dtype=float)
    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError("latlon must be n x 2")
    lat = np.deg2rad(x[:, 0])
    lon = np.deg2rad(x[:, 1])
    coslat = np.cos(lat)
    return EARTH_RADIUS_KM * np.column_stack(
        [coslat * np.cos(lon), coslat * np.sin(lon), np.sin(lat)]
    )


def percentile_or_none(values: Iterable[int | float], q: float) -> float | None:
    x = np.asarray(list(values), dtype=float)
    if len(x) == 0:
        return None
    return float(np.percentile(x, q))


def panel_record(
    class_name: str,
    fasta_path: Path,
    *,
    min_localities: int,
    min_endpoint_training_edges: int,
    neighbor_fraction: float,
) -> dict:
    stem = fasta_path.name[: -len(FASTA_SUFFIX)]
    if "-" not in stem:
        raise ValueError(f"cannot parse species/locus from {fasta_path.name}")
    species_key, locus = stem.rsplit("-", 1)
    occurrence_path = fasta_path.parent / f"{species_key}_occurfiltered.csv"
    ids = fasta_ids(fasta_path)

    base = {
        "class": class_name,
        "species": species_key.replace("-", " "),
        "species_key": species_key,
        "locus": locus,
        "fasta": str(fasta_path),
        "occurrence": str(occurrence_path),
        "fasta_sequences": len(ids),
        "unique_fasta_ids": len(set(ids)),
        "occurrence_file_exists": occurrence_path.exists(),
        "genetic_outcome_read": False,
        "sequence_characters_used": False,
    }
    if not occurrence_path.exists():
        return {
            **base,
            "occurrence_rows_total": None,
            "matched_records": 0,
            "matched_fraction_of_fasta": 0.0,
            "unique_localities": 0,
            "graph_k": None,
            "edge_count": 0,
            "min_endpoint_disjoint_training_edges": 0,
            "median_edge_length_km": None,
            "passes_min_localities": False,
            "passes_endpoint_safe_ibd": False,
            "precalibration_eligible": False,
        }

    latlon, occurrence_total = occurrence_coordinates(occurrence_path, set(ids))
    matched_ids = len(latlon)
    unique_latlon = np.unique(latlon, axis=0) if matched_ids else np.empty((0, 2))
    n_localities = int(len(unique_latlon))
    passes_localities = n_localities >= int(min_localities)

    graph_k = None
    edge_count = 0
    min_disjoint = 0
    median_edge = None
    if n_localities >= 3:
        ecef = latlon_to_ecef_km(latlon)
        geometry = prepare_density_scaled_genetic_geometry(
            ecef,
            neighbor_fraction=float(neighbor_fraction),
        )
        graph_k = int(geometry.graph_k)
        edge_count = int(geometry.n_edges)
        min_disjoint = int(geometry.min_endpoint_disjoint_training_edges)
        nodes = geometry.edge_nodes
        edge_length = np.linalg.norm(
            geometry.coordinates[nodes[:, 0]] - geometry.coordinates[nodes[:, 1]],
            axis=1,
        )
        median_edge = float(np.median(edge_length))

    endpoint_pass = min_disjoint >= int(min_endpoint_training_edges)
    matched_fraction = float(matched_ids / len(ids)) if len(ids) else 0.0
    return {
        **base,
        "occurrence_rows_total": int(occurrence_total),
        "matched_records": int(matched_ids),
        "matched_fraction_of_fasta": matched_fraction,
        "unique_localities": n_localities,
        "graph_k": graph_k,
        "edge_count": edge_count,
        "min_endpoint_disjoint_training_edges": min_disjoint,
        "median_edge_length_km": median_edge,
        "passes_min_localities": bool(passes_localities),
        "passes_endpoint_safe_ibd": bool(endpoint_pass),
        "precalibration_eligible": bool(passes_localities and endpoint_pass),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Response-blind species-locus geometry census for Decker phylogeography data."
    )
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-localities", type=int, default=12)
    parser.add_argument("--min-endpoint-training-edges", type=int, default=5)
    parser.add_argument("--neighbor-fraction", type=float, default=0.15)
    args = parser.parse_args()

    if args.min_localities < 3:
        raise ValueError("min-localities must be >= 3")
    if args.min_endpoint_training_edges < 3:
        raise ValueError("min-endpoint-training-edges must be >= 3")

    panels: list[dict] = []
    for class_name, directory in (
        ("Aves", args.root / "Aves_AlignOccur"),
        ("Chiroptera", args.root / "Chiroptera_AlignOccur"),
    ):
        if not directory.exists():
            raise FileNotFoundError(directory)
        for fasta_path in sorted(directory.glob(f"*{FASTA_SUFFIX}")):
            panels.append(
                panel_record(
                    class_name,
                    fasta_path,
                    min_localities=args.min_localities,
                    min_endpoint_training_edges=args.min_endpoint_training_edges,
                    neighbor_fraction=args.neighbor_fraction,
                )
            )

    by_class: dict[str, dict] = {}
    for class_name in ("Aves", "Chiroptera"):
        rows = [p for p in panels if p["class"] == class_name]
        eligible = [p for p in rows if p["precalibration_eligible"]]
        by_class[class_name] = {
            "panels": len(rows),
            "species": len({p["species"] for p in rows}),
            "occurrence_file_missing": sum(not p["occurrence_file_exists"] for p in rows),
            "passes_min_localities": sum(p["passes_min_localities"] for p in rows),
            "passes_endpoint_safe_ibd": sum(p["passes_endpoint_safe_ibd"] for p in rows),
            "precalibration_eligible": len(eligible),
            "localities_p10": percentile_or_none((p["unique_localities"] for p in rows), 10),
            "localities_median": percentile_or_none((p["unique_localities"] for p in rows), 50),
            "localities_p90": percentile_or_none((p["unique_localities"] for p in rows), 90),
        }

    eligible = [p for p in panels if p["precalibration_eligible"]]
    output = {
        "schema": "ttf_genetic_decker_geometry_census_v0.1",
        "status": "response_blind_geometry_census_only",
        "source_repository": "skdecker/PhylogeographicBreaks",
        "source_commit": SOURCE_COMMIT,
        "genetic_outcomes_opened": False,
        "sequence_characters_used": False,
        "monmonier_results_opened": False,
        "rules": {
            "coordinate_match": "phylogatr_id intersection of final FASTA headers and occurfiltered rows",
            "locality_collapse": "exact latitude-longitude pair only",
            "coordinate_geometry": "spherical ECEF km, radius 6371.0088",
            "neighbor_fraction": float(args.neighbor_fraction),
            "graph_rule": "density_scaled_k_after_exact_locality_collapse",
            "coarse_min_localities": int(args.min_localities),
            "min_endpoint_disjoint_training_edges": int(args.min_endpoint_training_edges),
            "precalibration_eligible": "coarse_min_localities AND endpoint_safe_ibd",
            "eligibility_is_not_final_qualification": True,
        },
        "summary": {
            "panels": len(panels),
            "species": len({p["species"] for p in panels}),
            "precalibration_eligible_panels": len(eligible),
            "precalibration_eligible_species": len({p["species"] for p in eligible}),
            "by_class": by_class,
        },
        "eligible_panels": [
            {
                "class": p["class"],
                "species": p["species"],
                "locus": p["locus"],
                "fasta_sequences": p["fasta_sequences"],
                "matched_records": p["matched_records"],
                "unique_localities": p["unique_localities"],
                "graph_k": p["graph_k"],
                "edge_count": p["edge_count"],
                "min_endpoint_disjoint_training_edges": p[
                    "min_endpoint_disjoint_training_edges"
                ],
                "median_edge_length_km": p["median_edge_length_km"],
            }
            for p in eligible
        ],
        "panels_ledger": panels,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
