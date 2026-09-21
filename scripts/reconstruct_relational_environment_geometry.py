#!/usr/bin/env python3
"""Reconstruct frozen fresh-1000 response-blind genetic geometry.

This script is intentionally restricted to aligned FASTA headers, occurrence
coordinates, the frozen candidate table, and the inherited graph rule. It never
decodes sequence lines or computes nucleotide identity, genetic distance, T_st,
or any ecology-genetics association.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.phylogatr_confirmatory import (
    canonical_latlon_for_geometry,
    coordinates_for_headers,
    fasta_headers_only,
    latlon_to_ecef_km,
    read_occurrence_rows,
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_directory(root: Path, raw_dir: str) -> Path:
    rel = Path(str(raw_dir))
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(f"unsafe candidate raw_dir: {raw_dir!r}")
    base = root.resolve()
    resolved = (base / rel).resolve()
    if resolved != base and base not in resolved.parents:
        raise ValueError(f"candidate raw_dir escapes root: {raw_dir!r}")
    return resolved


def reconstruct_candidate(root: Path, row: dict[str, str], *, neighbor_fraction: float):
    species = str(row["species"]).strip()
    raw_gene = str(row["raw_gene"])
    directory = safe_directory(root, str(row["raw_dir"]))
    fasta = directory / f"{raw_gene}.afa"
    occurrences = directory / "occurrences.txt"
    if not fasta.is_file():
        raise FileNotFoundError(f"{species}: missing aligned FASTA {fasta}")
    if not occurrences.is_file():
        raise FileNotFoundError(f"{species}: missing occurrences {occurrences}")

    # fasta_headers_only() reads the file in binary mode and decodes only ">"
    # header lines. Sequence bytes are discarded without interpretation.
    headers = fasta_headers_only(fasta)
    latlon = coordinates_for_headers(headers, read_occurrence_rows(occurrences))
    if len(latlon) == 0:
        raise RuntimeError(f"{species}: no header-occurrence coordinate matches")
    exact_unique = np.unique(latlon, axis=0)
    ecef = latlon_to_ecef_km(exact_unique)
    geometry = prepare_density_scaled_genetic_geometry(
        ecef, neighbor_fraction=float(neighbor_fraction)
    )
    canonical_latlon = canonical_latlon_for_geometry(exact_unique, geometry)
    return geometry, canonical_latlon


def verify_candidate(row: dict[str, str], geometry) -> None:
    expected = {
        "n_localities": int(row["n_localities"]),
        "graph_k": int(row["graph_k"]),
        "edges": int(row["edges"]),
        "min_endpoint_disjoint_training_edges": int(
            row["min_endpoint_disjoint_training_edges"]
        ),
    }
    observed = {
        "n_localities": int(geometry.n_localities),
        "graph_k": int(geometry.graph_k),
        "edges": int(geometry.n_edges),
        "min_endpoint_disjoint_training_edges": int(
            geometry.min_endpoint_disjoint_training_edges
        ),
    }
    if expected != observed:
        raise RuntimeError(
            f"{row['species']}: frozen geometry drift: expected={expected}, observed={observed}"
        )


def locality_rows(species: str, geometry, canonical_latlon: np.ndarray):
    xyz = np.asarray(geometry.coordinates, dtype=float)
    ll = np.asarray(canonical_latlon, dtype=float)
    if len(xyz) != len(ll):
        raise RuntimeError(f"{species}: locality coordinate alignment drift")
    for index, (latlon, point) in enumerate(zip(ll, xyz)):
        yield {
            "species": species,
            "locality_index": index,
            "latitude": repr(float(latlon[0])),
            "longitude": repr(float(latlon[1])),
            "x_km": repr(float(point[0])),
            "y_km": repr(float(point[1])),
            "z_km": repr(float(point[2])),
        }


def edge_rows(species: str, geometry):
    xyz = np.asarray(geometry.coordinates, dtype=float)
    nodes = np.asarray(geometry.edge_nodes, dtype=np.int64)
    for index, (left, right) in enumerate(nodes):
        midpoint = 0.5 * (xyz[int(left)] + xyz[int(right)])
        yield {
            "species": species,
            "edge_index": index,
            "node_left": int(left),
            "node_right": int(right),
            "mid_x_km": repr(float(midpoint[0])),
            "mid_y_km": repr(float(midpoint[1])),
            "mid_z_km": repr(float(midpoint[2])),
        }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--source-archive-sha256", required=True)
    ap.add_argument("--contract", type=Path, required=True)
    ap.add_argument("--neighbor-fraction", type=float, default=0.15)
    ap.add_argument("--output-localities", type=Path, required=True)
    ap.add_argument("--output-edges", type=Path, required=True)
    ap.add_argument("--output-receipt", type=Path, required=True)
    args = ap.parse_args()

    root = args.root.resolve()
    if not (root / "genes.txt").is_file() or not (root / "cite.txt").is_file():
        raise FileNotFoundError("root must contain genes.txt and cite.txt")
    if len(args.source_archive_sha256) != 64:
        raise ValueError("source archive SHA-256 must have 64 hex characters")
    if not 0 < float(args.neighbor_fraction) <= 1:
        raise ValueError("neighbor fraction must be in (0,1]")

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    if contract.get("schema") != "ttf_relational_environment_geometry_reconstruction_contract_v0.1":
        raise RuntimeError("unexpected geometry reconstruction contract")
    if contract["source_archive"]["sha256"] != args.source_archive_sha256:
        raise RuntimeError("source archive SHA drift from reconstruction contract")
    if not np.isclose(
        float(contract["locality_geometry"]["neighbor_fraction"]),
        float(args.neighbor_fraction),
        rtol=0.0,
        atol=0.0,
    ):
        raise RuntimeError("neighbor fraction drift from reconstruction contract")
    candidate_sha = sha256_path(args.candidates)
    if candidate_sha != contract["candidates"]["sha256"]:
        raise RuntimeError("candidate CSV SHA drift from reconstruction contract")

    with args.candidates.open(newline="", encoding="utf-8") as handle:
        candidates = list(csv.DictReader(handle))
    expected_species = int(contract["candidates"]["species"])
    if len(candidates) != expected_species:
        raise RuntimeError(
            f"expected exact {expected_species} candidates, found {len(candidates)}"
        )
    if len({row["species"] for row in candidates}) != len(candidates):
        raise RuntimeError("duplicate candidate species")

    locality_fields = [
        "species", "locality_index", "latitude", "longitude", "x_km", "y_km", "z_km"
    ]
    edge_fields = [
        "species", "edge_index", "node_left", "node_right",
        "mid_x_km", "mid_y_km", "mid_z_km"
    ]
    args.output_localities.parent.mkdir(parents=True, exist_ok=True)
    args.output_edges.parent.mkdir(parents=True, exist_ok=True)

    species_count = 0
    locality_count = 0
    edge_count = 0
    with (
        args.output_localities.open("w", newline="", encoding="utf-8") as local_handle,
        args.output_edges.open("w", newline="", encoding="utf-8") as edge_handle,
    ):
        lw = csv.DictWriter(local_handle, fieldnames=locality_fields, lineterminator="\n")
        ew = csv.DictWriter(edge_handle, fieldnames=edge_fields, lineterminator="\n")
        lw.writeheader()
        ew.writeheader()
        for row in candidates:
            species = str(row["species"])
            geometry, canonical_latlon = reconstruct_candidate(
                root, row, neighbor_fraction=float(args.neighbor_fraction)
            )
            verify_candidate(row, geometry)
            lr = list(locality_rows(species, geometry, canonical_latlon))
            er = list(edge_rows(species, geometry))
            lw.writerows(lr)
            ew.writerows(er)
            species_count += 1
            locality_count += len(lr)
            edge_count += len(er)

    aggregate = contract["frozen_aggregate_invariants"]
    observed_aggregate = {
        "species": species_count,
        "locality_rows": locality_count,
        "edge_rows": edge_count,
        "verification_failures": 0,
    }
    expected_aggregate = {
        key: int(aggregate[key])
        for key in ("species", "locality_rows", "edge_rows", "verification_failures")
    }
    if observed_aggregate != expected_aggregate:
        raise RuntimeError(
            f"aggregate frozen geometry drift: expected={expected_aggregate}, "
            f"observed={observed_aggregate}"
        )

    receipt = {
        "schema": "ttf_relational_environment_genetic_geometry_reconstruction_v0.2",
        "status": "PASS_RESPONSE_BLIND_FRESH_1000_GEOMETRY_RECONSTRUCTION",
        "contract_sha256": sha256_path(args.contract),
        "source_archive_sha256": args.source_archive_sha256,
        "candidate_csv_sha256": candidate_sha,
        "neighbor_fraction": float(args.neighbor_fraction),
        "species": species_count,
        "locality_rows": locality_count,
        "edge_rows": edge_count,
        "locality_csv_sha256": sha256_path(args.output_localities),
        "edge_csv_sha256": sha256_path(args.output_edges),
        "verification_failures": 0,
        "verification": (
            "Every species exactly reproduced the frozen candidate n_localities, "
            "graph_k, edge count, and minimum endpoint-disjoint training-edge count."
        ),
        "allowed_information": (
            "aligned FASTA headers + occurrence coordinates + frozen candidate "
            "taxonomy/graph rule only"
        ),
        "response_firewall": {
            "sequence_lines_decoded": False,
            "sequence_identity_opened": False,
            "pairwise_genetic_distances_opened": False,
            "T_st_computed": False,
            "beta_R_computed": False,
        },
    }
    args.output_receipt.parent.mkdir(parents=True, exist_ok=True)
    args.output_receipt.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
