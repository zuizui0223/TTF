#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

try:
    from scripts.reconstruct_relational_environment_geometry import (
        edge_rows,
        locality_rows,
        reconstruct_candidate,
        sha256_path,
        verify_candidate,
    )
except ModuleNotFoundError:
    from reconstruct_relational_environment_geometry import (
        edge_rows,
        locality_rows,
        reconstruct_candidate,
        sha256_path,
        verify_candidate,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--history-rule", type=Path, required=True)
    ap.add_argument("--neighbor-fraction", type=float, default=0.15)
    ap.add_argument("--output-localities", type=Path, required=True)
    ap.add_argument("--output-edges", type=Path, required=True)
    ap.add_argument("--output-receipt", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.history_rule.read_text())
    if rule.get("schema") != "ttf_relational_historical_climate_exposure_rule_v0.1":
        raise RuntimeError("unexpected Study-C history rule")
    geometry_rule = rule["independent_species_domain"]["geometry_eligibility"]
    if float(args.neighbor_fraction) != float(geometry_rule["neighbor_fraction"]):
        raise RuntimeError("Study-C neighbor fraction drift")

    root = args.root.resolve()
    if not (root / "genes.txt").is_file() or not (root / "cite.txt").is_file():
        raise FileNotFoundError("root must contain genes.txt and cite.txt")

    candidates = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    minimum = int(rule["independent_species_domain"]["minimum_candidate_species"])
    cap = int(geometry_rule["candidate_cap"])
    if not (minimum <= len(candidates) <= cap):
        raise RuntimeError("Study-C candidate count outside frozen 500-1000 range")
    if len({row["species"] for row in candidates}) != len(candidates):
        raise RuntimeError("duplicate Study-C candidate species")

    locality_fields = ["species", "locality_index", "latitude", "longitude", "x_km", "y_km", "z_km"]
    edge_fields = ["species", "edge_index", "node_left", "node_right", "mid_x_km", "mid_y_km", "mid_z_km"]
    args.output_localities.parent.mkdir(parents=True, exist_ok=True)
    args.output_edges.parent.mkdir(parents=True, exist_ok=True)

    n_localities = 0
    n_edges = 0
    with (
        args.output_localities.open("w", newline="", encoding="utf-8") as local_handle,
        args.output_edges.open("w", newline="", encoding="utf-8") as edge_handle,
    ):
        lw = csv.DictWriter(local_handle, fieldnames=locality_fields, lineterminator="\n")
        ew = csv.DictWriter(edge_handle, fieldnames=edge_fields, lineterminator="\n")
        lw.writeheader(); ew.writeheader()
        for row in candidates:
            species = str(row["species"])
            geometry, canonical_latlon = reconstruct_candidate(
                root, row, neighbor_fraction=float(args.neighbor_fraction)
            )
            verify_candidate(row, geometry)
            lr = list(locality_rows(species, geometry, canonical_latlon))
            er = list(edge_rows(species, geometry))
            lw.writerows(lr); ew.writerows(er)
            n_localities += len(lr); n_edges += len(er)

    receipt = {
        "schema": "ttf_relational_historical_genetic_geometry_reconstruction_v0.1",
        "status": "PASS_RESPONSE_BLIND_HISTORICAL_CANDIDATE_GEOMETRY_RECONSTRUCTION",
        "source_archive_sha256": rule["independent_species_domain"]["source_archive_sha256"],
        "candidate_csv_sha256": sha256_path(args.candidates),
        "history_rule_sha256": sha256_path(args.history_rule),
        "neighbor_fraction": float(args.neighbor_fraction),
        "species": len(candidates),
        "locality_rows": n_localities,
        "edge_rows": n_edges,
        "locality_csv_sha256": sha256_path(args.output_localities),
        "edge_csv_sha256": sha256_path(args.output_edges),
        "verification_failures": 0,
        "response_firewall": {
            "sequence_lines_decoded": False,
            "Study_C_sequence_identity_opened": False,
            "Study_C_pairwise_genetic_distances_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }
    args.output_receipt.parent.mkdir(parents=True, exist_ok=True)
    args.output_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
