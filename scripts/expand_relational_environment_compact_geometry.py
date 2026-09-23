#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Expand the frozen response-blind compact Study-B genetic geometry to the edge-midpoint CSV consumed by the opportunity stage."
    )
    ap.add_argument("--compact-npz", type=Path, required=True)
    ap.add_argument("--compact-summary", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--contract", type=Path, required=True)
    ap.add_argument("--output-edges", type=Path, required=True)
    ap.add_argument("--output-receipt", type=Path, required=True)
    args = ap.parse_args()

    summary = json.loads(args.compact_summary.read_text())
    if summary.get("schema") != "ttf_relational_environment_compact_geometry_v0.2":
        raise RuntimeError("unexpected compact geometry schema")
    if summary.get("status") != "PASS_RESPONSE_BLIND_FRESH_1000_COMPACT_GEOMETRY":
        raise RuntimeError("compact geometry did not pass frozen response-blind verification")
    if any(bool(v) for v in summary["response_firewall"].values()):
        raise RuntimeError("compact geometry response firewall is open")
    if sha256_path(args.compact_npz) != summary["compact_npz_sha256"]:
        raise RuntimeError("compact geometry NPZ SHA-256 drift")
    if args.compact_npz.stat().st_size != int(summary["compact_npz_bytes"]):
        raise RuntimeError("compact geometry NPZ byte-size drift")

    contract = json.loads(args.contract.read_text())
    if contract.get("schema") != "ttf_relational_environment_geometry_reconstruction_contract_v0.1":
        raise RuntimeError("unexpected geometry reconstruction contract")
    candidate_sha = sha256_path(args.candidates)
    if candidate_sha != summary["candidate_csv_sha256"]:
        raise RuntimeError("candidate CSV drift from compact geometry receipt")
    if candidate_sha != contract["candidates"]["sha256"]:
        raise RuntimeError("candidate CSV drift from geometry reconstruction contract")
    if summary["source_archive_sha256"] != contract["source_archive"]["sha256"]:
        raise RuntimeError("source archive SHA drift between compact and reconstruction contracts")
    if int(summary["source_archive_bytes"]) != int(contract["source_archive"]["size_bytes"]):
        raise RuntimeError("source archive byte-size drift between compact and reconstruction contracts")

    with args.candidates.open(newline="", encoding="utf-8") as handle:
        candidates = list(csv.DictReader(handle))
    candidate_species = [str(row["species"]) for row in candidates]
    if len(candidate_species) != int(summary["species"]):
        raise RuntimeError("candidate species count drift")
    if len(set(candidate_species)) != len(candidate_species):
        raise RuntimeError("duplicate candidate species")

    data = np.load(args.compact_npz, allow_pickle=False)
    required = {
        "species_order",
        "edge_offsets",
        "edge_midpoints_ecef_km",
        "n_localities",
        "n_edges",
    }
    missing = required - set(data.files)
    if missing:
        raise RuntimeError(f"compact geometry missing arrays: {sorted(missing)}")

    species = np.asarray(data["species_order"]).astype(str)
    offsets = np.asarray(data["edge_offsets"], dtype=np.int64)
    midpoints = np.asarray(data["edge_midpoints_ecef_km"], dtype=np.float64)
    n_localities = np.asarray(data["n_localities"], dtype=np.int64)
    n_edges = np.asarray(data["n_edges"], dtype=np.int64)

    if list(species) != candidate_species:
        raise RuntimeError("compact species order differs from frozen candidate CSV")
    if len(offsets) != len(species) + 1 or int(offsets[0]) != 0:
        raise RuntimeError("compact edge offsets shape drift")
    if not np.array_equal(np.diff(offsets), n_edges):
        raise RuntimeError("compact edge offsets disagree with n_edges")
    if midpoints.shape != (int(offsets[-1]), 3):
        raise RuntimeError("compact edge-midpoint shape drift")
    if midpoints.dtype != np.float64 or not np.isfinite(midpoints).all():
        raise RuntimeError("compact edge midpoints are not finite float64 ECEF-km")
    if int(np.sum(n_localities)) != int(summary["locality_rows"]):
        raise RuntimeError("compact locality count drift")
    if int(np.sum(n_edges)) != int(summary["edge_rows"]):
        raise RuntimeError("compact edge count drift")
    if int(offsets[-1]) != int(summary["edge_rows"]):
        raise RuntimeError("compact final edge offset drift")

    aggregate = contract["frozen_aggregate_invariants"]
    if (
        len(species) != int(aggregate["species"])
        or int(np.sum(n_localities)) != int(aggregate["locality_rows"])
        or int(offsets[-1]) != int(aggregate["edge_rows"])
        or int(aggregate["verification_failures"]) != 0
    ):
        raise RuntimeError("compact geometry aggregate invariants drift from frozen reconstruction contract")

    args.output_edges.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "edge_index", "mid_x_km", "mid_y_km", "mid_z_km"]
    with args.output_edges.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for i, name in enumerate(species):
            left = int(offsets[i])
            right = int(offsets[i + 1])
            for edge_index, xyz in enumerate(midpoints[left:right]):
                writer.writerow(
                    {
                        "species": name,
                        "edge_index": edge_index,
                        "mid_x_km": repr(float(xyz[0])),
                        "mid_y_km": repr(float(xyz[1])),
                        "mid_z_km": repr(float(xyz[2])),
                    }
                )

    receipt = {
        "schema": "ttf_relational_environment_genetic_geometry_reconstruction_v0.2",
        "status": "PASS_RESPONSE_BLIND_FRESH_1000_GEOMETRY_RECONSTRUCTION",
        "representation": "compact_edge_midpoints",
        "contract_sha256": sha256_path(args.contract),
        "source_archive_sha256": summary["source_archive_sha256"],
        "candidate_csv_sha256": candidate_sha,
        "neighbor_fraction": float(contract["locality_geometry"]["neighbor_fraction"]),
        "species": int(len(species)),
        "locality_rows": int(np.sum(n_localities)),
        "edge_rows": int(offsets[-1]),
        "locality_csv_sha256": None,
        "edge_csv_sha256": sha256_path(args.output_edges),
        "compact_npz_sha256": sha256_path(args.compact_npz),
        "compact_summary_sha256": sha256_path(args.compact_summary),
        "verification_failures": 0,
        "verification": (
            "The compact response-blind representation reproduces the exact frozen "
            "candidate species order, per-species edge counts, aggregate locality/edge "
            "counts, and float64 ECEF-km edge midpoints used by the opportunity stage."
        ),
        "allowed_information": (
            "response-blind frozen genetic edge midpoints and candidate taxonomy/graph "
            "metadata only; no sequence identity or genetic response"
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
    args.output_receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "representation": receipt["representation"],
                "species": receipt["species"],
                "locality_rows": receipt["locality_rows"],
                "edge_rows": receipt["edge_rows"],
                "edge_csv_sha256": receipt["edge_csv_sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
