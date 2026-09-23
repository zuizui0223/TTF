from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np


SCRIPT = Path("scripts/expand_relational_environment_compact_geometry.py")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fixture(tmp_path: Path, *, bad_edges: bool = False):
    candidates = tmp_path / "candidates.csv"
    with candidates.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["species", "n_localities", "edges"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"species": "sp_a", "n_localities": 3, "edges": 2})
        writer.writerow({"species": "sp_b", "n_localities": 4, "edges": 2 if bad_edges else 1})

    compact = tmp_path / "compact.npz"
    np.savez_compressed(
        compact,
        species_order=np.asarray(["sp_a", "sp_b"]),
        edge_offsets=np.asarray([0, 2, 3], dtype=np.int64),
        edge_midpoints_ecef_km=np.asarray(
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            dtype=np.float64,
        ),
        n_localities=np.asarray([3, 4], dtype=np.int64),
        n_edges=np.asarray([2, 1], dtype=np.int64),
    )

    source_sha = "a" * 64
    contract = tmp_path / "contract.json"
    contract.write_text(
        json.dumps(
            {
                "schema": "ttf_relational_environment_geometry_reconstruction_contract_v0.1",
                "source_archive": {"sha256": source_sha, "size_bytes": 123},
                "candidates": {"sha256": sha256(candidates), "species": 2},
                "locality_geometry": {"neighbor_fraction": 0.15},
                "frozen_aggregate_invariants": {
                    "species": 2,
                    "locality_rows": 7,
                    "edge_rows": 3,
                    "verification_failures": 0,
                },
            }
        )
        + "\n"
    )
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "schema": "ttf_relational_environment_compact_geometry_v0.2",
                "status": "PASS_RESPONSE_BLIND_FRESH_1000_COMPACT_GEOMETRY",
                "candidate_csv_sha256": sha256(candidates),
                "compact_npz_bytes": compact.stat().st_size,
                "compact_npz_sha256": sha256(compact),
                "edge_rows": 3,
                "locality_rows": 7,
                "source_archive_bytes": 123,
                "source_archive_sha256": source_sha,
                "species": 2,
                "verification_failures": 0,
                "response_firewall": {
                    "T_st_computed": False,
                    "beta_R_computed": False,
                    "pairwise_genetic_distances_opened": False,
                    "sequence_identity_opened": False,
                    "sequence_lines_decoded": False,
                },
            }
        )
        + "\n"
    )
    return candidates, compact, contract, summary


def run_fixture(tmp_path: Path, *, bad_edges: bool = False):
    candidates, compact, contract, summary = write_fixture(tmp_path, bad_edges=bad_edges)
    edges = tmp_path / "edges.csv"
    receipt = tmp_path / "receipt.json"
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--compact-npz",
            str(compact),
            "--compact-summary",
            str(summary),
            "--candidates",
            str(candidates),
            "--contract",
            str(contract),
            "--output-edges",
            str(edges),
            "--output-receipt",
            str(receipt),
        ],
        capture_output=True,
        text=True,
    )
    return done, edges, receipt


def test_compact_geometry_expands_exact_species_order_and_midpoints(tmp_path: Path) -> None:
    done, edges, receipt = run_fixture(tmp_path)
    assert done.returncode == 0, done.stderr

    with edges.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["species"] for row in rows] == ["sp_a", "sp_a", "sp_b"]
    assert [int(row["edge_index"]) for row in rows] == [0, 1, 0]
    assert [
        [float(row["mid_x_km"]), float(row["mid_y_km"]), float(row["mid_z_km"])]
        for row in rows
    ] == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]

    payload = json.loads(receipt.read_text())
    assert payload["schema"] == "ttf_relational_environment_genetic_geometry_reconstruction_v0.2"
    assert payload["status"] == "PASS_RESPONSE_BLIND_FRESH_1000_GEOMETRY_RECONSTRUCTION"
    assert payload["representation"] == "compact_edge_midpoints"
    assert payload["species"] == 2
    assert payload["locality_rows"] == 7
    assert payload["edge_rows"] == 3
    assert payload["locality_csv_sha256"] is None
    assert payload["edge_csv_sha256"] == sha256(edges)
    assert all(value is False for value in payload["response_firewall"].values())


def test_compact_geometry_rejects_per_species_edge_count_drift(tmp_path: Path) -> None:
    done, edges, receipt = run_fixture(tmp_path, bad_edges=True)
    assert done.returncode != 0
    assert "per-species edge counts drift" in done.stderr
    assert not edges.exists()
    assert not receipt.exists()
