#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from ttf.genetic_conditional_qualification import (
    make_conditional_order_worlds,
    prepare_conditional_order_qualification,
    score_conditional_order_world_batch,
)
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry


RULE_SCHEMA = "ttf_genetic_conditional_order_qualification_rule_v0.1"
GEOMETRY_SCHEMA = "ttf_genetic_conditional_response_blind_geometry_v0.1"


def sha256_path(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path: Path, schema: str) -> dict:
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def load_geometry(csv_path: Path, manifest: dict):
    if sha256_path(csv_path) != manifest["geometry_csv_sha256"]:
        raise RuntimeError("geometry CSV hash mismatch")
    if manifest.get("status") != "FROZEN_RESPONSE_BLIND_GEOMETRY_NO_NUCLEOTIDE_IDENTITY":
        raise RuntimeError("geometry manifest is not response-blind frozen geometry")
    firewall = manifest.get("outcome_firewall", {})
    if not firewall or any(value is not False for value in firewall.values()):
        raise RuntimeError("geometry manifest outcome firewall is not closed")
    if manifest.get("panel") != "development":
        raise RuntimeError("formal conditional qualification may use the development panel only")

    per_species: dict[str, list[dict[str, str]]] = {}
    with Path(csv_path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "species", "order", "locality_index", "x_km", "y_km", "z_km",
            "graph_k", "edge_count", "min_endpoint_disjoint_training_edges", "split",
        }
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("conditional geometry CSV schema mismatch")
        for row in reader:
            per_species.setdefault(str(row["species"]), []).append(row)

    if set(per_species) != set(manifest["panels"]):
        raise RuntimeError("geometry CSV species differ from manifest")
    geometries = {}
    taxonomy_order = {}
    for species in sorted(per_species):
        rows = sorted(per_species[species], key=lambda row: int(row["locality_index"]))
        if [int(row["locality_index"]) for row in rows] != list(range(len(rows))):
            raise RuntimeError(f"noncanonical locality indices for {species}")
        coords = np.asarray(
            [[float(row["x_km"]), float(row["y_km"]), float(row["z_km"])] for row in rows],
            dtype=float,
        )
        geometry = prepare_density_scaled_genetic_geometry(coords, neighbor_fraction=0.15)
        panel = manifest["panels"][species]
        if geometry.graph_k != int(panel["graph_k"]):
            raise RuntimeError(f"graph-k drift for {species}")
        if geometry.n_edges != int(panel["edge_count"]):
            raise RuntimeError(f"edge-count drift for {species}")
        if (
            geometry.min_endpoint_disjoint_training_edges
            != int(panel["min_endpoint_disjoint_training_edges"])
        ):
            raise RuntimeError(f"endpoint-support drift for {species}")
        geometries[species] = geometry
        orders = {str(row["order"]).strip() for row in rows}
        if len(orders) != 1:
            raise RuntimeError(f"order label drift within {species}")
        taxonomy_order[species] = orders.pop()

    train = tuple(map(str, manifest["train_species"]))
    evaluation = tuple(map(str, manifest["eval_species"]))
    if set(train) & set(evaluation):
        raise RuntimeError("development split overlap")
    if set(train) | set(evaluation) != set(geometries):
        raise RuntimeError("development split is not a geometry partition")
    return geometries, taxonomy_order, train, evaluation


def resolve_cell(rule: dict, cell_name: str) -> tuple[dict, str]:
    for row in rule["synthetic_worlds"]["private_null_cells"]:
        if row["cell"] == cell_name:
            return row, "private"
    positive = rule["synthetic_worlds"]["positive_control"]
    if positive["cell"] == cell_name:
        return positive, "same_order"
    raise RuntimeError(f"unknown frozen qualification cell: {cell_name}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry-csv", type=Path, required=True)
    parser.add_argument("--geometry-manifest", type=Path, required=True)
    parser.add_argument(
        "--rule",
        type=Path,
        default=Path("docs/supporting/genetic_conditional_order_qualification_rule_v0.1.json"),
    )
    parser.add_argument("--cell", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--edge-chunk-size", type=int, default=32)
    parser.add_argument("--train-chunk-size", type=int, default=4096)
    args = parser.parse_args()

    if args.start < 0 or args.count < 1:
        parser.error("--start must be >=0 and --count must be >=1")
    rule = load_json(args.rule, RULE_SCHEMA)
    if not str(rule.get("status", "")).startswith("FROZEN_BEFORE_ANY_CONDITIONAL"):
        raise RuntimeError("qualification rule is not frozen pre-outcome")
    if any(value is not False for value in rule["outcome_firewall"].values()):
        raise RuntimeError("qualification rule outcome firewall is not closed")
    manifest = load_json(args.geometry_manifest, GEOMETRY_SCHEMA)
    panel = rule["development_panel"]
    if sha256_path(args.geometry_csv) != panel["geometry_csv_sha256"]:
        raise RuntimeError("frozen development geometry CSV hash mismatch")
    if sha256_path(args.geometry_manifest) != panel["geometry_manifest_sha256"]:
        raise RuntimeError("frozen development geometry manifest hash mismatch")
    geometries, taxonomy_order, train, evaluation = load_geometry(
        args.geometry_csv, manifest
    )

    if len(geometries) != int(panel["species"]):
        raise RuntimeError("development species count drift")
    if len(train) != int(panel["train_species"]) or len(evaluation) != int(panel["eval_species"]):
        raise RuntimeError("development split count drift")

    cell, mode = resolve_cell(rule, args.cell)
    if args.start + args.count > int(cell["worlds"]):
        raise RuntimeError("requested shard exceeds frozen cell world count")
    estimator = rule["frozen_estimator"]
    design = prepare_conditional_order_qualification(
        geometries,
        taxonomy_order,
        train_species=train,
        eval_species=evaluation,
        bandwidth=float(estimator["field_bandwidth_km"]),
        support_radius=float(estimator["support_radius_km"]),
        minimum_target_coverage=float(estimator["minimum_target_coverage"]),
        minimum_source_species=int(estimator["minimum_source_species"]),
        edge_chunk_size=int(args.edge_chunk_size),
        train_chunk_size=int(args.train_chunk_size),
    )
    if len(design.eligible_eval_species) != int(estimator["eligible_same_order_eval_species"]):
        raise RuntimeError(
            f"response-blind supported target count drift: "
            f"{len(design.eligible_eval_species)} != {estimator['eligible_same_order_eval_species']}"
        )

    worlds_cfg = rule["synthetic_worlds"]
    worlds = make_conditional_order_worlds(
        design,
        cell=str(args.cell),
        residual_amplitude=float(cell["residual_amplitude"]),
        absolute_start=int(args.start),
        count=int(args.count),
        group_mode=mode,
        ibd_strength=float(worlds_cfg["ibd_strength"]),
        noise_sd=float(worlds_cfg["noise_sd"]),
        transition_width=float(worlds_cfg["transition_width"]),
        noise_dimensions=int(worlds_cfg["noise_dimensions"]),
        master_seed=int(worlds_cfg["master_seed"]),
    )
    scored = score_conditional_order_world_batch(
        design,
        worlds,
        cell=str(args.cell),
        absolute_start=int(args.start),
        bootstrap_resamples=int(estimator["bootstrap_resamples"]),
        master_seed=int(worlds_cfg["master_seed"]),
    )

    rows = []
    for offset in range(args.count):
        rows.append(
            {
                "absolute_replicate_index": int(args.start + offset),
                "statistic": float(scored.statistics[offset]),
                "p_value": float(scored.p_values[offset]),
                "reject": bool(scored.p_values[offset] <= float(estimator["alpha"])),
            }
        )
    payload = {
        "schema": "ttf_genetic_conditional_order_qualification_shard_v0.1",
        "status": "SYNTHETIC_QUALIFICATION_SHARD_COMPLETE",
        "cell": str(args.cell),
        "group_mode": mode,
        "residual_amplitude": float(cell["residual_amplitude"]),
        "start": int(args.start),
        "count": int(args.count),
        "n_eval_species": int(scored.n_eval_species),
        "rule_sha256": sha256_path(args.rule),
        "geometry_manifest_sha256": sha256_path(args.geometry_manifest),
        "geometry_csv_sha256": sha256_path(args.geometry_csv),
        "rows": rows,
        "outcome_firewall": {
            "development_sequence_identity_opened": False,
            "confirmatory_sequence_identity_opened": False,
            "conditional_empirical_result_opened": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "cell": args.cell,
                "start": args.start,
                "count": args.count,
                "rejections": sum(row["reject"] for row in rows),
                "mean_statistic": float(np.mean(scored.statistics)),
                "n_eval_species": scored.n_eval_species,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
