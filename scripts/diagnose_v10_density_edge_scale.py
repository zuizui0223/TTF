#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from ttf import knn_edges
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.private_geometry_control import pooled_affine_frame
from ttf.shared_geometry_observability import shared_transition_observability

DENSITIES = (20, 40, 60, 80, 100)
FIXED_K = (3, 3, 3, 3, 3)
FRACTION_K = (3, 6, 9, 12, 15)


def priority(seed: int, species: str, photo_id: str) -> str:
    return hashlib.sha256(f"{seed}|{species}|{photo_id}".encode("utf-8")).hexdigest()


def load_nested(path: Path, *, seed: int) -> dict[int, tuple[SpeciesGeometry, ...]]:
    grouped: dict[str, list[tuple[str, np.ndarray]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            species = str(row["species"])
            photo_id = str(row["source_photo_id"])
            coord = np.asarray([float(row["x_km"]), float(row["y_km"]), float(row["z_km"])], dtype=float)
            grouped[species].append((photo_id, coord))
    if len(grouped) != 250:
        raise RuntimeError(f"expected 250 species, found {len(grouped)}")
    ordered: dict[str, list[np.ndarray]] = {}
    for species, rows in grouped.items():
        if len(rows) != 100 or len({p for p, _ in rows}) != 100:
            raise RuntimeError(f"{species}: expected 100 unique records")
        rows.sort(key=lambda item: (priority(seed, species, item[0]), item[0]))
        ordered[species] = [coord for _, coord in rows]
    out: dict[int, tuple[SpeciesGeometry, ...]] = {}
    for n in DENSITIES:
        out[n] = tuple(
            SpeciesGeometry(species=s, coordinates=np.vstack(ordered[s][:n]))
            for s in sorted(ordered)
        )
    return out


def graph_metrics(geometries: tuple[SpeciesGeometry, ...], *, k: int) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    center, scale, _ = pooled_affine_frame(geometries)
    graphs: dict[str, np.ndarray] = {}
    species_median_raw: list[float] = []
    species_mean_raw: list[float] = []
    species_median_std: list[float] = []
    species_mean_std: list[float] = []
    edge_counts: list[int] = []
    for g in geometries:
        nodes = knn_edges(g.coordinates, k=k)
        graphs[g.species] = nodes
        a = g.coordinates[nodes[:, 0]]
        b = g.coordinates[nodes[:, 1]]
        raw = np.linalg.norm(a - b, axis=1)
        z = (g.coordinates - center) / scale
        std = np.linalg.norm(z[nodes[:, 0]] - z[nodes[:, 1]], axis=1)
        species_median_raw.append(float(np.median(raw)))
        species_mean_raw.append(float(np.mean(raw)))
        species_median_std.append(float(np.median(std)))
        species_mean_std.append(float(np.mean(std)))
        edge_counts.append(int(len(nodes)))
    sm = np.asarray(species_median_raw)
    ss = np.asarray(species_median_std)
    return graphs, {
        "equal_species_mean_edge_length_km": float(np.mean(species_mean_raw)),
        "equal_species_median_edge_length_km": float(np.mean(species_median_raw)),
        "species_median_edge_length_km": {
            "median": float(np.median(sm)),
            "q10": float(np.quantile(sm, 0.10)),
            "q90": float(np.quantile(sm, 0.90)),
        },
        "equal_species_mean_standardized_edge_length": float(np.mean(species_mean_std)),
        "equal_species_median_standardized_edge_length": float(np.mean(species_median_std)),
        "species_median_standardized_edge_length": {
            "median": float(np.median(ss)),
            "q10": float(np.quantile(ss, 0.10)),
            "q90": float(np.quantile(ss, 0.90)),
        },
        "equal_species_mean_edge_count": float(np.mean(edge_counts)),
    }


def log_slope(x: list[float], y: list[float]) -> float:
    lx = np.log(np.asarray(x, dtype=float))
    ly = np.log(np.asarray(y, dtype=float))
    if np.any(~np.isfinite(ly)) or np.any(np.asarray(y) <= 0):
        return float("nan")
    return float(np.polyfit(lx, ly, 1)[0])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--source-ledger", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    source = json.loads(args.source_ledger.read_text())
    if rule.get("status") != "frozen_before_v10_density_edge_scale_diagnostic_execution":
        raise RuntimeError("v0.10 rule drift")
    if source.get("schema") != "ttf_v09_rgfca_reserve_b_max_density_geometry_v0.1":
        raise RuntimeError("wrong source geometry")
    if source.get("empirical_trait_or_colour_fields_read") is not False:
        raise RuntimeError("empirical outcomes already opened")
    if source.get("output", {}).get("records") != 25000:
        raise RuntimeError("v0.9 source size drift")

    nested = load_nested(args.input, seed=20260921)
    rows: list[dict[str, object]] = []
    for graph_rule, ks in (("fixed_k3", FIXED_K), ("constant_neighbor_fraction", FRACTION_K)):
        for n, k in zip(DENSITIES, ks):
            geometries = nested[n]
            graphs, metrics = graph_metrics(geometries, k=k)
            audit = shared_transition_observability(
                geometries,
                graphs,
                transition_width=0.2,
                n_directions=2048,
                seed=20260922,
                nontrivial_edge_contrast=0.25,
            )
            nontrivial = np.asarray(list(audit.species_nontrivial_fraction.values()), dtype=float)
            rows.append({
                "graph_rule": graph_rule,
                "records_per_species": n,
                "k": k,
                "geometry_fingerprint_sha256": geometry_fingerprint(geometries),
                **metrics,
                "shared_observability": audit.summary(),
                "equal_species_mean_nontrivial_direction_fraction": float(np.mean(nontrivial)),
            })

    by_rule: dict[str, dict[str, object]] = {}
    for graph_rule in ("fixed_k3", "constant_neighbor_fraction"):
        sub = [r for r in rows if r["graph_rule"] == graph_rule]
        n = [float(r["records_per_species"]) for r in sub]
        edge = [float(r["equal_species_median_standardized_edge_length"]) for r in sub]
        contrast = [float(r["shared_observability"]["equal_species_mean_contrast"]["mean"]) for r in sub]
        by_rule[graph_rule] = {
            "loglog_slope_standardized_edge_length_vs_n": log_slope(n, edge),
            "loglog_slope_shared_mean_contrast_vs_n": log_slope(n, contrast),
            "n20_to_n100_edge_length_ratio": float(edge[-1] / edge[0]),
            "n20_to_n100_shared_contrast_ratio": float(contrast[-1] / contrast[0]),
            "monotone_nonincreasing_edge_scale": bool(all(b <= a + 1e-12 for a, b in zip(edge, edge[1:]))),
            "monotone_nonincreasing_shared_contrast": bool(all(b <= a + 1e-12 for a, b in zip(contrast, contrast[1:]))),
        }

    fixed = by_rule["fixed_k3"]
    scaled = by_rule["constant_neighbor_fraction"]
    mechanism_supported = bool(
        fixed["monotone_nonincreasing_edge_scale"]
        and fixed["monotone_nonincreasing_shared_contrast"]
        and float(fixed["n20_to_n100_edge_length_ratio"]) < 1.0
        and float(fixed["n20_to_n100_shared_contrast_ratio"]) < 1.0
    )
    fraction_rule_stabilizes = bool(
        abs(float(scaled["loglog_slope_standardized_edge_length_vs_n"]))
        < abs(float(fixed["loglog_slope_standardized_edge_length_vs_n"]))
        and abs(float(scaled["loglog_slope_shared_mean_contrast_vs_n"]))
        < abs(float(fixed["loglog_slope_shared_mean_contrast_vs_n"]))
    )

    payload = {
        "schema": "ttf_v10_density_edge_scale_diagnostic_v0.1",
        "status": "development_only_geometry_mechanism_audit",
        "rule": str(args.rule),
        "source_ledger": str(args.source_ledger),
        "trait_or_colour_outcomes_used": False,
        "synthetic_trait_worlds_used": False,
        "ttf_estimator_fit": False,
        "rows": rows,
        "trend_summary": by_rule,
        "predeclared_signature": {
            "fixed_k_edge_localization_supported": mechanism_supported,
            "constant_neighbor_fraction_stabilizes_geometry_signature": fraction_rule_stabilizes,
        },
        "claim_ready": False,
        "candidate_selectable": False,
        "interpretation_boundary": "This audit may explain failed-development behavior but cannot qualify k/n scaling, authorize an empirical analysis, or open a new confirmatory panel by itself."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"signature": payload["predeclared_signature"], "trend_summary": by_rule}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
