#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from ttf.relational_environment import (
    canonical_pca,
    frozen_grid_bounds,
    grid_centres,
    normalized_kde_grid,
    project_whiten,
    relation_pairs,
)
from ttf.relational_historical import (
    HISTORICAL_VARIABLES,
    LGM_TIME_INDEX,
    PRESENT_TIME_INDEX,
    historical_delta,
    logical_asset_basename,
    panel_order,
    panel_roles,
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sample_rasters(rows: list[dict[str, str]], paths: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    import rasterio

    datasets = [rasterio.open(path) for path in paths]
    try:
        coords = [(float(row["longitude"]), float(row["latitude"])) for row in rows]
        valid = np.ones(len(rows), dtype=bool)
        columns = []
        for ds in datasets:
            if ds.crs is None or not ds.crs.is_geographic:
                raise RuntimeError(f"raster is not geographic: {ds.name}")
            values = np.asarray([float(x[0]) for x in ds.sample(coords)], dtype=float)
            if ds.nodata is not None:
                valid &= ~np.isclose(values, float(ds.nodata), rtol=0.0, atol=0.0)
            valid &= np.isfinite(values)
            columns.append(values)
        return np.column_stack(columns), valid
    finally:
        for ds in datasets:
            ds.close()


def density_representation(by_species, names, values_key: str):
    pooled = []
    for name in names:
        pooled.extend(row[values_key] for row in by_species[name][:50])
    pooled = np.vstack(pooled)
    mean, sd, eigval, eigvec = canonical_pca(pooled)
    pooled_white = project_whiten(pooled, mean, sd, eigval, eigvec, axes=2)
    low, high = frozen_grid_bounds(pooled_white)
    centres = grid_centres(low, high, cells=50)

    density = np.empty((len(names), len(centres)), dtype=np.float64)
    for i, name in enumerate(names):
        x = np.vstack([row[values_key] for row in by_species[name]])
        white = project_whiten(x, mean, sd, eigval, eigvec, axes=2)
        density[i] = normalized_kde_grid(white, centres)
    return density, {
        "mean": mean,
        "sd": sd,
        "eigval": eigval,
        "eigvec": eigvec,
        "low": low,
        "high": high,
        "sample_rows": len(pooled),
    }


def qdict(x: np.ndarray) -> dict[str, float]:
    q = np.quantile(np.asarray(x, float), [0, .1, .25, .5, .75, .9, 1])
    return {k: float(v) for k, v in zip(
        ["min", "q10", "q25", "median", "q75", "q90", "max"], q
    )}


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    def ranks(a):
        a = np.asarray(a, float)
        order = np.argsort(a, kind="stable")
        out = np.empty(len(a), float)
        sx = a[order]
        starts = np.r_[0, 1 + np.flatnonzero(sx[1:] != sx[:-1])]
        stops = np.r_[starts[1:], len(a)]
        for start, stop in zip(starts, stops):
            out[order[start:stop]] = 0.5 * ((start + 1) + stop)
        return out
    a = ranks(x); b = ranks(y)
    a -= a.mean(); b -= b.mean()
    den = float(np.sqrt(np.dot(a, a) * np.dot(b, b)))
    return 0.0 if den <= np.finfo(float).eps else float(np.dot(a, b) / den)


def taxonomic_breadth_summary(
    names: list[str],
    metadata: dict[str, dict[str, str]],
    guardrail: dict[str, object],
) -> dict[str, object]:
    counts = Counter(
        (metadata[name].get("order") or "UNKNOWN").strip() or "UNKNOWN"
        for name in names
    )
    total = len(names)
    threshold = float(guardrail["order_fraction_threshold"])
    largest_max = float(guardrail["largest_single_order_fraction_max"])
    minimum_orders = int(guardrail["minimum_orders_at_or_above_fraction_threshold"])
    if total == 0:
        return {
            "species": 0,
            "orders": 0,
            "top_orders": {},
            "largest_order_fraction": None,
            "largest_single_order_fraction_max": largest_max,
            "order_fraction_threshold": threshold,
            "orders_at_or_above_fraction_threshold": [],
            "orders_at_or_above_fraction_threshold_count": 0,
            "minimum_orders_at_or_above_fraction_threshold": minimum_orders,
            "pass": False,
        }
    fractions = {name: count / total for name, count in counts.items()}
    qualifying = sorted(name for name, value in fractions.items() if value >= threshold)
    largest = max(fractions.values())
    return {
        "species": total,
        "orders": len(counts),
        "top_orders": dict(counts.most_common(12)),
        "largest_order_fraction": largest,
        "largest_single_order_fraction_max": largest_max,
        "order_fraction_threshold": threshold,
        "orders_at_or_above_fraction_threshold": qualifying,
        "orders_at_or_above_fraction_threshold_count": len(qualifying),
        "minimum_orders_at_or_above_fraction_threshold": minimum_orders,
        "pass": bool(largest <= largest_max and len(qualifying) >= minimum_orders),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--occurrences", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--asset-receipt", type=Path, required=True)
    ap.add_argument("--historical-asset", type=Path, action="append", required=True)
    ap.add_argument("--current-bio1", type=Path, required=True)
    ap.add_argument("--current-bio7", type=Path, required=True)
    ap.add_argument("--current-bio12", type=Path, required=True)
    ap.add_argument("--current-bio15", type=Path, required=True)
    ap.add_argument("--output-npz", type=Path, required=True)
    ap.add_argument("--output-summary", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_relational_historical_climate_exposure_rule_v0.1":
        raise RuntimeError("unexpected Study-C relation rule")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("Study C response firewall is open")

    receipt = json.loads(args.asset_receipt.read_text())
    if receipt.get("schema") != "ttf_relational_historical_climate_asset_receipt_v0.1":
        raise RuntimeError("unexpected historical asset receipt")
    if receipt.get("status") != "FROZEN_RESPONSE_BLIND_HISTORICAL_ASSETS":
        raise RuntimeError("historical assets are not frozen")

    required = list(rule["historical_climate"]["asset_contract"]["required_basenames"])
    asset_map = {path.name.casefold(): path for path in args.historical_asset}
    if len(asset_map) != 8:
        raise RuntimeError("exactly eight unique historical assets are required")
    receipt_map = {row["logical_basename"].casefold(): row for row in receipt["assets"]}
    ordered_hist = []
    for basename in required:
        key = basename.casefold()
        if key not in asset_map or key not in receipt_map:
            raise RuntimeError(f"missing historical asset: {basename}")
        path = asset_map[key]
        if sha256_path(path) != receipt_map[key]["sha256"]:
            raise RuntimeError(f"historical asset SHA drift: {basename}")
        ordered_hist.append(path)

    current_paths = [args.current_bio1, args.current_bio7, args.current_bio12, args.current_bio15]
    candidate_rows = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    metadata = {str(row["species"]): row for row in candidate_rows}
    if not (500 <= len(metadata) <= 1000):
        raise RuntimeError("Study-C candidate table must contain 500-1000 species")

    rows = list(csv.DictReader(args.occurrences.open(encoding="utf-8")))
    rows.sort(key=lambda row: (row["species"], int(row["priority_rank"]), int(row["source_key"])))
    if any(row["species"] not in metadata for row in rows):
        raise RuntimeError("occurrence species outside frozen Study-C candidates")

    hist_values, hist_valid = sample_rasters(rows, ordered_hist)
    current_values, current_valid = sample_rasters(rows, current_paths)
    valid = hist_valid & current_valid

    hist_col = {Path(path).name.casefold(): i for i, path in enumerate(ordered_hist)}
    past_idx = [hist_col[logical_asset_basename(v, LGM_TIME_INDEX).casefold()] for v in HISTORICAL_VARIABLES]
    present_idx = [hist_col[logical_asset_basename(v, PRESENT_TIME_INDEX).casefold()] for v in HISTORICAL_VARIABLES]
    delta = historical_delta(hist_values[:, past_idx], hist_values[:, present_idx])

    by_species = defaultdict(list)
    for i, (row, keep) in enumerate(zip(rows, valid)):
        if keep:
            by_species[row["species"]].append({
                "row": row,
                "delta": np.asarray(delta[i], float),
                "current": np.asarray(current_values[i], float),
            })

    minimum_occurrences = int(rule["external_occurrences"]["minimum_retained_occurrences_per_species"])
    admissible = sorted(
        name for name in metadata
        if len(by_species.get(name, ())) >= minimum_occurrences
    )
    minimum = int(rule["panel_assignment_after_historical_admissibility"]["minimum_admissible_species"])
    if len(admissible) < minimum:
        payload = {
            "schema": "ttf_relational_historical_climate_relation_design_v0.1",
            "status": "NOT_EVALUABLE_HISTORICAL_EXTERNAL_DATA",
            "candidate_species": len(metadata),
            "historical_environment_admissible_species": len(admissible),
            "minimum_required": minimum,
            "response_firewall": {
                "Study_C_sequence_identity_opened": False,
                "Study_C_pairwise_genetic_distances_opened": False,
                "Study_C_T_st_computed": False,
                "Study_C_beta_hist_computed": False,
            },
        }
        args.output_summary.parent.mkdir(parents=True, exist_ok=True)
        args.output_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(json.dumps(payload, sort_keys=True))
        return 0

    taxonomy_guard = taxonomic_breadth_summary(
        admissible,
        metadata,
        rule["taxonomic_breadth_guardrail"],
    )
    if not taxonomy_guard["pass"]:
        payload = {
            "schema": "ttf_relational_historical_climate_relation_design_v0.1",
            "status": str(rule["taxonomic_breadth_guardrail"]["if_failed"]),
            "candidate_species": len(metadata),
            "historical_environment_admissible_species": len(admissible),
            "minimum_required": minimum,
            "taxonomy": taxonomy_guard,
            "inputs": {
                "occurrence_csv_sha256": sha256_path(args.occurrences),
                "candidate_csv_sha256": sha256_path(args.candidates),
                "rule_sha256": sha256_path(args.rule),
                "asset_receipt_sha256": sha256_path(args.asset_receipt),
                "historical_assets_sha256": {p.name: sha256_path(p) for p in ordered_hist},
                "current_rasters_sha256": {p.name: sha256_path(p) for p in current_paths},
            },
            "response_firewall": {
                "Study_C_sequence_identity_opened": False,
                "Study_C_pairwise_genetic_distances_opened": False,
                "Study_C_T_st_computed": False,
                "Study_C_beta_hist_computed": False,
            },
        }
        args.output_summary.parent.mkdir(parents=True, exist_ok=True)
        args.output_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(json.dumps(payload, sort_keys=True))
        return 0

    hist_density, hist_repr = density_representation(by_species, admissible, "delta")
    current_density, current_repr = density_representation(by_species, admissible, "current")
    species_to_index = {name: i for i, name in enumerate(admissible)}
    source_sha = str(rule["independent_species_domain"]["source_archive_sha256"])
    ordered_panel = panel_order(admissible, source_sha)
    development = ordered_panel[:250]
    confirmatory = ordered_panel[250:]
    dev_source, dev_target = panel_roles(development, source_sha)
    con_source, con_target = panel_roles(confirmatory, source_sha)

    def make_pairs(source_names, target_names):
        source_idx = [species_to_index[x] for x in source_names]
        target_idx = [species_to_index[x] for x in target_names]
        hs, ht, rh = relation_pairs(hist_density, source_idx, target_idx)
        cs, ct, rc = relation_pairs(current_density, source_idx, target_idx)
        if not (np.array_equal(hs, cs) and np.array_equal(ht, ct)):
            raise RuntimeError("historical/current dyad geometry drift")
        return hs, ht, rh, rc

    ds, dt, drh, drc = make_pairs(dev_source, dev_target)
    cs, ct, crh, crc = make_pairs(con_source, con_target)
    occurrence_n = np.asarray([len(by_species[name]) for name in admissible], np.int64)
    class_name = np.asarray([metadata[name]["class"] for name in admissible], dtype="U64")
    order_name = np.asarray([metadata[name]["order"] for name in admissible], dtype="U64")
    family_name = np.asarray([metadata[name]["family"] for name in admissible], dtype="U96")

    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_npz,
        species_order=np.asarray(admissible, dtype="U160"),
        class_name=class_name,
        order_name=order_name,
        family_name=family_name,
        retained_environment_occurrences=occurrence_n,
        historical_environment_mean=hist_repr["mean"],
        historical_environment_sd=hist_repr["sd"],
        historical_pca_eigenvalues=hist_repr["eigval"],
        historical_pca_eigenvectors=hist_repr["eigvec"],
        historical_grid_low=hist_repr["low"],
        historical_grid_high=hist_repr["high"],
        current_environment_mean=current_repr["mean"],
        current_environment_sd=current_repr["sd"],
        current_pca_eigenvalues=current_repr["eigval"],
        current_pca_eigenvectors=current_repr["eigvec"],
        current_grid_low=current_repr["low"],
        current_grid_high=current_repr["high"],
        development_source_index=ds,
        development_target_index=dt,
        development_R_hist=drh,
        development_R_current=drc,
        confirmatory_source_index=cs,
        confirmatory_target_index=ct,
        confirmatory_R_hist=crh,
        confirmatory_R_current=crc,
    )

    all_hist = np.concatenate((drh, crh))
    all_current = np.concatenate((drc, crc))
    payload = {
        "schema": "ttf_relational_historical_climate_relation_design_v0.1",
        "status": "PASS_RESPONSE_BLIND_HISTORICAL_RELATION_DESIGN",
        "candidate_species": len(metadata),
        "historical_environment_admissible_species": len(admissible),
        "development_species": len(development),
        "development_sources": len(dev_source),
        "development_targets": len(dev_target),
        "development_dyads": int(len(drh)),
        "confirmatory_species": len(confirmatory),
        "confirmatory_sources": len(con_source),
        "confirmatory_targets": len(con_target),
        "confirmatory_dyads": int(len(crh)),
        "R_hist_quantiles_all_panel_dyads": qdict(all_hist),
        "R_current_quantiles_all_panel_dyads": qdict(all_current),
        "spearman_R_hist_vs_R_current": spearman(all_hist, all_current),
        "taxonomy": taxonomy_guard,
        "inputs": {
            "occurrence_csv_sha256": sha256_path(args.occurrences),
            "candidate_csv_sha256": sha256_path(args.candidates),
            "rule_sha256": sha256_path(args.rule),
            "asset_receipt_sha256": sha256_path(args.asset_receipt),
            "historical_assets_sha256": {p.name: sha256_path(p) for p in ordered_hist},
            "current_rasters_sha256": {p.name: sha256_path(p) for p in current_paths},
        },
        "design_npz_sha256": sha256_path(args.output_npz),
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_pairwise_genetic_distances_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
        "next_step": "Attach the frozen Study-C geographic-opportunity rule, then run development synthetic qualification. Genetic response remains closed.",
    }
    args.output_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "admissible": len(admissible),
        "development_dyads": len(drh),
        "confirmatory_dyads": len(crh),
        "spearman_R_hist_vs_R_current": payload["spearman_R_hist_vs_R_current"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
