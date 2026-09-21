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


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_key(tag: str, source_sha: str, species: str) -> tuple[str, str]:
    return hashlib.sha256(f"{tag}|{source_sha}|{species}".encode()).hexdigest(), species


def sample_rasters(rows: list[dict[str, str]], raster_paths: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    import rasterio

    datasets = [rasterio.open(path) for path in raster_paths]
    try:
        for dataset in datasets:
            if dataset.crs is None or not dataset.crs.is_geographic:
                raise RuntimeError(f"CHELSA raster is not geographic: {dataset.name}")
        coords = [(float(row["longitude"]), float(row["latitude"])) for row in rows]
        columns = []
        valid = np.ones(len(rows), dtype=bool)
        for dataset in datasets:
            values = np.asarray([float(sample[0]) for sample in dataset.sample(coords)], dtype=float)
            nodata = dataset.nodata
            if nodata is not None:
                valid &= ~np.isclose(values, float(nodata), rtol=0.0, atol=0.0)
            valid &= np.isfinite(values)
            columns.append(values)
        return np.column_stack(columns), valid
    finally:
        for dataset in datasets:
            dataset.close()


def panel_roles(names: list[str], source_sha: str, tag: str) -> tuple[list[str], list[str]]:
    ordered = sorted(names, key=lambda name: hash_key(tag, source_sha, name))
    n_source = len(ordered) // 2
    return sorted(ordered[:n_source]), sorted(ordered[n_source:])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--occurrences", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--freshness-rule", type=Path, required=True)
    ap.add_argument("--raster", type=Path, action="append", required=True)
    ap.add_argument("--output-npz", type=Path, required=True)
    ap.add_argument("--output-summary", type=Path, required=True)
    args = ap.parse_args()

    if len(args.raster) != 4:
        raise RuntimeError("exactly four frozen CHELSA rasters are required")
    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_relational_environment_relation_rule_v0.3":
        raise RuntimeError("unexpected environmental relation rule")
    freshness = json.loads(args.freshness_rule.read_text())
    if freshness.get("schema") != "ttf_relational_future_family_freshness_amendment_v0.1":
        raise RuntimeError("unexpected future-family freshness rule")
    if any(bool(v) for v in rule["response_firewall"].values()):
        raise RuntimeError("Study B response firewall is open")

    candidate_rows = list(csv.DictReader(args.candidates.open(encoding="utf-8")))
    metadata = {str(row["species"]): row for row in candidate_rows}
    if len(metadata) != 1000:
        raise RuntimeError("candidate metadata must contain exact fresh 1000")

    rows = list(csv.DictReader(args.occurrences.open(encoding="utf-8")))
    by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["species"] not in metadata:
            raise RuntimeError(f"occurrence species outside fresh 1000: {row['species']}")
        by_species[row["species"]].append(row)
    for values in by_species.values():
        values.sort(key=lambda row: (int(row["priority_rank"]), int(row["source_key"])))

    environment, valid = sample_rasters(rows, list(args.raster))
    valid_by_species: dict[str, list[tuple[dict[str, str], np.ndarray]]] = defaultdict(list)
    for row, env, keep in zip(rows, environment, valid):
        if keep:
            valid_by_species[row["species"]].append((row, np.asarray(env, dtype=float)))

    raw_admissible = sorted(name for name, values in valid_by_species.items() if len(values) >= 30)
    excluded = set(map(str, freshness["exact_overlap_with_study_B_candidates"]["union_species"]))
    if len(excluded) != int(freshness["exact_overlap_with_study_B_candidates"]["union_overlap_species"]):
        raise RuntimeError("freshness exclusion count drift")
    if not excluded.issubset(metadata):
        raise RuntimeError("freshness exclusion species outside frozen candidate 1000")
    excluded_admissible = sorted(set(raw_admissible) & excluded)
    admissible = sorted(name for name in raw_admissible if name not in excluded)
    minimum_after = int(rule["freshness_before_pca_and_panel"]["minimum_post_exclusion_environment_admissible_species"])
    if len(admissible) < minimum_after:
        payload = {
            "schema": "ttf_relational_environment_relation_design_v0.3",
            "status": "NOT_EVALUABLE_ENVIRONMENT_FRESHNESS",
            "fresh_candidate_species": 1000,
            "environment_admissible_before_freshness": len(raw_admissible),
            "freshness_excluded_total": len(excluded),
            "freshness_excluded_environment_admissible": len(excluded_admissible),
            "environment_admissible_species": len(admissible),
            "minimum_required_after_freshness": minimum_after,
            "response_firewall": {
                "study_B_sequence_identity_opened": False,
                "study_B_pairwise_genetic_distances_opened": False,
                "study_B_T_st_computed": False,
                "study_B_beta_R_computed": False,
            },
        }
        args.output_summary.parent.mkdir(parents=True, exist_ok=True)
        args.output_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(json.dumps(payload, sort_keys=True))
        # A frozen feasibility failure is a valid NOT_EVALUABLE program state,
        # not an infrastructure error. Return success so the workflow can
        # persist the gate receipt and stop downstream response-blind stages.
        return 0

    pca_sample_rows = []
    pca_sample_species = []
    for name in admissible:
        values = valid_by_species[name][:50]
        pca_sample_rows.extend(env for _, env in values)
        pca_sample_species.extend([name] * len(values))
    pca_sample = np.vstack(pca_sample_rows)
    mean, sd, eigval, eigvec = canonical_pca(pca_sample)
    pca_white = project_whiten(pca_sample, mean, sd, eigval, eigvec, axes=2)
    low, high = frozen_grid_bounds(pca_white)
    centres = grid_centres(low, high, cells=50)

    density = np.empty((len(admissible), len(centres)), dtype=np.float64)
    occurrence_n = np.empty(len(admissible), dtype=np.int64)
    species_to_index = {name: i for i, name in enumerate(admissible)}
    for name in admissible:
        env = np.vstack([value for _, value in valid_by_species[name]])
        white = project_whiten(env, mean, sd, eigval, eigvec, axes=2)
        density[species_to_index[name]] = normalized_kde_grid(white, centres)
        occurrence_n[species_to_index[name]] = len(env)

    source_sha = str(json.loads(Path("benchmarks/frozen/relational_fresh_candidate_census_v0.1.json").read_text())["source_archive_sha256"])
    panel_order = sorted(admissible, key=lambda name: hash_key("relational-env-panel-v0.1", source_sha, name))
    development = panel_order[:250]
    confirmatory = panel_order[250:]
    dev_source, dev_target = panel_roles(development, source_sha, "relational-env-role-v0.1")
    con_source, con_target = panel_roles(confirmatory, source_sha, "relational-env-role-v0.1")

    def make_pairs(source_names: list[str], target_names: list[str]):
        source_idx = [species_to_index[name] for name in source_names]
        target_idx = [species_to_index[name] for name in target_names]
        return relation_pairs(density, source_idx, target_idx)

    ds, dt, dr = make_pairs(dev_source, dev_target)
    cs, ct, cr = make_pairs(con_source, con_target)

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
        environment_mean=mean,
        environment_sd=sd,
        pca_eigenvalues=eigval,
        pca_eigenvectors=eigvec,
        grid_low=low,
        grid_high=high,
        density=density.astype(np.float32),
        development_source_index=ds,
        development_target_index=dt,
        development_R_env=dr,
        confirmatory_source_index=cs,
        confirmatory_target_index=ct,
        confirmatory_R_env=cr,
    )

    variance_fraction = eigval / eigval.sum()
    relation_all = np.concatenate((dr, cr))
    q = np.quantile(relation_all, [0, .1, .25, .5, .75, .9, 1])
    order_counts = Counter(metadata[name]["order"] for name in admissible)
    summary = {
        "schema": "ttf_relational_environment_relation_design_v0.3",
        "status": "PASS_RESPONSE_BLIND_ENVIRONMENT_RELATION_DESIGN",
        "fresh_candidate_species": 1000,
        "environment_admissible_before_freshness": len(raw_admissible),
        "freshness_excluded_total": len(excluded),
        "freshness_excluded_environment_admissible": len(excluded_admissible),
        "freshness_excluded_species_environment_admissible": excluded_admissible,
        "environment_admissible_species": len(admissible),
        "development_species": len(development),
        "development_sources": len(dev_source),
        "development_targets": len(dev_target),
        "development_dyads": int(len(dr)),
        "confirmatory_species": len(confirmatory),
        "confirmatory_sources": len(con_source),
        "confirmatory_targets": len(con_target),
        "confirmatory_dyads": int(len(cr)),
        "pca": {
            "sample_rows": len(pca_sample),
            "eigenvalues": eigval.tolist(),
            "variance_fraction": variance_fraction.tolist(),
            "mean": mean.tolist(),
            "sd": sd.tolist(),
            "eigenvectors_columns_are_PCs": eigvec.tolist(),
            "retained_axes": 2,
        },
        "R_env_quantiles_all_panel_dyads": {
            "min": float(q[0]), "q10": float(q[1]), "q25": float(q[2]),
            "median": float(q[3]), "q75": float(q[4]), "q90": float(q[5]), "max": float(q[6]),
        },
        "taxonomy": {
            "orders": len(order_counts),
            "top_orders": dict(order_counts.most_common(12)),
            "largest_order_fraction": max(order_counts.values()) / len(admissible),
        },
        "inputs": {
            "occurrence_csv_sha256": sha256_path(args.occurrences),
            "candidate_csv_sha256": sha256_path(args.candidates),
            "rule_sha256": sha256_path(args.rule),
            "freshness_rule_sha256": sha256_path(args.freshness_rule),
            "rasters_sha256": {path.name: sha256_path(path) for path in args.raster},
        },
        "design_npz_sha256": sha256_path(args.output_npz),
        "response_firewall": {
            "study_B_sequence_identity_opened": False,
            "study_B_pairwise_genetic_distances_opened": False,
            "study_B_T_st_computed": False,
            "study_B_beta_R_computed": False,
        },
        "next_step": "Attach frozen phylogatR geographic-opportunity controls to these exact panel dyads, then run the already frozen v0.2 development synthetic qualification at alpha=0.025. Do not open Study B character masks or nucleotide identity.",
    }
    args.output_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("status", "environment_admissible_species", "development_dyads", "confirmatory_dyads", "R_env_quantiles_all_panel_dyads")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
