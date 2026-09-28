#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from ttf.external_geographic_opportunity import directed_occurrence_coverage


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_occurrence_geometry(path: Path) -> dict[str, np.ndarray]:
    grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                lat = float(row["latitude"])
                lon = float(row["longitude"])
            except (KeyError, TypeError, ValueError):
                continue
            if (
                np.isfinite(lat)
                and np.isfinite(lon)
                and -90.0 <= lat <= 90.0
                and -180.0 <= lon <= 180.0
            ):
                grouped[str(row["species"])].append((lat, lon))
    return {
        species: np.asarray(values, dtype=float)
        for species, values in grouped.items()
        if values
    }


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return {str(row["species"]): row for row in csv.DictReader(handle)}


def qdict(values: np.ndarray) -> dict[str, float]:
    q = np.quantile(np.asarray(values, dtype=float), [0, .1, .25, .5, .75, .9, 1])
    return {
        key: float(value)
        for key, value in zip(
            ("min", "q10", "q25", "median", "q75", "q90", "max"),
            q,
        )
    }


def process_panel(
    panel: str,
    data: np.lib.npyio.NpzFile,
    species: np.ndarray,
    metadata: dict[str, dict[str, str]],
    occurrence_latlon: dict[str, np.ndarray],
    retained_environment_occurrences: np.ndarray,
    *,
    radius_km: float,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    source = np.asarray(data[f"{panel}_source_index"], dtype=np.int64)
    target = np.asarray(data[f"{panel}_target_index"], dtype=np.int64)
    r_hist = np.asarray(data[f"{panel}_R_hist"], dtype=float)
    r_current = np.asarray(data[f"{panel}_R_current"], dtype=float)
    if not (len(source) == len(target) == len(r_hist) == len(r_current)):
        raise RuntimeError("historical relation dyad arrays do not align")

    missing = sorted({
        str(species[int(index)])
        for index in np.concatenate((source, target))
        if str(species[int(index)]) not in occurrence_latlon
    })
    if missing:
        raise RuntimeError(
            f"missing response-blind occurrence geometry for {len(missing)} panel species"
        )

    coverage = np.empty(len(source), dtype=float)
    same_class = np.empty(len(source), dtype=np.int8)
    same_order = np.empty(len(source), dtype=np.int8)
    same_family = np.empty(len(source), dtype=np.int8)
    sample_ratio = np.empty(len(source), dtype=float)

    for i, (s, t) in enumerate(zip(source, target)):
        source_name = str(species[int(s)])
        target_name = str(species[int(t)])
        coverage[i] = directed_occurrence_coverage(
            occurrence_latlon[target_name],
            occurrence_latlon[source_name],
            radius_km=radius_km,
        )
        source_meta = metadata[source_name]
        target_meta = metadata[target_name]
        same_class[i] = int(source_meta.get("class", "") == target_meta.get("class", ""))
        same_order[i] = int(source_meta.get("order", "") == target_meta.get("order", ""))
        same_family[i] = int(source_meta.get("family", "") == target_meta.get("family", ""))
        sample_ratio[i] = (
            retained_environment_occurrences[int(s)]
            / retained_environment_occurrences[int(t)]
        )
        if (i + 1) % 5000 == 0:
            print(json.dumps({
                "panel": panel,
                "dyads_done": i + 1,
                "dyads_total": len(source),
            }, sort_keys=True))

    arrays = {
        f"{panel}_source_index": source,
        f"{panel}_target_index": target,
        f"{panel}_R_hist": r_hist,
        f"{panel}_R_current": r_current,
        f"{panel}_coverage": coverage,
        f"{panel}_same_class": same_class,
        f"{panel}_same_order": same_order,
        f"{panel}_same_family": same_family,
        f"{panel}_locality_count_ratio": sample_ratio,
    }
    summary = {
        "panel": panel,
        "dyads": int(len(source)),
        "source_species": int(len(np.unique(source))),
        "target_species": int(len(np.unique(target))),
        "geographic_coverage_quantiles": qdict(coverage),
        "same_class_fraction": float(same_class.mean()),
        "same_order_fraction": float(same_order.mean()),
        "same_family_fraction": float(same_family.mean()),
    }
    return arrays, summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--historical-design", type=Path, required=True)
    parser.add_argument("--occurrences", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument(
        "--rule",
        type=Path,
        default=Path("docs/supporting/ttf_q_c2_geographic_opportunity_v0.1.json"),
    )
    parser.add_argument(
        "--panel",
        choices=("development", "confirmatory", "both"),
        default="development",
    )
    parser.add_argument("--output-npz", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    args = parser.parse_args()

    rule = json.loads(args.rule.read_text())
    if rule.get("schema") != "ttf_q_c2_geographic_opportunity_v0.1":
        raise RuntimeError("unexpected TTF-Q C2 rule")
    if rule["program_boundary"]["genetic_response_opening_authorized"]:
        raise RuntimeError("TTF-Q C2 must remain response blind")
    radius = float(rule["geographic_opportunity"]["support_radius_km"])

    data = np.load(args.historical_design, allow_pickle=False)
    species = np.asarray(data["species_order"]).astype(str)
    retained = np.asarray(data["retained_environment_occurrences"], dtype=float)
    if retained.shape != species.shape or np.any(retained <= 0):
        raise RuntimeError("retained occurrence counts do not align with species")
    metadata = load_metadata(args.candidates)
    occurrence_latlon = load_occurrence_geometry(args.occurrences)

    panels = (
        ("development", "confirmatory")
        if args.panel == "both"
        else (args.panel,)
    )
    arrays: dict[str, np.ndarray] = {}
    summaries: dict[str, object] = {}
    for panel in panels:
        panel_arrays, panel_summary = process_panel(
            panel,
            data,
            species,
            metadata,
            occurrence_latlon,
            retained,
            radius_km=radius,
        )
        arrays.update(panel_arrays)
        summaries[panel] = panel_summary

    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output_npz, species_order=species, **arrays)
    payload = {
        "schema": "ttf_q_c2_external_opportunity_design_v0.1",
        "status": "CHARACTERIZED_EXTERNAL_GEOGRAPHIC_OPPORTUNITY",
        "historical_design_sha256": sha256_path(args.historical_design),
        "occurrence_csv_sha256": sha256_path(args.occurrences),
        "candidate_csv_sha256": sha256_path(args.candidates),
        "rule_sha256": sha256_path(args.rule),
        "geographic_opportunity": {
            "definition": (
                "directed fraction of retained target occurrence coordinates "
                "whose nearest retained source occurrence coordinate is within "
                f"{radius:g} km great-circle distance"
            ),
            "support_radius_km": radius,
            "hard_coverage_filter_applied": False,
            "role": "continuous response-blind spatial co-opportunity control",
        },
        "sampling_imbalance": {
            "definition": (
                "source/target retained environmental occurrence-count ratio, "
                "identical to the count information used in C1"
            ),
            "reason": "C2 differs from C1 by geography alone.",
        },
        "panels": summaries,
        "design_npz_sha256": sha256_path(args.output_npz),
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_pairwise_genetic_distances_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }
    args.output_summary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps({
        "status": payload["status"],
        "panels": summaries,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
