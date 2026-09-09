#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

TERMINAL = {
    "classified_four_state_morph",
    "not_evaluable_insufficient_flower_pixels",
    "not_evaluable_no_biological_palette_mass",
    "not_evaluable_ambiguous_palette_composition",
    "not_evaluable_roi_or_flip_gate",
    "image_acquisition_failed",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def keep_key(photo_id: object, seed: int) -> tuple[bytes, str]:
    value = str(photo_id)
    return hashlib.sha256(f"{int(seed)}\x1f{value}".encode("utf-8")).digest(), value


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--firewall-dir", type=Path, required=True)
    ap.add_argument("--rule", type=Path, required=True)
    ap.add_argument("--output-geometry", type=Path, required=True)
    ap.add_argument("--output-ledger", type=Path, required=True)
    args = ap.parse_args()

    rule = json.loads(args.rule.read_text())
    fw = json.loads((args.firewall_dir / "firewall_manifest.json").read_text())
    if rule.get("schema") != "ttf_v13_minimum_classifiability_geometry_rule_v0.1":
        raise RuntimeError("wrong v0.13 rule")
    if fw.get("status") != "actual_geometry_measurement_firewall_frozen_before_pixels":
        raise RuntimeError("wrong inherited classifiability firewall")
    if fw.get("colour_values_opened") is not False or fw.get("pixels_opened_at_firewall_freeze") is not False:
        raise RuntimeError("inherited firewall already opened forbidden outcomes")

    csvs = sorted(args.results_dir.glob("**/partition_s*_p*.csv"))
    receipts = sorted(args.results_dir.glob("**/partition_s*_p*.json"))
    if len(csvs) != 128 or len(receipts) != 128:
        raise RuntimeError(f"requires all 128 terminal partitions; got csv={len(csvs)} json={len(receipts)}")

    frames = []
    for path in csvs:
        # Read only the location-blind terminal state allowed by the inherited firewall.
        frame = pd.read_csv(path, usecols=["measurement_id", "measurement_status"], dtype=str).fillna("")
        frames.append(frame)
    terminal = pd.concat(frames, ignore_index=True)
    if len(terminal) != 25000 or terminal["measurement_id"].nunique() != 25000:
        raise RuntimeError("terminal coverage is not exactly 25,000 unique measurement IDs")
    if not set(terminal["measurement_status"]).issubset(TERMINAL):
        raise RuntimeError("unknown terminal measurement status")

    worker = pd.read_csv(args.firewall_dir / "worker_manifest.csv", dtype=str).fillna("")
    meta = pd.read_csv(
        args.firewall_dir / "metadata_join_key.csv",
        dtype={"measurement_id": str, "source_photo_id": str, "source_observation_id": str},
    )
    if len(worker) != 25000 or len(meta) != 25000:
        raise RuntimeError("firewall denominator drift")
    expected = set(worker["measurement_id"].astype(str))
    if set(terminal["measurement_id"].astype(str)) != expected or set(meta["measurement_id"].astype(str)) != expected:
        raise RuntimeError("terminal or metadata IDs differ from frozen worker denominator")

    joined = meta.merge(terminal, on="measurement_id", how="inner", validate="one_to_one")
    joined["measurement_evaluable"] = joined["measurement_status"].eq("classified_four_state_morph")
    support = joined.groupby("species", observed=True)["measurement_evaluable"].agg(["sum", "size"]).reset_index()
    support = support.rename(columns={"sum": "n_evaluable", "size": "n_frozen"})
    support["n_evaluable"] = support["n_evaluable"].astype(int)
    support["n_frozen"] = support["n_frozen"].astype(int)
    if len(support) != 250 or not (support["n_frozen"] == 100).all():
        raise RuntimeError("species support denominator drift")

    minimum = int(rule["geometry_gate"]["minimum_evaluable_records_per_species"])
    records_per_species = int(rule["geometry_gate"]["retained_records_per_species"])
    if minimum != 15 or records_per_species != 15:
        raise RuntimeError("v0.13 minimum geometry drift")
    minimum_observed = int(support["n_evaluable"].min())
    species_meeting = int((support["n_evaluable"] >= minimum).sum())
    pass_gate = bool(species_meeting == 250)

    # The v0.13 architecture is justified only by the already frozen aggregate v0.12 minimum.
    expected_predecessor_minimum = int(rule["adaptive_provenance"]["predecessor_minimum_observed_evaluable_records"])
    if minimum_observed != expected_predecessor_minimum:
        raise RuntimeError("terminal artifact set does not reproduce the frozen aggregate predecessor minimum")

    retained_sha = None
    retained_rows = 0
    if pass_gate:
        seed = int(rule["geometry_gate"]["selection_seed"])
        selected_parts = []
        for _, group in joined.loc[joined["measurement_evaluable"]].groupby("species", sort=True, observed=True):
            group = group.copy()
            order = sorted(group.index.tolist(), key=lambda idx: keep_key(group.loc[idx, "source_photo_id"], seed))
            selected_parts.append(group.loc[order[:records_per_species]].copy())
        selected = pd.concat(selected_parts, ignore_index=True)
        expected_rows = int(rule["actual_geometry_architecture"]["records"])
        if len(selected) != expected_rows or selected["species"].nunique() != 250:
            raise RuntimeError("retained v0.13 geometry size mismatch")
        counts = selected.groupby("species", observed=True).size()
        if not (counts == records_per_species).all():
            raise RuntimeError("retained v0.13 geometry is not species-balanced")
        keep = [
            "species",
            "inat_taxon_id",
            "x_km",
            "y_km",
            "z_km",
            "source_observation_id",
            "source_photo_id",
            "latitude",
            "longitude",
        ]
        selected = selected[keep].sort_values(["inat_taxon_id", "source_photo_id"], kind="mergesort").reset_index(drop=True)
        args.output_geometry.parent.mkdir(parents=True, exist_ok=True)
        selected.to_csv(args.output_geometry, index=False, lineterminator="\n")
        retained_sha = sha256_file(args.output_geometry)
        retained_rows = len(selected)
    elif args.output_geometry.exists():
        args.output_geometry.unlink()

    ledger = {
        "schema": "ttf_v13_minimum_classifiability_geometry_result_v0.1",
        "status": "minimum_classifiability_geometry_frozen" if pass_gate else "minimum_classifiability_geometry_failed",
        "classifiability_gate_pass": pass_gate,
        "frozen_species": 250,
        "frozen_records": 25000,
        "minimum_evaluable_records_per_species": minimum,
        "minimum_observed_evaluable_records": minimum_observed,
        "species_meeting_minimum": species_meeting,
        "species_specific_classifiability_counts_persisted": False,
        "new_image_acquisition_performed": False,
        "species_replacement_performed": False,
        "colour_vector_values_read_by_ttf": False,
        "palette_fractions_read_by_ttf": False,
        "pairwise_colour_distances_computed": False,
        "ttf_empirical_statistic_computed": False,
        "retained_geometry": {
            "path": str(args.output_geometry) if pass_gate else None,
            "sha256": retained_sha,
            "species": 250 if pass_gate else 0,
            "records": retained_rows,
            "records_per_species": records_per_species if pass_gate else None,
            "selection_seed": int(rule["geometry_gate"]["selection_seed"]),
            "k": int(rule["actual_geometry_architecture"]["k"]) if pass_gate else None,
            "bandwidth_km": float(rule["actual_geometry_architecture"]["bandwidth_km"]) if pass_gate else None,
        },
        "firewall_sha256": sha256_file(args.firewall_dir / "firewall_manifest.json"),
        "next_gate": "Run the separately authorized frozen v0.13 synthetic requalification on this exact retained geometry before any empirical colour TTF statistic." if pass_gate else "Stop empirical TTF.",
        "claim_boundary": "This result uses only evaluability plus pre-existing source/species/coordinate metadata. No colour-vector value, palette fraction, colour distance, edge turnover, empirical TTF statistic, shared-boundary map, or environmental association was inspected.",
    }
    args.output_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.output_ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "classifiability_gate_pass": pass_gate,
        "species_meeting_minimum": species_meeting,
        "minimum_observed_evaluable_records": minimum_observed,
        "retained_records": retained_rows,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
