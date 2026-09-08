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
    if rule.get("schema") != "ttf_v12_classifiability_actual_geometry_rule_v0.1":
        raise RuntimeError("wrong v0.12 rule")
    if fw.get("status") != "actual_geometry_measurement_firewall_frozen_before_pixels":
        raise RuntimeError("wrong v0.12 firewall")
    if fw.get("colour_values_opened") is not False or fw.get("pixels_opened_at_firewall_freeze") is not False:
        raise RuntimeError("v0.12 firewall already opened forbidden outcomes")

    csvs = sorted(args.results_dir.glob("**/partition_s*_p*.csv"))
    receipts = sorted(args.results_dir.glob("**/partition_s*_p*.json"))
    if len(csvs) != 128 or len(receipts) != 128:
        raise RuntimeError(f"requires all 128 terminal partitions; got csv={len(csvs)} json={len(receipts)}")
    frames = []
    for path in csvs:
        # Intentionally read only the two fields permitted by the v0.12 firewall.
        frame = pd.read_csv(path, usecols=["measurement_id", "measurement_status"], dtype=str).fillna("")
        frames.append(frame)
    terminal = pd.concat(frames, ignore_index=True)
    if len(terminal) != 25000 or terminal["measurement_id"].nunique() != 25000:
        raise RuntimeError("terminal coverage is not exactly 25,000 unique measurement IDs")
    if not set(terminal["measurement_status"]).issubset(TERMINAL):
        raise RuntimeError("unknown terminal measurement status")

    worker = pd.read_csv(args.firewall_dir / "worker_manifest.csv", dtype=str).fillna("")
    meta = pd.read_csv(args.firewall_dir / "metadata_join_key.csv", dtype={"measurement_id": str, "source_photo_id": str, "source_observation_id": str})
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
    minimum = int(rule["classifiability_gate"]["minimum_evaluable_records_per_species"])
    pass_gate = bool((support["n_evaluable"] >= minimum).all())
    failing = support.loc[support["n_evaluable"] < minimum, ["species", "n_evaluable", "n_frozen"]].sort_values(["n_evaluable", "species"]).to_dict(orient="records")

    retained_sha = None
    retained_rows = 0
    if pass_gate:
        seed = int(rule["classifiability_gate"]["selection_seed"])
        selected_parts = []
        for species, group in joined.loc[joined["measurement_evaluable"]].groupby("species", sort=True, observed=True):
            group = group.copy()
            order = sorted(group.index.tolist(), key=lambda idx: keep_key(group.loc[idx, "source_photo_id"], seed))
            selected_parts.append(group.loc[order[:40]].copy())
        selected = pd.concat(selected_parts, ignore_index=True)
        if len(selected) != 10000 or selected["species"].nunique() != 250:
            raise RuntimeError("retained actual geometry is not exactly 250 x 40")
        counts = selected.groupby("species", observed=True).size()
        if not (counts == 40).all():
            raise RuntimeError("retained geometry does not have exactly 40 rows per species")
        # Do not persist measurement status or any colour-derived measurement fields.
        keep = [
            "species", "inat_taxon_id", "x_km", "y_km", "z_km",
            "source_observation_id", "source_photo_id", "latitude", "longitude",
        ]
        selected = selected[keep].sort_values(["inat_taxon_id", "source_photo_id"], kind="mergesort").reset_index(drop=True)
        args.output_geometry.parent.mkdir(parents=True, exist_ok=True)
        selected.to_csv(args.output_geometry, index=False, lineterminator="\n")
        retained_sha = sha256_file(args.output_geometry)
        retained_rows = len(selected)
    elif args.output_geometry.exists():
        args.output_geometry.unlink()

    status_counts = terminal["measurement_status"].value_counts().to_dict()
    ledger = {
        "schema": "ttf_v12_classifiability_actual_geometry_result_v0.1",
        "status": "actual_geometry_frozen_after_classifiability_gate" if pass_gate else "not_evaluable_actual_geometry_classifiability_gate_failed",
        "classifiability_gate_pass": pass_gate,
        "frozen_species": 250,
        "frozen_records": 25000,
        "minimum_evaluable_records_per_species": minimum,
        "minimum_observed_evaluable_records": int(support["n_evaluable"].min()),
        "median_evaluable_records": float(support["n_evaluable"].median()),
        "species_meeting_minimum": int((support["n_evaluable"] >= minimum).sum()),
        "failing_species": failing,
        "terminal_status_counts": {str(k): int(v) for k, v in status_counts.items()},
        "measurement_evaluable_terminal_status": "classified_four_state_morph",
        "no_failure_replacement": True,
        "additional_photos_after_evaluability_known": False,
        "colour_vector_values_read_by_ttf": False,
        "pairwise_colour_distances_computed": False,
        "ttf_empirical_statistic_computed": False,
        "retained_geometry": {
            "path": str(args.output_geometry) if pass_gate else None,
            "sha256": retained_sha,
            "species": 250 if pass_gate else 0,
            "records": retained_rows,
            "records_per_species": 40 if pass_gate else None,
            "selection_seed": int(rule["classifiability_gate"]["selection_seed"]),
            "k": int(rule["actual_geometry_architecture"]["k"]) if pass_gate else None,
            "bandwidth_km": float(rule["actual_geometry_architecture"]["bandwidth_km"]) if pass_gate else None,
        },
        "firewall_sha256": sha256_file(args.firewall_dir / "firewall_manifest.json"),
        "next_gate": "Run the frozen v0.12 synthetic requalification on this exact retained geometry before any empirical colour TTF statistic." if pass_gate else "Stop empirical TTF without replacement, target relaxation, or tuning.",
        "claim_boundary": "This result contains photo-level evaluability only. No colour-vector value, colour distance, edge turnover, TTF transfer statistic, shared-boundary map, or environment association was inspected.",
    }
    args.output_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.output_ledger.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"classifiability_gate_pass": pass_gate, "species_meeting_minimum": ledger["species_meeting_minimum"], "minimum_evaluable": ledger["minimum_observed_evaluable_records"], "retained_records": retained_rows}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
