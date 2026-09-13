#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ttf.genetic_geometry import prepare_genetic_sampling_geometry
from ttf.geometry import SpeciesGeometry, geometry_fingerprint
from ttf.phylogatr_character_mask import (
    CharacterMaskError,
    canonical_mask_sha256,
    edge_mask_support,
    masks_by_frozen_locality,
    read_canonical_mask_alignment,
)
from ttf.phylogatr_confirmatory import (
    fasta_header_sha256,
    fasta_headers_only,
    read_occurrence_rows,
    sha256_path,
)


def _load_json(path: Path, schema: str) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _assert_firewall(rule: dict) -> None:
    firewall = rule.get("outcome_firewall")
    if not isinstance(firewall, dict) or not firewall:
        raise RuntimeError("phase-2 outcome firewall missing")
    if any(value is not False for value in firewall.values()):
        raise RuntimeError("phase-2 outcome firewall is open")


def _read_phase1_csv(path: Path) -> tuple[list[str], dict[str, list[dict[str, str]]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or ())
        required = {
            "species",
            "locality_index",
            "latitude",
            "longitude",
            "x_km",
            "y_km",
            "z_km",
            "graph_k",
        }
        if not required.issubset(fieldnames):
            raise RuntimeError(f"phase-1 geometry CSV missing required fields: {sorted(required-set(fieldnames))}")
        grouped: dict[str, list[dict[str, str]]] = {}
        for row in reader:
            grouped.setdefault(str(row["species"]), []).append(dict(row))
    for species, rows in grouped.items():
        rows.sort(key=lambda row: int(row["locality_index"]))
        if [int(row["locality_index"]) for row in rows] != list(range(len(rows))):
            raise RuntimeError(f"noncanonical locality indices in phase-1 CSV for {species}")
    return fieldnames, grouped


def _frozen_arrays(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray, int]:
    latlon = np.asarray(
        [[float(row["latitude"]), float(row["longitude"])] for row in rows],
        dtype=float,
    )
    xyz = np.asarray(
        [[float(row["x_km"]), float(row["y_km"]), float(row["z_km"])] for row in rows],
        dtype=float,
    )
    graph_k_values = {int(row["graph_k"]) for row in rows}
    if len(graph_k_values) != 1:
        raise RuntimeError("graph_k varies within a frozen species geometry")
    if not np.isfinite(latlon).all() or not np.isfinite(xyz).all():
        raise RuntimeError("non-finite phase-1 geometry")
    return latlon, xyz, int(next(iter(graph_k_values)))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Freeze the response-blind canonical-valid character-mask survivor geometry."
    )
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--phase1-manifest", type=Path, required=True)
    ap.add_argument("--phase1-csv", type=Path, required=True)
    ap.add_argument("--phase2-rule", type=Path, required=True)
    ap.add_argument("--output-csv", type=Path, required=True)
    ap.add_argument("--output-manifest", type=Path, required=True)
    args = ap.parse_args()

    phase1 = _load_json(
        args.phase1_manifest,
        "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1",
    )
    rule = _load_json(args.phase2_rule, "ttf_genetic_phylogatr_phase2_mask_rule_v0.1")
    _assert_firewall(rule)
    required_status = str(rule["source_integrity"]["require_phase1_status"])
    if phase1.get("status") != required_status:
        raise RuntimeError(
            f"phase-2 mask opening requires phase-1 status {required_status!r}, got {phase1.get('status')!r}"
        )
    if phase1.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("phase-1 manifest indicates nucleotide identity was opened")
    if phase1.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("phase-1 manifest indicates genetic distances were opened")
    if sha256_path(args.phase1_csv) != phase1["geometry_csv_sha256"]:
        raise RuntimeError("phase-1 geometry CSV SHA256 drift")

    root = args.root.resolve()
    if sha256_path(root / "cite.txt") != phase1["provenance"]["cite_sha256"]:
        raise RuntimeError("cite.txt drift since phase-1 freeze")
    if sha256_path(root / "genes.txt") != phase1["provenance"]["genes_sha256"]:
        raise RuntimeError("genes.txt drift since phase-1 freeze")

    fieldnames, csv_by_species = _read_phase1_csv(args.phase1_csv)
    panel_meta: dict[str, dict] = {
        str(name): dict(meta) for name, meta in phase1["selected_panels"].items()
    }
    if set(csv_by_species) != set(panel_meta):
        raise RuntimeError("phase-1 CSV species and manifest selected_panels disagree")

    fraction = float(rule["edge_admissibility"]["minimum_comparable_fraction"])
    if fraction != 0.50:
        raise RuntimeError("phase-2 comparable-site fraction drift")

    ledger: dict[str, dict] = {}
    survivors: list[str] = []
    survivor_geometry: dict[str, SpeciesGeometry] = {}

    for species in sorted(panel_meta):
        meta = panel_meta[species]
        rows = csv_by_species[species]
        fasta_path = root / str(meta["aligned_fasta_relative"])
        occurrence_path = root / str(meta["occurrence_relative"])
        if not fasta_path.is_file() or not occurrence_path.is_file():
            raise RuntimeError(f"phase-2 source file missing for {species}")

        # Source integrity is a global STOP condition, not a biological negative.
        if sha256_path(occurrence_path) != str(meta["occurrence_sha256"]):
            raise RuntimeError(f"occurrence source drift since phase-1 freeze for {species}")
        header_only = fasta_headers_only(fasta_path)
        if fasta_header_sha256(header_only) != str(meta["aligned_header_sha256"]):
            raise RuntimeError(f"aligned FASTA header drift since phase-1 freeze for {species}")

        latlon, xyz, graph_k = _frozen_arrays(rows)
        if len(latlon) != int(meta["unique_localities"]):
            raise RuntimeError(f"phase-1 locality-count drift for {species}")
        geometry = prepare_genetic_sampling_geometry(xyz, k=graph_k)
        if geometry.n_localities != len(latlon):
            raise RuntimeError(f"phase-2 geometry locality drift for {species}")
        if geometry.graph_k != graph_k:
            raise RuntimeError(f"phase-2 graph-k drift for {species}")
        if geometry.n_edges != int(meta["edge_count"]):
            raise RuntimeError(f"phase-2 edge-count drift for {species}")
        if not np.array_equal(np.asarray(geometry.coordinates), xyz):
            raise RuntimeError(f"phase-2 canonical coordinate order drift for {species}")

        try:
            alignment = read_canonical_mask_alignment(fasta_path)
        except CharacterMaskError as exc:
            ledger[species] = {
                "status": "NOT_EVALUABLE_CHARACTER_SUPPORT",
                "reason": str(exc),
                "survives": False,
                "nucleotide_identity_persisted": False,
            }
            continue

        # Header hash is deliberately rechecked after mask parsing too.
        if fasta_header_sha256(alignment.headers) != str(meta["aligned_header_sha256"]):
            raise RuntimeError(f"aligned FASTA header/order drift during phase-2 parse for {species}")
        occurrence_rows = read_occurrence_rows(occurrence_path)
        grouped_masks = masks_by_frozen_locality(alignment, occurrence_rows, latlon)
        support = edge_mask_support(
            grouped_masks,
            geometry.edge_nodes,
            alignment_length=alignment.alignment_length,
            minimum_comparable_fraction=fraction,
        )
        invalid_edges = int(np.count_nonzero(~support.valid_edges))
        survives = bool(support.all_edges_valid)
        ledger[species] = {
            "status": "ADMISSIBLE_CHARACTER_SUPPORT" if survives else "NOT_EVALUABLE_CHARACTER_SUPPORT",
            "reason": None if survives else "one_or_more_frozen_edges_lack_a_valid_cross_locality_sequence_pair",
            "survives": survives,
            "alignment_length": int(alignment.alignment_length),
            "minimum_comparable_columns": int(support.minimum_comparable_columns),
            "aligned_records": int(len(alignment.headers)),
            "localities_with_one_or_more_masks": int(len(grouped_masks)),
            "frozen_edges": int(len(support.valid_edges)),
            "valid_edges": int(np.count_nonzero(support.valid_edges)),
            "invalid_edges": invalid_edges,
            "minimum_best_comparable_columns": int(np.min(support.best_comparable_columns)),
            "canonical_mask_sha256": canonical_mask_sha256(alignment),
            "nucleotide_identity_persisted": False,
        }
        if survives:
            survivors.append(species)
            survivor_geometry[species] = SpeciesGeometry(species=species, coordinates=xyz)

    survivor_set = set(survivors)
    phase1_train = tuple(map(str, phase1["split"]["train_species"]))
    phase1_eval = tuple(map(str, phase1["split"]["eval_species"]))
    train = tuple(name for name in phase1_train if name in survivor_set)
    evaluation = tuple(name for name in phase1_eval if name in survivor_set)
    if set(train) & set(evaluation):
        raise RuntimeError("phase-2 inherited train/evaluation overlap")
    if set(train) | set(evaluation) != survivor_set:
        raise RuntimeError("phase-2 survivor set is not inherited from phase-1 split")

    min_species = int(rule["survivor_panel"]["minimum_species"])
    passed = len(survivors) >= min_species and len(train) > 0 and len(evaluation) > 0
    status = "PASS_TO_SYNTHETIC_GATE" if passed else "NOT_EVALUABLE_PHASE2_PANEL_TOO_SMALL"

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for species in sorted(survivor_set):
            for row in csv_by_species[species]:
                writer.writerow(row)

    manifest = {
        "schema": "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1",
        "status": status,
        "phase1": {
            "dataset_digest_sha256": phase1["dataset_digest_sha256"],
            "geometry_fingerprint_sha256": phase1["geometry_fingerprint_sha256"],
            "geometry_csv_sha256": phase1["geometry_csv_sha256"],
            "species_count": int(phase1["census"]["final_species"]),
            "train_count": int(phase1["split"]["train_count"]),
            "eval_count": int(phase1["split"]["eval_count"]),
        },
        "rule_sha256": sha256_path(args.phase2_rule),
        "source_integrity": {
            "cite_sha256": sha256_path(root / "cite.txt"),
            "genes_sha256": sha256_path(root / "genes.txt"),
            "phase1_csv_sha256_verified": True,
            "per_species_occurrence_sha256_verified": True,
            "per_species_aligned_header_sha256_verified": True,
        },
        "mask_contract": {
            "canonical_valid_characters": ["A", "C", "G", "T", "a", "c", "g", "t"],
            "minimum_comparable_fraction": fraction,
            "nucleotide_identity_persisted": False,
            "pairwise_nucleotide_differences_computed": False,
        },
        "species": {
            "phase1": int(len(panel_meta)),
            "survivors": int(len(survivors)),
            "failed_character_support": int(len(panel_meta) - len(survivors)),
            "minimum_required": min_species,
        },
        "split": {
            "inherit_phase1_without_resplitting": True,
            "train_species": list(train),
            "eval_species": list(evaluation),
            "train_count": len(train),
            "eval_count": len(evaluation),
        },
        "geometry_csv_sha256": sha256_path(args.output_csv),
        "geometry_fingerprint_sha256": (
            geometry_fingerprint(
                [survivor_geometry[name] for name in sorted(survivor_geometry)]
            )
            if survivor_geometry
            else None
        ),
        "species_ledger": ledger,
        "character_mask_opened": True,
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "qualification_claim_made": False,
        "next_gate": (
            "repeat_full_dataset_specific_profiled_private_Gate_D_on_phase2_survivor_geometry"
            if passed
            else "STOP_DO_NOT_OPEN_NUCLEOTIDE_IDENTITIES"
        ),
    }
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": status,
                "phase1_species": len(panel_meta),
                "survivors": len(survivors),
                "failed_character_support": len(panel_meta) - len(survivors),
                "train": len(train),
                "eval": len(evaluation),
                "geometry_fingerprint_sha256": manifest["geometry_fingerprint_sha256"],
                "nucleotide_identity_opened": False,
                "pairwise_genetic_distances_opened": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
