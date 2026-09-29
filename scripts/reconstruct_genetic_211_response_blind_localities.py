#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import tempfile
import zipfile
from pathlib import Path

from ttf.phylogatr_confirmatory import (
    apply_species_cap,
    choose_one_panel_per_species,
    phase1_dataset_digest,
    scan_phylogatr_phase1,
    sha256_path,
)
from ttf.phylogatr_source_projection import project_genes_file_in_place


EXPECTED_ARCHIVE_SHA256 = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"


def _safe_extract_zip(archive: Path, destination: Path) -> None:
    with zipfile.ZipFile(archive) as handle:
        for info in handle.infolist():
            path = Path(info.filename)
            if path.is_absolute() or ".." in path.parts:
                raise RuntimeError(f"unsafe archive member: {info.filename!r}")
        handle.extractall(destination)


def _find_root(extracted: Path) -> Path:
    candidates = []
    if (extracted / "genes.txt").is_file() and (extracted / "cite.txt").is_file():
        candidates.append(extracted)
    candidates.extend(
        path.parent
        for path in extracted.rglob("genes.txt")
        if (path.parent / "cite.txt").is_file()
    )
    unique = sorted(set(candidates))
    if len(unique) != 1:
        raise RuntimeError(f"expected one phylogatR root, found {len(unique)}")
    return unique[0]


def _load_json(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--authorization", type=Path, required=True)
    ap.add_argument("--protocol", type=Path, required=True)
    ap.add_argument("--execution-rule", type=Path, required=True)
    ap.add_argument("--decker-exclusion", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    archive = args.archive.resolve()
    if sha256_path(archive) != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError("exact phylogatR archive SHA256 mismatch")

    authorization = _load_json(args.authorization)
    if authorization.get("schema") != "ttf_genetic_phylogatr_phase4_identity_opening_authorization_v0.1":
        raise RuntimeError("unexpected Phase-4 authorization schema")
    survivor_block = authorization.get("species", {})
    survivor_names = sorted(
        set(map(str, survivor_block.get("train_species", [])))
        | set(map(str, survivor_block.get("eval_species", [])))
    )
    if len(survivor_names) != 211 or int(survivor_block.get("survivors", -1)) != 211:
        raise RuntimeError("authorization does not define the exact 211 survivors")

    protocol = _load_json(args.protocol)
    execution = _load_json(args.execution_rule)
    exclusion = _load_json(args.decker_exclusion)
    aliases = tuple(protocol["marker_contract"]["normalized_aliases"])
    phase1 = execution["phase1"]
    excluded = tuple(map(str, exclusion["species"]))

    with tempfile.TemporaryDirectory(prefix="ttf_palearctic_phase1_") as td:
        extracted = Path(td)
        _safe_extract_zip(archive, extracted)
        root = _find_root(extracted)

        projection = project_genes_file_in_place(root)
        if projection["projected_genes_sha256"] != "3a09b31e466985e1120cef94f5805f5bb23ebc76b4f3b650f3e9188239415e0e":
            raise RuntimeError("Animalia projected genes SHA256 drift")

        scan = scan_phylogatr_phase1(
            root,
            aliases=aliases,
            excluded_species=excluded,
            min_localities=int(phase1["minimum_unique_localities"]),
            min_endpoint_training_edges=int(
                phase1["minimum_endpoint_disjoint_ibd_training_edges"]
            ),
            neighbor_fraction=float(phase1["neighbor_fraction"]),
        )
        one_panel = choose_one_panel_per_species(scan.candidates)
        dataset_digest, _ = phase1_dataset_digest(root, one_panel)
        if dataset_digest != authorization["dataset_digest_sha256"]:
            raise RuntimeError("response-blind Phase-1 dataset digest drift")

        selected = apply_species_cap(
            one_panel,
            dataset_digest=dataset_digest,
            maximum_species=int(phase1["maximum_species"]),
        )
        selected_by_species = {panel.species: panel for panel in selected}
        if len(selected_by_species) != 250:
            raise RuntimeError("response-blind Phase-1 panel did not reconstruct to 250 species")
        if not set(survivor_names).issubset(selected_by_species):
            missing = sorted(set(survivor_names) - set(selected_by_species))
            raise RuntimeError(f"authorized survivors missing from reconstructed Phase 1: {missing}")

        species = {}
        for name in survivor_names:
            panel = selected_by_species[name]
            species[name] = {
                "taxonomy": {
                    "kingdom": panel.kingdom,
                    "phylum": panel.phylum,
                    "class": panel.class_name,
                    "order": panel.order,
                    "family": panel.family,
                    "genus": panel.genus,
                },
                "selected_gene": panel.raw_gene,
                "localities": [
                    [float(lat), float(lon)]
                    for lat, lon in panel.canonical_latlon
                ],
            }

    payload = {
        "schema": "ttf_genetic_211_response_blind_selected_panel_localities_v0.1",
        "status": "RECONSTRUCTED_FROM_EXACT_ARCHIVE_WITHOUT_GENETIC_RESPONSE",
        "source_archive_sha256": EXPECTED_ARCHIVE_SHA256,
        "dataset_digest_sha256": authorization["dataset_digest_sha256"],
        "species_count": len(species),
        "insecta_species": sum(
            row["taxonomy"]["class"] == "Insecta" for row in species.values()
        ),
        "lepidoptera_species": sum(
            row["taxonomy"]["class"] == "Insecta"
            and row["taxonomy"]["order"] == "Lepidoptera"
            for row in species.values()
        ),
        "species": species,
        "outcome_firewall": {
            "nucleotide_identity_used": False,
            "pairwise_genetic_distances_used": False,
            "species_level_phase4_scores_used": False,
        },
        "authority": "subgroup_design_input_only",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "species": payload["species_count"],
        "insecta": payload["insecta_species"],
        "lepidoptera": payload["lepidoptera_species"],
        "genetic_response_used": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
