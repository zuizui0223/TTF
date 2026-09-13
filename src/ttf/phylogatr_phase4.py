from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry
from .genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from .geometry import SpeciesGeometry, geometry_fingerprint


@dataclass(frozen=True)
class PhylogatrPhase4Context:
    geometries: dict[str, GeneticSamplingGeometry]
    frozen_latlon: dict[str, np.ndarray]
    train_species: tuple[str, ...]
    eval_species: tuple[str, ...]
    phase1: dict
    phase2: dict
    phase3_rule: dict
    phase3_authorization: dict
    references: dict
    qualification: dict
    phase4_rule: dict
    authorization: dict


def _load_json(path: Path, schema: str) -> dict:
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _assert_false_mapping(payload: Mapping[str, object], *, label: str) -> None:
    if not payload:
        raise RuntimeError(f"missing {label}")
    bad = {str(key): value for key, value in payload.items() if value is not False}
    if bad:
        raise RuntimeError(f"{label} is open: {bad}")


def _load_latlon_by_species(path: Path) -> dict[str, np.ndarray]:
    grouped: dict[str, list[tuple[int, float, float]]] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"species", "locality_index", "latitude", "longitude"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(f"Phase-2 geometry CSV lacks lat/lon columns: {sorted(missing)}")
        for row in reader:
            species = str(row["species"])
            index = int(row["locality_index"])
            lat = float(row["latitude"])
            lon = float(row["longitude"])
            if not np.isfinite(lat) or not np.isfinite(lon):
                raise RuntimeError(f"non-finite frozen lat/lon for {species}")
            grouped.setdefault(species, []).append((index, lat, lon))
    out: dict[str, np.ndarray] = {}
    for species in sorted(grouped):
        rows = sorted(grouped[species], key=lambda value: value[0])
        if [row[0] for row in rows] != list(range(len(rows))):
            raise RuntimeError(f"noncanonical frozen locality indices for {species}")
        out[species] = np.asarray([[row[1], row[2]] for row in rows], dtype=float)
    return out


def verify_phase4_authorized_code_files(
    authorization: Mapping[str, object], *, repo_root: Path = Path(".")
) -> None:
    frozen = authorization.get("frozen_code_sha256")
    if not isinstance(frozen, dict) or not frozen:
        raise RuntimeError("Phase-4 authorization lacks frozen_code_sha256")
    for relative, expected in sorted(frozen.items()):
        path = Path(repo_root) / str(relative)
        if not path.is_file():
            raise RuntimeError(f"authorized Phase-4 code file missing: {relative}")
        if sha256_path(path) != str(expected):
            raise RuntimeError(f"authorized Phase-4 code file drift: {relative}")


def load_phylogatr_phase4_context(
    geometry_csv: Path,
    phase1_manifest_path: Path,
    phase2_manifest_path: Path,
    phase3_rule_path: Path,
    phase3_authorization_path: Path,
    references_path: Path,
    qualification_path: Path,
    phase4_rule_path: Path,
    phase4_authorization_path: Path,
    opening_state_path: Path,
    *,
    verify_code: bool = True,
    repo_root: Path = Path("."),
) -> PhylogatrPhase4Context:
    phase1 = _load_json(
        phase1_manifest_path, "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1"
    )
    phase2 = _load_json(
        phase2_manifest_path, "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1"
    )
    phase3_rule = _load_json(
        phase3_rule_path, "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1"
    )
    phase3_auth = _load_json(
        phase3_authorization_path,
        "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1",
    )
    references = _load_json(
        references_path, "ttf_genetic_phylogatr_phase3_references_v0.1"
    )
    qualification = _load_json(
        qualification_path, "ttf_genetic_phylogatr_phase3_qualification_v0.1"
    )
    phase4_rule = _load_json(
        phase4_rule_path, "ttf_genetic_phylogatr_phase4_response_rule_v0.1"
    )
    authorization = _load_json(
        phase4_authorization_path,
        "ttf_genetic_phylogatr_phase4_identity_opening_authorization_v0.1",
    )
    opening_state = _load_json(
        opening_state_path, "ttf_genetic_empirical_opening_state_v0.1"
    )

    if authorization.get("status") != "AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING":
        raise RuntimeError("fresh Phase-4 identity opening is not authorized")
    _assert_false_mapping(
        phase4_rule["outcome_firewall_at_rule_freeze"],
        label="Phase-4 rule freeze firewall",
    )
    if opening_state.get("development_panel", {}).get("opening_decision") != "PERMANENTLY_CLOSED_UNDER_V0_2_DESIGN":
        raise RuntimeError("Decker development opening state drift")
    _assert_false_mapping(opening_state["global_firewall"], label="canonical genetic opening-state firewall")

    hashes = {
        "phase1_manifest_sha256": sha256_path(phase1_manifest_path),
        "phase2_manifest_sha256": sha256_path(phase2_manifest_path),
        "phase3_rule_sha256": sha256_path(phase3_rule_path),
        "phase3_authorization_sha256": sha256_path(phase3_authorization_path),
        "phase3_references_sha256": sha256_path(references_path),
        "phase3_qualification_sha256": sha256_path(qualification_path),
        "phase4_rule_sha256": sha256_path(phase4_rule_path),
        "opening_state_sha256": sha256_path(opening_state_path),
        "geometry_csv_sha256": sha256_path(geometry_csv),
    }
    for key, observed in hashes.items():
        if observed != authorization[key]:
            raise RuntimeError(f"Phase-4 authorization source drift: {key}")

    if phase3_auth.get("phase1_manifest_sha256") != hashes["phase1_manifest_sha256"]:
        raise RuntimeError("Phase-3/Phase-4 Phase-1 manifest drift")
    if phase3_auth.get("phase2_manifest_sha256") != hashes["phase2_manifest_sha256"]:
        raise RuntimeError("Phase-3/Phase-4 Phase-2 manifest drift")
    if phase3_auth.get("phase3_rule_sha256") != hashes["phase3_rule_sha256"]:
        raise RuntimeError("Phase-3/Phase-4 rule drift")
    if phase3_auth.get("geometry_csv_sha256") != hashes["geometry_csv_sha256"]:
        raise RuntimeError("Phase-3/Phase-4 geometry CSV drift")

    if qualification.get("status") != "PASS" or qualification.get("passed") is not True:
        raise RuntimeError("fresh Phase-3 qualification did not PASS")
    if qualification.get("phase4_identity_opening_eligible") is not True:
        raise RuntimeError("fresh Phase-3 qualification does not authorize Phase 4")
    if qualification.get("type1_gate", {}).get("pass") is not True:
        raise RuntimeError("fresh Phase-3 Type-I gate did not PASS")
    if qualification.get("power_gate", {}).get("pass") is not True:
        raise RuntimeError("fresh Phase-3 power gate did not PASS")
    if qualification.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("fresh identity was already opened before Phase-4 authorization")
    if qualification.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("fresh genetic distances were already opened before Phase 4")
    if qualification.get("confirmatory_ttf_statistic_opened") is not False:
        raise RuntimeError("fresh empirical TTF was already opened before Phase 4")

    fingerprint = str(phase3_auth["geometry_fingerprint_sha256"])
    if references.get("status") != "ordered_private_reference_families_complete":
        raise RuntimeError("fresh Phase-3 reference families are incomplete")
    if references.get("geometry_fingerprint_sha256") != fingerprint:
        raise RuntimeError("fresh Phase-3 references use a different geometry")
    if qualification.get("geometry_fingerprint_sha256") != fingerprint:
        raise RuntimeError("fresh Phase-3 qualification uses a different geometry")
    if references.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("fresh Phase-3 references indicate identity opening")
    if references.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("fresh Phase-3 references indicate distance opening")

    table = load_frozen_genetic_geometry_csv(
        geometry_csv,
        expected_sha256=phase3_auth["geometry_csv_sha256"],
        expected_neighbor_fraction=float(phase3_rule["geometry_contract"]["neighbor_fraction"]),
    )
    reconstructed = geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=table.geometries[name].coordinates)
            for name in table.species
        ]
    )
    if reconstructed != fingerprint:
        raise RuntimeError("fresh Phase-4 reconstructed geometry fingerprint drift")

    train = tuple(map(str, phase3_auth["species"]["train_species"]))
    evaluation = tuple(map(str, phase3_auth["species"]["eval_species"]))
    if set(train) | set(evaluation) != set(table.species) or set(train) & set(evaluation):
        raise RuntimeError("fresh Phase-4 split/geometry drift")
    if set(phase2["split"]["train_species"]) != set(train):
        raise RuntimeError("fresh Phase-4 training split differs from Phase 2")
    if set(phase2["split"]["eval_species"]) != set(evaluation):
        raise RuntimeError("fresh Phase-4 evaluation split differs from Phase 2")

    selected = set(map(str, phase1.get("selected_panels", {}).keys()))
    if not set(table.species).issubset(selected):
        raise RuntimeError("fresh Phase-4 survivor species are absent from Phase-1 selected panels")

    latlon = _load_latlon_by_species(geometry_csv)
    if set(latlon) != set(table.species):
        raise RuntimeError("fresh Phase-4 lat/lon species set drift")
    for species in table.species:
        if len(latlon[species]) != table.geometries[species].n_localities:
            raise RuntimeError(f"fresh Phase-4 lat/lon locality-count drift for {species}")

    if verify_code:
        verify_phase4_authorized_code_files(authorization, repo_root=repo_root)

    return PhylogatrPhase4Context(
        geometries=table.geometries,
        frozen_latlon=latlon,
        train_species=train,
        eval_species=evaluation,
        phase1=phase1,
        phase2=phase2,
        phase3_rule=phase3_rule,
        phase3_authorization=phase3_auth,
        references=references,
        qualification=qualification,
        phase4_rule=phase4_rule,
        authorization=authorization,
    )


__all__ = [
    "PhylogatrPhase4Context",
    "load_phylogatr_phase4_context",
    "verify_phase4_authorized_code_files",
]
