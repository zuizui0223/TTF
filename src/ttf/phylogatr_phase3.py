from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Mapping

from .genetic_geometry import GeneticSamplingGeometry
from .genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from .geometry import SpeciesGeometry, geometry_fingerprint


_PHASE2_SCHEMA = "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1"
_RULE_SCHEMA = "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1"
_AUTH_SCHEMA = "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1"


@dataclass(frozen=True)
class PhylogatrPhase3Context:
    geometries: dict[str, GeneticSamplingGeometry]
    train_species: tuple[str, ...]
    eval_species: tuple[str, ...]
    master_seed: int
    phase2: dict
    rule: dict
    authorization: dict


def derive_phase3_master_seed(dataset_digest_sha256: str, geometry_fingerprint_sha256: str) -> int:
    payload = (
        "ttf|phylogatr|phase3|"
        + str(dataset_digest_sha256)
        + "|"
        + str(geometry_fingerprint_sha256)
    ).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


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


def verify_authorized_code_files(authorization: Mapping[str, object], *, repo_root: Path = Path(".")) -> None:
    frozen = authorization.get("frozen_code_sha256")
    if not isinstance(frozen, dict) or not frozen:
        raise RuntimeError("phase-3 authorization lacks frozen_code_sha256")
    for relative, expected in sorted(frozen.items()):
        path = Path(repo_root) / str(relative)
        if not path.is_file():
            raise RuntimeError(f"authorized code file missing: {relative}")
        observed = sha256_path(path)
        if observed != str(expected):
            raise RuntimeError(f"authorized code file drift: {relative}")


def load_phylogatr_phase3_context(
    geometry_csv: Path,
    phase2_manifest_path: Path,
    phase3_rule_path: Path,
    authorization_path: Path,
    *,
    verify_code: bool = True,
    repo_root: Path = Path("."),
) -> PhylogatrPhase3Context:
    phase2 = _load_json(phase2_manifest_path, _PHASE2_SCHEMA)
    rule = _load_json(phase3_rule_path, _RULE_SCHEMA)
    authorization = _load_json(authorization_path, _AUTH_SCHEMA)

    if phase2.get("status") != rule["geometry_contract"]["required_phase2_status"]:
        raise RuntimeError("fresh phase-3 requires a PASS_TO_SYNTHETIC_GATE phase-2 manifest")
    if authorization.get("status") != "authorize_frozen_fresh_phylogatr_phase3_gate_d":
        raise RuntimeError("fresh phase-3 Gate-D is not authorized")
    if phase2.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("fresh nucleotide identity is already open")
    if phase2.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("fresh pairwise genetic distances are already open")
    if phase2.get("confirmatory_ttf_statistic_opened") is not False:
        raise RuntimeError("fresh empirical TTF statistic is already open")
    _assert_false_mapping(rule["outcome_firewall"], label="phase-3 rule outcome firewall")
    _assert_false_mapping(authorization["outcome_firewall"], label="phase-3 authorization outcome firewall")

    if sha256_path(phase2_manifest_path) != authorization["phase2_manifest_sha256"]:
        raise RuntimeError("phase-2 manifest SHA256 drift after authorization")
    if sha256_path(phase3_rule_path) != authorization["phase3_rule_sha256"]:
        raise RuntimeError("phase-3 rule SHA256 drift after authorization")
    if sha256_path(geometry_csv) != authorization["geometry_csv_sha256"]:
        raise RuntimeError("phase-2 survivor geometry CSV SHA256 drift after authorization")
    if authorization["geometry_csv_sha256"] != phase2["geometry_csv_sha256"]:
        raise RuntimeError("authorization/phase-2 geometry CSV digest mismatch")

    table = load_frozen_genetic_geometry_csv(
        geometry_csv,
        expected_sha256=phase2["geometry_csv_sha256"],
        expected_neighbor_fraction=float(rule["geometry_contract"]["neighbor_fraction"]),
    )
    species = set(table.species)
    train = tuple(map(str, phase2["split"]["train_species"]))
    evaluation = tuple(map(str, phase2["split"]["eval_species"]))
    if not train or not evaluation or set(train) & set(evaluation):
        raise RuntimeError("phase-2 split is empty or overlapping")
    if set(train) | set(evaluation) != species:
        raise RuntimeError("phase-2 split does not exactly cover survivor geometry")
    if len(species) != int(phase2["species"]["survivors"]):
        raise RuntimeError("phase-2 survivor count disagrees with geometry CSV")
    if len(species) < int(rule["geometry_contract"]["minimum_surviving_species"]):
        raise RuntimeError("phase-2 survivor panel is below the frozen phase-3 minimum")

    fingerprint = geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=table.geometries[name].coordinates)
            for name in sorted(table.geometries)
        ]
    )
    if fingerprint != phase2["geometry_fingerprint_sha256"]:
        raise RuntimeError("reconstructed fresh phase-2 geometry fingerprint drift")
    if fingerprint != authorization["geometry_fingerprint_sha256"]:
        raise RuntimeError("authorization geometry fingerprint drift")

    dataset_digest = str(phase2["phase1"]["dataset_digest_sha256"])
    master_seed = derive_phase3_master_seed(dataset_digest, fingerprint)
    if master_seed != int(authorization["master_seed"]):
        raise RuntimeError("fresh phase-3 deterministic master seed drift")

    if verify_code:
        verify_authorized_code_files(authorization, repo_root=repo_root)

    return PhylogatrPhase3Context(
        geometries=table.geometries,
        train_species=train,
        eval_species=evaluation,
        master_seed=master_seed,
        phase2=phase2,
        rule=rule,
        authorization=authorization,
    )


__all__ = [
    "PhylogatrPhase3Context",
    "derive_phase3_master_seed",
    "load_phylogatr_phase3_context",
    "verify_authorized_code_files",
]
