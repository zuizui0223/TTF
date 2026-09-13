#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from ttf.geometry import SpeciesGeometry, geometry_fingerprint


PHASE4_CODE_PATHS = (
    "src/ttf/phylogatr_phase4.py",
    "src/ttf/phylogatr_empirical.py",
    "src/ttf/genetic_empirical_score.py",
    "src/ttf/phylogatr_character_mask.py",
    "src/ttf/phylogatr_confirmatory.py",
    "src/ttf/genetic_self_detectability.py",
    "src/ttf/private_null_inference.py",
    "src/ttf/profiled_private_null.py",
    "scripts/authorize_phylogatr_phase4_identity_opening.py",
    "scripts/run_phylogatr_phase4_empirical_test.py",
)


def _load(path: Path, schema: str) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _closed(payload: dict, label: str) -> None:
    if not payload or any(value is not False for value in payload.values()):
        raise RuntimeError(f"{label} is not fully closed")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Authorize exact fresh phylogatR nucleotide-identity opening after Phase-3 qualification."
    )
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase1-manifest", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--phase3-authorization", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--qualification", type=Path, required=True)
    ap.add_argument("--self-rule", type=Path, required=True)
    ap.add_argument("--self-references", type=Path, required=True)
    ap.add_argument("--self-qualification", type=Path, required=True)
    ap.add_argument("--phase4-rule", type=Path, required=True)
    ap.add_argument("--opening-state", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    phase1 = _load(
        args.phase1_manifest, "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1"
    )
    phase2 = _load(
        args.phase2_manifest, "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1"
    )
    phase3_rule = _load(
        args.phase3_rule, "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1"
    )
    phase3_auth = _load(
        args.phase3_authorization, "ttf_genetic_phylogatr_phase3_gate_d_authorization_v0.1"
    )
    refs = _load(args.references, "ttf_genetic_phylogatr_phase3_references_v0.1")
    qualification = _load(
        args.qualification, "ttf_genetic_phylogatr_phase3_qualification_v0.1"
    )
    self_rule = _load(
        args.self_rule, "ttf_genetic_phylogatr_phase3_self_detectability_rule_v0.1"
    )
    self_refs = _load(
        args.self_references, "ttf_genetic_phylogatr_phase3_self_references_v0.1"
    )
    self_qualification = _load(
        args.self_qualification, "ttf_genetic_phylogatr_phase3_self_qualification_v0.1"
    )
    phase4_rule = _load(
        args.phase4_rule, "ttf_genetic_phylogatr_phase4_response_rule_v0.1"
    )
    opening = _load(args.opening_state, "ttf_genetic_empirical_opening_state_v0.1")

    _closed(phase4_rule["outcome_firewall_at_rule_freeze"], "Phase-4 rule firewall")
    _closed(self_rule["outcome_firewall"], "fresh self-rule firewall")
    _closed(opening["global_firewall"], "canonical opening-state firewall")
    if opening.get("development_panel", {}).get("opening_decision") != "PERMANENTLY_CLOSED_UNDER_V0_2_DESIGN":
        raise RuntimeError("Decker development outcome opening state drift")

    if qualification.get("status") != "PASS" or qualification.get("passed") is not True:
        raise RuntimeError("fresh Phase-3 Gate-D did not PASS")
    if qualification.get("phase4_identity_opening_eligible") is not True:
        raise RuntimeError("fresh Phase-3 receipt does not authorize Phase 4")
    if qualification.get("type1_gate", {}).get("pass") is not True:
        raise RuntimeError("fresh Phase-3 Type-I gate did not PASS")
    if qualification.get("power_gate", {}).get("pass") is not True:
        raise RuntimeError("fresh Phase-3 power gate did not PASS")
    for key in (
        "confirmatory_sequence_identity_opened",
        "confirmatory_pairwise_genetic_distances_opened",
        "confirmatory_ttf_statistic_opened",
    ):
        if qualification.get(key) is not False:
            raise RuntimeError(f"fresh empirical firewall already open before authorization: {key}")

    if refs.get("status") != "ordered_private_reference_families_complete":
        raise RuntimeError("fresh Phase-3 references are incomplete")
    if self_refs.get("status") != "complete_independent_null_reference":
        raise RuntimeError("fresh self references are incomplete")
    if self_qualification.get("status") not in {"PASS", "SELF_DETECTABILITY_NOT_QUALIFIED"}:
        raise RuntimeError("fresh self qualification is incomplete")
    if self_qualification.get("passed") not in {True, False}:
        raise RuntimeError("fresh self qualification lacks a final boolean decision")
    for payload, label in ((refs, "references"), (self_refs, "self references"), (self_qualification, "self qualification")):
        if payload.get("confirmatory_sequence_identity_opened") is not False:
            raise RuntimeError(f"fresh identity was opened in {label}")
        if payload.get("confirmatory_pairwise_genetic_distances_opened") is not False:
            raise RuntimeError(f"fresh genetic distances were opened in {label}")

    hashes = {
        "phase1_manifest_sha256": sha256_path(args.phase1_manifest),
        "phase2_manifest_sha256": sha256_path(args.phase2_manifest),
        "phase3_rule_sha256": sha256_path(args.phase3_rule),
        "phase3_authorization_sha256": sha256_path(args.phase3_authorization),
        "phase3_references_sha256": sha256_path(args.references),
        "phase3_qualification_sha256": sha256_path(args.qualification),
        "phase3_self_rule_sha256": sha256_path(args.self_rule),
        "phase3_self_references_sha256": sha256_path(args.self_references),
        "phase3_self_qualification_sha256": sha256_path(args.self_qualification),
        "phase4_rule_sha256": sha256_path(args.phase4_rule),
        "opening_state_sha256": sha256_path(args.opening_state),
        "geometry_csv_sha256": sha256_path(args.geometry),
    }
    if phase3_auth.get("phase1_manifest_sha256") != hashes["phase1_manifest_sha256"]:
        raise RuntimeError("Phase-1 manifest differs from Phase-3 authorization")
    if phase3_auth.get("phase2_manifest_sha256") != hashes["phase2_manifest_sha256"]:
        raise RuntimeError("Phase-2 manifest differs from Phase-3 authorization")
    if phase3_auth.get("phase3_rule_sha256") != hashes["phase3_rule_sha256"]:
        raise RuntimeError("Phase-3 rule differs from Phase-3 authorization")
    if phase3_auth.get("geometry_csv_sha256") != hashes["geometry_csv_sha256"]:
        raise RuntimeError("geometry CSV differs from Phase-3 authorization")

    fingerprint = str(phase3_auth["geometry_fingerprint_sha256"])
    if refs.get("geometry_fingerprint_sha256") != fingerprint:
        raise RuntimeError("Phase-3 reference geometry drift")
    if qualification.get("geometry_fingerprint_sha256") != fingerprint:
        raise RuntimeError("Phase-3 qualification geometry drift")
    if self_refs.get("geometry_fingerprint_sha256") != fingerprint:
        raise RuntimeError("fresh self-reference geometry drift")
    if self_qualification.get("geometry_fingerprint_sha256") != fingerprint:
        raise RuntimeError("fresh self-qualification geometry drift")

    table = load_frozen_genetic_geometry_csv(
        args.geometry,
        expected_sha256=hashes["geometry_csv_sha256"],
        expected_neighbor_fraction=float(phase3_rule["geometry_contract"]["neighbor_fraction"]),
    )
    reconstructed = geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=table.geometries[name].coordinates)
            for name in table.species
        ]
    )
    if reconstructed != fingerprint:
        raise RuntimeError("fresh Phase-4 geometry fingerprint drift")
    if set(phase3_auth["species"]["train_species"]) | set(phase3_auth["species"]["eval_species"]) != set(table.species):
        raise RuntimeError("fresh Phase-4 split does not cover geometry")
    if set(phase3_auth["species"]["train_species"]) & set(phase3_auth["species"]["eval_species"]):
        raise RuntimeError("fresh Phase-4 split overlaps")

    repo_root = args.repo_root.resolve()
    phase3_frozen = phase3_auth.get("frozen_code_sha256")
    if not isinstance(phase3_frozen, dict) or not phase3_frozen:
        raise RuntimeError("Phase-3 authorization lacks frozen core code hashes")
    for relative, expected in sorted(phase3_frozen.items()):
        path = repo_root / str(relative)
        if not path.is_file() or sha256_path(path) != str(expected):
            raise RuntimeError(f"qualified Phase-3 core code drift before Phase 4: {relative}")

    frozen_code: dict[str, str] = dict(phase3_frozen)
    for relative in PHASE4_CODE_PATHS:
        path = repo_root / relative
        if not path.is_file():
            raise RuntimeError(f"Phase-4 critical code file missing: {relative}")
        frozen_code[relative] = sha256_path(path)

    out = {
        "schema": "ttf_genetic_phylogatr_phase4_identity_opening_authorization_v0.1",
        "status": "AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING",
        **hashes,
        "dataset_digest_sha256": phase3_auth["dataset_digest_sha256"],
        "geometry_fingerprint_sha256": fingerprint,
        "species": phase3_auth["species"],
        "phase3_gate_d": {
            "status": qualification["status"],
            "passed": True,
            "type1_gate_pass": True,
            "power_gate_pass": True,
        },
        "fresh_self_detectability": {
            "status": self_qualification["status"],
            "passed": bool(self_qualification["passed"]),
            "negative_interpretation_if_not_passed": "NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING",
        },
        "frozen_code_sha256": frozen_code,
        "opening_scope": {
            "allowed": "Only the selected aligned FASTA nucleotide identities needed to compute the frozen Phase-4 response on exact survivor edges.",
            "forbidden": [
                "Decker empirical genetic outcomes",
                "marker switching",
                "edge rewiring",
                "sequence deletion based on divergence",
                "changing the comparable-site threshold, p-distance definition, IBD model, bandwidth, alpha, or reference family",
            ],
        },
        "outcome_firewall_before_execution": {
            "confirmatory_sequence_identity_opened": False,
            "confirmatory_pairwise_genetic_distances_opened": False,
            "confirmatory_ttf_statistic_opened": False,
            "decker_empirical_genetic_outcomes_opened": False,
        },
        "claim_boundary": "This authorization opens only the exact fresh confirmatory Phase-4 response path. The empirical conclusion exists only after the authorized runner completes successfully.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": out["status"],
                "geometry_fingerprint_sha256": fingerprint,
                "species": len(table.species),
                "self_method_qualified": bool(self_qualification["passed"]),
                "confirmatory_sequence_identity_opened": False,
                "decker_empirical_genetic_outcomes_opened": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
