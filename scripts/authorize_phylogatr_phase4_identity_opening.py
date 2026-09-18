#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from ttf.genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from ttf.geometry import SpeciesGeometry, geometry_fingerprint


FRAGILITY_EXECUTION_RULE_GIT_BLOB_SHA = "731928fef0435921e7d4966e67e04e1f3435d493"
PHASE1_PROVENANCE_RECOVERY = Path("docs/supporting/genetic_phylogatr_phase1_provenance_recovery_v0.1.json")

PHASE4_CODE_PATHS = (
    "src/ttf/phylogatr_phase4.py",
    "src/ttf/phylogatr_empirical.py",
    "src/ttf/genetic_empirical_score.py",
    "src/ttf/phylogatr_compact_empirical.py",
    "src/ttf/phylogatr_compact_execution.py",
    "src/ttf/phylogatr_compact_ibd.py",
    "src/ttf/phylogatr_compact_self_detectability.py",
    "src/ttf/phylogatr_character_mask.py",
    "src/ttf/phylogatr_confirmatory.py",
    "src/ttf/genetic_self_detectability.py",
    "src/ttf/private_null_inference.py",
    "src/ttf/profiled_private_null.py",
    "src/ttf/phylogatr_fragility_execution.py",
    "src/ttf/phylogatr_fragility_formal_anchor.py",
    "scripts/run_phylogatr_phase3_self_reference_shard.py",
    "scripts/run_phylogatr_phase3_self_evaluation_shard.py",
    "scripts/plan_phylogatr_gate_d_fragility_synthetic.py",
    "scripts/run_phylogatr_gate_d_fragility_reference_shard.py",
    "scripts/aggregate_phylogatr_gate_d_fragility_references.py",
    "scripts/run_phylogatr_gate_d_fragility_observed_shard.py",
    "scripts/aggregate_phylogatr_gate_d_fragility_curve.py",
    "scripts/authorize_phylogatr_phase4_identity_opening.py",
    "scripts/run_phylogatr_phase4_empirical_test.py",
    "scripts/run_phylogatr_projected_phase4_empirical_test.py",
    "src/ttf/phylogatr_source_projection.py",
    "docs/supporting/genetic_phylogatr_source_projection_rule_v0.1.json",
    "docs/supporting/genetic_phylogatr_phase1_provenance_recovery_v0.1.json",
    "docs/supporting/genetic_phylogatr_phase4_compact_execution_v0.1.json",
)


def _load(path: Path, schema: str) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _closed(payload: dict, label: str) -> None:
    if not payload or any(value is not False for value in payload.values()):
        raise RuntimeError(f"{label} is not fully closed")


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _validate_phase1_provenance_recovery(
    recovery: dict,
    phase2: dict,
    phase3_auth: dict,
    *,
    phase2_manifest_sha256: str,
    phase3_authorization_sha256: str,
) -> str:
    if recovery.get("status") != "FROZEN_RESPONSE_BLIND_PHASE1_PROVENANCE_RECOVERY":
        raise RuntimeError("Phase-1 provenance recovery is not frozen")
    _closed(recovery.get("outcome_firewall", {}), "Phase-1 provenance recovery firewall")

    authority = recovery.get("authority")
    if not isinstance(authority, dict):
        raise RuntimeError("Phase-1 provenance recovery authority block missing")
    for key in (
        "can_reconstruct_missing_phase1_manifest_bytes",
        "can_change_phase1_or_phase2_selection",
        "can_change_phase3_decision",
        "can_change_thresholds_or_estimands",
        "can_open_sequence_identity",
        "can_open_pairwise_genetic_distances",
        "can_open_empirical_ttf",
        "can_authorize_phase4_by_itself",
    ):
        if authority.get(key) is not False:
            raise RuntimeError(f"Phase-1 provenance recovery gained forbidden authority: {key}")

    manifest_sha = str(recovery.get("phase1_manifest_sha256", ""))
    if manifest_sha != str(phase3_auth.get("phase1_manifest_sha256", "")):
        raise RuntimeError("Phase-1 recovery manifest hash differs from Phase-3 authorization")
    if recovery.get("phase2_manifest_sha256") != phase2_manifest_sha256:
        raise RuntimeError("Phase-1 recovery Phase-2 manifest provenance drift")
    if recovery.get("phase3_authorization_sha256") != phase3_authorization_sha256:
        raise RuntimeError("Phase-1 recovery Phase-3 authorization provenance drift")
    if recovery.get("phase1_summary") != phase2.get("phase1"):
        raise RuntimeError("Phase-1 recovery summary differs from exact Phase-2 witness")
    if recovery.get("source_integrity") != phase2.get("source_integrity"):
        raise RuntimeError("Phase-1 recovery source-integrity witness drift")
    if recovery["phase1_summary"].get("dataset_digest_sha256") != phase3_auth.get("dataset_digest_sha256"):
        raise RuntimeError("Phase-1 recovery dataset digest differs from Phase-3 authorization")
    return manifest_sha


def _validate_fragility_completion(
    curve: dict,
    execution_rule: dict,
    *,
    phase3_rule_sha256: str,
    formal_qualification_sha256: str,
    execution_rule_sha256: str,
    expected_fingerprint: str,
    expected_dataset_digest: str,
) -> list[float]:
    policy = execution_rule.get("phase4_completion_policy")
    if not isinstance(policy, dict):
        raise RuntimeError("fragility execution rule lacks Phase-4 completion policy")
    if policy.get("completed_curve_required_before_identity_opening") is not True:
        raise RuntimeError("fragility completion is not frozen as a Phase-4 prerequisite")
    if policy.get("diagnostic_metrics_do_not_gate_phase4") is not True:
        raise RuntimeError("fragility diagnostic metrics unexpectedly gate Phase 4")
    if policy.get("only_formal_phase3_pass_controls_scientific_gate") is not True:
        raise RuntimeError("formal Phase-3 is not the sole scientific gate")
    if curve.get("status") != policy.get("required_curve_status"):
        raise RuntimeError("fragility diagnostic curve is incomplete")
    if curve.get("schema") != policy.get("required_curve_schema"):
        raise RuntimeError("fragility diagnostic curve schema drift")

    if curve.get("phase3_rule_sha256") != phase3_rule_sha256:
        raise RuntimeError("fragility curve formal Phase-3 rule provenance drift")
    if curve.get("formal_qualification_sha256") != formal_qualification_sha256:
        raise RuntimeError("fragility curve formal qualification provenance drift")
    if curve.get("execution_rule_sha256") != execution_rule_sha256:
        raise RuntimeError("fragility curve execution-rule provenance drift")
    if curve.get("full_geometry_fingerprint_sha256") != expected_fingerprint:
        raise RuntimeError("fragility curve full-geometry fingerprint drift")
    if curve.get("dataset_digest_sha256") != expected_dataset_digest:
        raise RuntimeError("fragility curve dataset digest drift")
    if curve.get("formal_full_geometry_status") != "PASS":
        raise RuntimeError("fragility curve is not anchored to a formal Phase-3 PASS")

    _closed(curve.get("authority_firewall", {}), "fragility curve authority firewall")
    for key in (
        "confirmatory_sequence_identity_opened",
        "confirmatory_pairwise_genetic_distances_opened",
        "confirmatory_ttf_statistic_opened",
        "formal_gate_d_decision_made_by_this_curve",
        "phase4_identity_opening_authorized_by_this_curve",
    ):
        if curve.get(key) is not False:
            raise RuntimeError(f"fragility curve firewall/authority drift: {key}")
    if curve.get("diagnostic_completion_is_procedural_prerequisite_only") is not True:
        raise RuntimeError("fragility completion is not marked procedural-only")

    expected_retentions = [1.0] + [
        float(value)
        for value in execution_rule["thinned_geometry_policy"]["run_retention_fractions"]
    ]
    declared = [float(value) for value in curve.get("expected_retention_fractions", [])]
    if declared != expected_retentions:
        raise RuntimeError("fragility curve declared retention levels drift")
    rows = curve.get("curve")
    if not isinstance(rows, list):
        raise RuntimeError("fragility curve rows missing")
    observed = [float(row["retention_fraction"]) for row in rows]
    if observed != expected_retentions:
        raise RuntimeError("fragility curve is incomplete or reordered")

    anchor = rows[0]
    if anchor.get("source") != "formal_phase3_qualification_receipt":
        raise RuntimeError("fragility retention=1.0 is not the formal Phase-3 anchor")
    if anchor.get("formal_status") != "PASS":
        raise RuntimeError("fragility formal anchor is not PASS")
    if anchor.get("formal_gate_d_decision_authority") is not True:
        raise RuntimeError("fragility formal anchor lost formal decision authority")
    if anchor.get("diagnostic_decision_authority") is not False:
        raise RuntimeError("fragility formal anchor gained diagnostic decision authority")

    allowed_thinned_status = {
        "DIAGNOSTIC_SYNTHETIC_SUMMARY",
        "STRUCTURAL_SUPPORT_BELOW_FORMAL_MINIMUM",
    }
    for row in rows[1:]:
        if row.get("status") not in allowed_thinned_status:
            raise RuntimeError("fragility thinned level is incomplete")
        if row.get("formal_gate_d_decision_authority") is not False:
            raise RuntimeError("fragility thinned level gained formal decision authority")
        if row.get("phase4_identity_opening_authority") is not False:
            raise RuntimeError("fragility thinned level gained Phase-4 authority")
    return observed


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Authorize exact fresh phylogatR nucleotide-identity opening after Phase-3 qualification and completed response-blind fragility diagnostics."
    )
    ap.add_argument("--geometry", type=Path, required=True)
    phase1_source = ap.add_mutually_exclusive_group(required=True)
    phase1_source.add_argument("--phase1-manifest", type=Path)
    phase1_source.add_argument("--use-frozen-phase1-provenance-recovery", action="store_true")
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--phase3-authorization", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--qualification", type=Path, required=True)
    ap.add_argument("--self-rule", type=Path, required=True)
    ap.add_argument("--self-references", type=Path, required=True)
    ap.add_argument("--self-qualification", type=Path, required=True)
    ap.add_argument("--fragility-execution-rule", type=Path, required=True)
    ap.add_argument("--fragility-curve", type=Path, required=True)
    ap.add_argument("--phase4-rule", type=Path, required=True)
    ap.add_argument("--opening-state", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    repo_root = args.repo_root.resolve()
    phase1 = None
    phase1_recovery = None
    phase1_recovery_path = repo_root / PHASE1_PROVENANCE_RECOVERY
    if args.phase1_manifest is not None:
        phase1 = _load(
            args.phase1_manifest, "ttf_genetic_phylogatr_confirmatory_phase1_geometry_v0.1"
        )
    else:
        phase1_recovery = _load(
            phase1_recovery_path, "ttf_genetic_phylogatr_phase1_provenance_recovery_v0.1"
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
    fragility_execution_rule = _load(
        args.fragility_execution_rule,
        "ttf_genetic_phylogatr_gate_d_fragility_execution_rule_v0.1",
    )
    fragility_curve = _load(
        args.fragility_curve,
        "ttf_genetic_phylogatr_gate_d_fragility_curve_v0.1",
    )
    phase4_rule = _load(
        args.phase4_rule, "ttf_genetic_phylogatr_phase4_response_rule_v0.1"
    )
    opening = _load(args.opening_state, "ttf_genetic_empirical_opening_state_v0.1")

    if _git_blob_sha1(args.fragility_execution_rule) != FRAGILITY_EXECUTION_RULE_GIT_BLOB_SHA:
        raise RuntimeError("frozen fragility execution rule Git blob SHA drift before Phase 4")

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
        "phase2_manifest_sha256": sha256_path(args.phase2_manifest),
        "phase3_rule_sha256": sha256_path(args.phase3_rule),
        "phase3_authorization_sha256": sha256_path(args.phase3_authorization),
        "phase3_references_sha256": sha256_path(args.references),
        "phase3_qualification_sha256": sha256_path(args.qualification),
        "phase3_self_rule_sha256": sha256_path(args.self_rule),
        "phase3_self_references_sha256": sha256_path(args.self_references),
        "phase3_self_qualification_sha256": sha256_path(args.self_qualification),
        "fragility_execution_rule_sha256": sha256_path(args.fragility_execution_rule),
        "fragility_curve_sha256": sha256_path(args.fragility_curve),
        "phase4_rule_sha256": sha256_path(args.phase4_rule),
        "opening_state_sha256": sha256_path(args.opening_state),
        "geometry_csv_sha256": sha256_path(args.geometry),
    }
    if phase1 is not None:
        hashes["phase1_manifest_sha256"] = sha256_path(args.phase1_manifest)
        phase1_provenance_mode = "exact_manifest"
    else:
        hashes["phase1_manifest_sha256"] = _validate_phase1_provenance_recovery(
            phase1_recovery,
            phase2,
            phase3_auth,
            phase2_manifest_sha256=hashes["phase2_manifest_sha256"],
            phase3_authorization_sha256=hashes["phase3_authorization_sha256"],
        )
        hashes["phase1_provenance_recovery_sha256"] = sha256_path(phase1_recovery_path)
        phase1_provenance_mode = "frozen_downstream_hash_witness"
    if phase3_auth.get("phase1_manifest_sha256") != hashes["phase1_manifest_sha256"]:
        raise RuntimeError("Phase-1 manifest differs from Phase-3 authorization")
    if phase3_auth.get("phase2_manifest_sha256") != hashes["phase2_manifest_sha256"]:
        raise RuntimeError("Phase-2 manifest differs from Phase-3 authorization")
    if phase3_auth.get("phase3_rule_sha256") != hashes["phase3_rule_sha256"]:
        raise RuntimeError("Phase-3 rule differs from Phase-3 authorization")
    if phase3_auth.get("geometry_csv_sha256") != hashes["geometry_csv_sha256"]:
        raise RuntimeError("geometry CSV differs from Phase-3 authorization")

    fingerprint = str(phase3_auth["geometry_fingerprint_sha256"])
    dataset_digest = str(phase3_auth["dataset_digest_sha256"])
    fragility_retentions = _validate_fragility_completion(
        fragility_curve,
        fragility_execution_rule,
        phase3_rule_sha256=hashes["phase3_rule_sha256"],
        formal_qualification_sha256=hashes["phase3_qualification_sha256"],
        execution_rule_sha256=hashes["fragility_execution_rule_sha256"],
        expected_fingerprint=fingerprint,
        expected_dataset_digest=dataset_digest,
    )

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
        "dataset_digest_sha256": dataset_digest,
        "geometry_fingerprint_sha256": fingerprint,
        "species": phase3_auth["species"],
        "phase1_provenance": {
            "mode": phase1_provenance_mode,
            "manifest_sha256": hashes["phase1_manifest_sha256"],
            "recovery_sha256": hashes.get("phase1_provenance_recovery_sha256"),
            "recovery_changes_scientific_gate": False,
        },
        "phase3_gate_d": {
            "status": qualification["status"],
            "passed": True,
            "type1_gate_pass": True,
            "power_gate_pass": True,
        },
        "fragility_diagnostic": {
            "status": fragility_curve["status"],
            "completion_required_before_opening": True,
            "completion_is_procedural_prerequisite_only": True,
            "diagnostic_metrics_gate_phase4": False,
            "curve_sha256": hashes["fragility_curve_sha256"],
            "execution_rule_sha256": hashes["fragility_execution_rule_sha256"],
            "retention_levels": fragility_retentions,
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
                "using fragility diagnostic metrics to alter, rescue, or reinterpret the formal Phase-3 decision",
            ],
        },
        "outcome_firewall_before_execution": {
            "confirmatory_sequence_identity_opened": False,
            "confirmatory_pairwise_genetic_distances_opened": False,
            "confirmatory_ttf_statistic_opened": False,
            "decker_empirical_genetic_outcomes_opened": False,
        },
        "claim_boundary": (
            "This authorization opens only the exact fresh confirmatory Phase-4 response path. "
            "Formal Phase-3 PASS is the sole scientific opening gate; the completed fragility curve "
            "is a response-blind procedural prerequisite and has no rescue or PASS/FAIL authority. "
            "The empirical conclusion exists only after the authorized runner completes successfully."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": out["status"],
                "geometry_fingerprint_sha256": fingerprint,
                "species": len(table.species),
                "fragility_curve_complete": True,
                "fragility_metrics_gate_phase4": False,
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
