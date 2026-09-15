from __future__ import annotations

import hashlib
import json
from pathlib import Path


RECEIPT = Path("benchmarks/frozen/genetic_phylogatr_response_blind_method_closure_v0.1.json")


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def test_response_blind_method_closure_is_complete_and_empirical_firewall_remains_closed() -> None:
    payload = json.loads(RECEIPT.read_text())

    assert payload["schema"] == "ttf_genetic_phylogatr_response_blind_method_closure_v0.1"
    assert payload["status"] == "RESPONSE_BLIND_METHOD_PATH_COMPLETE_WAITING_FOR_AUTHENTICATED_ARCHIVE"
    assert payload["source_main_sha_after_method_completion"] == "a0450dbf81c6ed71dbd8c446ee05892dcb02bd6d"

    opening = payload["canonical_empirical_opening_state"]
    assert opening["status"] == "EMPIRICAL_OPENING_CLOSED_SINGLE_EXTERNAL_BLOCKER_FRESH_PHYLOGATR_ARCHIVE"
    assert opening["unchanged_by_this_receipt"] is True
    assert _git_blob_sha1(Path(opening["path"])) == opening["git_blob_sha"]

    assert payload["external_blocker"]["status"] == "WAITING_FOR_UNTOUCHED_AUTHENTICATED_PHYLOGATR_ARCHIVE"
    assert payload["external_blocker"]["same_purpose_source_hunting_closed"] is True
    assert all(value is False for value in payload["global_firewall"].values())


def test_response_blind_method_closure_pins_all_new_method_surfaces() -> None:
    payload = json.loads(RECEIPT.read_text())
    lineage = payload["method_lineage"]

    geometry = lineage["fragility_geometry_diagnostic"]
    assert geometry["pull_request"] == 31
    assert geometry["merge_commit"] == "cd73c4537263ace177e58e644c6c9d3f149c35f0"
    assert _git_blob_sha1(Path(geometry["rule"])) == geometry["rule_git_blob_sha"]

    execution = lineage["fragility_synthetic_execution"]
    assert execution["pull_request"] == 32
    assert execution["merge_commit"] == "bbfcffce1554d3ff6bf5c43ca9d8b1014ff772de"
    assert _git_blob_sha1(Path(execution["rule"])) == execution["rule_git_blob_sha"]

    phase2 = lineage["phase2_archive_intake"]
    assert phase2["pull_request"] == 33
    assert phase2["merge_commit"] == "a0450dbf81c6ed71dbd8c446ee05892dcb02bd6d"
    assert phase2["full_ci_conclusion"] == "success"
    assert _git_blob_sha1(Path(phase2["script"])) == phase2["script_git_blob_sha"]

    handoff = payload["execution_handoff"]
    assert _git_blob_sha1(Path(handoff["runbook"])) == handoff["runbook_git_blob_sha"]
    assert handoff["canonical_source_mode"] == "one untouched authenticated phylogatR portal archive"
    assert handoff["phase1_archive_direct"] is True
    assert handoff["phase2_same_archive_direct"] is True
    assert handoff["manual_archive_extraction_required"] is False


def test_response_blind_method_closure_cannot_promote_fragility_to_a_scientific_gate() -> None:
    payload = json.loads(RECEIPT.read_text())
    contract = payload["scientific_contract_unchanged"]
    handoff = payload["execution_handoff"]

    assert contract["primary_estimand"] == "place_beyond_ibd"
    assert contract["formal_cross_species_gate"] == "retention=1.0 exact Phase-2 survivor geometry only"
    assert contract["type1_wilson95_upper_ceiling"] == 0.10
    assert contract["shared_A2_wilson95_lower_floor"] == 0.80
    assert contract["bandwidth_km"] == 500.0
    assert contract["alpha"] == 0.05
    assert contract["post_result_rescue_forbidden"] is True
    assert contract["diagnostic_values_can_change_formal_gate"] is False
    assert contract["diagnostic_values_can_rescue_not_evaluable"] is False
    assert contract["diagnostic_values_can_authorize_phase4"] is False
    assert contract["failed_geometry_mask_or_qualification_is_biological_null"] is False

    assert handoff["phase3_formal_decision_is_retention_1_only"] is True
    assert handoff["phase4_requires_formal_phase3_pass"] is True
    assert handoff["phase4_requires_completed_fragility_curve_procedurally"] is True
    assert handoff["fragility_metrics_gate_phase4_scientifically"] is False
    assert handoff["phase4_opening_state_baseline"] == "benchmarks/frozen/genetic_empirical_opening_state_v0.1.json"
