#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def load(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


def sha256_path(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob_sha1(path: str | Path) -> str:
    data = Path(path).read_bytes()
    header = f"blob {len(data)}\0".encode("utf-8")
    return hashlib.sha1(header + data).hexdigest()


def assert_closed_firewall(payload: dict, key: str = "response_firewall") -> None:
    fw = payload.get(key)
    if not isinstance(fw, dict) or not fw:
        raise RuntimeError(f"missing {key}")
    if any(value is not False for value in fw.values()):
        raise RuntimeError(f"open response firewall in {key}: {fw}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    paths = {
        "program": Path("docs/supporting/relational_program_v0.3.json"),
        "transition": Path("docs/supporting/relational_program_transition_rule_v0.1.json"),
        "freshness": Path("docs/supporting/relational_future_family_freshness_amendment_v0.1.json"),
        "s3_receipt": Path("benchmarks/frozen/relational_host_resource_empirical_repaired_receipt_v0.2.json"),
        "s3_reproduction": Path("benchmarks/frozen/relational_host_resource_empirical_reproduction_v0.2.json"),
        "s3_overlap": Path("benchmarks/frozen/relational_prior_host_test_overlap_audit_v0.1.json"),
        "c_s1_exclusion": Path("benchmarks/frozen/relational_prior_S1_species_exclusion_v0.1.json"),
        "c_s2_exclusion": Path("benchmarks/frozen/relational_prior_S2_species_exclusion_v0.1.json"),
        "b_relation": Path("docs/supporting/relational_environment_relation_rule_v0.3.json"),
        "b_transport": Path("benchmarks/frozen/relational_environment_transport_execution_v0.7.json"),
        "b_transport_audit": Path("benchmarks/frozen/relational_environment_transport_partition_audit_v0.7.json"),
        "b_local_geometry": Path("benchmarks/frozen/relational_environment_geometry_reconstruction_local_v0.2.json"),
        "b_preretry": Path("benchmarks/frozen/relational_environment_transport_preretry_receipt_v0.7.json"),
        "b_retry_rule_v02": Path("docs/supporting/relational_environment_transport_retry_rule_v0.2.json"),
        "b_recovery_v02": Path("benchmarks/frozen/relational_environment_request_error_recovery_v0.2.json"),
        "b_final_recovery_audit": Path("benchmarks/frozen/relational_environment_final_recovery_audit_contract_v0.3.json"),
        "b_final_producer_promotion": Path("docs/supporting/relational_environment_final_producer_promotion_rule_v0.2.json"),
        "b_final_producer_contract": Path("benchmarks/frozen/relational_environment_final_producer_contract_v0.2.json"),
        "b_unbound_long_repair": Path("benchmarks/frozen/relational_environment_cutover_long_repair_v0.1.json"),
        "b_producer": Path("benchmarks/frozen/relational_environment_relation_producer_v0.1.json"),
        "b_opportunity": Path("docs/supporting/relational_environment_opportunity_rule_v0.2.json"),
        "b_qualification": Path("docs/supporting/relational_environment_qualification_rule_v0.2.json"),
        "b_mask": Path("docs/supporting/relational_environment_character_mask_rule_v0.1.json"),
        "b_empirical": Path("docs/supporting/relational_environment_empirical_opening_rule_v0.1.json"),
        "c_relation": Path("docs/supporting/relational_historical_climate_exposure_rule_v0.1.json"),
        "c_opportunity": Path("docs/supporting/relational_historical_climate_opportunity_rule_v0.1.json"),
        "c_qualification": Path("docs/supporting/relational_historical_climate_qualification_rule_v0.1.json"),
        "c_mask": Path("docs/supporting/relational_historical_climate_character_mask_rule_v0.1.json"),
        "c_empirical": Path("docs/supporting/relational_historical_climate_empirical_opening_rule_v0.1.json"),
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing v0.3 contract inputs: {missing}")

    p = {name: load(path) for name, path in paths.items()}
    program = p["program"]
    if program.get("schema") != "ttf_relational_program_v0.3":
        raise RuntimeError("unexpected relational program schema")
    if program.get("status") != "FROZEN_FINITE_PROSPECTIVE_FAMILY_BEFORE_STUDY_B_FORMAL_QUALIFICATION":
        raise RuntimeError("unexpected relational program status")

    family = program["prospective_family"]
    if family["slots_fixed_before_first_v03_qualification"] != ["B", "C"]:
        raise RuntimeError("v0.3 prospective slots drift")
    allocation = family["allocation"]
    alpha_b = float(allocation["B_environmental_niche_similarity"])
    alpha_c = float(allocation["C_historical_climate_exposure_similarity"])
    if abs(float(family["familywise_alpha"]) - 0.05) > 1e-15:
        raise RuntimeError("familywise alpha drift")
    if abs(alpha_b - 0.025) > 1e-15 or abs(alpha_c - 0.025) > 1e-15:
        raise RuntimeError("slot alpha drift")
    if abs(alpha_b + alpha_c - 0.05) > 1e-15:
        raise RuntimeError("Bonferroni allocation no longer sums to 0.05")
    if family["no_alpha_recycling"] is not True or family["no_replacement_slots"] is not True:
        raise RuntimeError("finite-family stopping firewall drift")
    if int(family["maximum_confirmatory_openings"]) != 2:
        raise RuntimeError("unexpected number of confirmatory openings")

    prior = {
        row["label"]: row for row in program["completed_predictor_tests_reported_as_prior_evidence"]
    }
    s3 = prior["S3 relational larval host-resource similarity"]
    receipt = p["s3_receipt"]
    if receipt.get("schema") != "ttf_relational_host_resource_empirical_repaired_receipt_v0.2":
        raise RuntimeError("S3 receipt schema drift")
    if s3["receipt"] != str(paths["s3_receipt"]):
        raise RuntimeError("program does not point to repaired S3 receipt")
    if s3["reproduction_audit"] != str(paths["s3_reproduction"]):
        raise RuntimeError("program does not point to exact S3 reproduction audit")
    if abs(float(s3["p_one_sided"]) - float(receipt["primary"]["p_value_one_sided"])) > 1e-15:
        raise RuntimeError("S3 p-value drift between program and repaired receipt")
    if s3["result"] != receipt["decision"]:
        raise RuntimeError("S3 decision drift")
    if receipt["one_shot"]["post_result_retuning_allowed"] is not False:
        raise RuntimeError("S3 post-result retuning firewall drift")
    if receipt["one_shot"]["result_selection_rerun_allowed"] is not False:
        raise RuntimeError("S3 result-selection firewall drift")
    if int(receipt["empirical_design"]["species_union"]) != 399:
        raise RuntimeError("repaired S3 empirical species count drift")
    if int(receipt["empirical_design"]["directed_source_target_dyads"]) != 11696:
        raise RuntimeError("repaired S3 dyad count drift")

    assert_closed_firewall(program)
    transition = p["transition"]
    if program.get("transition_rule") != str(paths["transition"]):
        raise RuntimeError("B-to-C transition-rule pointer drift")
    if transition.get("schema") != "ttf_relational_program_transition_rule_v0.1":
        raise RuntimeError("unexpected B-to-C transition schema")
    if transition.get("status") != "FROZEN_BEFORE_STUDY_B_EMPIRICAL_RESULT_AND_BEFORE_ANY_STUDY_C_EXTERNAL_OR_GENETIC_OPENING":
        raise RuntimeError("B-to-C transition rule was not frozen pre-result")
    states = transition["allowed_B_states"]
    if states["DETECTED_B"]["C_open_authorized"] is not False:
        raise RuntimeError("DETECTED_B must forbid Study C")
    for state, entry in (
        ("STUDY_B_ENVIRONMENT_RELATIONAL_NULL_WITH_QUALIFIED_POWER", "PROCEED_TO_C_AFTER_B_NULL"),
        ("NOT_EVALUABLE_B", "PROCEED_TO_C_AFTER_B_NOT_EVALUABLE"),
    ):
        if states[state]["C_open_authorized"] is not True or states[state]["C_entry_state"] != entry:
            raise RuntimeError(f"Study-C transition drift for {state}")
    assert_closed_firewall(transition)
    for key in ("b_relation", "b_transport", "b_transport_audit", "b_preretry", "b_retry_rule_v02", "b_recovery_v02", "b_final_recovery_audit", "b_final_producer_promotion", "b_final_producer_contract", "b_unbound_long_repair", "b_producer", "b_qualification", "b_mask", "b_empirical", "c_relation", "c_opportunity", "c_qualification", "c_mask", "c_empirical"):
        assert_closed_firewall(p[key])

    transport = p["b_transport"]
    if transport.get("schema") != "ttf_relational_environment_transport_execution_v0.7":
        raise RuntimeError("Study B transport execution schema drift")
    if transport.get("status") != "FROZEN_FINAL_AUTHORITATIVE_COMPOSITE_TRANSPORT_BEFORE_RELATION_RESULT":
        raise RuntimeError("Study B transport execution is not frozen final authoritative v0.7")
    base = transport["frozen_base_component"]
    cutover = transport["cutover_component"]
    short_repair = transport["short_timeout_repair_component"]
    long_repair = transport["long_timeout_repair_component"]
    final_partition = transport["final_partition"]
    if int(base["workflow_run_id"]) != 35691973637 or int(base["accepted_species_count"]) != 165:
        raise RuntimeError("Study B frozen base transport drift")
    if int(cutover["workflow_run_id"]) != 35712027961:
        raise RuntimeError("Study B cutover transport run drift")
    if cutover["workflow_name"] != "relational-environment-transport-cutover":
        raise RuntimeError("Study B cutover workflow drift")
    if int(cutover["accepted_success_batch_count"]) != 207 or int(cutover["accepted_species_count"]) != 827:
        raise RuntimeError("Study B cutover-success partition drift")
    if int(short_repair["workflow_run_id"]) != 35735717833:
        raise RuntimeError("Study B short repair run drift")
    if list(short_repair["accepted_repair_indices"]) != [2, 3, 4, 5, 7]:
        raise RuntimeError("Study B short repair accepted-index drift")
    if int(short_repair["accepted_species_count"]) != 5:
        raise RuntimeError("Study B short repair species count drift")
    if list(short_repair["ignored_repair_indices"]) != [0, 1, 6]:
        raise RuntimeError("Study B short repair ignored-index drift")
    if int(long_repair["workflow_run_id"]) != 35749326437:
        raise RuntimeError("Study B long repair run drift")
    if list(long_repair["repair_indices"]) != [0, 1, 2]:
        raise RuntimeError("Study B long repair index drift")
    if list(long_repair["replaces_short_repair_indices"]) != [0, 1, 6]:
        raise RuntimeError("Study B long/short replacement mapping drift")
    if int(long_repair["species_count"]) != 3 or int(long_repair["timeout_minutes"]) != 240:
        raise RuntimeError("Study B long repair contract drift")
    if int(final_partition["exact_species"]) != 1000:
        raise RuntimeError("Study B final transport species count drift")
    if (
        int(final_partition["base_species"]) != 165
        or int(final_partition["cutover_success_species"]) != 827
        or int(final_partition["short_repair_species"]) != 5
        or int(final_partition["long_repair_species"]) != 3
    ):
        raise RuntimeError("Study B final transport split drift")
    if final_partition["overlap_allowed"] is not False or final_partition["backfill_allowed"] is not False:
        raise RuntimeError("Study B final transport overlap/backfill firewall drift")
    if int(transport["acceptance_gate"]["exact_species_ledgers"]) != 1000:
        raise RuntimeError("Study B transport exact-species gate drift")
    if int(transport["acceptance_gate"]["request_error_count_must_equal"]) != 0:
        raise RuntimeError("Study B transport zero-error gate drift")
    if 35740715490 not in set(map(int, transport["obsolete_or_ignored_runs"])):
        raise RuntimeError("obsolete v0.6 relation producer run is not explicitly ignored")

    audit = p["b_transport_audit"]
    if audit.get("schema") != "ttf_relational_environment_transport_partition_audit_v0.7":
        raise RuntimeError("Study B v0.7 partition-audit schema drift")
    if audit.get("status") != "PASS_EXACT_1000_SPECIES_DISJOINT_COMPOSITE_PARTITION":
        raise RuntimeError("Study B v0.7 partition audit did not pass")
    if audit.get("transport_execution") != str(paths["b_transport"]):
        raise RuntimeError("Study B v0.7 partition audit transport pointer drift")
    if int(audit["union_species_count"]) != 1000:
        raise RuntimeError("Study B v0.7 partition-audit union drift")
    if audit["overlap_candidate_indices"] or audit["missing_candidate_indices"] or audit["species_name_mismatches"]:
        raise RuntimeError("Study B v0.7 partition audit found overlap/missing/name drift")
    checks = audit["checks"]
    for key in (
        "candidate_count",
        "base_species_165",
        "cutover_success_species_827",
        "short_repair_species_5",
        "long_repair_species_3",
        "overlap_zero",
        "missing_zero",
        "union_1000",
        "species_names_match",
    ):
        if checks.get(key) is not True:
            raise RuntimeError(f"Study B v0.7 partition audit failed: {key}")

    preretry = p["b_preretry"]
    if preretry.get("schema") != "ttf_relational_environment_transport_preretry_receipt_v0.7":
        raise RuntimeError("Study B pre-retry transport receipt schema drift")
    if preretry.get("status") != "FROZEN_RESPONSE_BLIND_PRE_RETRY_TRANSPORT_STATE":
        raise RuntimeError("Study B pre-retry transport receipt status drift")
    if preretry.get("transport_execution") != str(paths["b_transport"]):
        raise RuntimeError("Study B pre-retry transport pointer drift")
    if int(preretry["species"]) != 1000 or int(preretry["selected_artifacts"]) != 248:
        raise RuntimeError("Study B pre-retry exact-universe drift")
    expected_counts = {
        "PASS_OCCURRENCE_GEOMETRY": 209,
        "FAIL_OCCURRENCE_GEOMETRY": 72,
        "REJECTED_GBIF_TAXON_MATCH": 15,
        "REQUEST_ERROR": 704,
    }
    if preretry["status_counts"] != expected_counts:
        raise RuntimeError("Study B pre-retry status-count drift")
    if int(preretry["request_error_count"]) != 704:
        raise RuntimeError("Study B pre-retry REQUEST_ERROR count drift")
    if preretry["request_error_species_sorted_sha256"] != "9902e6d65c951cb752c06a2ce18fc7ae96fede0ab3e40089624402c185cdd034":
        raise RuntimeError("Study B pre-retry species digest drift")

    retry_v02 = p["b_retry_rule_v02"]
    if retry_v02.get("schema") != "ttf_relational_environment_transport_retry_rule_v0.2":
        raise RuntimeError("Study B final retry-rule schema drift")
    if retry_v02.get("status") != "FROZEN_FINAL_TECHNICAL_RETRY_AMENDMENT_BEFORE_RELATION_RESULT":
        raise RuntimeError("Study B final retry rule is not frozen pre-result")
    if retry_v02.get("supersedes") != transport.get("retry_rule"):
        raise RuntimeError("Study B final retry rule does not supersede the v0.7 retry rule")
    if int(retry_v02["prior_maximum_retry_rounds"]) != 3:
        raise RuntimeError("Study B prior retry-round count drift")
    if int(retry_v02["maximum_retry_rounds_total"]) != 4 or int(retry_v02["final_round"]) != 4:
        raise RuntimeError("Study B final retry-round ceiling drift")
    if retry_v02.get("no_fifth_round_authorized") is not True:
        raise RuntimeError("Study B fifth retry round is not explicitly forbidden")

    recovery_v02 = p["b_recovery_v02"]
    if recovery_v02.get("schema") != "ttf_relational_environment_request_error_recovery_v0.2":
        raise RuntimeError("Study B final recovery schema drift")
    if recovery_v02.get("status") != "FROZEN_FINAL_TRANSPORT_ROUND_BEFORE_RELATION_RESULT":
        raise RuntimeError("Study B final recovery is not frozen pre-result")
    if recovery_v02.get("retry_rule") != str(paths["b_retry_rule_v02"]):
        raise RuntimeError("Study B final recovery retry-rule pointer drift")
    if int(recovery_v02["retry_round"]) != 4:
        raise RuntimeError("Study B final recovery round drift")
    if int(recovery_v02["source_fallback_run_id"]) != 35809546248:
        raise RuntimeError("Study B round-3 fallback run binding drift")
    if int(recovery_v02["source_audit_run_id"]) != 35824707898:
        raise RuntimeError("Study B round-3 audit run binding drift")
    if int(recovery_v02["source_audit_artifact_id"]) != 10734318651:
        raise RuntimeError("Study B round-3 audit artifact binding drift")
    if recovery_v02["source_audit_artifact_digest"] != "sha256:3a91136cbebd1518acce1a56112c300730397db8a82ecb74a07bd2d441abd800":
        raise RuntimeError("Study B round-3 audit artifact digest drift")
    if int(recovery_v02["source_original_request_error_count"]) != 704:
        raise RuntimeError("Study B round-3 source universe drift")
    if int(recovery_v02["source_unresolved_request_error_count"]) != 556:
        raise RuntimeError("Study B final recovery unresolved-species count drift")
    if int(recovery_v02["batch_size"]) != 3 or int(recovery_v02["batch_count"]) != 186:
        raise RuntimeError("Study B final recovery batching drift")
    if int(recovery_v02["max_parallel"]) != 4:
        raise RuntimeError("Study B final recovery max-parallel drift")
    if recovery_v02["transport_implementation"]["scientific_query_change"] is not False:
        raise RuntimeError("Study B final recovery changed scientific query")
    workflow_path = recovery_v02["transport_implementation"]["workflow"]
    if workflow_path != ".github/workflows/relational-environment-request-error-recovery-v02.yml":
        raise RuntimeError("Study B final recovery workflow pointer drift")
    if git_blob_sha1(workflow_path) != recovery_v02["transport_implementation"]["workflow_git_blob"]:
        raise RuntimeError("Study B final recovery workflow blob drift")
    for path, expected in recovery_v02["corrected_transport_core_git_blobs"].items():
        if git_blob_sha1(path) != expected:
            raise RuntimeError(f"Study B corrected transport-core blob drift: {path}")
    final_audit = recovery_v02["final_round_audit"]
    if int(final_audit["expected_species"]) != 556:
        raise RuntimeError("Study B final-round audit universe drift")
    if int(final_audit["required_request_error_count_for_relation_binding"]) != 0:
        raise RuntimeError("Study B final-round zero-error binding gate drift")
    if final_audit.get("no_fifth_round_authorized") is not True:
        raise RuntimeError("Study B final-round fifth-retry firewall drift")

    audit_only = p["b_final_recovery_audit"]
    if audit_only.get("schema") != "ttf_relational_environment_final_recovery_audit_contract_v0.3":
        raise RuntimeError("Study B audit-only final recovery schema drift")
    if audit_only.get("status") != "FROZEN_AUDIT_ONLY_BEFORE_FINAL_RECOVERY_RESULT":
        raise RuntimeError("Study B audit-only recovery verifier was not frozen pre-result")
    if int(audit_only["recovery_run_id"]) != 35827478403:
        raise RuntimeError("Study B audit-only recovery-run binding drift")
    if int(audit_only["source_round3_audit_run_id"]) != 35824707898:
        raise RuntimeError("Study B audit-only round-3 binding drift")
    if int(audit_only["expected_recovery_batches"]) != 186 or int(audit_only["expected_species"]) != 556:
        raise RuntimeError("Study B audit-only universe drift")
    if audit_only.get("network_retrieval_allowed") is not False:
        raise RuntimeError("Study B audit-only verifier permits network acquisition")
    if audit_only.get("no_fifth_round_authorized") is not True:
        raise RuntimeError("Study B audit-only verifier permits a fifth retry")
    if audit_only["workflow"] != ".github/workflows/relational-environment-final-recovery-audit-v03.yml":
        raise RuntimeError("Study B audit-only workflow pointer drift")
    if git_blob_sha1(audit_only["workflow"]) != audit_only["workflow_git_blob"]:
        raise RuntimeError("Study B audit-only workflow Git blob drift")

    final_promotion = p["b_final_producer_promotion"]
    if final_promotion.get("schema") != "ttf_relational_environment_final_producer_promotion_rule_v0.2":
        raise RuntimeError("Study B final producer promotion schema drift")
    if final_promotion.get("status") != "FROZEN_BEFORE_FINAL_ROUND_RESULT_AND_BEFORE_RELATION_RESULT":
        raise RuntimeError("Study B final producer promotion was not frozen pre-result")
    if int(final_promotion["canonical_relation_producer_run_id"]) != 35795821824:
        raise RuntimeError("Study B canonical producer binding drift")
    if int(final_promotion["round3_fallback_run_id"]) != 35809546248:
        raise RuntimeError("Study B round-3 producer-promotion binding drift")
    if int(final_promotion["round3_audit_run_id"]) != 35824707898:
        raise RuntimeError("Study B round-3 audit promotion binding drift")
    if int(final_promotion["final_recovery_run_id"]) != 35827478403:
        raise RuntimeError("Study B final recovery producer-promotion binding drift")
    if int(final_promotion["final_recovery_audit_workflow_run_id"]) != 35828444286:
        raise RuntimeError("Study B audit-only producer-promotion binding drift")
    if final_promotion["final_recovery_audit_artifact_name"] != "relational-environment-request-error-recovery-audit-v0.3":
        raise RuntimeError("Study B authoritative final audit artifact drift")
    if final_promotion.get("no_fifth_round_authorized") is not True:
        raise RuntimeError("Study B producer promotion permits a fifth retry")
    merge_rule = final_promotion["final_merge_rule"]
    if (
        int(merge_rule["pre_retry_request_error_species"]) != 704
        or int(merge_rule["round3_resolved_species"]) != 148
        or int(merge_rule["round3_unresolved_species"]) != 556
        or int(merge_rule["round4_species_exactly"]) != 556
    ):
        raise RuntimeError("Study B final producer merge universe drift")
    if merge_rule["round4_may_replace_only_round3_request_error"] is not True:
        raise RuntimeError("Study B final producer round-4 replacement scope drift")
    if merge_rule["overwrite_round3_non_request_error_forbidden"] is not True:
        raise RuntimeError("Study B final producer may overwrite resolved round-3 outcomes")
    if int(merge_rule["final_request_error_required_for_relation_binding"]) != 0:
        raise RuntimeError("Study B final producer zero-error relation-binding gate drift")

    final_producer_contract = p["b_final_producer_contract"]
    if final_producer_contract.get("schema") != "ttf_relational_environment_final_producer_contract_v0.2":
        raise RuntimeError("Study B final producer implementation schema drift")
    if final_producer_contract.get("status") != "FROZEN_BEFORE_FINAL_ROUND_RESULT_AND_BEFORE_RELATION_RESULT":
        raise RuntimeError("Study B final producer implementation was not frozen pre-result")
    if final_producer_contract.get("promotion_rule") != str(paths["b_final_producer_promotion"]):
        raise RuntimeError("Study B final producer promotion-rule pointer drift")
    if final_producer_contract.get("final_recovery_contract") != str(paths["b_recovery_v02"]):
        raise RuntimeError("Study B final producer recovery-contract pointer drift")
    if int(final_producer_contract["final_recovery_run_id"]) != 35827478403:
        raise RuntimeError("Study B final producer recovery-run drift")
    if final_producer_contract.get("final_recovery_audit_contract") != str(paths["b_final_recovery_audit"]):
        raise RuntimeError("Study B final producer audit-contract pointer drift")
    if int(final_producer_contract["final_recovery_audit_workflow_run_id"]) != 35828444286:
        raise RuntimeError("Study B final producer audit-run drift")
    if final_producer_contract["final_recovery_audit_artifact_name"] != "relational-environment-request-error-recovery-audit-v0.3":
        raise RuntimeError("Study B final producer audit-artifact drift")
    if final_producer_contract.get("trigger_must_not_exist_before_zero_error_audit") is not True:
        raise RuntimeError("Study B final producer trigger timing firewall drift")
    final_merge = final_producer_contract["merge_contract"]
    if (
        int(final_merge["authoritative_base_species"]) != 1000
        or int(final_merge["pre_retry_request_error_species"]) != 704
        or int(final_merge["round3_artifacts"]) != 235
        or int(final_merge["round3_resolved_species"]) != 148
        or int(final_merge["round3_unresolved_species"]) != 556
        or int(final_merge["round4_artifacts"]) != 186
        or int(final_merge["round4_species"]) != 556
        or int(final_merge["final_request_error_required"]) != 0
    ):
        raise RuntimeError("Study B final producer merge contract drift")
    if final_merge["overwrite_round3_resolved_forbidden"] is not True:
        raise RuntimeError("Study B final producer round-3 resolution firewall drift")
    if final_merge["no_fifth_round_authorized"] is not True:
        raise RuntimeError("Study B final producer merge permits a fifth retry")
    downstream = final_producer_contract["downstream_preflight"]
    if set(downstream["relation_binding_accepts_occurrence_schemas"]) != {
        "ttf_relational_environment_occurrence_acquisition_v0.2",
        "ttf_relational_environment_occurrence_acquisition_v0.3",
    }:
        raise RuntimeError("Study B relation-binding occurrence schema compatibility drift")
    if downstream["final_occurrence_filename"] != "occurrence_ledger_v0.3.json":
        raise RuntimeError("Study B final occurrence filename drift")
    if downstream["final_transport_receipt_filename"] != "final_transport_merge_receipt_v0.2.json":
        raise RuntimeError("Study B final transport receipt filename drift")
    if downstream["final_transport_receipt_schema"] != "ttf_relational_environment_final_transport_merge_v0.2":
        raise RuntimeError("Study B final transport receipt schema drift")
    if downstream["qualification_final_receipt_status"] != "PASS_ZERO_REQUEST_ERROR_FINAL_MERGE":
        raise RuntimeError("Study B final qualification transport status drift")
    if int(downstream["final_retry_rounds_completed"]) != 4:
        raise RuntimeError("Study B final downstream retry-round drift")
    if downstream["no_fifth_round_authorized"] is not True:
        raise RuntimeError("Study B downstream preflight permits fifth retry")

    relation_contract = final_producer_contract["relation_contract"]
    if relation_contract["relation_rule"] != str(paths["b_relation"]):
        raise RuntimeError("Study B final producer relation-rule pointer drift")
    if relation_contract["freshness_rule"] != str(paths["freshness"]):
        raise RuntimeError("Study B final producer freshness-rule pointer drift")
    if relation_contract["relation_artifact_name"] != "relational-environment-relation-v0.3":
        raise RuntimeError("Study B final producer relation artifact drift")
    if relation_contract["relation_result_may_be_constructed_only_after_final_request_error_zero"] is not True:
        raise RuntimeError("Study B final producer may construct relation before zero-error gate")
    for path, expected in final_producer_contract["git_blobs"].items():
        if git_blob_sha1(path) != expected:
            raise RuntimeError(f"Study B final producer Git blob drift: {path}")

    unbound = p["b_unbound_long_repair"]
    if unbound.get("schema") != "ttf_relational_environment_cutover_long_repair_v0.1":
        raise RuntimeError("unbound Study B technical experiment schema drift")
    if unbound.get("status") != "OBSOLETE_UNBOUND_TECHNICAL_EXPERIMENT":
        raise RuntimeError("unbound two-species repair experiment reactivated")
    serialized_transport = json.dumps(transport, sort_keys=True)
    serialized_audit = json.dumps(audit, sort_keys=True)
    if "relational_environment_cutover_long_repair_v0.1.json" in serialized_transport:
        raise RuntimeError("authoritative v0.7 references obsolete two-species experiment")
    if "relational_environment_cutover_long_repair_v0.1.json" in serialized_audit:
        raise RuntimeError("v0.7 partition audit references obsolete two-species experiment")
    if any(bool(v) for v in unbound["response_firewall"].values()):
        raise RuntimeError("obsolete two-species experiment response firewall is open")

    producer = p["b_producer"]
    if producer.get("schema") != "ttf_relational_environment_relation_producer_v0.1":
        raise RuntimeError("Study B producer schema drift")
    if producer.get("status") != "FROZEN_RELATION_PRODUCER_BEFORE_RESULT":
        raise RuntimeError("Study B producer is not frozen pre-result")
    if producer.get("transport_execution") != str(paths["b_transport"]):
        raise RuntimeError("Study B producer transport pointer drift")
    if int(producer["workflow_run_id"]) <= 0:
        raise RuntimeError("Study B final producer run is invalid")
    if len(str(producer["workflow_head_sha"])) != 40:
        raise RuntimeError("Study B final producer head SHA is invalid")
    allowed_producer_workflows = {
        "relational-environment-transport-assemble-v07",
        "relational-environment-final-recovery-producer",
    }
    if producer["workflow_name"] not in allowed_producer_workflows:
        raise RuntimeError("Study B producer workflow drift")
    if producer["workflow_name"] == "relational-environment-transport-assemble-v07":
        if int(producer["workflow_run_id"]) != int(final_promotion["canonical_relation_producer_run_id"]):
            raise RuntimeError("Study B pre-promotion canonical producer run drift")
    else:
        if producer.get("producer_rule") != str(paths["b_final_producer_promotion"]):
            raise RuntimeError("Study B final-recovery producer-rule pointer drift")
        if int(producer.get("final_transport_recovery_run_id", -1)) != 35827478403:
            raise RuntimeError("Study B final-recovery producer transport binding drift")
        if int(final_promotion["canonical_relation_producer_run_id"]) not in set(
            map(int, producer.get("supersedes_obsolete_producer_runs", []))
        ):
            raise RuntimeError("Study B final-recovery producer did not supersede canonical producer")
    if int(producer["workflow_run_id"]) in set(map(int, transport["obsolete_or_ignored_runs"])):
        raise RuntimeError("Study B final producer is listed obsolete")
    if producer["relation_result_seen"] is not False or producer["genetic_response_used"] is not False:
        raise RuntimeError("Study B producer was not frozen before relation/genetic response")

    local_geometry = p["b_local_geometry"]
    if local_geometry.get("schema") != "ttf_relational_environment_geometry_reconstruction_local_v0.2":
        raise RuntimeError("Study B local geometry receipt schema drift")
    if local_geometry.get("status") != "PASS_RESPONSE_BLIND_FRESH_1000_GEOMETRY_RECONSTRUCTION":
        raise RuntimeError("Study B local geometry reconstruction did not pass")
    if local_geometry["source_archive"]["sha256"] != "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5":
        raise RuntimeError("Study B recovered source archive SHA drift")
    if int(local_geometry["source_archive"]["size_bytes"]) != 274988692:
        raise RuntimeError("Study B recovered source archive size drift")
    if (
        int(local_geometry["species"]) != 1000
        or int(local_geometry["locality_rows"]) != 28355
        or int(local_geometry["edge_rows"]) != 167828
        or int(local_geometry["verification_failures"]) != 0
    ):
        raise RuntimeError("Study B local geometry aggregate invariant drift")
    if local_geometry["relation_result_seen"] is not False or local_geometry["genetic_response_used"] is not False:
        raise RuntimeError("Study B local geometry was not frozen pre-result")
    if not all(v is False for v in local_geometry["response_firewall"].values()):
        raise RuntimeError("Study B local geometry response firewall opened")

    b = program["slot_B"]
    expected_b_chain = [
        str(paths["b_relation"]),
        str(paths["b_opportunity"]),
        str(paths["b_qualification"]),
        str(paths["b_mask"]),
        str(paths["b_empirical"]),
    ]
    if b["execution_chain"] != expected_b_chain:
        raise RuntimeError("Study B execution-chain drift")
    if float(b["alpha_one_sided"]) != 0.025:
        raise RuntimeError("Study B alpha drift")
    if float(p["b_qualification"]["inference"]["alpha"]) != 0.025:
        raise RuntimeError("Study B qualification alpha drift")
    if float(p["b_qualification"]["gate_numeric"]["p_value_cutoff"]) != 0.025:
        raise RuntimeError("Study B synthetic p cutoff drift")
    if float(p["b_qualification"]["gate_numeric"]["private_type1_wilson95_upper_max"]) != 0.05:
        raise RuntimeError("Study B Type-I Wilson gate drift")
    if float(p["b_qualification"]["gate_numeric"]["relational_power_wilson95_lower_min"]) != 0.80:
        raise RuntimeError("Study B power gate drift")
    if float(p["b_mask"]["required_survivor_requalification"]["same_alpha"]) != 0.025:
        raise RuntimeError("Study B mask requalification alpha drift")
    if float(p["b_empirical"]["relational_model"]["alpha_one_sided"]) != 0.025:
        raise RuntimeError("Study B empirical alpha drift")
    if p["b_empirical"]["relational_model"]["predictors_in_order"][0] != "z_R_env":
        raise RuntimeError("Study B primary relation drift")

    c = program["slot_C"]
    expected_c_chain = [
        str(paths["c_relation"]),
        str(paths["c_opportunity"]),
        str(paths["c_qualification"]),
        str(paths["c_mask"]),
        str(paths["c_empirical"]),
    ]
    if c["execution_chain"] != expected_c_chain:
        raise RuntimeError("Study C execution-chain drift")
    if float(c["alpha_one_sided"]) != 0.025:
        raise RuntimeError("Study C alpha drift")
    if c.get("opportunity_rule") != str(paths["c_opportunity"]):
        raise RuntimeError("Study C opportunity-rule pointer drift")
    c_opp = p["c_opportunity"]
    if c_opp.get("schema") != "ttf_relational_historical_climate_opportunity_rule_v0.1":
        raise RuntimeError("Study C opportunity schema drift")
    inherited = p["b_opportunity"]["directed_geographic_opportunity"]
    bound = c_opp["directed_geographic_opportunity"]
    for key in ("support_radius_km", "minimum_target_coverage", "minimum_source_species_per_target"):
        if bound[key] != inherited[key]:
            raise RuntimeError(f"Study C inherited opportunity constant drift: {key}")
    b_struct = p["b_opportunity"]["structural_gates_before_synthetic_qualification"]
    c_struct = c_opp["structural_gates_before_synthetic_qualification"]
    for key in ("minimum_supported_target_species", "minimum_supported_source_species", "minimum_supported_directed_dyads"):
        if c_struct[key] != b_struct[key]:
            raise RuntimeError(f"Study C inherited structural gate drift: {key}")
    breadth_keys = (
        "largest_single_order_fraction_max",
        "order_fraction_threshold",
        "minimum_orders_at_or_above_fraction_threshold",
    )
    b_relation_breadth = p["b_relation"]["taxonomic_breadth_guardrail"]
    c_relation_breadth = p["c_relation"]["taxonomic_breadth_guardrail"]
    b_opportunity_breadth = b_struct["taxonomic_breadth"]
    c_opportunity_breadth = c_struct["taxonomic_breadth"]
    expected_breadth = {
        "largest_single_order_fraction_max": 0.50,
        "order_fraction_threshold": 0.05,
        "minimum_orders_at_or_above_fraction_threshold": 4,
    }
    for key in breadth_keys:
        expected = expected_breadth[key]
        for label, contract in (
            ("Study B relation", b_relation_breadth),
            ("Study B opportunity", b_opportunity_breadth),
            ("Study C relation", c_relation_breadth),
            ("Study C opportunity", c_opportunity_breadth),
        ):
            if contract[key] != expected:
                raise RuntimeError(f"{label} taxonomic breadth drift: {key}")
    required_groups = {
        "development supported source species",
        "development supported target species",
        "confirmatory supported source species",
        "confirmatory supported target species",
    }
    for label, contract in (
        ("Study B opportunity", b_opportunity_breadth),
        ("Study C opportunity", c_opportunity_breadth),
    ):
        if set(contract["required_groups"]) != required_groups:
            raise RuntimeError(f"{label} supported taxonomic groups drift")
    if float(p["c_qualification"]["inference"]["alpha"]) != 0.025:
        raise RuntimeError("Study C qualification alpha drift")
    if float(p["c_qualification"]["gate_numeric"]["private_type1_wilson95_upper_max"]) != 0.05:
        raise RuntimeError("Study C Type-I Wilson gate drift")
    if float(p["c_qualification"]["gate_numeric"]["historical_power_wilson95_lower_min"]) != 0.80:
        raise RuntimeError("Study C power gate drift")
    if float(p["c_mask"]["required_survivor_requalification"]["same_alpha"]) != 0.025:
        raise RuntimeError("Study C mask requalification alpha drift")
    if float(p["c_empirical"]["relational_model"]["alpha_one_sided"]) != 0.025:
        raise RuntimeError("Study C empirical alpha drift")
    if p["c_empirical"]["relational_model"]["predictors_in_order"][:2] != ["z_R_hist", "z_R_current"]:
        raise RuntimeError("Study C primary/current-environment relation drift")

    c_source_contract = p["c_relation"]["independent_species_domain"]["prior_universe_reproduction"]
    c_s1 = p["c_s1_exclusion"]
    c_s2 = p["c_s2_exclusion"]
    if c_s1.get("schema") != "ttf_relational_prior_S1_species_exclusion_v0.1":
        raise RuntimeError("Study C S1 compact exclusion schema drift")
    if c_s2.get("schema") != "ttf_relational_prior_S2_species_exclusion_v0.1":
        raise RuntimeError("Study C S2 compact exclusion schema drift")
    if int(c_s1["species_count"]) != int(c_source_contract["S1_butterfly_trait"]["expected_species"]):
        raise RuntimeError("Study C S1 compact exclusion count drift")
    if c_s1["source_design_sha256"] != c_source_contract["S1_butterfly_trait"]["expected_design_sha256"]:
        raise RuntimeError("Study C S1 compact source-design hash drift")
    if c_s1["species_list_sha256"] != c_source_contract["S1_butterfly_trait"]["expected_species_list_sha256"]:
        raise RuntimeError("Study C S1 compact species digest drift")
    if int(c_s2["species_count"]) != int(c_source_contract["S2_host_resource_geography"]["expected_species"]):
        raise RuntimeError("Study C S2 compact exclusion count drift")
    if c_s2["source_full_census_sha256"] != c_source_contract["S2_host_resource_geography"]["expected_full_census_sha256"]:
        raise RuntimeError("Study C S2 compact source-census hash drift")
    if c_s2["species_list_sha256"] != c_source_contract["S2_host_resource_geography"]["expected_selected_species_sha256"]:
        raise RuntimeError("Study C S2 compact species digest drift")
    for key in ("c_s1_exclusion", "c_s2_exclusion"):
        assert_closed_firewall(p[key])

    freshness = p["freshness"]
    if int(freshness["exact_overlap_with_study_B_candidates"]["union_overlap_species"]) != 70:
        raise RuntimeError("Study B frozen prior-overlap count drift")
    if p["b_relation"]["freshness_before_pca_and_panel"]["backfill"] is not False:
        raise RuntimeError("Study B backfill firewall drift")
    if int(p["c_relation"]["independent_species_domain"]["minimum_candidate_species"]) != 500:
        raise RuntimeError("Study C candidate minimum drift")
    if p["c_relation"]["independent_species_domain"]["backfill_after_any_C_external_data_or_qualification_result"] is not False:
        raise RuntimeError("Study C backfill firewall drift")

    terminal = set(program["terminal_states"])
    expected_terminal = {
        "DETECTED_B",
        "DETECTED_C_AFTER_B_NULL",
        "DETECTED_C_AFTER_B_NOT_EVALUABLE",
        "CLOSED_TWO_NULLS",
        "CLOSED_PARTIAL_NOT_EVALUABLE",
        "CLOSED_NO_EVALUABLE_TEST",
    }
    if terminal != expected_terminal:
        raise RuntimeError("terminal-state set drift")

    payload = {
        "schema": "ttf_relational_program_v03_closure_audit_v0.1",
        "status": "PASS_FROZEN_V03_CONTRACT_CLOSURE",
        "program_status": program["status"],
        "prospective_slots": ["B", "C"],
        "familywise_alpha": 0.05,
        "slot_alpha": {"B": 0.025, "C": 0.025},
        "current_S3": {
            "decision": receipt["decision"],
            "p_one_sided": receipt["primary"]["p_value_one_sided"],
            "species_union": receipt["empirical_design"]["species_union"],
            "dyads": receipt["empirical_design"]["directed_source_target_dyads"],
        },
        "study_B_full_downstream_contract_frozen": True,
        "study_B_final_transport_retry_round": 4,
        "study_B_no_fifth_transport_retry": True,
        "study_B_final_relation_producer_prefrozen": True,
        "study_B_final_transport_audit_only_prefrozen": True,
        "study_C_full_downstream_contract_frozen": True,
        "all_future_genetic_response_firewalls_closed": True,
        "inputs_sha256": {name: sha256_path(path) for name, path in paths.items()},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: payload[k] for k in (
        "status", "prospective_slots", "familywise_alpha",
        "study_B_full_downstream_contract_frozen",
        "study_C_full_downstream_contract_frozen",
    )}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
