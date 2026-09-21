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
        "freshness": Path("docs/supporting/relational_future_family_freshness_amendment_v0.1.json"),
        "s3_receipt": Path("benchmarks/frozen/relational_host_resource_empirical_repaired_receipt_v0.2.json"),
        "s3_reproduction": Path("benchmarks/frozen/relational_host_resource_empirical_reproduction_v0.2.json"),
        "s3_overlap": Path("benchmarks/frozen/relational_prior_host_test_overlap_audit_v0.1.json"),
        "b_relation": Path("docs/supporting/relational_environment_relation_rule_v0.3.json"),
        "b_opportunity": Path("docs/supporting/relational_environment_opportunity_rule_v0.2.json"),
        "b_qualification": Path("docs/supporting/relational_environment_qualification_rule_v0.2.json"),
        "b_mask": Path("docs/supporting/relational_environment_character_mask_rule_v0.1.json"),
        "b_empirical": Path("docs/supporting/relational_environment_empirical_opening_rule_v0.1.json"),
        "c_relation": Path("docs/supporting/relational_historical_climate_exposure_rule_v0.1.json"),
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
    for key in ("b_relation", "b_qualification", "b_mask", "b_empirical", "c_relation", "c_qualification", "c_mask", "c_empirical"):
        assert_closed_firewall(p[key])

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
        str(paths["c_qualification"]),
        str(paths["c_mask"]),
        str(paths["c_empirical"]),
    ]
    if c["execution_chain"] != expected_c_chain:
        raise RuntimeError("Study C execution-chain drift")
    if float(c["alpha_one_sided"]) != 0.025:
        raise RuntimeError("Study C alpha drift")
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
