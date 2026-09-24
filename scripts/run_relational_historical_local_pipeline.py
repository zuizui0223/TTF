#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


SOURCE_SHA = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
SOURCE_SIZE = 274_988_692


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob_sha_path(path: Path) -> str:
    data = path.read_bytes()
    h = hashlib.sha1()
    h.update(f"blob {len(data)}\\0".encode())
    h.update(data)
    return h.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def run(*args: str) -> None:
    subprocess.run([sys.executable, *args], check=True)


def verify_implementation_binding(path: Path) -> dict:
    payload = load(path)
    if payload.get("schema") != "ttf_relational_historical_implementation_binding_v0.1":
        raise RuntimeError("unexpected Study-C implementation binding")
    if payload.get("status") != "FROZEN_BEFORE_STUDY_C_RELATION_RESULT":
        raise RuntimeError("Study-C implementation binding is not frozen")
    if payload.get("relation_result_seen") is not False or payload.get("genetic_response_used") is not False:
        raise RuntimeError("Study-C implementation binding firewall is open")
    if any(bool(v) for v in payload["response_firewall"].values()):
        raise RuntimeError("Study-C implementation binding response firewall is open")
    for rel, expected in payload["git_blobs"].items():
        file_path = Path(rel)
        if not file_path.is_file():
            raise FileNotFoundError(f"missing Study-C bound implementation file: {rel}")
        actual = git_blob_sha_path(file_path)
        if actual != expected:
            raise RuntimeError(f"Study-C bound implementation drift: {rel}")
    return payload


def verify_local_execution_rule(path: Path) -> dict:
    payload = load(path)
    if payload.get("schema") != "ttf_relational_historical_local_execution_rule_v0.1":
        raise RuntimeError("unexpected Study-C local execution rule")
    if payload.get("status") != "FROZEN_RESPONSE_BLIND_EXECUTION_EQUIVALENCE_BEFORE_STUDY_C_RELATION_RESULT":
        raise RuntimeError("Study-C local execution rule is not frozen")
    if payload.get("relation_result_seen") is not False or payload.get("genetic_response_used") is not False:
        raise RuntimeError("Study-C local execution rule firewall is open")
    if any(bool(v) for v in payload["response_firewall"].values()):
        raise RuntimeError("Study-C local execution response firewall is open")
    expected = str(payload.get("executor_git_blob", ""))
    actual = git_blob_sha_path(Path(__file__))
    if actual != expected:
        raise RuntimeError("Study-C local executor blob drift")
    return payload


def verify_source_archive(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size != SOURCE_SIZE:
        raise RuntimeError("exact phylogatR archive size drift")
    if sha256_path(path) != SOURCE_SHA:
        raise RuntimeError("exact phylogatR archive SHA-256 drift")


def verify_occurrence_binding(binding_path: Path, occurrence_csv: Path, occurrence_ledger: Path) -> dict:
    binding = load(binding_path)
    if binding.get("schema") != "ttf_relational_historical_occurrence_artifact_binding_v0.2":
        raise RuntimeError("unexpected bounded Study-C occurrence binding")
    if binding.get("status") != "FROZEN_FINAL_BOUNDED_RESPONSE_BLIND_OCCURRENCE_ARTIFACT":
        raise RuntimeError("bounded Study-C occurrence binding is not frozen")
    if int(binding.get("workflow_run_id", -1)) != 35941577015:
        raise RuntimeError("Study-C occurrence run identity drift")
    if binding.get("workflow_head_sha") != "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab":
        raise RuntimeError("Study-C occurrence run head drift")
    if binding.get("genetic_response_used") is not False:
        raise RuntimeError("Study-C occurrence binding genetic firewall is open")
    if any(bool(v) for v in binding["response_firewall"].values()):
        raise RuntimeError("Study-C occurrence binding response firewall is open")
    if sha256_path(occurrence_csv) != binding["occurrence_csv_sha256"]:
        raise RuntimeError("Study-C occurrence CSV SHA drift")
    if sha256_path(occurrence_ledger) != binding["occurrence_ledger_sha256"]:
        raise RuntimeError("Study-C occurrence ledger SHA drift")
    return binding


def final_state(output_dir: Path) -> dict:
    transition = load(Path("benchmarks/frozen/relational_environment_program_transition_v0.3.json"))
    b_state = transition["B_state"]
    empirical = output_dir / "relational_historical_empirical_result_v0.1.json"
    status_paths = {
        "occurrence": output_dir / "occurrence_ledger_v0.2.json",
        "external": output_dir / "relational_historical_external_data_status_v0.1.json",
        "relation": output_dir / "relational_historical_design_v0.1.json",
        "opportunity": output_dir / "relational_historical_opportunity_v0.1.json",
        "development_qualification": output_dir / "relational_historical_qualification_v0.1.json",
        "survivor_qualification": output_dir / "relational_historical_survivor_qualification_v0.1.json",
    }
    c_state = None
    c_decision = None
    c_p = None
    if empirical.is_file():
        result = load(empirical)
        c_decision = result["decision"]
        c_p = result["primary"]["p_value_one_sided"]
        if c_decision == "STUDY_C_HISTORICAL_RELATIONAL_POSITIVE":
            c_state = "POSITIVE"
        elif c_decision == "STUDY_C_HISTORICAL_RELATIONAL_NULL_WITH_QUALIFIED_POWER":
            c_state = "QUALIFIED_NULL"
        else:
            raise RuntimeError(f"unexpected Study-C decision {c_decision}")
    else:
        statuses = {
            key: load(path).get("status")
            for key, path in status_paths.items()
            if path.is_file()
        }
        not_eval = any(
            isinstance(value, str)
            and (
                value.startswith("NOT_EVALUABLE_")
                or value == "NOT_EVALUABLE_C_CHARACTER_MASK_OR_SURVIVOR_GEOMETRY"
            )
            for value in statuses.values()
        )
        c_state = "NOT_EVALUABLE" if not_eval else "INCOMPLETE_TECHNICAL_EXECUTION"

    if c_state == "POSITIVE":
        final = (
            "DETECTED_C_AFTER_B_NULL"
            if b_state == "STUDY_B_ENVIRONMENT_RELATIONAL_NULL_WITH_QUALIFIED_POWER"
            else "DETECTED_C_AFTER_B_NOT_EVALUABLE"
        )
    elif c_state == "QUALIFIED_NULL":
        final = (
            "CLOSED_TWO_NULLS"
            if b_state == "STUDY_B_ENVIRONMENT_RELATIONAL_NULL_WITH_QUALIFIED_POWER"
            else "CLOSED_PARTIAL_NOT_EVALUABLE"
        )
    elif c_state == "NOT_EVALUABLE":
        final = (
            "CLOSED_PARTIAL_NOT_EVALUABLE"
            if b_state == "STUDY_B_ENVIRONMENT_RELATIONAL_NULL_WITH_QUALIFIED_POWER"
            else "CLOSED_NO_EVALUABLE_TEST"
        )
    else:
        final = "INCOMPLETE_TECHNICAL_EXECUTION"
    payload = {
        "schema": "ttf_relational_program_final_state_v0.3",
        "B_state": b_state,
        "C_state": c_state,
        "C_decision": c_decision,
        "C_p_one_sided": c_p,
        "final_state": final,
        "no_additional_predictor_authorized": True,
    }
    (output_dir / "relational_program_final_state_v0.3.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    return payload


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-archive", type=Path, required=True)
    ap.add_argument("--occurrences", type=Path, required=True)
    ap.add_argument("--occurrence-ledger", type=Path, required=True)
    ap.add_argument(
        "--occurrence-binding",
        type=Path,
        default=Path("benchmarks/frozen/relational_historical_occurrence_artifact_binding_v0.2.json"),
    )
    ap.add_argument(
        "--implementation-binding",
        type=Path,
        default=Path("benchmarks/frozen/relational_historical_implementation_binding_v0.1.json"),
    )
    ap.add_argument(
        "--local-execution-rule",
        type=Path,
        default=Path("docs/supporting/relational_historical_local_execution_rule_v0.1.json"),
    )
    ap.add_argument("--historical-asset", type=Path, action="append", required=True)
    ap.add_argument("--historical-url", action="append", required=True)
    ap.add_argument("--current-bio1", type=Path, required=True)
    ap.add_argument("--current-bio7", type=Path, required=True)
    ap.add_argument("--current-bio12", type=Path, required=True)
    ap.add_argument("--current-bio15", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    if len(args.historical_asset) != 8 or len(args.historical_url) != 8:
        raise RuntimeError("Study-C local executor requires exact eight historical assets and URLs")

    verify_local_execution_rule(args.local_execution_rule)
    verify_implementation_binding(args.implementation_binding)
    verify_source_archive(args.source_archive)
    binding = verify_occurrence_binding(
        args.occurrence_binding, args.occurrences, args.occurrence_ledger
    )

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.occurrence_ledger, out / "occurrence_ledger_v0.2.json")
    occurrence = load(args.occurrence_ledger)
    if occurrence.get("schema") != "ttf_relational_historical_occurrence_acquisition_v0.2":
        raise RuntimeError("unexpected Study-C occurrence ledger")
    if occurrence.get("status") != binding.get("occurrence_status"):
        raise RuntimeError("Study-C occurrence status drift from binding")
    if occurrence.get("status") != "PASS_TO_HISTORICAL_ASSET_EXTRACTION":
        print(json.dumps(final_state(out), sort_keys=True))
        return 0

    asset_receipt = out / "historical_asset_receipt_v0.1.json"
    cmd = [
        "scripts/freeze_relational_historical_asset_receipt.py",
        "--rule", "docs/supporting/relational_historical_climate_exposure_rule_v0.1.json",
    ]
    for path in args.historical_asset:
        cmd += ["--asset", str(path)]
    for url in args.historical_url:
        cmd += ["--resolved-url", str(url)]
    cmd += ["--output", str(asset_receipt)]
    run(*cmd)

    design = out / "relational_historical_design_v0.1.npz"
    design_summary = out / "relational_historical_design_v0.1.json"
    run(
        "scripts/build_relational_historical_relation.py",
        "--occurrences", str(args.occurrences),
        "--candidates", "benchmarks/frozen/relational_historical_candidates_v0.1.csv",
        "--rule", "docs/supporting/relational_historical_climate_exposure_rule_v0.1.json",
        "--asset-receipt", str(asset_receipt),
        *sum([["--historical-asset", str(p)] for p in args.historical_asset], []),
        "--current-bio1", str(args.current_bio1),
        "--current-bio7", str(args.current_bio7),
        "--current-bio12", str(args.current_bio12),
        "--current-bio15", str(args.current_bio15),
        "--output-npz", str(design),
        "--output-summary", str(design_summary),
    )
    if load(design_summary).get("status") != "PASS_RESPONSE_BLIND_HISTORICAL_RELATION_DESIGN":
        print(json.dumps(final_state(out), sort_keys=True))
        return 0

    source_root = out / "exact-source"
    if source_root.exists():
        shutil.rmtree(source_root)
    source_root.mkdir(parents=True)
    with zipfile.ZipFile(args.source_archive) as archive:
        archive.extractall(source_root)
    genes = next(source_root.rglob("genes.txt"))
    root = genes.parent
    if not (root / "cite.txt").is_file():
        raise RuntimeError("exact source extraction lacks cite.txt")

    candidate_out = out / "candidate-geometry"
    candidate_out.mkdir(parents=True, exist_ok=True)
    localities = candidate_out / "relational_historical_localities_v0.1.csv"
    edges = candidate_out / "relational_historical_edges_v0.1.csv"
    geometry_receipt = candidate_out / "relational_historical_geometry_v0.1.json"
    run(
        "scripts/reconstruct_relational_historical_geometry.py",
        "--root", str(root),
        "--candidates", "benchmarks/frozen/relational_historical_candidates_v0.1.csv",
        "--history-rule", "docs/supporting/relational_historical_climate_exposure_rule_v0.1.json",
        "--neighbor-fraction", "0.15",
        "--output-localities", str(localities),
        "--output-edges", str(edges),
        "--output-receipt", str(geometry_receipt),
    )

    independent = load(Path("benchmarks/frozen/relational_historical_geometry_local_reconstruction_v0.1.json"))
    if sha256_path(localities) != independent["geometry"]["localities_csv_sha256"]:
        raise RuntimeError("Study-C locality geometry SHA drift")
    if sha256_path(edges) != independent["geometry"]["edges_csv_sha256"]:
        raise RuntimeError("Study-C edge geometry SHA drift")

    opportunity = out / "relational_historical_opportunity_v0.1.npz"
    opportunity_summary = out / "relational_historical_opportunity_v0.1.json"
    run(
        "scripts/attach_relational_historical_opportunity.py",
        "--historical-design", str(design),
        "--edges", str(edges),
        "--candidates", "benchmarks/frozen/relational_historical_candidates_v0.1.csv",
        "--opportunity-rule", "docs/supporting/relational_historical_climate_opportunity_rule_v0.1.json",
        "--output-npz", str(opportunity),
        "--output-summary", str(opportunity_summary),
    )
    if load(opportunity_summary).get("status") != "PASS_TO_HISTORICAL_DEVELOPMENT_SYNTHETIC_QUALIFICATION":
        print(json.dumps(final_state(out), sort_keys=True))
        return 0

    qualification = out / "relational_historical_qualification_v0.1.json"
    run(
        "scripts/run_relational_historical_climate_qualification.py",
        "--opportunity-design", str(opportunity),
        "--opportunity-summary", str(opportunity_summary),
        "--rule", "docs/supporting/relational_historical_climate_qualification_rule_v0.1.json",
        "--output", str(qualification),
    )
    q = load(qualification)
    if not (
        q.get("status") == "PASS_TO_HISTORICAL_CONFIRMATORY_CHARACTER_MASK_PREPARATION"
        and q.get("gates", {}).get("overall_pass") is True
    ):
        print(json.dumps(final_state(out), sort_keys=True))
        return 0

    mask_result = out / "relational_historical_character_mask_result_v0.1.json"
    survivor_design = out / "relational_historical_survivor_design_v0.1.npz"
    survivor_summary = out / "relational_historical_survivor_summary_v0.1.json"
    run(
        "scripts/freeze_relational_historical_character_mask.py",
        "--root", str(root),
        "--opportunity-design", str(opportunity),
        "--opportunity-summary", str(opportunity_summary),
        "--development-qualification", str(qualification),
        "--candidates", "benchmarks/frozen/relational_historical_candidates_v0.1.csv",
        "--mask-rule", "docs/supporting/relational_historical_climate_character_mask_rule_v0.1.json",
        "--output-mask-result", str(mask_result),
        "--output-survivor-design", str(survivor_design),
        "--output-survivor-summary", str(survivor_summary),
    )
    if not survivor_design.is_file():
        print(json.dumps(final_state(out), sort_keys=True))
        return 0

    survivor_qualification = out / "relational_historical_survivor_qualification_v0.1.json"
    run(
        "scripts/run_relational_historical_climate_qualification.py",
        "--opportunity-design", str(survivor_design),
        "--opportunity-summary", str(survivor_summary),
        "--rule", "docs/supporting/relational_historical_climate_qualification_rule_v0.1.json",
        "--qualification-stage", "confirmatory_survivor",
        "--output", str(survivor_qualification),
    )
    sq = load(survivor_qualification)
    if not (
        sq.get("status") == "PASS_TO_HISTORICAL_EMPIRICAL_IDENTITY_OPENING_PREPARATION"
        and sq.get("gates", {}).get("overall_pass") is True
    ):
        print(json.dumps(final_state(out), sort_keys=True))
        return 0

    authorization = out / "relational_historical_identity_opening_authorization_v0.1.json"
    run(
        "scripts/freeze_relational_historical_empirical_authorization.py",
        "--root", str(root),
        "--source-archive-sha256", SOURCE_SHA,
        "--survivor-design", str(survivor_design),
        "--survivor-summary", str(survivor_summary),
        "--mask-result", str(mask_result),
        "--survivor-qualification", str(survivor_qualification),
        "--candidates", "benchmarks/frozen/relational_historical_candidates_v0.1.csv",
        "--opening-rule", "docs/supporting/relational_historical_climate_empirical_opening_rule_v0.1.json",
        "--output", str(authorization),
    )

    empirical = out / "relational_historical_empirical_result_v0.1.json"
    run(
        "scripts/run_relational_historical_empirical.py",
        "--root", str(root),
        "--source-archive-sha256", SOURCE_SHA,
        "--survivor-design", str(survivor_design),
        "--survivor-summary", str(survivor_summary),
        "--mask-result", str(mask_result),
        "--survivor-qualification", str(survivor_qualification),
        "--candidates", "benchmarks/frozen/relational_historical_candidates_v0.1.csv",
        "--opening-rule", "docs/supporting/relational_historical_climate_empirical_opening_rule_v0.1.json",
        "--authorization", str(authorization),
        "--output", str(empirical),
    )

    print(json.dumps(final_state(out), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
