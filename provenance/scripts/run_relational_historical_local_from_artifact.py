#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import zipfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_RUN_ID = 35941577015
EXPECTED_RUN_HEAD = "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab"


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def verify_final_artifact(binding_path: Path, artifact_zip: Path) -> dict:
    binding = load(binding_path)
    if binding.get("schema") != "ttf_relational_historical_occurrence_artifact_binding_v0.2":
        raise RuntimeError("unexpected Study-C occurrence artifact binding")
    if binding.get("status") != "FROZEN_FINAL_BOUNDED_RESPONSE_BLIND_OCCURRENCE_ARTIFACT":
        raise RuntimeError("Study-C occurrence artifact binding is not frozen")
    if int(binding.get("workflow_run_id", -1)) != EXPECTED_RUN_ID:
        raise RuntimeError("Study-C occurrence workflow run drift")
    if binding.get("workflow_head_sha") != EXPECTED_RUN_HEAD:
        raise RuntimeError("Study-C occurrence workflow head drift")
    if binding.get("artifact_name") != "relational-historical-occurrence-final-v0.2":
        raise RuntimeError("Study-C occurrence artifact name drift")
    if binding.get("genetic_response_used") is not False:
        raise RuntimeError("Study-C occurrence binding genetic firewall is open")
    if any(bool(v) for v in binding["response_firewall"].values()):
        raise RuntimeError("Study-C occurrence binding response firewall is open")

    expected = str(binding.get("artifact_digest", ""))
    if not expected.startswith("sha256:"):
        raise RuntimeError("Study-C occurrence binding lacks artifact SHA-256 digest")
    actual = sha256_path(artifact_zip)
    if actual != expected.removeprefix("sha256:"):
        raise RuntimeError("Study-C final artifact ZIP SHA-256 drift")
    return binding


def locate_exact_files(root: Path) -> tuple[Path, Path]:
    csvs = list(root.rglob("occurrences_v0.2.csv"))
    ledgers = list(root.rglob("occurrence_ledger_v0.2.json"))
    if len(csvs) != 1 or len(ledgers) != 1:
        raise RuntimeError(
            f"Study-C final artifact requires exactly one occurrence CSV and ledger; "
            f"found csv={len(csvs)} ledger={len(ledgers)}"
        )
    return csvs[0], ledgers[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--final-artifact-zip", type=Path, required=True)
    ap.add_argument("--occurrence-binding", type=Path, required=True)
    ap.add_argument("--source-archive", type=Path, required=True)
    ap.add_argument("--chelsa-staging-manifest", type=Path, required=True)
    ap.add_argument("--historical-asset", type=Path, action="append", required=True)
    ap.add_argument("--historical-url", action="append", required=True)
    ap.add_argument("--current-bio1", type=Path, required=True)
    ap.add_argument("--current-bio7", type=Path, required=True)
    ap.add_argument("--current-bio12", type=Path, required=True)
    ap.add_argument("--current-bio15", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    artifact_zip = args.final_artifact_zip.resolve()
    binding_path = args.occurrence_binding.resolve()
    work_dir = args.work_dir.resolve()
    output_dir = args.output_dir.resolve()

    binding = verify_final_artifact(binding_path, artifact_zip)

    if work_dir.exists():
        import shutil
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)
    with zipfile.ZipFile(artifact_zip) as archive:
        archive.extractall(work_dir)

    occurrence_csv, occurrence_ledger = locate_exact_files(work_dir)
    if sha256_path(occurrence_csv) != binding["occurrence_csv_sha256"]:
        raise RuntimeError("Study-C occurrence CSV SHA drift from binding")
    if sha256_path(occurrence_ledger) != binding["occurrence_ledger_sha256"]:
        raise RuntimeError("Study-C occurrence ledger SHA drift from binding")

    ledger = load(occurrence_ledger)
    if ledger.get("schema") != "ttf_relational_historical_occurrence_acquisition_v0.2":
        raise RuntimeError("unexpected Study-C final occurrence ledger")
    if ledger.get("status") != binding.get("occurrence_status"):
        raise RuntimeError("Study-C final occurrence status drift from binding")
    if int(ledger.get("species", -1)) != 1000:
        raise RuntimeError("Study-C final occurrence species universe drift")
    if int(ledger.get("maximum_retry_rounds", -1)) != 1:
        raise RuntimeError("Study-C final occurrence retry-count drift")
    if ledger.get("additional_retry_authorized") is not False:
        raise RuntimeError("Study-C final occurrence authorizes forbidden extra retry")
    if ledger.get("scientific_query_change") is not False:
        raise RuntimeError("Study-C final occurrence changed the scientific query")
    if ledger.get("genetic_response_used") is not False:
        raise RuntimeError("Study-C final occurrence ledger genetic firewall is open")
    if any(bool(v) for v in ledger["response_firewall"].values()):
        raise RuntimeError("Study-C final occurrence ledger response firewall is open")

    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts/run_relational_historical_local_pipeline.py"),
        "--source-archive", str(args.source_archive.resolve()),
        "--occurrences", str(occurrence_csv),
        "--occurrence-ledger", str(occurrence_ledger),
        "--occurrence-binding", str(binding_path),
        "--chelsa-staging-manifest", str(args.chelsa_staging_manifest.resolve()),
    ]
    for path in args.historical_asset:
        cmd += ["--historical-asset", str(path.resolve())]
    for url in args.historical_url:
        cmd += ["--historical-url", str(url)]
    cmd += [
        "--current-bio1", str(args.current_bio1.resolve()),
        "--current-bio7", str(args.current_bio7.resolve()),
        "--current-bio12", str(args.current_bio12.resolve()),
        "--current-bio15", str(args.current_bio15.resolve()),
        "--output-dir", str(output_dir),
    ]

    print(json.dumps({
        "status": "PASS_FINAL_OCCURRENCE_ARTIFACT_HANDOFF_TO_FROZEN_LOCAL_STUDY_C",
        "workflow_run_id": binding["workflow_run_id"],
        "artifact_id": binding["artifact_id"],
        "occurrence_status": binding["occurrence_status"],
        "species_passing_ge_30": binding["species_passing_ge_30"],
        "final_request_error_count": binding["final_request_error_count"],
    }, sort_keys=True))
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
