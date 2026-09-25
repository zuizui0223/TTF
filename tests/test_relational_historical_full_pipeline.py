from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_historical_qualification_runner_has_survivor_stage():
    text = (ROOT / "scripts/run_relational_historical_climate_qualification.py").read_text()
    assert '--qualification-stage' in text
    assert 'choices=("development", "confirmatory_survivor")' in text
    assert '"PASS_TO_HISTORICAL_EMPIRICAL_IDENTITY_OPENING_PREPARATION"' in text
    assert '"NOT_EVALUABLE_C_CHARACTER_MASK_OR_SURVIVOR_GEOMETRY"' in text
    assert '"qualification_stage": args.qualification_stage' in text


def test_historical_workflow_is_gated_by_B_transition_and_closes_program():
    text = (ROOT / ".github/workflows/relational-historical-prequalification.yml").read_text()
    assert "Study-C opening forbidden" in text
    assert "--qualification-stage confirmatory_survivor" in text
    assert "freeze_relational_historical_empirical_authorization.py" in text
    assert "run_relational_historical_empirical.py" in text
    assert '"DETECTED_C_AFTER_B_NULL"' in text
    assert '"DETECTED_C_AFTER_B_NOT_EVALUABLE"' in text
    assert '"CLOSED_TWO_NULLS"' in text
    assert '"CLOSED_PARTIAL_NOT_EVALUABLE"' in text
    assert '"CLOSED_NO_EVALUABLE_TEST"' in text


def test_bounded_occurrence_binding_and_full_c_handoff_are_frozen():
    binder = (ROOT / ".github/workflows/relational-historical-bind-bounded-occurrence.yml").read_text()
    downstream = (ROOT / ".github/workflows/relational-historical-from-bounded-binding-v02.yml").read_text()

    assert "relational-historical-occurrence-bounded-v02" in binder
    assert "relational-historical-occurrence-final-v0.2" in binder
    assert "ttf_relational_historical_occurrence_acquisition_v0.2" in binder
    assert "FROZEN_FINAL_BOUNDED_RESPONSE_BLIND_OCCURRENCE_ARTIFACT" in binder

    assert "relational_historical_occurrence_artifact_binding_v0.2.json" in downstream
    assert "reconstruct_relational_historical_geometry.py" in downstream
    assert "run_relational_historical_climate_qualification.py" in downstream
    assert "--qualification-stage confirmatory_survivor" in downstream
    assert "freeze_relational_historical_empirical_authorization.py" in downstream
    assert "run_relational_historical_empirical.py" in downstream
    assert '"DETECTED_C_AFTER_B_NOT_EVALUABLE"' in downstream
    assert '"CLOSED_PARTIAL_NOT_EVALUABLE"' in downstream
    assert '"CLOSED_NO_EVALUABLE_TEST"' in downstream
    assert '"no_additional_predictor_authorized":True' in downstream


def test_historical_geometry_materialization_is_hash_bound():
    import json
    rule = json.loads(
        (ROOT / "docs/supporting/relational_historical_geometry_materialization_rule_v0.1.json").read_text()
    )
    receipt = json.loads(
        (ROOT / "benchmarks/frozen/relational_historical_geometry_reproduction_v0.1.json").read_text()
    )
    independent = json.loads(
        (ROOT / "benchmarks/frozen/relational_historical_geometry_local_reconstruction_v0.1.json").read_text()
    )
    reproducer = (ROOT / "scripts/reconstruct_relational_historical_geometry_snapshot.py").read_text()

    assert rule["status"] == "FROZEN_BEFORE_STUDY_C_RELATION_RESULT"
    assert rule["qualification_without_exact_geometry_allowed"] is False
    assert rule["genetic_opening_without_exact_geometry_allowed"] is False
    assert rule["alternate_geometry_source_allowed"] is False

    exact = rule["exact_geometry"]
    assert exact["species"] == receipt["species"] == independent["geometry"]["species"] == 1000
    assert exact["locality_rows"] == receipt["locality_rows"] == independent["geometry"]["locality_rows"] == 28735
    assert exact["edge_rows"] == receipt["edge_rows"] == independent["geometry"]["edge_rows"] == 143642
    assert exact["locality_csv_sha256"] == receipt["locality_csv_sha256"] == independent["geometry"]["localities_csv_sha256"] == "cbc09656be79b38ea10b38473c029a242023d1e179a2835a531fddf0deebd4a5"
    assert exact["edge_csv_sha256"] == receipt["edge_csv_sha256"] == independent["geometry"]["edges_csv_sha256"] == "f142467873eb5aec08230e436ff1fc92ac0cbfc8b9a1c1797e9ea1d31c5b98f2"

    assert "snapshot_base64_file_sha256" in reproducer
    assert "snapshot_gzip_sha256" in reproducer
    assert "snapshot_json_sha256" in reproducer
    assert "historical reconstructed locality CSV SHA drift" in reproducer
    assert "historical reconstructed edge CSV SHA drift" in reproducer


def test_bounded_handoff_requires_exact_source_transport_before_relation_but_delays_extraction():
    text = (ROOT / ".github/workflows/relational-historical-from-bounded-binding-v02.yml").read_text()
    cache_gate = text.index("Hold Study C unless exact source archive is runner-readable")
    hist_assets = text.index("Resolve and freeze exact CHELSA-TraCE21k assets")
    current_assets = text.index("Download frozen current CHELSA V2.1 nuisance layers")
    relation = text.index("Build frozen historical/current relation")
    extraction = text.index("Require and extract exact archive")
    geometry = text.index("Reconstruct frozen Study-C response-blind geometry")
    opportunity = text.index("Attach frozen historical opportunity geometry")

    assert cache_gate < hist_assets < current_assets < relation < extraction < geometry < opportunity
    assert "relation_construction_before_geometry_materialization_allowed" not in text
    # Transport readiness is checked before any relation is opened. Archive
    # extraction and exact genetic geometry remain delayed until the frozen
    # historical/current relation itself has passed.


def test_historical_archive_hash_clerical_repair_is_explicit_and_downstream_is_canonical():
    import json

    canonical = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
    legacy = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce61a5"
    workflow = (
        ROOT / ".github/workflows/relational-historical-from-bounded-binding-v02.yml"
    ).read_text()
    correction = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_archive_hash_provenance_correction_v0.1.json"
        ).read_text()
    )
    local = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_geometry_local_reconstruction_v0.1.json"
        ).read_text()
    )
    reproduction = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_geometry_reproduction_v0.1.json"
        ).read_text()
    )
    materialization = json.loads(
        (
            ROOT
            / "docs/supporting/relational_historical_geometry_materialization_rule_v0.1.json"
        ).read_text()
    )

    assert len(canonical) == 64
    assert len(legacy) != 64
    assert legacy not in workflow
    assert canonical in workflow
    assert materialization["source_archive_sha256"] == canonical

    # Preserve the original response-blind receipts as immutable audit history.
    # Their truncated literal is repaired only through this explicit provenance map.
    assert local["source_archive"]["sha256"] == legacy
    assert reproduction["source_archive"]["sha256"] == canonical
    assert correction["legacy_malformed_literal"] == legacy
    assert correction["canonical_source_archive_sha256"] == canonical
    assert correction["affected_immutable_receipts"] == [
        "benchmarks/frozen/relational_historical_geometry_local_reconstruction_v0.1.json",
    ]
    assert correction["canonical_receipts_already_correct"] == [
        "benchmarks/frozen/relational_historical_geometry_reproduction_v0.1.json",
    ]
    assert correction["relation_result_seen"] is False
    assert correction["genetic_response_used"] is False
    assert correction["scientific_rule_change"] is False
    assert not (
        ROOT
        / "benchmarks/frozen/relational_historical_geometry_snapshot_v0.1/part00.txt"
    ).exists()


def test_historical_downstream_is_reusable_and_checks_out_study_branch():
    import json

    downstream = (
        ROOT / ".github/workflows/relational-historical-from-bounded-binding-v02.yml"
    ).read_text()
    receipt = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_downstream_activation_repair_v0.1.json"
        ).read_text()
    )

    binder = (
        ROOT / ".github/workflows/relational-historical-bind-bounded-occurrence.yml"
    ).read_text()

    assert "workflow_call:" in downstream
    assert "ref: relational-ttf-v01" in downstream
    assert "needs: bind" in binder
    assert (
        "uses: zuizui0223/TTF/.github/workflows/relational-historical-from-bounded-binding-v02.yml@"
        "17923b1719d9a3ebce842c959e4547378ac28338"
    ) in binder
    assert receipt["status"] == "FROZEN_RESPONSE_BLIND_ORCHESTRATION_REPAIR_BEFORE_STUDY_C_RELATION_RESULT"
    assert receipt["scientific_rule_change"] is False
    assert receipt["predictor_change"] is False
    assert receipt["alpha_change"] is False
    assert receipt["called_workflow_commit_sha"] == "17923b1719d9a3ebce842c959e4547378ac28338"
    assert receipt["default_branch_binder_commit_sha"] == "0b5bded97d7e02db86f59fbd89d3114b19a651b2"
    assert receipt["bot_push_is_not_activation_dependency"] is True
    assert receipt["relation_result_seen"] is False
    assert receipt["genetic_response_used"] is False
    assert all(v is False for v in receipt["response_firewall"].values())


def test_historical_implementation_binding_is_enforced_before_binding_open():
    import json

    receipt = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_implementation_binding_v0.1.json"
        ).read_text()
    )
    downstream = (
        ROOT / ".github/workflows/relational-historical-from-bounded-binding-v02.yml"
    ).read_text()

    assert receipt["status"] == "FROZEN_BEFORE_STUDY_C_RELATION_RESULT"
    assert len(receipt["git_blobs"]) == 31
    assert "scripts/build_relational_historical_relation.py" in receipt["git_blobs"]
    assert "scripts/run_relational_historical_climate_qualification.py" in receipt["git_blobs"]
    assert "scripts/run_relational_historical_empirical.py" in receipt["git_blobs"]
    assert "src/ttf/relational_dyadic.py" in receipt["git_blobs"]
    assert "Verify frozen Study-C implementation binding" in downstream
    assert 'git","rev-parse",f"HEAD:{path}"' in downstream
    assert receipt["relation_result_seen"] is False
    assert receipt["genetic_response_used"] is False
    assert all(v is False for v in receipt["response_firewall"].values())


def test_historical_source_transport_holds_before_relation_without_closing_c():
    import json

    rule = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_source_transport_hold_rule_v0.1.json"
        ).read_text()
    )
    downstream = (
        ROOT / ".github/workflows/relational-historical-from-bounded-binding-v02.yml"
    ).read_text()

    cache_gate = downstream.index("Hold Study C unless exact source archive is runner-readable")
    occurrence = downstream.index("Verify exact bounded occurrence artifact")
    historical_assets = downstream.index("Resolve and freeze exact CHELSA-TraCE21k assets")
    relation = downstream.index("Build frozen historical/current relation")
    assert cache_gate < occurrence < historical_assets < relation

    assert "HOLD_EXACT_SOURCE_ARCHIVE_CACHE_UNAVAILABLE" in downstream
    assert "HOLD_EXACT_SOURCE_ARCHIVE_CACHE_MISMATCH" in downstream
    assert 'c_state="HOLD_TECHNICAL_SOURCE_TRANSPORT"' in downstream
    assert 'final="HOLD_C_SOURCE_TRANSPORT"' in downstream
    terminal_block = downstream.split('terminal={',1)[1].split('}',1)[0]
    assert '"HOLD_C_SOURCE_TRANSPORT"' not in terminal_block

    assert rule["terminal"] is False
    assert rule["biological_null"] is False
    assert rule["consumes_C_slot"] is False
    assert rule["alpha_consumed"] is False
    assert rule["relation_result_seen"] is False
    assert rule["genetic_response_used"] is False
    assert all(v is False for v in rule["response_firewall"].values())


def test_historical_cache_resume_cannot_select_a_new_run_or_reopen_terminal_c():
    import json

    resume = (
        ROOT / ".github/workflows/relational-historical-resume-after-source-cache.yml"
    ).read_text()
    receipt = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_downstream_activation_repair_v0.1.json"
        ).read_text()
    )

    assert '"seed frozen phylogatr archive cache"' in resume
    assert 'b["workflow_run_id"]==35941577015' in resume
    assert 'b["workflow_head_sha"]=="ffacbd51d58689a4b18f7a2cb920f5a1385a74ab"' in resume
    assert '"DETECTED_C_AFTER_B_NOT_EVALUABLE"' in resume
    assert '"CLOSED_PARTIAL_NOT_EVALUABLE"' in resume
    assert (
        "uses: zuizui0223/TTF/.github/workflows/relational-historical-from-bounded-binding-v02.yml@"
        "17923b1719d9a3ebce842c959e4547378ac28338"
    ) in resume
    assert receipt["resume_result_selection_allowed"] is False
    assert receipt["default_branch_resume_commit_sha"] == "b80a248242646e69c1ed060049b09b6e453ce6c1"


def test_local_study_c_executor_preserves_one_shot_order_and_exact_inputs():
    import json

    text = (ROOT / "scripts/run_relational_historical_local_pipeline.py").read_text()
    rule = json.loads(
        (
            ROOT
            / "docs/supporting/relational_historical_local_execution_rule_v0.1.json"
        ).read_text()
    )

    verify_impl = text.index("verify_implementation_binding(args.implementation_binding)")
    verify_source = text.index("verify_source_archive(args.source_archive)")
    verify_binding = text.index("binding = verify_occurrence_binding(")
    relation = text.index('"scripts/build_relational_historical_relation.py"')
    extract = text.index("with zipfile.ZipFile(args.source_archive) as archive:")
    development_q = text.index('"scripts/run_relational_historical_climate_qualification.py"')
    mask = text.index('"scripts/freeze_relational_historical_character_mask.py"')
    survivor_q = text.index('"--qualification-stage", "confirmatory_survivor"')
    authorization = text.index('"scripts/freeze_relational_historical_empirical_authorization.py"')
    empirical = text.index('"scripts/run_relational_historical_empirical.py"')

    assert verify_impl < verify_source < verify_binding < relation < extract
    assert relation < development_q < mask < survivor_q < authorization < empirical
    assert "35941577015" in text
    assert "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab" in text
    assert "274_988_692" in text
    assert "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5" in text

    assert rule["status"] == "FROZEN_RESPONSE_BLIND_EXECUTION_EQUIVALENCE_BEFORE_STUDY_C_RELATION_RESULT"
    assert rule["result_selection_rerun_allowed"] is False
    assert rule["alternate_predictor_allowed"] is False
    assert rule["source_repacking_allowed"] is False
    assert rule["executor_git_blob"]
    assert rule["source_bundle_compatible"] is True
    assert "git_blob_sha_path" in text
    assert "git rev-parse" not in text
    assert 'env["PYTHONPATH"]' in text
    assert "cwd=REPO_ROOT" in text
    assert rule["editable_install_required"] is False
    assert rule["chelsa_staging_manifest_required"] is True
    assert "verify_chelsa_staging_manifest" in text
    assert rule["chelsa_staging_manifest_sha256"] == "7d82833996006fc1e99b86acb4ea7e80cde29c7f500df739ef48cb9ce145d450"
    assert rule["chelsa_staging_manifest_source_commit"] == "4ea89a703bce2ed88c715cb7cfcca94d1e3b27e4"
    assert "Study-C CHELSA staging manifest SHA-256 drift" in text
    assert rule["relation_result_seen"] is False
    assert rule["genetic_response_used"] is False
    assert all(v is False for v in rule["response_firewall"].values())
    compile(text, "scripts/run_relational_historical_local_pipeline.py", "exec")


def test_local_study_c_final_state_matches_current_v03_transition(tmp_path):
    import importlib.util
    import json

    script = ROOT / "scripts/run_relational_historical_local_pipeline.py"
    spec = importlib.util.spec_from_file_location("hist_local_pipeline", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    # Current frozen transition is B NOT_EVALUABLE, so C positive/null/gate
    # failure have three distinct licensed outcomes.
    empirical = tmp_path / "relational_historical_empirical_result_v0.1.json"
    empirical.write_text(
        json.dumps(
            {
                "decision": "STUDY_C_HISTORICAL_RELATIONAL_POSITIVE",
                "primary": {"p_value_one_sided": 0.01},
            }
        )
    )
    assert module.final_state(tmp_path)["final_state"] == "DETECTED_C_AFTER_B_NOT_EVALUABLE"

    empirical.write_text(
        json.dumps(
            {
                "decision": "STUDY_C_HISTORICAL_RELATIONAL_NULL_WITH_QUALIFIED_POWER",
                "primary": {"p_value_one_sided": 0.40},
            }
        )
    )
    assert module.final_state(tmp_path)["final_state"] == "CLOSED_PARTIAL_NOT_EVALUABLE"

    empirical.unlink()
    (tmp_path / "relational_historical_qualification_v0.1.json").write_text(
        json.dumps({"status": "NOT_EVALUABLE_HISTORICAL_SYNTHETIC_QUALIFICATION"})
    )
    result = module.final_state(tmp_path)
    assert result["C_state"] == "NOT_EVALUABLE"
    assert result["final_state"] == "CLOSED_NO_EVALUABLE_TEST"
    assert result["no_additional_predictor_authorized"] is True


def test_local_study_c_blob_binding_matches_without_git_metadata():
    import importlib.util
    import json

    script = ROOT / "scripts/run_relational_historical_local_pipeline.py"
    spec = importlib.util.spec_from_file_location("hist_local_blob_binding", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    binding = json.loads(
        (ROOT / "benchmarks/frozen/relational_historical_implementation_binding_v0.1.json").read_text()
    )
    for rel, expected in binding["git_blobs"].items():
        assert module.git_blob_sha_path(ROOT / rel) == expected

    rule = json.loads(
        (ROOT / "docs/supporting/relational_historical_local_execution_rule_v0.1.json").read_text()
    )
    assert module.git_blob_sha_path(script) == rule["executor_git_blob"]


def test_local_study_c_subprocess_imports_from_source_bundle_without_install():
    import os
    import subprocess
    import sys

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/build_relational_historical_relation.py",
            "--help",
        ],
        cwd=ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "--occurrences" in completed.stdout


def test_local_study_c_chelsa_staging_manifest_detects_asset_drift(tmp_path):
    import hashlib
    import importlib.util
    import json
    import pytest

    script = ROOT / "scripts/run_relational_historical_local_pipeline.py"
    spec = importlib.util.spec_from_file_location("hist_local_chelsa", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    staging_rule = json.loads(
        (ROOT / "docs/supporting/relational_historical_chelsa_asset_staging_rule_v0.1.json").read_text()
    )
    hist = []
    hist_paths = []
    hist_urls = []
    for i, name in enumerate(staging_rule["historical"]["exact_logical_filenames"]):
        path = tmp_path / name
        path.write_bytes(f"historical-{i}".encode())
        url = f"https://example.invalid/{name}"
        hist_paths.append(path)
        hist_urls.append(url)
        hist.append(
            {
                "filename": name,
                "resolved_url": url,
                "size_bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    current = []
    current_paths = []
    for i, name in enumerate(staging_rule["current"]["exact_filenames"]):
        path = tmp_path / name
        path.write_bytes(f"current-{i}".encode())
        current_paths.append(path)
        current.append(
            {
                "filename": name,
                "resolved_url": f"https://example.invalid/{name}",
                "size_bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "ttf_relational_historical_chelsa_asset_staging_manifest_v0.1",
                "status": "PASS_EXACT_RESPONSE_BLIND_CHELSA_ASSET_STAGING",
                "historical": hist,
                "current": current,
                "relation_result_seen": False,
                "genetic_response_used": False,
                "response_firewall": {
                    "Study_C_sequence_identity_opened": False,
                    "Study_C_pairwise_genetic_distances_opened": False,
                    "Study_C_T_st_computed": False,
                    "Study_C_beta_hist_computed": False,
                },
            }
        )
    )

    module.verify_chelsa_staging_manifest(manifest, hist_paths, hist_urls, current_paths)
    hist_paths[0].write_bytes(b"drift")
    with pytest.raises(RuntimeError, match="size drift|SHA drift"):
        module.verify_chelsa_staging_manifest(manifest, hist_paths, hist_urls, current_paths)


def test_chelsa_staging_and_chunk_transport_are_bound_response_blind():
    import json

    staging_rule = json.loads(
        (
            ROOT
            / "docs/supporting/relational_historical_chelsa_asset_staging_rule_v0.1.json"
        ).read_text()
    )
    chunk_rule = json.loads(
        (
            ROOT
            / "docs/supporting/relational_historical_chelsa_artifact_chunk_transport_v0.1.json"
        ).read_text()
    )
    staging_workflow = (
        ROOT / ".github/workflows/relational-c-chelsa-asset-staging.yml"
    ).read_text()
    chunk_workflow = (
        ROOT / ".github/workflows/relational-c-chelsa-chunk-transport.yml"
    ).read_text()

    assert staging_rule["status"] == "FROZEN_RESPONSE_BLIND_ASSET_TRANSPORT_BEFORE_STUDY_C_RELATION_RESULT"
    assert staging_rule["alternate_asset_selection_allowed"] is False
    assert staging_rule["relation_result_seen"] is False
    assert staging_rule["genetic_response_used"] is False
    assert len(staging_rule["historical"]["exact_logical_filenames"]) == 8
    assert len(staging_rule["current"]["exact_filenames"]) == 4
    assert all(v is False for v in staging_rule["response_firewall"].values())

    assert chunk_rule["status"] == "FROZEN_TRANSPORT_ONLY_FROM_EXACT_STAGED_ARTIFACTS"
    assert chunk_rule["source_run"]["workflow_run_id"] == 35959224240
    assert chunk_rule["source_run"]["workflow_head_sha"] == "4ea89a703bce2ed88c715cb7cfcca94d1e3b27e4"
    assert chunk_rule["source_artifacts"]["historical"]["artifact_id"] == 10791389409
    assert chunk_rule["source_artifacts"]["historical"]["digest_sha256"] == "0a797c609d26a185ba08d209a71203a584fa60f1045c4a8c335934372cee30b2"
    assert chunk_rule["source_artifacts"]["historical"]["expected_chunks"] == 24
    assert chunk_rule["source_artifacts"]["current"]["artifact_id"] == 10791867875
    assert chunk_rule["source_artifacts"]["current"]["digest_sha256"] == "724743d8666d32bc8717e52200a798369665c2f1965ba33021d3e5c68c5eff4c"
    assert chunk_rule["source_artifacts"]["current"]["expected_chunks"] == 6
    assert chunk_rule["chunk_bytes"] == 200000000
    assert chunk_rule["scientific_asset_change"] is False
    assert chunk_rule["alternate_asset_selection_allowed"] is False
    assert chunk_rule["relation_result_seen"] is False
    assert chunk_rule["genetic_response_used"] is False

    assert "PASS_EXACT_RESPONSE_BLIND_CHELSA_ASSET_STAGING" in staging_workflow
    assert "10791389409" in chunk_workflow
    assert "10791867875" in chunk_workflow
    assert "0a797c609d26a185ba08d209a71203a584fa60f1045c4a8c335934372cee30b2" in chunk_workflow
    assert "724743d8666d32bc8717e52200a798369665c2f1965ba33021d3e5c68c5eff4c" in chunk_workflow
    assert "split -b 200000000" in chunk_workflow


def test_local_final_occurrence_artifact_handoff_is_digest_bound_and_run_fixed():
    import json

    wrapper = (
        ROOT / "scripts/run_relational_historical_local_from_artifact.py"
    ).read_text()
    rule = json.loads(
        (
            ROOT
            / "docs/supporting/relational_historical_local_execution_rule_v0.1.json"
        ).read_text()
    )
    source_bundle = (
        ROOT / ".github/workflows/relational-c-local-source-bundle.yml"
    ).read_text()

    assert "EXPECTED_RUN_ID = 35941577015" in wrapper
    assert 'EXPECTED_RUN_HEAD = "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab"' in wrapper
    assert 'binding.get("artifact_name") != "relational-historical-occurrence-final-v0.2"' in wrapper
    assert 'binding.get("artifact_digest", "")' in wrapper
    assert "Study-C final artifact ZIP SHA-256 drift" in wrapper
    assert 'binding["occurrence_csv_sha256"]' in wrapper
    assert 'binding["occurrence_ledger_sha256"]' in wrapper
    assert 'ledger.get("maximum_retry_rounds", -1)' in wrapper
    assert 'ledger.get("additional_retry_authorized") is not False' in wrapper
    assert "scripts/run_relational_historical_local_pipeline.py" in wrapper

    assert rule["final_artifact_handoff_executor"] == "scripts/run_relational_historical_local_from_artifact.py"
    assert rule["final_artifact_handoff_executor_git_blob"] == "487de64b3e4abe411344edeb992a6294252c191c"
    assert rule["final_artifact_alternative_run_allowed"] is False
    assert rule["final_artifact_repacking_allowed"] is False
    assert "artifact_digest" in rule["final_artifact_handoff_contract"]
    assert "run_relational_historical_local_from_artifact.py" in source_bundle
    compile(wrapper, "scripts/run_relational_historical_local_from_artifact.py", "exec")


def test_local_final_artifact_handoff_runtime_rejects_digest_and_inner_hash_drift(tmp_path):
    import hashlib
    import importlib.util
    import json
    import zipfile

    script = ROOT / "scripts/run_relational_historical_local_from_artifact.py"
    spec = importlib.util.spec_from_file_location("hist_local_artifact_handoff", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    payload_dir = tmp_path / "payload"
    payload_dir.mkdir()
    occ = payload_dir / "occurrences_v0.2.csv"
    ledger = payload_dir / "occurrence_ledger_v0.2.json"
    occ.write_text("species,source_key,latitude,longitude,priority_rank,priority_sha256\n")
    ledger.write_text(
        json.dumps(
            {
                "schema": "ttf_relational_historical_occurrence_acquisition_v0.2",
                "status": "PASS_TO_HISTORICAL_ASSET_EXTRACTION",
                "species": 1000,
                "maximum_retry_rounds": 1,
                "additional_retry_authorized": False,
                "scientific_query_change": False,
                "genetic_response_used": False,
                "response_firewall": {
                    "Study_C_sequence_identity_opened": False,
                    "Study_C_pairwise_genetic_distances_opened": False,
                    "Study_C_T_st_computed": False,
                    "Study_C_beta_hist_computed": False,
                },
            }
        )
    )

    artifact = tmp_path / "final.zip"
    with zipfile.ZipFile(artifact, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.write(occ, arcname="occurrences_v0.2.csv")
        z.write(ledger, arcname="occurrence_ledger_v0.2.json")

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    binding = tmp_path / "binding.json"
    binding_payload = {
        "schema": "ttf_relational_historical_occurrence_artifact_binding_v0.2",
        "status": "FROZEN_FINAL_BOUNDED_RESPONSE_BLIND_OCCURRENCE_ARTIFACT",
        "workflow_run_id": 35941577015,
        "workflow_head_sha": "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab",
        "artifact_id": 123,
        "artifact_name": "relational-historical-occurrence-final-v0.2",
        "artifact_digest": "sha256:" + sha(artifact),
        "occurrence_csv_sha256": sha(occ),
        "occurrence_ledger_sha256": sha(ledger),
        "occurrence_status": "PASS_TO_HISTORICAL_ASSET_EXTRACTION",
        "species_passing_ge_30": 999,
        "final_request_error_count": 0,
        "genetic_response_used": False,
        "response_firewall": {
            "Study_C_sequence_identity_opened": False,
            "Study_C_pairwise_genetic_distances_opened": False,
            "Study_C_T_st_computed": False,
            "Study_C_beta_hist_computed": False,
        },
    }
    binding.write_text(json.dumps(binding_payload))

    checked = module.verify_final_artifact(binding, artifact)
    assert checked["workflow_run_id"] == 35941577015

    extracted = tmp_path / "extracted"
    extracted.mkdir()
    with zipfile.ZipFile(artifact) as z:
        z.extractall(extracted)
    found_occ, found_ledger = module.locate_exact_files(extracted)
    assert sha(found_occ) == binding_payload["occurrence_csv_sha256"]
    assert sha(found_ledger) == binding_payload["occurrence_ledger_sha256"]

    binding_payload["artifact_digest"] = "sha256:" + "0" * 64
    binding.write_text(json.dumps(binding_payload))
    import pytest
    with pytest.raises(RuntimeError, match="artifact ZIP SHA-256 drift"):
        module.verify_final_artifact(binding, artifact)

    binding_payload["artifact_digest"] = "sha256:" + sha(artifact)
    binding_payload["workflow_run_id"] = 35941577016
    binding.write_text(json.dumps(binding_payload))
    with pytest.raises(RuntimeError, match="workflow run drift"):
        module.verify_final_artifact(binding, artifact)

    binding_payload["workflow_run_id"] = 35941577015
    binding_payload["occurrence_csv_sha256"] = "0" * 64
    binding.write_text(json.dumps(binding_payload))
    checked = module.verify_final_artifact(binding, artifact)
    assert checked["occurrence_csv_sha256"] == "0" * 64
    assert sha(found_occ) != checked["occurrence_csv_sha256"]


def test_local_study_c_readiness_receipt_is_response_blind_and_single_run_bound():
    import json

    p = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_local_readiness_v0.1.json"
        ).read_text()
    )

    assert p["status"] == "READY_AWAITING_FROZEN_OCCURRENCE_FINAL_ONLY"
    assert p["occurrence_gate"]["workflow_run_id"] == 35941577015
    assert p["occurrence_gate"]["workflow_head_sha"] == "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab"
    assert p["occurrence_gate"]["alternative_run_allowed"] is False
    assert p["occurrence_gate"]["final_artifact_seen_at_freeze"] is False
    assert p["occurrence_gate"]["binding_seen_at_freeze"] is False

    assert p["exact_phylogatr_source"]["sha256"] == "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
    assert p["exact_phylogatr_source"]["size_bytes"] == 274988692
    assert p["exact_phylogatr_source"]["source_repacking_allowed"] is False

    assert p["chelsa"]["staging_manifest_sha256"] == "7d82833996006fc1e99b86acb4ea7e80cde29c7f500df739ef48cb9ce145d450"
    assert p["chelsa"]["historical_asset_count"] == 8
    assert p["chelsa"]["current_asset_count"] == 4
    assert p["chelsa"]["alternate_asset_selection_allowed"] is False

    assert p["local_execution_bundle"]["source_bundle_commit_sha"] == "8c5fd223a518e19453c8273069fd8b4d72a1882c"
    assert p["local_execution_bundle"]["executor_git_blob_sha1"] == "0f91a4ce29c9d0a1454673da2e2c979365eea86e"
    assert p["local_execution_bundle"]["final_artifact_handoff_git_blob_sha1"] == "487de64b3e4abe411344edeb992a6294252c191c"
    assert p["local_execution_bundle"]["implementation_bound_files"] == 31

    assert p["validation"]["closure_status"] == "PASS_FROZEN_V03_CONTRACT_CLOSURE"
    assert p["validation"]["focused_tests_passed"] == 41
    assert p["relation_result_seen"] is False
    assert p["genetic_response_used"] is False
    assert all(v is False for v in p["response_firewall"].values())


def test_default_branch_binder_packages_exact_local_final_handoff():
    import json

    binder = (
        ROOT / ".github/workflows/relational-historical-bind-bounded-occurrence.yml"
    ).read_text()
    rule = json.loads(
        (
            ROOT
            / "docs/supporting/relational_historical_local_final_handoff_transport_v0.1.json"
        ).read_text()
    )
    receipt = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_downstream_activation_repair_v0.1.json"
        ).read_text()
    )

    assert "outputs:" in binder
    assert "artifact_id: ${{ steps.frozen.outputs.artifact_id }}" in binder
    assert "artifact_digest: ${{ steps.frozen.outputs.artifact_digest }}" in binder
    assert "package_local_handoff:" in binder
    assert "needs: bind" in binder
    assert 'gh api "/repos/$GITHUB_REPOSITORY/actions/artifacts/$ARTIFACT_ID/zip"' in binder
    assert 'binding["artifact_digest"]==os.environ["EXPECTED_DIGEST"]==digest' in binder
    assert "relational-historical-occurrence-final-v0.2.zip" in binder
    assert "relational-c-local-final-handoff-v0.1" in binder
    assert "compression-level: 0" in binder

    assert rule["status"] == "FROZEN_TRANSPORT_ONLY_BEFORE_FINAL_OCCURRENCE_RESULT"
    assert rule["source_run"]["workflow_run_id"] == 35941577015
    assert rule["source_run"]["workflow_head_sha"] == "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab"
    assert rule["result_selection_allowed"] is False
    assert rule["alternate_run_allowed"] is False
    assert rule["repack_inner_artifact_allowed"] is False
    assert rule["relation_result_seen"] is False
    assert rule["genetic_response_used"] is False
    assert all(v is False for v in rule["response_firewall"].values())

    assert receipt["default_branch_binder_commit_sha"] == "0b5bded97d7e02db86f59fbd89d3114b19a651b2"
    assert receipt["default_branch_binder_file_git_blob"] == "9fd8a23d2d24463f5758eb0847c5d7e2e33ed153"
    assert receipt["local_final_handoff_output_artifact"] == "relational-c-local-final-handoff-v0.1"
    assert receipt["local_final_handoff_inner_artifact_repacked"] is False
    assert receipt["local_final_handoff_result_selection_allowed"] is False


def test_default_branch_binder_is_valid_yaml_and_has_required_jobs():
    import yaml

    path = ROOT / ".github/workflows/relational-historical-bind-bounded-occurrence.yml"
    payload = yaml.safe_load(path.read_text())

    assert isinstance(payload, dict)
    assert payload["name"] == "relational-historical-bind-bounded-occurrence"
    jobs = payload["jobs"]
    assert set(["bind", "package_local_handoff", "execute_c"]).issubset(jobs)
    assert jobs["package_local_handoff"]["needs"] == "bind"
    assert jobs["execute_c"]["needs"] == "bind"
    assert jobs["execute_c"]["uses"].endswith(
        "@17923b1719d9a3ebce842c959e4547378ac28338"
    )
    assert jobs["bind"]["outputs"]["artifact_id"] == "${{ steps.frozen.outputs.artifact_id }}"
    assert jobs["bind"]["outputs"]["artifact_digest"] == "${{ steps.frozen.outputs.artifact_digest }}"


def test_retry_import_recovery_is_exact_and_prequery():
    import json

    receipt = json.loads(
        (
            ROOT
            / "benchmarks/frozen/relational_historical_retry_import_recovery_v0.1.json"
        ).read_text()
    )
    workflow = (
        ROOT / ".github/workflows/relational-historical-retry-import-recovery-v01.yml"
    ).read_text()

    assert receipt["status"] == "FROZEN_TECHNICAL_RECOVERY_BEFORE_ANY_RETRY_QUERY"
    assert receipt["original_execution"]["workflow_run_id"] == 35941577015
    assert receipt["original_execution"]["workflow_head_sha"] == "ffacbd51d58689a4b18f7a2cb920f5a1385a74ab"
    assert receipt["original_execution"]["request_error_count"] == 95
    assert receipt["original_execution"]["retry_batch_count"] == 24
    assert receipt["original_execution"]["retry_jobs_failure"] == 24
    assert receipt["original_execution"]["retry_output_artifacts"] == 0
    assert receipt["failure_class"]["bounded_fetch_reached"] is False
    assert receipt["failure_class"]["gbif_retry_query_consumed"] is False
    assert receipt["recovery_execution"]["exact_retry_plan_artifact_id"] == 10878157732
    assert receipt["recovery_execution"]["request_error_species"] == 95
    assert receipt["recovery_execution"]["retry_batches"] == 24
    assert receipt["recovery_execution"]["batch_size"] == 4
    assert receipt["recovery_execution"]["max_parallel"] == 2
    assert receipt["recovery_execution"]["per_species_wall_timeout_seconds"] == 2700
    assert receipt["recovery_execution"]["retry_round"] == 1
    assert receipt["recovery_execution"]["retry_script_changed"] is False
    assert receipt["recovery_execution"]["scientific_query_change"] is False
    assert receipt["recovery_execution"]["additional_retry_round_authorized"] is False
    assert receipt["result_selection_allowed"] is False
    assert receipt["relation_result_seen"] is False
    assert receipt["genetic_response_used"] is False
    assert all(v is False for v in receipt["response_firewall"].values())

    assert "ref: ffacbd51d58689a4b18f7a2cb920f5a1385a74ab" in workflow
    assert 'PYTHONPATH="$PWD" python scripts/retry_relational_historical_occurrences_bounded.py' in workflow
    assert "--batch-size 4" in workflow
    assert "--per-species-timeout-seconds 2700" in workflow
    assert "max-parallel: 2" in workflow
    assert "10878157732" in workflow
    assert "b4291a54514837f944f1cb9dbf14d320fd52293d24aa33b26d36ac063b7c7103" in workflow
    assert "scripts/finalize_relational_historical_occurrence_bounded.py" in workflow
    assert "relational-historical-occurrence-final-v0.2" in workflow
    assert "relational-c-local-final-handoff-v0.1" in workflow
