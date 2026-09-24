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
    assert receipt["default_branch_binder_commit_sha"] == "46c2c9fc9a7f2933dadf629fca30eb559a57d43e"
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
    assert rule["final_artifact_alternative_run_allowed"] is False
    assert rule["final_artifact_repacking_allowed"] is False
    assert "artifact_digest" in rule["final_artifact_handoff_contract"]
    assert "run_relational_historical_local_from_artifact.py" in source_bundle
    compile(wrapper, "scripts/run_relational_historical_local_from_artifact.py", "exec")
