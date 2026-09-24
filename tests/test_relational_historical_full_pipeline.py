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


def test_bounded_handoff_delays_archive_until_relation_passes():
    text = (ROOT / ".github/workflows/relational-historical-from-bounded-binding-v02.yml").read_text()
    hist_assets = text.index("Resolve and freeze exact CHELSA-TraCE21k assets")
    current_assets = text.index("Download frozen current CHELSA V2.1 nuisance layers")
    relation = text.index("Build frozen historical/current relation")
    archive = text.index("Restore exact frozen phylogatR archive")
    geometry = text.index("Reconstruct frozen Study-C response-blind geometry")
    opportunity = text.index("Attach frozen historical opportunity geometry")

    assert hist_assets < current_assets < relation < archive < geometry < opportunity
    assert "relation_construction_before_geometry_materialization_allowed" not in text
    # Execution ordering, rather than a changed scientific rule, ensures the
    # response-blind relation can be evaluated before exact genetic geometry is needed.


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
        "d9e2c51793312b2d23504908ad0470549939ba7a"
    ) in binder
    assert receipt["status"] == "FROZEN_RESPONSE_BLIND_ORCHESTRATION_REPAIR_BEFORE_STUDY_C_RELATION_RESULT"
    assert receipt["scientific_rule_change"] is False
    assert receipt["predictor_change"] is False
    assert receipt["alpha_change"] is False
    assert receipt["called_workflow_commit_sha"] == "d9e2c51793312b2d23504908ad0470549939ba7a"
    assert receipt["default_branch_binder_commit_sha"] == "01dbcc466c46a7655ee8c8fda6735783f9116f34"
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
