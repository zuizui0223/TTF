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
