from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_postgate_transport_completion_contract_is_uniform_and_nonselective():
    contract = json.loads(
        (
            ROOT
            / "docs/exploratory/butterfly_climate_release_postgate_transport_completion_v0.1.json"
        ).read_text()
    )
    assert contract["status"] == (
        "FROZEN_POSTGATE_TECHNICAL_COMPLETION_AFTER_NOT_EVALUABLE_BEFORE_ANY_CLIMATE_RESULT"
    )
    assert contract["triggering_execution"]["run_id"] == 36241187187
    assert contract["triggering_execution"]["qualified_species"] == 10
    assert contract["triggering_execution"]["minimum_required"] == 12
    assert contract["triggering_execution"]["climate_job_conclusion"] == "skipped"

    census = contract["transport_census"]
    assert census == {
        "panel_species": 32,
        "complete": 17,
        "partial_transport": 14,
        "rejected_gbif_taxon_match": 1,
        "partial_missing_300_record_windows_total": 42,
    }

    population = contract["recovery_population"]
    assert population["species_count"] == 14
    assert len(population["species"]) == 14
    assert len(set(population["species"])) == 14
    assert population["complete_species_refetched"] is False
    assert population["rejected_taxon_species_replaced_or_recovered"] is False
    assert population["species_replacement_authorized"] is False

    technical = contract["technical_completion_rule"]
    assert technical["completed_windows_immutable"] is True
    assert technical["only_missing_windows_may_be_queried"] is True
    assert technical["gbif_filters_change"] is False
    assert technical["period_change"] is False
    assert technical["deterministic_window_offsets_change"] is False
    assert technical["maximum_windows_per_species_change"] is False
    assert technical["ordinal_window_size_records"] == 300
    assert technical["recovery_transport_chunk_size_records"] == 300
    assert technical["request_timeout_cap_seconds"] == 600
    assert technical["species_total_deadline_seconds_per_pass"] == 3600
    assert technical["maximum_passes_per_job"] == 2
    assert technical["max_parallel_species"] == 1

    science = contract["unchanged_scientific_contract"]
    assert science["independent_panel_species_sha256"] == (
        "bbaaa28cf361fb02d558988bdd1a3f948a243b79c8d3bb1709e908edb1f0f204"
    )
    assert science["protocol_schema"] == (
        "ttf_butterfly_climate_release_independent_test_v0.2.1"
    )
    assert science["quality_threshold_change"] is False
    assert science["minimum_preclimate_qualified_species"] == 12
    assert science["climate_crossfit_change"] is False
    assert science["minimum_climate_informative_species"] == 12
    assert science["primary_effect_score_change"] is False
    assert science["primary_predictor_change"] is False
    assert science["primary_control_change"] is False
    assert science["primary_inference_change"] is False
    assert science["primary_alternative"] == "negative"
    assert science["alpha"] == 0.05
    assert science["permutation_iterations"] == 9999


def test_postgate_recovery_workflow_uses_exact_authoritative_state_and_gate():
    workflow = (
        ROOT
        / ".github/workflows/butterfly-climate-release-postgate-recovery-v01.yml"
    ).read_text()

    assert "run-id: 36241187187" in workflow
    assert (
        "key: butterfly-climate-release-gbif-v01-${{ matrix.id }}-36241187187"
        in workflow
    )
    assert "max-parallel: 1" in workflow
    assert "--transport-chunk-size 300" in workflow
    assert "--request-seconds 600" in workflow
    assert "--species-seconds 3600" in workflow
    assert "--maximum-pages 6" in workflow
    assert (
        "--protocol-json docs/exploratory/"
        "butterfly_climate_release_independent_test_v0.2.1.json"
        in workflow
    )
    assert (
        "butterfly_climate_release_postgate_transport_completion_v0.1.json"
        in workflow
    )
    assert "completed_window_manifest_before.json" in workflow
    assert "completed_window_manifest_after.json" in workflow
