from __future__ import annotations

import json
from pathlib import Path


RESULT = Path(
    "benchmarks/frozen/genetic_phylogatr_phase4_empirical_handoff_v0.1.json"
)


def test_phase4_empirical_result_is_frozen_without_post_result_retuning() -> None:
    result = json.loads(RESULT.read_text())

    assert result["schema"] == "ttf_genetic_phylogatr_phase4_empirical_handoff_v0.1"
    assert result["status"] == "EMPIRICAL_PHASE4_COMPLETE"
    assert result["source_main_sha"] == "55ff72a4273d2119a361c3978aa3bea873ed8243"

    source = result["exact_source_archive"]
    assert source["size_bytes"] == 274988692
    assert source["sha256"] == "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"

    phase1 = result["phase1_regeneration"]
    assert phase1["status"] == "EXACT_PHASE1_MANIFEST_REGENERATED_RESPONSE_BLIND"
    assert phase1["phase1_manifest_sha256"] == "daeedabbeb99fa304568f0f2b2e1c46e0dd6221fbac02b0491faf2435ef7b7b3"
    assert phase1["phase1_geometry_csv_sha256"] == "d6f0b1f360a77de279263145386e7956af0d43c9386048f319117522f933cb74"
    assert phase1["selected_panels_count"] == 250
    assert phase1["sequence_identity_opened_during_regeneration"] is False

    authorization = result["phase4_authorization"]
    assert authorization["status"] == "AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING"
    assert authorization["phase1_provenance_mode"] == "exact_manifest"
    assert authorization["authorization_json_sha256"] == "a66edba0a6378f2fa2e24eba73aa5d48f6a58905420e13b53be802f603c3a0e1"
    assert authorization["survivor_species"] == 211

    empirical = result["empirical_result"]
    assert empirical["result_json_sha256"] == "f5d19fa50c7c18cbf2a110c0afb9c70dcd543c44b77521c015938013843e72fb"
    primary = empirical["primary_place_beyond_ibd"]
    assert primary["statistic"] == 0.03256548471362711
    assert primary["profiled_private_p_value"] == 0.7392607392607392
    assert primary["positive"] is False
    assert primary["selected_configurations"] == ["A3", "A5"]

    self_diag = empirical["within_species_self_diagnostic"]
    assert self_diag["synthetic_method_qualified"] is True
    assert self_diag["p_value"] == 0.008991008991008992
    assert self_diag["positive"] is True

    secondary = empirical["secondary_total_genetic_transfer"]
    assert secondary["status"] == "DESCRIPTIVE_ONLY"
    assert secondary["inferentially_qualified"] is False
    assert secondary["used_for_primary_decision"] is False

    assert empirical["decision"] == "LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN"

    opened = result["outcome_state"]
    assert opened["confirmatory_sequence_identity_opened"] is True
    assert opened["confirmatory_pairwise_genetic_distances_opened"] is True
    assert opened["confirmatory_ttf_statistic_opened"] is True
    assert opened["decker_empirical_genetic_outcomes_opened"] is False
    assert opened["serialized_sequence_identity"] is False
    assert opened["serialized_edge_genetic_distance_vectors"] is False

    one_shot = result["one_shot_execution"]
    assert one_shot["frozen_empirical_test_completed"] is True
    assert one_shot["result_selection_rerun_allowed"] is False
    assert one_shot["post_result_retuning_allowed"] is False
