from __future__ import annotations

import json
from pathlib import Path


RESULT=Path("benchmarks/frozen/lepidoptera_species_transferability_exploratory_summary_v0.1.json")


def test_species_transferability_result_is_exploratory_and_cannot_rewrite_primary():
    payload=json.loads(RESULT.read_text())
    assert payload["schema"]=="ttf_lepidoptera_species_transferability_exploratory_summary_v0.1"
    assert payload["status"]=="EXPLORATORY_SPECIES_PROPERTY_WEAK"
    assert payload["parent_primary_decision"]=="NO_DETECTED_POSITIVE_TRAIT_SIMILARITY_GRADIENT"
    assert payload["analysis"]["pairs"]==2428
    assert payload["analysis"]["source_species_with_supported_pairs"]==106
    assert payload["analysis"]["target_species"]==110
    variance=payload["variance_decomposition"]
    assert variance["source_fraction"]<0.01
    assert variance["target_fraction"]<0.02
    assert variance["residual_fraction"]>0.98
    assert payload["pair_cross_validation"]["rmse_improvement"]<=0
    assert payload["claim_boundary"]["confirmatory"] is False
    assert payload["claim_boundary"]["trait_predictors_tested"] is False
    assert payload["claim_boundary"]["primary_decision_rewritten"] is False
    assert payload["serialized_sequence_identity"] is False
    assert payload["serialized_edge_genetic_distance_vectors"] is False
