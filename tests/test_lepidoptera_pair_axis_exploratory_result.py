from __future__ import annotations

import json
from pathlib import Path

RESULT=Path("benchmarks/frozen/lepidoptera_pair_axis_exploratory_v0.1.json")

def test_pair_axis_exploratory_is_post_primary_and_nonconfirmatory():
    payload=json.loads(RESULT.read_text())
    assert payload["schema"]=="ttf_lepidoptera_pair_axis_exploratory_v0.1"
    assert payload["status"]=="POST_PRIMARY_EXPLORATORY_ONLY"
    assert payload["parent_primary_decision"]=="NO_DETECTED_POSITIVE_TRAIT_SIMILARITY_GRADIENT"
    assert payload["pairs"]==2428
    assert payload["source_species"]==106
    assert payload["target_species"]==110
    for name,row in payload["models"].items():
        assert row["rmse_improvement"]<=0.0, name
    assert payload["claim_boundary"]["confirmatory"] is False
    assert payload["claim_boundary"]["p_values_computed"] is False
    assert payload["claim_boundary"]["primary_result_rewritten"] is False
    assert payload["serialized_sequence_identity"] is False
    assert payload["serialized_edge_genetic_distance_vectors"] is False
    assert payload["post_result_retuning_allowed"] is False
