from __future__ import annotations

import json
from pathlib import Path


RESULT = Path("benchmarks/frozen/lepidoptera_within_order_descriptive_baseline_v0.1.json")


def test_lepidoptera_descriptive_baseline_is_noninferential_and_does_not_rewrite_primary() -> None:
    payload = json.loads(RESULT.read_text())
    assert payload["schema"] == "ttf_lepidoptera_within_order_descriptive_baseline_v0.1"
    assert payload["status"] == "DESCRIPTIVE_ONLY_AFTER_PRIMARY_DECISION"
    assert payload["n_eval_species"] == 110
    assert payload["statistic"] == 0.051674313053575054
    assert payload["inferential_claim_authorized"] is False
    assert payload["primary_decision_unchanged"]["decision"] == "NO_DETECTED_POSITIVE_TRAIT_SIMILARITY_GRADIENT"
    assert payload["primary_decision_unchanged"]["primary_envelope_p_value"] == 0.89
    assert payload["serialized_sequence_identity"] is False
    assert payload["serialized_edge_genetic_distance_vectors"] is False
    assert payload["post_result_retuning_allowed"] is False
