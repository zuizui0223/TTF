from __future__ import annotations

import json
from pathlib import Path


RESULT=Path("benchmarks/frozen/lepidoptera_host_resource_empirical_primary_result_v0.1.json")


def test_host_resource_one_shot_result_is_terminal_and_not_retunable():
    payload=json.loads(RESULT.read_text())
    assert payload["schema"]=="ttf_lepidoptera_host_resource_empirical_primary_result_v0.1"
    assert payload["status"]=="PRIMARY_ONE_SHOT_DECISION_COMPLETE"
    primary=payload["primary_host_resource_gradient"]
    assert primary["statistic"]==-0.00704481173406459
    assert primary["envelope_p_value"]==0.698
    assert primary["positive"] is False
    assert primary["finite_target_correlations"]==201
    assert primary["target_correlation_summary"]["positive_targets"]==101
    assert primary["target_correlation_summary"]["negative_targets"]==100
    assert payload["decision"]=="NO_DETECTED_POSITIVE_HOST_RESOURCE_GEOGRAPHY_GRADIENT"
    assert payload["qualification_context"]["host_resource_positive_power_wilson95_lower"]>0.80
    assert payload["outcome_state"]["sequence_identity_opened"] is True
    assert payload["outcome_state"]["serialized_sequence_identity"] is False
    assert payload["outcome_state"]["serialized_edge_genetic_distance_vectors"] is False
    assert payload["post_result_retuning_allowed"] is False
