from __future__ import annotations

import json
from pathlib import Path


RESULT=Path("benchmarks/frozen/lepidoptera_pair_repeatability_v01_qualification_result.json")


def test_pair_repeatability_v01_is_retired_without_empirical_opening():
    p=json.loads(RESULT.read_text())
    assert p["schema"]=="ttf_lepidoptera_pair_repeatability_v01_qualification_result"
    assert p["status"]=="FAIL_NOT_EVALUABLE_POWER"
    assert p["gates"]["private_type1_pass"] is True
    assert p["gates"]["geometry_confounded_trap_type1_pass"] is True
    assert p["gates"]["host_breadth_confounded_trap_type1_pass"] is True
    assert p["gates"]["host_resource_gradient_trap_type1_pass"] is True
    assert p["gates"]["relational_positive_power_pass"] is False
    assert p["decision"]=="RETIRE_V0_1_KEEP_EMPIRICAL_SPLIT_HALF_STATISTIC_CLOSED"
    assert p["outcome_firewall"]["empirical_split_half_repeatability_computed"] is False
