from __future__ import annotations

import json
from pathlib import Path

RESULT = Path("benchmarks/frozen/genetic_conditional_geometry_matched_development_screen_v0.3.json")

def test_v03_screen_stops_taxonomic_order_route_without_opening_identity() -> None:
    r = json.loads(RESULT.read_text())
    assert r["schema"] == "ttf_genetic_conditional_geometry_matched_development_screen_v0.3"
    assert r["status"] == "STOP_TAXONOMIC_ORDER_ROUTE_V03_SURVIVOR_POWER_FAIL"
    dev = r["development_geometry"]
    surv = r["confirmatory_survivor_geometry"]
    assert dev["eligible_eval_species"] == 52
    assert surv["eligible_eval_species"] == 63
    assert dev["private_screen_pass"] is True
    assert dev["positive_screen_pass"] is True
    assert surv["private_screen_pass"] is True
    assert surv["positive_screen_pass"] is False
    assert surv["cells"]["same_order_A2"]["rejections"] == 19
    assert surv["positive_requirement_rejections"] == 35
    decision = r["decision"]
    assert decision["type1_problem_from_v02_resolved_in_screen"] is True
    assert decision["survivor_positive_control_power_adequate"] is False
    assert decision["formal_v03_qualification_authorized"] is False
    assert decision["confirmatory_identity_opening_authorized"] is False
    assert decision["taxonomic_order_route_status"] == "STOP_AFTER_V03_SCREEN"
    assert all(value is False for value in r["outcome_firewall"].values())
