from __future__ import annotations

import json
from pathlib import Path


RULE=Path("docs/supporting/lepidoptera_postprimary_self_detectability_rule_v0.1.json")


def test_postprimary_self_anchor_cannot_rewrite_cross_species_primaries():
    rule=json.loads(RULE.read_text())
    assert rule["schema"]=="ttf_lepidoptera_postprimary_self_detectability_rule_v0.1"
    assert rule["status"]=="FROZEN_BEFORE_INSECT_PANEL_SELF_STATISTICS"
    assert rule["method_inheritance"]["changes_to_self_estimator"]=="none"
    assert rule["panels"]["butterfly_trait_panel"]["cross_species_primary_known"] is True
    assert rule["panels"]["host_resource_panel"]["cross_species_primary_known"] is True
    assert rule["panels"]["butterfly_trait_panel"]["self_statistic_known"] is False
    assert rule["panels"]["host_resource_panel"]["self_statistic_known"] is False
    assert rule["empirical_rule"]["no_species_dropping"] is True
    assert rule["empirical_rule"]["no_parameter_retuning"] is True
