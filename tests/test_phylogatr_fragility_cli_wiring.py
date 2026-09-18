from pathlib import Path
import json
import pytest

SCRIPTS = (
    Path("scripts/run_phylogatr_gate_d_fragility_reference_shard.py"),
    Path("scripts/run_phylogatr_gate_d_fragility_observed_shard.py"),
)

@pytest.mark.parametrize("script", SCRIPTS)
def test_fragility_cli_uses_exact_compact_execution(script: Path) -> None:
    source = script.read_text()
    assert "prepare_phylogatr_compact_ttf_design" in source
    assert "prepare_phylogatr_compact_cached_transfer" in source
    assert "score_phylogatr_compact_world_batch" in source
    assert "prepare_genetic_ttf_design(" not in source
    assert "prepare_genetic_cached_transfer(" not in source
    assert "score_genetic_world_batch(" not in source

def test_fragility_compact_execution_rule_is_frozen_before_fragility_results() -> None:
    payload = json.loads(
        Path("docs/supporting/genetic_phylogatr_fragility_compact_execution_v0.1.json").read_text()
    )
    assert payload["schema"] == "ttf_genetic_phylogatr_fragility_compact_execution_v0.1"
    assert payload["status"] == (
        "FROZEN_AFTER_GATE_D_AND_SELF_PASS_BEFORE_FRAGILITY_SYNTHETIC_RESULT_OR_IDENTITY_OPENING"
    )
    assert payload["trigger"]["formal_gate_d_status"] == "PASS"
    assert payload["trigger"]["self_detectability_status"] == "PASS"
    assert payload["trigger"]["fragility_synthetic_result_seen_before_amendment"] is False
    assert all(value is False for value in payload["firewall"].values())
