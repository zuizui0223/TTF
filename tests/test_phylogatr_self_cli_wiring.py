from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "script",
    [
        Path("scripts/run_phylogatr_phase3_self_reference_shard.py"),
        Path("scripts/run_phylogatr_phase3_self_evaluation_shard.py"),
    ],
)
def test_fresh_self_cli_uses_compact_design_and_scorer(script: Path) -> None:
    source = script.read_text()

    assert "prepare_phylogatr_compact_ttf_design" in source
    assert "score_phylogatr_compact_self_world_batch" in source
    assert "prepare_genetic_ttf_design(" not in source
    assert "score_genetic_self_world_batch(" not in source
