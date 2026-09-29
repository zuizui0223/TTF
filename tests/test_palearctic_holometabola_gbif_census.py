import ast
import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_palearctic_gbif_census_is_response_blind():
    for path in (
        ROOT/"scripts/census_palearctic_holometabola_gbif.py",
        ROOT/"scripts/aggregate_palearctic_holometabola_gbif.py",
    ):
        source=path.read_text()
        ast.parse(source)
        assert "species_scores" not in source
        assert "pairwise_genetic_concordance_opened" in source
    contract=json.loads(
        (ROOT/"docs/supporting/genetic_palearctic_holometabola_lgm_program_v0.1.json").read_text()
    )
    assert contract["lgm_sdm"]["minimum_thinned_occurrences_for_sdm"]==50
    assert contract["lgm_sdm"]["post_panel_failure_policy"]["resplit_after_occurrence_or_sdm_failure"] is False
    assert contract["lgm_sdm"]["minimum_species_with_usable_sdm"]==30
