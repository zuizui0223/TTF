import ast
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_occurrence_aggregator_has_closed_genetic_response_firewall():
    source=(ROOT/"scripts/aggregate_palearctic_holometabola_occurrences.py").read_text()
    ast.parse(source)
    assert '"pairwise_genetic_concordance_opened":False' in source
    assert '"lgm_relation_computed":False' in source
    assert "species_scores" not in source
