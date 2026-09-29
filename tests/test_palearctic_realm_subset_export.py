import ast
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_palearctic_subset_export_has_no_genetic_inputs():
    source=(ROOT/"scripts/export_palearctic_realm_subset.py").read_text()
    ast.parse(source)
    assert '"REALM"' in source
    assert '"Palearctic"' in source
    assert "species_scores" not in source
    assert "genetic_phase4" not in source
    assert '"panel_membership_opened":False' in source
