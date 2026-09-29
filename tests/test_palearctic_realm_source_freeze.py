import ast
import json
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_palearctic_source_freezer_is_response_independent():
    source=(ROOT/"scripts/freeze_palearctic_realm_source.py").read_text()
    tree=ast.parse(source)
    assert "genetic_pair_response_opened" in source
    assert '"REALM"' in source
    assert '"Palearctic"' in source
    assert "species_scores" not in source
    assert "phase4_empirical_result" not in source


def test_palearctic_contract_binds_ecoregions2017_before_response():
    contract=json.loads(
        (ROOT/"docs/supporting/genetic_palearctic_holometabola_lgm_program_v0.1.json").read_text()
    )
    realm=contract["panel_rule"]["realm_source"]
    assert realm["dataset"]=="Ecoregions2017"
    assert realm["realm_field"]=="REALM"
    assert realm["realm_value"]=="Palearctic"
    assert contract["ttf_q_qualification"]["must_precede_pairwise_genetic_response_construction"] is True
