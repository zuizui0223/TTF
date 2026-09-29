import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCREEN = runpy.run_path(
    str(ROOT / "scripts/screen_palearctic_lepidoptera_response_blind.py")
)


def test_polygon_boundary_counts_as_inside_and_hole_counts_as_outside():
    polygon = {
        "type": "Polygon",
        "coordinates": [[
            [0,0],[10,0],[10,10],[0,10],[0,0]
        ],[
            [4,4],[6,4],[6,6],[4,6],[4,4]
        ]],
    }
    contains = SCREEN["geometry_contains"]
    assert contains(1,1,polygon)
    assert contains(0,5,polygon)
    assert not contains(5,5,polygon)
    assert not contains(11,5,polygon)


def test_split_is_deterministic_species_only():
    split = SCREEN["_split"]
    species = ["Species c", "Species a", "Species b", "Species d"]
    first = split(species, "abc")
    second = split(list(reversed(species)), "abc")
    assert first == second
    assert set(first[0]) | set(first[1]) == set(species)
    assert set(first[0]).isdisjoint(first[1])


def test_frozen_palearctic_lgm_contracts_keep_genetic_response_closed():
    for path in (
        ROOT / "docs/supporting/genetic_palearctic_lepidoptera_eligibility_v0.1.json",
        ROOT / "docs/supporting/genetic_palearctic_lgm_refugia_predictor_v0.1.json",
        ROOT / "docs/supporting/genetic_palearctic_lgm_opening_rule_v0.1.json",
    ):
        payload = json.loads(path.read_text())
        text = json.dumps(payload).lower()
        assert '"species_level_phase4_scores_used": false' in text or '"species_level_genetic_scores_used": false' in text
