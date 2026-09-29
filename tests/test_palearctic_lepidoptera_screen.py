import runpy
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
NS=runpy.run_path(str(ROOT/"scripts/screen_palearctic_lepidoptera_subpanel.py"))


def test_core_fraction_is_deterministic_and_boundary_inclusive():
    rule={
        "latitude_min":30.0,"latitude_max":75.0,
        "longitude_min":-25.0,"longitude_max":180.0,
    }
    coords=[(30.0,-25.0),(75.0,180.0),(29.9,10.0),(50.0,0.0)]
    assert NS["core_fraction"](coords,rule)==0.75


def test_locus_normalization_matches_frozen_coi_aliases():
    f=NS["normalized_locus"]
    assert f("Acleris-bergmanniana-COI","Acleris bergmanniana")=="COI"
    assert f("X-Y-MT-CO1","X Y")=="MT CO1"


def test_screen_source_never_imports_genetic_response_modules():
    source=(ROOT/"scripts/screen_palearctic_lepidoptera_subpanel.py").read_text()
    forbidden=("genetic_empirical","phase4_empirical","species_scores")
    assert not any(token in source for token in forbidden)
