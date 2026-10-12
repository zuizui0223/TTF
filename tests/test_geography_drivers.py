"""Synthetic geography only; no real nucleotide identity or genetic distances."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import pytest
from audit_locality_geography_drivers import analyze_geographic_influence

def toy():
    species = [{"species":"ToyOne","localities":"5","edges":"6",
                "unscaled_slope":"0.4","locality_deletions_attempted":"5","locality_deletions_reversing_sign":"1"}]
    xyz = [(0,0,0),(1,0,0),(2,0,0),(3,0,0),(40,0,0)]
    loc=[dict(species="ToyOne",locality_index=str(i),
              x_km=str(x),y_km=str(y),z_km=str(z))
         for i,(x,y,z) in enumerate(xyz)]
    pairs=[(0,1),(1,2),(2,3),(3,4),(0,2),(1,3)]
    edges=[dict(species="ToyOne",edge_index=str(i),node_left=str(a),node_right=str(b))
           for i,(a,b) in enumerate(pairs)]
    dels=[dict(species="ToyOne",locality_index=str(j),status="defined",
               slope="-0.2" if j==4 else "0.4") for j in range(5)]
    return species,dels,loc,edges


def test_isolated_flip_site_ranks_farther_than_other_sites():
    report=analyze_geographic_influence(*toy())
    assert report["original_species"]==1
    assert report["original_sign_flipping_species"]==1
    assert report["site_deletions_reversing_sign"]==1
    nn=report["within_species_geographic_comparisons"]["nearest_neighbor_chord_km"]
    assert nn["species_with_higher_flip_site_rank"]==1
    assert nn["flipping_site_mean_rank"]>nn["other_site_mean_rank"]
    assert "species" not in report


def test_out_of_frozen_geometry_fails_closed():
    s,d,l,e=toy()
    e[0]["node_left"]="10"
    with pytest.raises(ValueError,match="nonexistent frozen locality"):
        analyze_geographic_influence(s,d,l,e)


def test_zero_baseline_is_not_a_directional_flip_or_invalid_species():
    s,d,l,e=toy()
    s[0]["unscaled_slope"]="0"
    s[0]["locality_deletions_reversing_sign"]="0"
    for row in d: row["slope"]="0"
    report=analyze_geographic_influence(s,d,l,e)
    assert report["original_species"]==1
    assert report["original_zero_species"]==1
    assert report["original_sign_flipping_species"]==0
    assert report["site_deletions_reversing_sign"]==0


def test_missing_deletion_record_fails_closed():
    s,d,l,e=toy()
    with pytest.raises(ValueError,match="frozen locality"):
        analyze_geographic_influence(s,d[:-1],l,e)



def test_average_tie_rank_matches_original_byte_parsing_logic():
    from audit_locality_geography_drivers import rank_avg
    import numpy as np
    x=rank_avg([1.,1.,2.])
    np.testing.assert_allclose(x,[0.5,0.5,1.0])
