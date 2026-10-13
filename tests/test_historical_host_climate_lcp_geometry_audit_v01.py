"""Synthetic arithmetic checks for the post-outcome CHELSA LCP geography audit."""
import math
import pytest

from scripts.audit_historical_host_climate_lcp_geometry import audit_rows


def artificial_edges():
    distances=(1.,2.,4.,8.)
    h=(0.,2.,4.,8.)
    s=(0.,4.,8.,16.)
    rows=[]
    for i,(d,dh,ds) in enumerate(zip(distances,h,s)):
        rows.append(dict(species="Artificial species",edge_index=str(i),
            node_left=str(i), node_right=str(i+1),geographic_chord_km=str(d),
            R_host_0BP_km=str(1.*d if dh>0 else 0),
            R_host_LGM_km=str(1.*d+dh if dh>0 else 0),
            R_self_0BP_km=str(2.*d if ds>0 else 0),
            R_self_LGM_km=str(2.*d+ds if ds>0 else 0),
            M_host_lcp_km=str(dh),M_self_lcp_km=str(ds),
            M_host_minus_self_lcp_km=str(dh-ds)))
    return rows


def test_geometry_normalization_does_not_invent_genetic_support():
    obj=audit_rows(artificial_edges(), expected_species="Artificial species", expected_edges=4)
    assert obj["status"] == "EXPLORATORY_CLIMATE_ONLY_GEO_NORMALIZATION"
    assert obj["counts"]["host_LGM_cost_higher"] == 3
    assert obj["counts"]["host_minus_self_negative"] == 3
    assert obj["counts"]["same_cell_zero_path_edges"] == 1
    assert obj["pairwise_km_normalized"]["host_mean_delta_per_original_chord_km"] == pytest.approx(.75)
    assert obj["pairwise_km_normalized"]["self_mean_delta_per_original_chord_km"] == pytest.approx(1.5)
    assert obj["diagnostics"]["no_genetic_effect_or_p_value"] is True
    assert obj["diagnostics"]["no_paleo_land_ice_mask"] is True
    assert math.isfinite(obj["geometry_association"]["pearson_original_chord_vs_self_delta"])


def test_duplicate_undirected_edges_and_bad_chord_fail_closed():
    rows=artificial_edges()
    rows[-1]["node_left"]="1";rows[-1]["node_right"]="0"
    with pytest.raises(ValueError,match="duplicate"):
        audit_rows(rows,expected_species="Artificial species",expected_edges=4)
    rows=artificial_edges()
    rows[1]["geographic_chord_km"]="0"
    with pytest.raises(ValueError,match="positive"):
        audit_rows(rows,expected_species="Artificial species",expected_edges=4)


def test_arithmetic_and_species_source_are_not_optional():
    rows=artificial_edges()
    rows[0]["species"]="Another species"
    with pytest.raises(ValueError,match="source species"):
        audit_rows(rows,expected_species="Artificial species",expected_edges=4)
    rows=artificial_edges()
    rows[0]["M_host_minus_self_lcp_km"]="99"
    with pytest.raises(ValueError,match="delta"):
        audit_rows(rows,expected_species="Artificial species",expected_edges=4)
