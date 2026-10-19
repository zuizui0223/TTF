import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest

from scripts import analyze_historical_host_same_site_design as mod


def candidate(sp, hosts='1.0', family='F'):
    return {'species':sp,'family':family,'host_species':str(len(hosts.split(';'))),
            'host_families':'1','accepted_host_ids':hosts,'accepted_host_names':'host'}


def loc(sp, idx, lat, lon):
    return {'species':sp,'locality_index':str(idx),'latitude':str(lat),'longitude':str(lon)}


def test_numeric_host_identity_not_text_identity():
    assert mod.host_ids(' 1.0 ;2;2.00 ') == {'1','2'}


def test_fractional_or_invalid_host_identity_rejected():
    for value in ['1.5', 'NaN', '', 'Infinity', '0', '-1']:
        with pytest.raises(ValueError):
            mod.host_ids(value)


def test_common_coordinates_are_numeric_without_tolerance_matching():
    c=[candidate('A'),candidate('B','2')]
    l=[loc('A',0,'10.0','20.00'),loc('A',1,10.1,20.1),
       loc('B',0,'10.00','20.0'),loc('B',1,10.10001,20.1)]
    result=mod.census(c,l,(19,9,21,11),minimum_localities=2)
    assert result['pairs'][0]['shared_recorded_coordinates']==1
    assert result['pairs'][0]['host_jaccard']==0


def test_duplicate_recorded_coordinate_within_species_rejected():
    with pytest.raises(ValueError,match='duplicate coordinate'):
        mod.census([candidate('A')],[loc('A',0,10,20),loc('A',1,'10.0','20.0')],(19,9,21,11),1)


def test_unrecognized_species_rejected():
    with pytest.raises(ValueError,match='candidate'):
        mod.census([candidate('A')],[loc('B',0,10,20)],(19,9,21,11),1)


def test_boundary_inclusive_and_all_species_reported():
    c=[candidate('A'),candidate('B')]
    l=[loc('A',0,9,19),loc('A',1,11,21),loc('B',0,11.00001,21)]
    r=mod.census(c,l,(19,9,21,11),2)
    assert len(r['screen'])==2
    assert [x['species'] for x in r['panel']]==['A']
    assert r['screen'][1]['localities_inside']==0


def test_row_order_invariance():
    c=[candidate('A'),candidate('B','2')]
    l=[loc('A',0,10,20),loc('A',1,10.1,20.1),loc('B',0,10,20),loc('B',1,10.1,20.1)]
    assert mod.census(c,l,(19,9,21,11),2)==mod.census(c[::-1],l[::-1],(19,9,21,11),2)


def test_shared_site_queries_keep_both_species_original_locality_indices():
    c=[candidate('A'),candidate('B','1.00')]
    l=[loc('A',10,10,20),loc('A',20,10.1,20.1),loc('B',4,10.1,20.1),loc('B',8,10,20)]
    r=mod.census(c,l,(19,9,21,11),2)
    q=mod.shared_queries(r,minimum_shared=2)
    assert len(q)==1
    assert (q[0]['source_node_left'],q[0]['source_node_right'])==(10,20)
    assert (q[0]['target_node_left'],q[0]['target_node_right'])==(8,4)
    assert r['pairs'][0]['host_jaccard']==1.0


def edge(idx, g, h0,h1,s0,s1):
    return {'edge_index':idx,'geographic_chord_km':g,'R_host_0BP_km':h0,
            'R_host_LGM_km':h1,'R_self_0BP_km':s0,'R_self_LGM_km':s1}


def test_same_pixel_zero_is_separated_not_imputed():
    d=[edge(0,1,0,0,0,0),edge(1,2,1,3,1,5),edge(2,4,1,4,1,7),edge(3,8,1,5,1,9)]
    r=mod.shape_diagnostic(d)
    assert r['all_edges']==4 and r['positive_path_edges']==3
    assert r['same_pixel_edge_ids']==[0]
    assert r['spearman_delta_raw']==pytest.approx(1)
    assert r['spearman_delta_per_chord_km']==pytest.approx(1)


def test_amplitude_alone_does_not_create_shape_difference():
    d=[edge(i,i+1,1,1+h,1,1+3*h) for i,h in enumerate([2,3,7,4,6])]
    r=mod.shape_diagnostic(d)
    assert r['spearman_delta_raw']==pytest.approx(1)
    assert r['spearman_delta_per_chord_km']==pytest.approx(1)


def test_source_mismatch_stops(tmp_path):
    p=tmp_path/'x';p.write_bytes(b'changed')
    with pytest.raises(ValueError,match='SHA'):
        mod.verify_sha(p,'0'*64)


def test_constancy_and_nonfinite_are_not_zero_correlations():
    assert mod.rank_correlation([1,1,1],[1,2,3]) is None
    with pytest.raises(ValueError):
        mod.rank_correlation([1,2,np.nan],[1,2,3])
