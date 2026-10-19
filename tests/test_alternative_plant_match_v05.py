import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from alternative_plant_match_v05 import (genus_of,stable_sample_indices,
    climate_signature,distance_to_signature,rank_alternatives,
    group_alternatives,observation_signature,equal_reference_indices)

def test_genus_exclusion_of_all_known_hosts():
    taxa={'a':{'name':'Quercus robur'},'b':{'name':'Quercus rubra'},
          'c':{'name':'Acer rubrum'},'d':{'name':'Betula lenta'}}
    x={'center':[0,0,0,0],'spread':[1,1,1,1],'n':60}
    res=rank_alternatives(x,{'a'},taxa,{k:x for k in taxa})
    assert [k for _,k in res]==['c','d']

def test_plant_moment_distance_zero_for_identical():
    x=climate_signature(np.ones((45,4)))
    assert distance_to_signature(x,x)==pytest.approx(0)
    assert climate_signature(np.ones((29,4))) is None
    assert climate_signature(np.full((30,4),np.nan)) is None

def test_observation_sampling_and_groups_are_deterministic():
    obs=[{'gbif_id':i} for i in range(400)]
    assert stable_sample_indices(obs)==stable_sample_indices(obs)
    assert len(stable_sample_indices(obs))==64
    xs=group_alternatives([str(i) for i in range(50)],3)
    assert len(xs)==8 and all(len(x)==3 for x in xs)
    assert len(set(sum(xs,[])))==24

def test_alternative_exclusion_does_not_rename_hosts():
    taxa={'123':{'name':'Ulmus minor subsp. minor'},
          '456':{'name':'Ulmus americana'},
          '789':{'name':'Picea abies'}}
    sig={'center':[.1]*4,'spread':[.2]*4,'n':64}
    assert [x[1] for x in rank_alternatives(sig,{'123'},taxa,{k:sig for k in taxa})]==['789']

def test_observation_digest_rejects_duplicate_GBIF_record():
    obs=[{'gbif_id':3},{'gbif_id':4}]
    assert observation_signature(obs)==observation_signature(obs)
    with pytest.raises(ValueError):observation_signature(obs+obs)

def test_equal_reference_counts_and_draw_synchronization():
    a,b,c=equal_reference_indices(80,140,120,'Species test',4,2)
    assert len(a)==len(b)==len(c)==60
    assert len(set(a))==len(a)
    d=equal_reference_indices(80,140,120,'Species test',4,2)
    assert all(np.array_equal(x,y) for x,y in zip((a,b,c),d))

@pytest.mark.parametrize('a,b,c',[(29,200,200),(200,29,200),(200,200,59)])
def test_ineligible_clouds_fail_closed(a,b,c):
    with pytest.raises(ValueError):equal_reference_indices(a,b,c,'A b',0,0)

def test_raw_plant_checkpoint_gets_frozen_taxon_identity_before_thinning():
    # Regression: directly calling filter_host_points(plants[id],...) once failed
    # KeyError('wcvp_id') because raw v0.2 observations do not carry WCVP IDs.
    from alternative_plant_match_v05 import tag_donor_reference
    source={'400000':[{'gbif_id':'123','latitude':45.0,'longitude':10.0,
                        'uncertainty_m':100}]}
    labeled=tag_donor_reference('400000',source)
    assert labeled[0]['wcvp_id']=='400000'
    assert labeled[0]['gbif_id']=='123'
    assert source['400000'][0].get('wcvp_id') is None
