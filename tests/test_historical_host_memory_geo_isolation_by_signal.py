"""Synthetic geography, coefficients and locality-deletion statuses only."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))

import numpy as np
import pandas as pd
import pytest
from audit_historical_host_memory_geo_isolation_by_signal import rank_average, stratified_geographic_fragility


def fixture(n_species=8):
    species=[]; deletions=[]; localities=[];edges=[]
    xyz=[0,1,2,3,4,5,6,20]
    for i in range(n_species):
        sp=f'Sp{i:02d}'
        beta=(-1)**i*(i+1)*.075
        slope=beta*.005
        flip=i%2==0
        species.append({
            'species':sp,'edges':28,'localities':8,'standardized_beta':beta,
            'semipartial_correlation':beta*np.sqrt(.5),
            'unique_fraction':.5,'unscaled_slope':slope,
            'locality_deletions_reversing_sign':int(flip),
            'locality_deletions_defined':8,'locality_deletions_undefined':0,
        })
        for j,x in enumerate(xyz):
            localities.append({'species':sp,'locality_index':j,'x_km':x,'y_km':0.,'z_km':0.})
            deletions.append({'species':sp,'locality_index':j,'status':'defined',
                              'slope': -slope if flip and j==7 else slope})
        ctr=0
        for a in range(8):
            for b in range(a+1,8):
                edges.append({'species':sp,'edge_index':ctr,'node_left':a,'node_right':b});ctr+=1
    return tuple(pd.DataFrame(t) for t in [species,deletions,localities,edges])


def test_within_species_rank_ties_are_averaged():
    np.testing.assert_allclose(rank_average([1.,1.,2.]),[.5,.5,1.])


def test_aggregate_stratifies_previously_defined_beta_and_no_extra_test():
    r=stratified_geographic_fragility(*fixture())
    assert r['status']=='EXPLORATORY_POST_OUTCOME_SAMPLING_GEOGRAPHY'
    assert r['original_confirmatory_status']=='NOT_EVALUABLE_EMPIRICAL_RESPONSE'
    assert r['n_frozen_survivors']==8 and r['n_nonzero_defined']==8
    assert r['flip_event_total']==4
    assert sum(x['species'] for x in r['semipartial_quartile_geography'])==8
    assert all(x['species']==2 and x['flip_events']==1 for x in r['semipartial_quartile_geography'])
    assert all(x['mean_flip_nearest_rank']>x['mean_other_nearest_rank']
               for x in r['semipartial_quartile_geography'])
    assert all(x['paired_median_rank_difference']>0 for x in r['semipartial_quartile_geography'])
    assert 'Sp00' not in str(r)
    assert r['new_p_values']==0 and r['new_hypothesis_tests']==0


def test_swapped_frozen_coordinates_and_drifted_edges_fail_closed():
    s,d,l,e=fixture()
    bad=l.drop(l.index[0])
    with pytest.raises(ValueError,match='frozen locality index'):
        stratified_geographic_fragility(s,d,bad,e)
    bad=e.drop(e.index[0])
    with pytest.raises(ValueError,match='original frozen edge'):
        stratified_geographic_fragility(s,d,l,bad)
    bad=e.copy();bad.loc[0,'node_left']=100
    with pytest.raises(ValueError,match='nonexistent frozen locality'):
        stratified_geographic_fragility(s,d,l,bad)


def test_wrong_reversal_tally_and_undefined_values_are_not_silently_zeroed():
    s,d,l,e=fixture()
    s.loc[0,'locality_deletions_reversing_sign']=0
    with pytest.raises(ValueError,match='source reversal count mismatch'):
        stratified_geographic_fragility(s,d,l,e)
    s,d,l,e=fixture()
    d.loc[0,'status']='undefined';d.loc[0,'slope']=np.nan
    s.loc[0,'locality_deletions_undefined']=1
    s.loc[0,'locality_deletions_defined']=7
    # The first record was not a sign-reversal site, so the original flip tally holds.
    out=stratified_geographic_fragility(s,d,l,e)
    assert out['undefined_deletion_variants']==1
    assert out['defined_site_deletions']==63


def test_original_four_invariants_stay_separate_from_exploratory_quartiles():
    s,d,l,e=fixture()
    # Replace first species with a zero-p-distance invariant, without deleting it.
    s.loc[0,['standardized_beta','semipartial_correlation']]=np.nan
    s.loc[0,'unscaled_slope']=0.
    s.loc[0,'locality_deletions_reversing_sign']=0
    d.loc[d.species=='Sp00','slope']=0.
    r=stratified_geographic_fragility(s,d,l,e)
    assert r['n_frozen_survivors']==8
    assert r['n_invariant']==1
    assert r['n_nonzero_defined']==7
    assert sum(q['species'] for q in r['semipartial_quartile_geography'])==7
