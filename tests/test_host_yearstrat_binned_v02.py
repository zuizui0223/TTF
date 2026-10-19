import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import pytest
from acquire_hosts_yearstrat_v02 import BINS,make_query,acquire_binned

DATASET = '50c9509d-22c7-4a22-a47d-8c48425ef4a7'

def rec(key, ident, year):
    return {
        'key':ident, 'occurrenceID':f'https://www.inaturalist.org/observations/{ident}',
        'datasetKey':DATASET, 'basisOfRecord':'HUMAN_OBSERVATION',
        'occurrenceStatus':'PRESENT','speciesKey':key,
        'taxonKey':key,'acceptedTaxonKey':key,'issues':[],
        'decimalLatitude':45.5,'decimalLongitude':10.25,
        'coordinateUncertaintyInMeters':500,'eventDate':f'{year}-06-16'
    }

def plant_match():
    return {'matchType':'EXACT','rank':'SPECIES','kingdom':'Plantae',
            'canonicalName':'Quercus robur','confidence':99,'usageKey':10}

def test_year_bins_partition_original_period():
    assert BINS == ((2010,2015),(2016,2020),(2021,2026))
    assert [y for start,end in BINS for y in range(start,end+1)] == list(range(2010,2027))

def test_no_deep_offsets_and_same_data_quality_filters():
    for start,end in BINS:
        q=make_query(10,start,end)
        assert (q['offset'],q['limit'])==(0,300)
        assert q['year']==f'{start},{end}'
        assert q['datasetKey']==DATASET
        assert q['coordinateUncertaintyInMeters']=='0,1000'
        assert q['basisOfRecord']=='HUMAN_OBSERVATION'
        assert q['occurrenceStatus']=='PRESENT'
        assert q['eventDate'].split(',')[0]==f'{start}-01-01'
    assert make_query(10,2021,2026)['eventDate'].endswith('2026-10-18')

def test_shallow_acquisition_preserves_individual_record_provenance_without_spatial_filter():
    calls=[]
    def fetch(session,url,params,ledger):
        calls.append((url,params))
        if 'species/match' in url:return plant_match()
        year=int(params['year'].split(',')[0])
        return {'count':1,'offset':0,'results':[rec(10,year,year)]}
    r=acquire_binned({'wcvp_id':'12','name':'Quercus robur','rank_hint':'SPECIES'},fetch=fetch)
    assert r['status']=='OBSERVATIONS_ACQUIRED'
    assert r['quality_records']==3
    assert {p['gbif_id'] for p in r['points']}=={2010,2016,2021}
    assert len(calls)==4
    assert all(p['offset']==0 for url,p in calls if 'occurrence/search' in url)
    assert not r['spatial_filter_applied'] and not r['native_range_filter_applied']
    assert not r['climate_or_genetic_response_used']
    assert r['scheme']=='THREE_NONOVERLAPPING_YEAR_BINS_SHALLOW_OFFSET_ZERO_V02'

def test_duplicated_gbif_observation_is_counted_once():
    def fetch(session,url,params,ledger):
        if 'species/match' in url:return plant_match()
        year=int(params['year'].split(',')[0])
        return {'count':1,'offset':0,'results':[rec(10,5,year)]}
    r=acquire_binned({'wcvp_id':'12','name':'Quercus robur','rank_hint':'SPECIES'},fetch=fetch)
    assert r['quality_records']==1
    assert r['rejection_counts']['duplicate_gbif_id']==2

def test_section_name_is_not_silently_promoted_to_species():
    result=acquire_binned({'wcvp_id':'12','name':'Taraxacum sect. Taraxacum',
          'rank_hint':'REQUIRES_TAXONOMIC_REVIEW'}, fetch=lambda *args:pytest.fail('wrong rank lookup'))
    assert result['status']=='TAXONOMIC_RANK_REVIEW'
    assert result['quality_records'] is None

def test_query_matching_zero_is_not_confused_with_species_absence():
    def fetch(session,url,params,ledger):
        if 'species/match' in url:return plant_match()
        return {'count':0,'offset':0,'results':[]}
    r=acquire_binned({'wcvp_id':'12','name':'Quercus robur','rank_hint':'SPECIES'},fetch=fetch)
    assert r['status']=='OBSERVATIONS_ACQUIRED'
    assert r['quality_records']==0 and r['points']==[]

def test_transport_error_never_claims_observed_zero():
    def fail(*args):raise OSError('fixture network unavailable')
    r=acquire_binned({'wcvp_id':'12','name':'Quercus robur','rank_hint':'SPECIES'},fetch=fail)
    assert r['status']=='REQUEST_ERROR'
    assert r['quality_records'] is None and r['points']==[]
