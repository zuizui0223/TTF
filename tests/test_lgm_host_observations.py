import importlib.util
import json
import gzip
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))

def module():
    p=Path(__file__).resolve().parents[1]/'scripts'/'acquire_lgm_host_observations.py'
    assert p.exists(), 'host acquisition module is not implemented'
    spec=importlib.util.spec_from_file_location('hosts',p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def match():
    return dict(matchType='EXACT',rank='SPECIES',kingdom='Plantae',canonicalName='Picea abies',confidence=99,usageKey=123)

def test_exact_plant_key():
    m=module(); assert m.plant_key(match(),'Picea abies')==123

def test_synonym_accepted_key_is_explicit():
    m=module(); q=match();q['acceptedUsageKey']=456
    assert m.plant_key(q,'Picea abies')==456

@pytest.mark.parametrize('field,value',[('rank','GENUS'),('kingdom','Animalia'),('matchType','FUZZY'),('confidence',94),('canonicalName','Picea sitchensis')])
def test_wrong_plant_match_rejected(field,value):
    m=module();q=match();q[field]=value
    assert m.plant_key(q,'Picea abies') is None

def test_manifest_reuses_shared_host():
    m=module();rows=[dict(insect_species=s,accepted_host_name='Picea abies',wcvp_id_as_frozen=i) for s,i in [('A a','123.0'),('B b','1.23e2')]]
    assert m.host_map(rows)=={'123':'Picea abies'}

@pytest.mark.parametrize('bad',['nan','inf','-1','0','1.5',''])
def test_invalid_wcvp_ids(bad):
    m=module()
    with pytest.raises(ValueError):m.host_map([dict(insect_species='A a',accepted_host_name='Picea abies',wcvp_id_as_frozen=bad)])

def test_conflicting_manifest_rejected():
    m=module();r=[dict(insect_species='A a',accepted_host_name=n,wcvp_id_as_frozen='123') for n in ['Picea abies','Picea sitchensis']]
    with pytest.raises(ValueError):m.host_map(r)

def test_checkpoint_roundtrip_and_binding(tmp_path):
    m=module();v={'host_id':'123','status':'RETRIEVAL_COMPLETE','points':[],'ledger':[]}
    m.save_result(tmp_path,'123',{'x':'one'},v)
    assert m.load_result(tmp_path,'123',{'x':'one'})==v
    with pytest.raises(ValueError):m.load_result(tmp_path,'123',{'x':'two'})

def test_corrupt_checkpoint_rejected(tmp_path):
    m=module();v={'host_id':'123','status':'REQUEST_ERROR','points':[],'ledger':[]}
    m.save_result(tmp_path,'123',{},v)
    path=tmp_path/'123.json.gz';q=json.loads(gzip.decompress(path.read_bytes()));q['result']['status']='RETRIEVAL_COMPLETE';path.write_bytes(gzip.compress(json.dumps(q).encode(),mtime=0))
    with pytest.raises(ValueError):m.load_result(tmp_path,'123',{})

def test_completed_failures_not_retried(tmp_path):
    m=module();calls=[]
    def f(k,n):calls.append(k);return dict(host_id=k,status='REQUEST_ERROR',points=[],ledger=[])
    mapping={'123':'Picea abies','456':'Picea sitchensis'}
    result=m.collect(mapping,tmp_path,{},f,workers=1,max_seconds=10)
    result2=m.collect(mapping,tmp_path,{},f,workers=1,max_seconds=10)
    assert calls==['123','456']; assert result==result2

def test_record_query_filter_keeps_same_insect_criteria():
    m=module();q=m.observation_query(123)
    assert q['datasetKey']==m.base.DATASET
    assert q['eventDate']=='2010-01-01,2026-10-18'
    assert q['coordinateUncertaintyInMeters']=='0,1000'
    assert q['basisOfRecord']=='HUMAN_OBSERVATION'
    assert q['taxonKey']==123

def test_clean_host_record_not_claimed_spatially_separated():
    m=module();raw=dict(datasetKey=m.base.DATASET,basisOfRecord='HUMAN_OBSERVATION',occurrenceStatus='PRESENT',speciesKey=123,
    decimalLatitude=45,decimalLongitude=10,coordinateUncertaintyInMeters=10,eventDate='2026-01-02',occurrenceID='https://www.inaturalist.org/observations/9',key=42,license='CC0_1_0')
    row,reason=m.base.clean_record(raw,123)
    assert reason=='retained' and 'nearest_evaluation_km' not in row

def test_empty_success_is_not_request_error(monkeypatch):
    m=module()
    responses=iter([match(),{'count':0}])
    monkeypatch.setattr(m.base,'get_json',lambda *a:next(responses))
    result=m.acquire_host('123','Picea abies')
    assert result['status']=='RETRIEVAL_COMPLETE'
    assert result['quality_unique_records']==0 and result['points']==[]

def test_taxonomy_failure_does_not_query_occurrences(monkeypatch):
    m=module();q=match();q['kingdom']='Animalia';calls=[]
    def fetch(*a):calls.append(a[1]);return q
    monkeypatch.setattr(m.base,'get_json',fetch)
    result=m.acquire_host('123','Picea abies')
    assert result['status']=='TAXONOMY_NOT_EXACT' and len(calls)==1

def test_transport_failure_is_not_zero_observations(monkeypatch):
    m=module()
    def fetch(*args): raise RuntimeError('test transport failure')
    monkeypatch.setattr(m.base,'get_json',fetch)
    result=m.acquire_host('123','Picea abies')
    assert result['status']=='REQUEST_ERROR' and 'quality_unique_records' not in result
