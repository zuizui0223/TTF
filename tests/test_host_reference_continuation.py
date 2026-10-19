from pathlib import Path
import importlib.util, sys, json
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
SPEC=importlib.util.find_spec('acquire_lgm_host_references')

def module():
    assert SPEC is not None, 'host-reference continuation implementation is missing'
    import acquire_lgm_host_references as m
    return m

def match(name='Quercus robur',rank='SPECIES',key=10):
    return dict(matchType='EXACT',rank=rank,kingdom='Plantae',canonicalName=name,
                confidence=99,usageKey=key)

def record(key=10):
    return dict(datasetKey='50c9509d-22c7-4a22-a47d-8c48425ef4a7',basisOfRecord='HUMAN_OBSERVATION',
        occurrenceStatus='PRESENT',speciesKey=key,taxonKey=key,acceptedTaxonKey=key,issues=[],
        decimalLatitude=45.0,decimalLongitude=10.0,coordinateUncertaintyInMeters=10,
        eventDate='2020-06-20',occurrenceID='https://www.inaturalist.org/observations/123',key=123,
        establishmentMeans='INTRODUCED')

def test_plan_deduplicates_shared_host_without_dropping_insects():
    m=module();rows=[dict(insect_species=s,accepted_host_name='Quercus robur',wcvp_id_as_frozen=i)
        for s,i in [('A a','4e+05'),('B b','400000.0')]]
    p=m.host_plan(rows)
    assert len(p)==1 and p[0]['wcvp_id']=='400000'
    assert p[0]['insect_species']==['A a','B b']

@pytest.mark.parametrize('bad',['nan','1.5','0','-1','inf',''])
def test_bad_wcvp_id_rejected(bad):
    m=module()
    with pytest.raises(ValueError): m.normalize_id(bad)

def test_conflicting_id_names_rejected():
    m=module();rows=[dict(insect_species='A a',accepted_host_name=n,wcvp_id_as_frozen='10')
        for n in ['Quercus robur','Quercus rubra']]
    with pytest.raises(ValueError):m.host_plan(rows)

def test_subspecies_is_preserved_and_section_not_a_species():
    m=module()
    assert m.request_rank('Ulmus minor subsp. minor')=='SUBSPECIES'
    assert m.request_rank('Taraxacum sect. Taraxacum')=='REQUIRES_TAXONOMIC_REVIEW'
    assert m.request_rank('Quercus robur')=='SPECIES'

def test_exact_plant_match_accepted():
    m=module(); assert m.plant_key(match(),'Quercus robur','SPECIES')==10

@pytest.mark.parametrize('change',[{'kingdom':'Animalia'},{'matchType':'FUZZY'},{'rank':'GENUS'},
                                    {'canonicalName':'Quercus rubra'},{'confidence':90}])
def test_inappropriate_taxonomy_rejected(change):
    m=module();x=match();x.update(change)
    assert m.plant_key(x,'Quercus robur','SPECIES') is None

def test_subspecies_not_replaced_by_parent_species():
    m=module();name='Ulmus minor subsp. minor'
    assert m.plant_key(match(name,'SPECIES'),name,'SUBSPECIES') is None
    assert m.plant_key(match(name,'SUBSPECIES',20),name,'SUBSPECIES')==20
    r=record(10)
    assert m.clean_host_record(r,20,'SUBSPECIES')[0] is None
    r.update(taxonKey=20,acceptedTaxonKey=20)
    item,reason=m.clean_host_record(r,20,'SUBSPECIES')
    assert item is not None and item['species_key']==10

def test_plant_observation_status_retained_not_called_native():
    m=module();item,reason=m.clean_host_record(record(),10,'SPECIES')
    assert item['establishment_means']=='INTRODUCED'
    assert 'nearest_evaluation_km' not in item

def test_wrong_species_record_rejected():
    m=module();assert m.clean_host_record(record(11),10,'SPECIES')[0] is None

def test_acquisition_uses_fixed_windows_and_saves_raw_points():
    m=module();calls=[]
    def fetch(session,url,params,ledger):
        calls.append((url,params))
        if 'species/match' in url:return match()
        if params['limit']==0:return {'count':1}
        return {'count':1,'offset':0,'results':[record()]}
    result=m.acquire_host({'wcvp_id':'100','name':'Quercus robur','rank_hint':'SPECIES'},fetch=fetch)
    assert result['status']=='OBSERVATIONS_ACQUIRED'
    assert result['quality_records']==1 and len(result['points'])==1
    assert calls[0][1]['rank']=='SPECIES'
    assert result['spatial_filter_applied'] is False
    assert result['native_range_filter_applied'] is False

def test_section_is_review_not_zero_observations():
    m=module()
    def forbidden(*args):raise AssertionError('must not reinterpret section as species')
    result=m.acquire_host({'wcvp_id':'200','name':'Taraxacum sect. Taraxacum',
                           'rank_hint':'REQUIRES_TAXONOMIC_REVIEW'},fetch=forbidden)
    assert result['status']=='TAXONOMIC_RANK_REVIEW'
    assert result['quality_records'] is None

def test_request_error_remains_unresolved():
    m=module()
    def fail(*args):raise OSError('offline fixture')
    r=m.acquire_host({'wcvp_id':'100','name':'Quercus robur','rank_hint':'SPECIES'},fetch=fail)
    assert r['status']=='REQUEST_ERROR' and r['quality_records'] is None

def test_checkpoint_identity_and_tamper(tmp_path):
    m=module();r={'wcvp_id':'10','name':'Quercus robur','status':'REQUEST_ERROR','points':[],
                'quality_records':None}
    m.save_checkpoint(tmp_path,'10',{'contract':'a'},r)
    assert m.load_checkpoint(tmp_path,'10',{'contract':'a'})==r
    with pytest.raises(ValueError):m.load_checkpoint(tmp_path,'10',{'contract':'b'})
    f=tmp_path/'10.json';v=json.loads(f.read_text());v['result']['name']='Wrong';f.write_text(json.dumps(v))
    with pytest.raises(ValueError):m.load_checkpoint(tmp_path,'10',{'contract':'a'})

def test_completed_error_is_not_retried(tmp_path):
    m=module();p=[{'wcvp_id':'10','name':'Quercus robur','rank_hint':'SPECIES'}]
    r={'wcvp_id':'10','name':'Quercus robur','status':'REQUEST_ERROR','points':[], 'quality_records':None}
    m.save_checkpoint(tmp_path,'10',{'x':1},r)
    def forbidden(_):raise AssertionError('completed request failure must not be silently re-requested')
    completed,pending=m.run_queue(p,tmp_path,{'x':1},acquire=forbidden,max_seconds=1)
    assert completed==[r] and pending==[]

def test_observation_count_reconciled(tmp_path):
    m=module();r={'wcvp_id':'10','name':'Quercus robur','status':'OBSERVATIONS_ACQUIRED',
                'quality_records':2,'points':[record()]}
    with pytest.raises(ValueError):m.save_checkpoint(tmp_path,'10',{},r)
