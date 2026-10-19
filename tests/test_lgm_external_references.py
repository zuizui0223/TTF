import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0, str(Path(__file__).parents[1]/'scripts'))
from acquire_lgm_external_references import window_offsets, unit_vectors, nearest_km, clean_record, thin_records, accepted_key

DATASET='50c9509d-22c7-4a22-a47d-8c48425ef4a7'
def record(**changes):
    d=dict(key=10, speciesKey=100, datasetKey=DATASET, basisOfRecord='HUMAN_OBSERVATION',
           occurrenceStatus='PRESENT',decimalLatitude=40,decimalLongitude=-75,
           coordinateUncertaintyInMeters=100,eventDate='2020-05-05',
           occurrenceID='https://www.inaturalist.org/observations/10',
           license='http://creativecommons.org/licenses/by/4.0/legalcode',issues=[])
    d.update(changes); return d

def test_window_complete_small():
    assert window_offsets(0)==[]
    assert window_offsets(200)==[0]
    assert window_offsets(601)==[0,300,600]

def test_windows_fixed_nonoverlap_and_limit():
    assert window_offsets(901)==[0,300,601]
    assert window_offsets(100001)==[0,49850,99700]

@pytest.mark.parametrize('bad', [-1, -20])
def test_bad_count(bad):
    with pytest.raises(ValueError): window_offsets(bad)

def test_spherical_dateline():
    d=nearest_km(np.array([[0,179.9],[0,0]]),np.array([[0,-179.9]]))
    assert d[0]==pytest.approx(22.239,abs=.01)
    assert d[1]>19000

def test_spherical_identity_and_norm():
    a=np.array([[40,-75],[90,130],[-90,-180]])
    assert np.allclose(np.linalg.norm(unit_vectors(a),axis=1),1)
    assert np.allclose(nearest_km(a,a),0)

def test_quality_retains_required_fields_only():
    d,reason=clean_record(record(recordedBy='unneeded personal attribution'),100)
    assert reason=='retained'
    assert d['gbif_id']==10 and d['latitude']==40
    assert 'recordedBy' not in d

@pytest.mark.parametrize('change',[
    {'basisOfRecord':'MATERIAL_SAMPLE'}, {'datasetKey':'other'},
    {'speciesKey':999}, {'occurrenceStatus':'ABSENT'},
    {'decimalLatitude':91}, {'decimalLongitude':float('nan')},
    {'coordinateUncertaintyInMeters':None}, {'coordinateUncertaintyInMeters':1001},
    {'coordinateUncertaintyInMeters':-1}, {'eventDate':'2009-12-31'},
    {'eventDate':'2026-10-19'}, {'eventDate':None},
    {'issues':['ZERO_COORDINATE']}, {'occurrenceID':'https://example.org/10'},
    {'decimalLatitude':0,'decimalLongitude':0}, {'eventDate':'2026-13-01'},
])
def test_reject_ineligible(change):
    assert clean_record(record(**change),100)[0] is None

def test_clearance_uncertainty_and_dedup():
    def row(i,lat,near,unc=.1):
        return dict(gbif_id=i,occurrence_id=f'https://www.inaturalist.org/observations/{i}',
                    latitude=lat,longitude=0,nearest_evaluation_km=near,uncertainty_m=unc*1000)
    rows=[row(1,0,50.05),row(2,1,50.2),row(3,1,60),row(4,1.01,100),row(5,2,100)]
    kept=thin_records(rows,50)
    assert len(kept)==2
    assert all(x['nearest_evaluation_km']-x['uncertainty_m']/1000>=50 for x in kept)
    assert len(thin_records(rows,50,cap=1))==1

def test_thinning_is_order_invariant():
    rows=[dict(gbif_id=i,occurrence_id=f'https://www.inaturalist.org/observations/{i}',latitude=i/10,
               longitude=0,nearest_evaluation_km=100,uncertainty_m=100) for i in range(20)]
    a=thin_records(rows,50); b=thin_records(list(reversed(rows)),50)
    assert a==b
    d=nearest_km(np.array([[a[0]['latitude'],0]]),np.array([[a[-1]['latitude'],0]]))
    assert d[0]>=10

def test_taxonomy_exact_only():
    m=dict(usageKey=123,rank='SPECIES',class_='Insecta',matchType='EXACT',confidence=100,
           canonicalName='Bungalotis erythus')
    m['class']=m.pop('class_')
    assert accepted_key(m,'Bungalotis erythus')==123
    assert accepted_key(dict(m,acceptedUsageKey=456),'Bungalotis erythus')==456
    assert accepted_key(dict(m,matchType='FUZZY'),'Bungalotis erythus') is None
    assert accepted_key(dict(m,rank='GENUS'),'Bungalotis erythus') is None
    assert accepted_key(dict(m,canonicalName='Wrong name'),'Bungalotis erythus') is None
