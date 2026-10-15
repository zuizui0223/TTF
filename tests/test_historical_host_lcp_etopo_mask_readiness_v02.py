"""Only synthetic ETOPO tiles; no actual climate, topo or genetic data used."""
import sys
from pathlib import Path
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_historical_host_lcp_etopo_mask_readiness_v02 import (
    audit_native_topography, FixedSourceError, verify_source_pin,
)


def tile(tmp_path,crs='EPSG:9518',nodata=False):
    dst=tmp_path/'fixture.tif'
    a=np.full((10,10),120,dtype=np.int16)
    a[4,3]=-50
    a[5,4]=-140
    a[5,5]=0
    if nodata:a[5,6]=-9999
    with rasterio.open(dst,'w',driver='GTiff',width=10,height=10,count=1,
      dtype='int16',crs=crs,transform=from_origin(-86,12,0.1,0.1),nodata=-9999) as ds:
        ds.write(a,1)
    return dst


def test_compound_9518_coordinate_grid_and_thresholds(tmp_path):
    z=audit_native_topography(tile(tmp_path),[(-85.7,11.7),(-85.4,11.4)])
    w=z['native_etopo_window']
    assert z['native_source_crs']=='EPSG:9518'
    assert w['valid_elevation_cells']>0
    assert w['modern_land_cells']+w['current_submerged_cells']==w['valid_elevation_cells']
    assert w['static_minus120m_land_cells']==w['modern_land_cells']+w['newly_emerged_under_static_minus120m_cells']
    assert w['newly_emerged_under_static_minus120m_cells']>0
    assert z['endpoint_sampling']['n_original_localities']==2
    assert z['scientific_boundary']['not_a_verified_paleo_coastline'] is True
    assert z['scientific_boundary']['original_climate_lcp_unchanged'] is True
    assert z['scientific_boundary']['COI_genetic_effect_used'] is False


def test_plain_4326_rejected_to_avoid_datum_drift(tmp_path):
    with pytest.raises(FixedSourceError,match='9518'):
        audit_native_topography(tile(tmp_path,crs='EPSG:4326'),[(-85.7,11.7),(-85.4,11.4)])


def test_wrong_tile_or_uncovered_roi_fails_closed(tmp_path):
    with pytest.raises(FixedSourceError,match='ROI outside'):
        audit_native_topography(tile(tmp_path),[(-85.2,11.2),(-85.1,11.1)])


def test_static_sea_level_is_unchangeable(tmp_path):
    with pytest.raises(ValueError,match='fixed original'):
        audit_native_topography(tile(tmp_path),[(-85.7,11.7),(-85.4,11.4)],threshold_m=-100)
    with pytest.raises(ValueError,match='fixed original'):
        audit_native_topography(tile(tmp_path),[(-85.7,11.7),(-85.4,11.4)],padding=0)


def test_nodata_is_not_counted_as_submerged(tmp_path):
    z=audit_native_topography(tile(tmp_path,nodata=True),[(-85.7,11.7),(-85.4,11.4)])
    w=z['native_etopo_window']
    assert w['valid_elevation_cells']==w['cells']-1
    assert w['current_submerged_cells']>=2


def test_missing_hash_or_wrong_file_refuses_computation(tmp_path):
    fake=tile(tmp_path)
    with pytest.raises((FixedSourceError,RuntimeError),match='source SHA256|source identity'):
        verify_source_pin(fake,{'exact_sha256':'0'*64,'exact_bytes':fake.stat().st_size,
                                'canonical_file_name':fake.name,
                                'header':{'compound_crs':'EPSG:9518','bounds':[-86,11,-85,12],'shape':[10,10]}},
                          source_only=True)
