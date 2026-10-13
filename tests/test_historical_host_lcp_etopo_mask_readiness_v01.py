"""Synthetic raster-test only: ETOPO bedrock/surface mask is not 21ka coastlines."""
from pathlib import Path
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from scripts.audit_historical_host_lcp_etopo_mask_readiness import audit_topography_mask, FixedLcpGeometryError


def artificial_tif(tmp_path):
    raster=tmp_path/'etopo_mock.tif'
    arr=np.full((10,10),120,dtype=np.int16)
    arr[4,3]=-50
    arr[5,4]=-140
    arr[5,5]=0
    with rasterio.open(raster,'w',driver='GTiff',width=10,height=10,
                       count=1,dtype='int16',crs='EPSG:4326',
                       transform=from_origin(-86,12,0.1,0.1),nodata=-32767) as ds:
        ds.write(arr,1)
    return raster


def test_static_topography_and_eustatic_sensitivity_are_distinct(tmp_path):
    path=artificial_tif(tmp_path)
    sites=[(-85.7,11.7),(-85.4,11.4)]
    obj=audit_topography_mask(path,sites,0.12,threshold_m=-120)
    c=obj['native_etopo_window']
    assert c['cells']>0 and c['modern_land_cells'] > 0
    assert c['current_submerged_cells']>=2
    assert c['newly_emerged_under_static_minus120m_cells']>=1
    assert c['static_minus120m_land_cells'] >= c['modern_land_cells']
    assert obj['scientific_boundary']['not_a_verified_paleo_coastline'] is True
    assert obj['scientific_boundary']['original_climate_lcp_unchanged'] is True
    assert obj['scientific_boundary']['COI_genetic_effect_used'] is False


def test_invalid_grid_reference_or_insufficient_extent_rejected(tmp_path):
    path=artificial_tif(tmp_path)
    with pytest.raises(FixedLcpGeometryError,match='ROI outside'):
        audit_topography_mask(path,[(-85.1,11.1),(-85.2,11.2)],0.12,threshold_m=-120)
    with pytest.raises((FixedLcpGeometryError,ValueError)):
        audit_topography_mask(path,[(-85.7,11.7),(-85.4,11.4)],0.0,threshold_m=-120)
    with pytest.raises((FixedLcpGeometryError,ValueError)):
        audit_topography_mask(path,[(-85.7,11.7),(-85.4,11.4)],0.12,threshold_m=-100)


def test_all_endpoints_are_assessed_without_snapping(tmp_path):
    path=artificial_tif(tmp_path)
    obj=audit_topography_mask(path,[(-85.65,11.65),(-85.55,11.55)],0.12,threshold_m=-120)
    assert obj['endpoint_sampling']['n_original_localities']==2
    assert obj['endpoint_sampling']['n_finite_elevations']==2
    assert obj['endpoint_sampling']['n_land_above_present_sea_level'] <= 2


def test_exact_static_eustatic_rule_must_be_hash_bound(tmp_path):
    from scripts.audit_historical_host_lcp_etopo_mask_readiness import verify_contract
    contract=Path(__file__).resolve().parents[1]/'docs/supporting/historical_host_lcp_etopo_static_mask_contract_v0.1.json'
    assert verify_contract(contract) is True
    bad=tmp_path/'forged.json';bad.write_bytes(contract.read_bytes().replace(b'-120',b'-100'))
    with pytest.raises(RuntimeError,match='frozen ETOPO mask contract'):
        verify_contract(bad)
