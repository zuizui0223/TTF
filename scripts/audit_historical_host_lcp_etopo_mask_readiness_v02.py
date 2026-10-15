#!/usr/bin/env python3
"""Corrected exact-source ETOPO topography screen — not a 21 ka coastline.

v0.1's N00W090 tile covers -15..0 latitude, whereas the frozen 26
Bungalotis localities lie 10.455..11.05 N. ETOPO GeoTIFF CRS EPSG:9518
means horizontal WGS84 degrees PLUS vertical EGM2008 metres. This v0.2
only corrects the source and CRS handling. Scientific and sea-level
thresholds remain identical. No nucleotide, p-distance, or COI response.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window, from_bounds

SOURCE_SHA = '2e0f33efebcbe6cc50f3a6dc7c13c78f46e04f648b5f6f3a5ca180c8a5c2f35c'
SOURCE_BYTES = 23_311_599
SOURCE_BASENAME = 'ETOPO_2022_v1_15s_N15W090_surface.tif'
SOURCE_URL = ('https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/15s/'
              '15s_surface_elev_gtif/' + SOURCE_BASENAME)
CONTRACT_SHA = 'af9a7c3aa0b6c9708d9fc26e13a1ca52ac8181f3254a4d4d92618792dd9431a6'
GEOMETRY_SHA = '037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9'
SPECIES = 'Bungalotis erythus'
PADDING = 0.12
SEA_LEVEL = -120


class FixedSourceError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20), b''):
            h.update(chunk)
    return h.hexdigest()


def require_native_header(ds, *, enforce_full_tile:bool=False) -> dict:
    if not ds.crs or ds.crs.to_epsg()!=9518 or ds.count!=1:
        raise FixedSourceError('ETOPO source must be single-band EPSG:9518, with WGS84 horizontal and EGM2008 vertical coordinates')
    t=ds.transform
    if abs(t.b)>1e-14 or abs(t.d)>1e-14 or t.a<=0 or t.e>=0:
        raise FixedSourceError('ETOPO source must use original north-up geographic grid')
    out=dict(driver=ds.driver, crs=str(ds.crs), epsg=ds.crs.to_epsg(),
             bounds=[float(z) for z in ds.bounds],shape=[int(ds.height),int(ds.width)],
             pixel_size=[float(t.a),float(t.e)],nodata=ds.nodata)
    if enforce_full_tile:
        if (out['driver']!='GTiff' or out['shape']!=[3600,3600]
                or any(abs(a-b)>1e-8 for a,b in zip(out['bounds'],[-90,0,-75,15]))
                or abs(t.a-1/240)>1e-10 or abs(t.e+1/240)>1e-10):
            raise FixedSourceError('ETOPO corrected N15W090 tile native footprint, resolution or format mismatch')
    return out


def verify_source_pin(source:Path, pin:dict, *, source_only:bool=False) -> dict:
    """Verify exact source *before* reading any ETOPO elevation pixel.

    source_only is reserved for isolated verification tests; it does not relax
    digest, size or header checks. In the actual CLI the contract is also pinned.
    """
    p=Path(source)
    if sha256(p)!=SOURCE_SHA or p.stat().st_size!=SOURCE_BYTES:
        raise FixedSourceError('ETOPO exact source SHA256 identity mismatch')
    if p.name!=SOURCE_BASENAME:
        raise FixedSourceError('ETOPO exact source filename mismatch')
    if (pin.get('schema')!='ttf_etopo_N15W090_source_pin_v0.2'
            or pin.get('exact_sha256')!=SOURCE_SHA
            or pin.get('exact_bytes')!=SOURCE_BYTES
            or pin.get('canonical_file_name')!=SOURCE_BASENAME
            or pin.get('frozen_contract_sha256')!=CONTRACT_SHA
            or pin.get('source_stage')!='SOURCE_BYTES_AND_HEADER_ONLY_NO_ELEVATION_VALUES'):
        raise FixedSourceError('ETOPO source pin provenance does not match independently recorded metadata')
    with rasterio.open(p) as ds:
        return require_native_header(ds,enforce_full_tile=True)


def audit_native_topography(source:Path, lonlat:list[tuple[float,float]],
                            padding:float=PADDING, *, threshold_m:int=SEA_LEVEL) -> dict:
    sites=np.asarray(lonlat,dtype=float)
    if (sites.ndim!=2 or sites.shape[1]!=2 or len(sites)<2 or not np.isfinite(sites).all()
            or padding!=PADDING or threshold_m!=SEA_LEVEL):
        raise ValueError('fixed original sampling geometry, ROI padding or -120m threshold changed')
    if np.any(sites[:,0]<-180) or np.any(sites[:,0]>180) or np.any(sites[:,1]<-90) or np.any(sites[:,1]>90):
        raise ValueError('invalid geographic endpoint coordinate')
    with rasterio.open(source) as ds:
        meta=require_native_header(ds)
        lon,lat=sites[:,0],sites[:,1]
        bbox=[float(lon.min()-padding),float(lat.min()-padding),float(lon.max()+padding),float(lat.max()+padding)]
        win=from_bounds(*bbox,transform=ds.transform)
        r0=int(math.floor(win.row_off));c0=int(math.floor(win.col_off))
        r1=int(math.ceil(win.row_off+win.height));c1=int(math.ceil(win.col_off+win.width))
        if r0<0 or c0<0 or r1>ds.height or c1>ds.width:
            raise FixedSourceError('ROI outside fixed native ETOPO tile; no replacement or relocation')
        rows,cols=r1-r0,c1-c0
        if rows*cols<=0 or rows*cols>150_000:
            raise FixedSourceError('ETOPO native window exceeds frozen safety budget')
        src=ds.read(1,window=Window(c0,r0,cols,rows),masked=True)
        z=np.asarray(np.ma.filled(src.astype('float64'),np.nan))
        valid=np.isfinite(z)
        if not np.any(valid):
            raise FixedSourceError('ETOPO has no valid elevation cells in frozen ROI')
        original_samples=list(ds.sample([tuple(row) for row in sites],masked=True))
        e=np.asarray([float(s[0]) if not np.ma.is_masked(s[0]) else np.nan for s in original_samples])
        if not np.isfinite(e).all():
            raise FixedSourceError('Original geographic endpoints cannot be sampled on verified ETOPO')
        land=valid&(z>0)
        lowered=valid&(z>SEA_LEVEL)
        additional=valid&(z>SEA_LEVEL)&(z<=0)
        return {
            'schema':'ttf_hhm_etopo_N15W090_static_topographic_mask_v0.2',
            'status':'EXPLORATORY_ETOPO_PRESENT_LAND_AND_ILLUSTRATIVE_MINUS120M_ONLY',
            'native_source_crs':meta['crs'],
            'native_source_epsg':meta['epsg'],
            'source_file_sha256':sha256(source),
            'source_file_bytes':int(Path(source).stat().st_size),
            'original_bbox':[float(lon.min()),float(lat.min()),float(lon.max()),float(lat.max())],
            'native_etopo_window':{
                'rows':rows,'cols':cols,'cells':rows*cols,
                'valid_elevation_cells':int(valid.sum()),
                'modern_land_cells':int(land.sum()),
                'current_submerged_cells':int((valid&~land).sum()),
                'newly_emerged_under_static_minus120m_cells':int(additional.sum()),
                'static_minus120m_land_cells':int(lowered.sum()),
                'modern_land_fraction_of_valid':float(land.sum()/valid.sum()),
                'static_minus120m_land_fraction_of_valid':float(lowered.sum()/valid.sum()),
                'row_offset':r0,'col_offset':c0,
            },
            'endpoint_sampling':{
                'n_original_localities':len(e),
                'n_finite_elevations':int(np.isfinite(e).sum()),
                'n_land_above_present_sea_level':int((e>0).sum()),
                'n_land_under_static_minus120m_sensitivity':int((e>SEA_LEVEL).sum()),
                'minimum_elevation_m':float(e.min()),
                'maximum_elevation_m':float(e.max())
            },
            'scientific_boundary':{
                'not_a_verified_paleo_coastline':True,
                'ignores_GIA_tectonics_eustatic_history_ice':True,
                'ETOPO_present_topography_not_LGM_reconstruction':True,
                'not_resampled_to_original_CHELSA_grid':True,
                'source_not_connected_to_original_Dijkstra_raster_yet':True,
                'original_climate_lcp_unchanged':True,
                'COI_genetic_effect_used':False,
                'original_259_primary':'NOT_EVALUABLE_EMPIRICAL_RESPONSE',
                'no_new_p_value_or_genetic_retest':True,
            }
        }


def main():
    a=argparse.ArgumentParser()
    for name in ('contract','source_pin','etopo_15s_surface_tif','frozen_original_localities','output'):
        a.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    args=a.parse_args()
    if sha256(args.contract)!=CONTRACT_SHA:
        raise FixedSourceError('corrected frozen v0.2 ETOPO contract SHA256 drift')
    c=json.loads(args.contract.read_text())
    if (c.get('schema')!='ttf_historical_host_lcp_etopo_static_mask_contract_v0.2'
        or c['comparison']['fixed_threshold_m']!=SEA_LEVEL
        or c['original_geometry']['fixed_bbox_buffer_each_side_deg']!=PADDING
        or c['original_geometry']['species']!=SPECIES):
        raise FixedSourceError('corrected ETOPO scientific contract changed')
    pin=json.loads(args.source_pin.read_text())
    metadata=verify_source_pin(args.etopo_15s_surface_tif,pin)
    if sha256(args.frozen_original_localities)!=GEOMETRY_SHA:
        raise FixedSourceError('original exact response-blind locality geometry mismatch')
    with args.frozen_original_localities.open(newline='',encoding='utf8') as f:
        rows=[r for r in csv.DictReader(f) if r['species']==SPECIES]
    if len(rows)!=26 or sorted(int(r['locality_index']) for r in rows)!=list(range(26)):
        raise FixedSourceError('original selected sample localities changed')
    rows.sort(key=lambda x:int(x['locality_index']))
    sites=[(float(r['longitude']),float(r['latitude'])) for r in rows]
    extent=(min(l for l,_ in sites),max(l for l,_ in sites),min(a for _,a in sites),max(a for _,a in sites))
    if extent!=(-85.51,-85.2804,10.455,11.05):
        raise FixedSourceError('original frozen Bungalotis geographic extent changed')
    result=audit_native_topography(args.etopo_15s_surface_tif,sites)
    result['source_pin_sha256']=sha256(args.source_pin)
    result['frozen_contract_sha256']=CONTRACT_SHA
    result['original_locality_geometry_sha256']=GEOMETRY_SHA
    result['source_metadata_verified_before_topography_read']=True
    result['original_source_metadata']=metadata
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'source_file_sha256':SOURCE_SHA,
        'native_etopo_window':result['native_etopo_window'],
        'endpoint_sampling':result['endpoint_sampling']},sort_keys=True))

if __name__=='__main__': main()
