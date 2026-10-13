#!/usr/bin/env python3
"""Source-backed ETOPO land-mask *readiness*, not a LGM coastline model.

Modern topography is unchanged; an illustrative uniform -120 m eustatic sea
level is not a glacio-isostatically adjusted regional palaeoshoreline. The
original 64-edge climate LCP result and 259-species genetic terminal state
are never recomputed or altered by this independent geography-only screen.
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

ORIGINAL_GEO_SHA = '037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9'
SPECIES = 'Bungalotis erythus'
SOURCE_URL = ('https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2022/data/15s/'
              '15s_surface_elev_gtif/ETOPO_2022_v1_15s_N00W090_surface.tif')
SEA_LEVEL_THRESHOLD_M = -120
PADDING_DEG = 0.12
FROZEN_CONTRACT_SHA256 = '026432857596123d1628f81be969cfb3d6bc2f5f8193a17f3916c22f841edb8a'


class FixedLcpGeometryError(RuntimeError):
    pass


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda:f.read(1<<20),b''):
            h.update(data)
    return h.hexdigest()


def verify_contract(path:Path)->bool:
    if sha256(path)!=FROZEN_CONTRACT_SHA256:
        raise RuntimeError('frozen ETOPO mask contract SHA256 mismatch')
    c=json.loads(path.read_text())
    if (c['status']!='FROZEN_BEFORE_FIRST_ETOPO_MASK_OBSERVATION_POST_CLIMATE_LCP_EXPOSURE_NOT_GENETIC_PREREGISTRATION'
        or c['original_geometry']['species']!=SPECIES
        or c['comparison']['fixed_sensitivity_sea_level_m']!=SEA_LEVEL_THRESHOLD_M
        or c['original_geometry']['fixed_bbox_buffer_each_side_deg']!=PADDING_DEG):
        raise RuntimeError('frozen ETOPO mask contract scientific scope changed')
    return True


def audit_topography_mask(path:Path, lonlat:list[tuple[float,float]], padding:float=PADDING_DEG,
                           *,threshold_m:int=SEA_LEVEL_THRESHOLD_M)->dict:
    coords=np.asarray(lonlat,dtype=float)
    if (coords.ndim!=2 or coords.shape[1]!=2 or len(coords)<2 or not np.isfinite(coords).all() or
        padding!=PADDING_DEG or threshold_m!=SEA_LEVEL_THRESHOLD_M):
        raise ValueError('fixed original ROI padding/sea-level sensitivity drift')
    with rasterio.open(path) as ds:
        if not ds.crs or ds.crs.to_epsg()!=4326 or ds.count!=1:
            raise FixedLcpGeometryError('ETOPO source not one-band geographic EPSG:4326')
        t=ds.transform
        if t.b!=0 or t.d!=0 or t.a<=0 or t.e>=0:
            raise FixedLcpGeometryError('ETOPO raster not north-up axis-aligned')
        lon=coords[:,0];lat=coords[:,1]
        if np.any(lon<-180) or np.any(lon>180) or np.any(lat<-90) or np.any(lat>90):
            raise ValueError('geographic coordinates outside globe')
        bb=[float(lon.min()-padding),float(lat.min()-padding),float(lon.max()+padding),float(lat.max()+padding)]
        f=from_bounds(*bb,transform=t)
        r0=int(math.floor(f.row_off));c0=int(math.floor(f.col_off))
        r1=int(math.ceil(f.row_off+f.height));c1=int(math.ceil(f.col_off+f.width))
        if r0<0 or c0<0 or r1>ds.height or c1>ds.width:
            raise FixedLcpGeometryError('ROI outside fixed ETOPO tile; no substitution')
        nr=r1-r0;nc=c1-c0
        if nr*nc<=0 or nr*nc>150_000:
            raise FixedLcpGeometryError('ETOPO grid cap or empty ROI')
        arr=ds.read(1,window=Window(c0,r0,nc,nr),masked=True)
        elev=np.asarray(np.ma.filled(arr.astype('float64'),np.nan),dtype=float)
        good=np.isfinite(elev)
        if not np.any(good):
            raise FixedLcpGeometryError('ETOPO ROI has no valid elevation')
        samples=list(ds.sample([tuple(x) for x in coords],masked=True))
        endpoint=np.asarray([float(x[0]) if not np.ma.is_masked(x[0]) else np.nan for x in samples])
        if not np.isfinite(endpoint).all():
            raise FixedLcpGeometryError('Original locality not represented by fixed ETOPO source')
        current=good & (elev>0)
        illustrative=good & (elev>SEA_LEVEL_THRESHOLD_M)
        newly=good & (elev>SEA_LEVEL_THRESHOLD_M) & (elev<=0)
        return {
            'schema':'ttf_etopo_static_eustatic_mask_readiness_v0.1',
            'status':'EXPLORATORY_MODERN_BATHYMETRY_MASK_ONLY_NOT_VALIDATED_LGM_COASTLINE',
            'source_dataset':'NOAA ETOPO 2022 v1, 15-arcsec ice-surface elevation GeoTIFF (tropical ROI)',
            'source_file_sha256':sha256(path),
            'source_file_bytes':path.stat().st_size,
            'source_url':SOURCE_URL,
            'geographic_roi_buffer_deg':padding,
            'illustrative_static_eustatic_sea_level_m':SEA_LEVEL_THRESHOLD_M,
            'original_geometry':{'species':SPECIES,'source_sha256':ORIGINAL_GEO_SHA},
            'original_bbox':[float(lon.min()),float(lat.min()),float(lon.max()),float(lat.max())],
            'native_etopo_window':{
                'rows':nr,'cols':nc,'cells':nr*nc,'valid_elevation_cells':int(good.sum()),
                'modern_land_cells':int(current.sum()),
                'current_submerged_cells':int((good&~current).sum()),
                'newly_emerged_under_static_minus120m_cells':int(newly.sum()),
                'static_minus120m_land_cells':int(illustrative.sum()),
                'modern_land_fraction_of_valid':float(current.sum()/good.sum()),
                'static_minus120m_land_fraction_of_valid':float(illustrative.sum()/good.sum())
            },
            'endpoint_sampling':{
                'n_original_localities':len(coords),
                'n_finite_elevations':int(np.isfinite(endpoint).sum()),
                'n_land_above_present_sea_level':int((endpoint>0).sum()),
                'n_land_under_static_minus120m_sensitivity':int((endpoint>SEA_LEVEL_THRESHOLD_M).sum()),
                'minimum_elevation_m':float(np.min(endpoint)),
                'maximum_elevation_m':float(np.max(endpoint)),
            },
            'scientific_boundary':{
                'not_a_verified_paleo_coastline':True,
                'ignores_local_GIA_and_tectonics_and_paleoshoreline_history':True,
                'ETOPO_present_topography_not_a_LGM_reconstruction':True,
                'not_resampled_to_original_CHELSA_101x57_pixel_grid':True,
                'source_not_connected_to_original_Dijkstra_raster_yet':True,
                'original_climate_lcp_unchanged':True,
                'COI_genetic_effect_used':False,
                'original_259_primary':'NOT_EVALUABLE_EMPIRICAL_RESPONSE',
                'no_new_p_value_or_genetic_retest':True,
                'source_licence_and_citation_must_be_respected':True
            }
        }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--etopo-15s-surface-tif',type=Path,required=True)
    ap.add_argument('--frozen-original-localities',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    verify_contract(a.contract)
    if sha256(a.frozen_original_localities)!=ORIGINAL_GEO_SHA:
        raise RuntimeError('original frozen response-blind geometry digest drift')
    with a.frozen_original_localities.open(newline='',encoding='utf8') as f:
        s=[r for r in csv.DictReader(f) if r['species']==SPECIES]
    if len(s)!=26 or sorted(int(r['locality_index']) for r in s)!=list(range(26)):
        raise RuntimeError('original selected 26 localities changed')
    s.sort(key=lambda x:int(x['locality_index']))
    sites=[(float(r['longitude']),float(r['latitude'])) for r in s]
    if (min(l for l,_ in sites),max(l for l,_ in sites),min(a for _,a in sites),max(a for _,a in sites))!=(-85.51,-85.2804,10.455,11.05):
        raise RuntimeError('exact original Bungalotis geographic extent drift')
    result=audit_topography_mask(a.etopo_15s_surface_tif,sites)
    result['frozen_contract_sha256']=FROZEN_CONTRACT_SHA256
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':result['status'],'native_etopo_window':result['native_etopo_window'],
                      'endpoint_sampling':result['endpoint_sampling']},sort_keys=True))


if __name__=='__main__': main()
