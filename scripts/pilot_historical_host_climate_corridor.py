#!/usr/bin/env python3
"""CHELSA climate-only 4-neighbor LCP pilot; NEVER opens genetic responses.

Exploratory: a modern host-climate analogue projection, NOT true ancient plant
occupancy and NOT a palaeocoastline/ice-sheet movement mask. The fixed selected
original insect's 26 localities/64 edges are only geography, not new genetic data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

SPECIES = 'Bungalotis erythus'
BASE_SOURCE_SHA={
    'candidates':'c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b',
    'localities':'037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9',
    'edges':'ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9',
    'host_pairs':'0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9',
    'native_units':'c731906315f7452f83ce302c1d36c75ad94afc24ce7e24d461360312244ac558',
    'wgsrpd_support':'d0fc12f635ec56a442dd37dec06cf3aa06684a00b14f1bb1f2db2f6fe3c4a21f',
    'frozen_predictors':'70b748d588d9ab10309ace080f85b8053c6e14aa96d7a84288724c1000272a5b',
    'manifest':'7d82833996006fc1e99b86acb4ea7e80cde29c7f500df739ef48cb9ce145d450',
}
TIMES={'0BP':20,'LGM':-190}
VARS=('bio01','bio07','bio12','bio15')


def verify_pilot_source_bindings(rule):
    """Validate exact bytes through their *published* field names.

    This function repairs only a 2026-10-13 field-name mapping error,
    not the pilot's scientific parameters, frozen sources, or estimand.
    """
    declared=rule.get('source_sha256')
    if not isinstance(declared,dict):
        raise RuntimeError('frozen source binding dictionary missing')
    mapping={
        'candidates':'candidates',
        'localities':'localities',
        'edges':'edges',
        'host_pairs':'host_pairs',
        'native_units':'native_units',
        'wgsrpd_support':'wgsrpd_support',
        'frozen_predictors':'frozen_edge_predictors',
        'manifest':'raster_manifest',
    }
    if any(declared.get(public_key)!=BASE_SOURCE_SHA[internal_key]
           for internal_key,public_key in mapping.items()):
        raise RuntimeError('frozen pilot source binding drift')
    if declared.get('original_eight_raster_staged_zip')!='0a797c609d26a185ba08d209a71203a584fa60f1045c4a8c335934372cee30b2':
        raise RuntimeError('frozen pilot external CHELSA source binding drift')
    return True


class ClimateCorridorNotEvaluable(RuntimeError):
    """Terminal failure for this fixed exploratory climate-only pilot."""


def sha256_file(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for piece in iter(lambda:handle.read(1<<20),b''):
            h.update(piece)
    return h.hexdigest()


def load_csv(path:Path):
    with path.open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))


def haversine_km(lon1,lat1,lon2,lat2):
    p1,p2=np.deg2rad([lat1,lat2]);dphi=p2-p1; dl=np.deg2rad(lon2-lon1)
    a=np.sin(dphi/2)**2+np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    return float(6371.0088*2*np.arctan2(np.sqrt(a),np.sqrt(max(0.,1.-a))))


def native_bounding_window(transform,lon,lat,padding,*,shape,maximum_cells):
    """Native pixel alignment with an outward-rounded fixed 0.12-degree ROI."""
    from rasterio.windows import from_bounds, Window
    lo=np.asarray(lon,dtype=float);la=np.asarray(lat,dtype=float)
    if (lo.ndim!=1 or la.shape!=lo.shape or len(lo)<2
            or not np.all(np.isfinite(lo)) or not np.all(np.isfinite(la))
            or padding!=0.12 or maximum_cells<1):
        raise ValueError('invalid fixed ROI coordinates, padding or cell cap')
    if transform.b!=0 or transform.d!=0 or transform.a<=0 or transform.e>=0:
        raise ClimateCorridorNotEvaluable('non-north-up CHELSA geographic raster')
    full=from_bounds(float(lo.min()-padding),float(la.min()-padding),float(lo.max()+padding),float(la.max()+padding),transform)
    c0=int(math.floor(full.col_off));r0=int(math.floor(full.row_off))
    c1=int(math.ceil(full.col_off+full.width));r1=int(math.ceil(full.row_off+full.height))
    if c0<0 or r0<0 or c1>shape[1] or r1>shape[0]:
        raise ClimateCorridorNotEvaluable('ROI outside original CHELSA native grid')
    w=Window(c0,r0,c1-c0,r1-r0)
    if w.width*w.height>maximum_cells:
        raise ClimateCorridorNotEvaluable('native ROI cell cap exceeded; no adaptive downsampling')
    return w


def analogue_cell_cost(climate,cloud):
    """1 + native frozen nearest-four-var Euclidean climate-analogue distance."""
    values=np.asarray(climate,dtype=float)
    reference=np.asarray(cloud,dtype=float)
    if values.ndim!=3 or reference.ndim!=2 or values.shape[2]!=reference.shape[1] or len(reference)==0 or not np.isfinite(reference).all():
        raise ValueError('climate and reference cloud dimension/finite mismatch')
    result=np.full(values.shape[:2],np.inf,dtype=float)
    valid=np.all(np.isfinite(values),axis=2)
    if np.any(valid):
        distance,_=cKDTree(reference).query(values[valid],k=1,workers=1)
        result[valid]=1.0+np.asarray(distance,dtype=float)
    if np.any(result[np.isfinite(result)]<=0):
        raise RuntimeError('unexpected nonpositive analogue cost')
    return result


def _make_grid_graph(cost, transform):
    z=np.asarray(cost,dtype=float)
    if z.ndim!=2 or min(z.shape)<2 or np.any(np.isnan(z)) or np.any(z[np.isfinite(z)]<=0):
        raise ValueError('raster must be 2D with positive finite or inf-blocked costs')
    if not (transform.b==0 and transform.d==0 and transform.a>0 and transform.e<0):
        raise ValueError('rotated or non-geographic grid')
    nrow,ncol=z.shape
    valid=np.isfinite(z)
    # Native E/W lengths depend on each row's geodetic latitude.
    from rasterio.transform import xy
    rows=np.arange(nrow)
    lat=np.asarray([xy(transform,int(i),0,offset='center')[1] for i in rows])
    lon0=xy(transform,0,0,offset='center')[0]
    lon1=xy(transform,0,1,offset='center')[0]
    horizontal=np.asarray([haversine_km(lon0,l,lon1,l) for l in lat],dtype=float)
    vertical=np.asarray([haversine_km(lon0,float(lat[i]),lon0,float(lat[i+1])) for i in range(nrow-1)],dtype=float)
    if np.any(horizontal<=0) or np.any(vertical<=0):
        raise ValueError('invalid native CHELSA cell-center kilometer distances')
    hh=valid[:,:-1] & valid[:,1:]
    ri,ci=np.where(hh)
    src_h=ri*ncol+ci; dst_h=src_h+1
    w_h=horizontal[ri]*(z[ri,ci]+z[ri,ci+1])/2
    vv=valid[:-1,:] & valid[1:,:]
    rv,cv=np.where(vv)
    src_v=rv*ncol+cv;dst_v=src_v+ncol
    w_v=vertical[rv]*(z[rv,cv]+z[rv+1,cv])/2
    start=np.r_[src_h,src_v];end=np.r_[dst_h,dst_v];weight=np.r_[w_h,w_v]
    if not np.all(np.isfinite(weight)) or np.any(weight<=0):
        raise RuntimeError('invalid positive graph edge weight')
    u=np.r_[start,end];v=np.r_[end,start];w=np.r_[weight,weight]
    graph=csr_matrix((w,(u,v)),shape=(nrow*ncol,nrow*ncol))
    return graph


def cost_distance_between_pixels(cost, transform, pixel_nodes, edge_nodes,*,forbid_self_edges=False):
    """Unmodified native raster 4-neighbor shortest path for exact original edges."""
    cost=np.asarray(cost,dtype=float)
    if cost.ndim!=2 or np.any(np.isnan(cost)) or np.any(cost[np.isfinite(cost)]<=0):
        raise ValueError('all traversable cells must have strictly positive costs')
    nodes=np.asarray(pixel_nodes)
    pairs=np.asarray(edge_nodes)
    if nodes.ndim!=2 or nodes.shape[1]!=2 or not np.issubdtype(nodes.dtype,np.integer):
        raise ValueError('pixel nodes need integer row/col')
    if pairs.ndim!=2 or pairs.shape[1]!=2 or not np.issubdtype(pairs.dtype,np.integer):
        raise ValueError('edge nodes must index exact frozen sample points')
    if (np.any(nodes[:,0]<0) or np.any(nodes[:,0]>=cost.shape[0]) or
        np.any(nodes[:,1]<0) or np.any(nodes[:,1]>=cost.shape[1]) or
        np.any(pairs<0) or np.any(pairs>=len(nodes))):
        raise ValueError('frozen endpoint or edge index outside native raster')
    if forbid_self_edges and np.any(pairs[:,0]==pairs[:,1]):
        raise ValueError('frozen self edges disallowed')
    if np.any(~np.isfinite(cost[nodes[:,0],nodes[:,1]])):
        raise ClimateCorridorNotEvaluable('native raster endpoint is invalid, no snapping')
    flat=nodes[:,0]*cost.shape[1]+nodes[:,1]
    graph=_make_grid_graph(cost,transform)
    unique,reverse=np.unique(flat,return_inverse=True)
    dis=dijkstra(graph,directed=False,indices=unique,return_predecessors=False)
    # Intermediates with no active neighbors will have only 0 at self.
    result=dis[reverse[pairs[:,0]],flat[pairs[:,1]]]
    if np.any(~np.isfinite(result)):
        raise ClimateCorridorNotEvaluable('unreachable frozen original edge under climate-only native grid')
    return np.asarray(result,dtype=float)


def check_assets_and_read_native_window(assets,lon,lat,padding=0.12,max_cells=15000):
    from rasterio.windows import transform as win_transform
    with ExitStack() as stack:
        sources={}
        for time in TIMES:
            for v in VARS:
                path=assets[(v,TIMES[time])]
                sources[(time,v)]=stack.enter_context(__import__('rasterio').open(path))
        first=sources[('0BP','bio01')]
        if not first.crs or first.crs.to_epsg()!=4326:
            raise ClimateCorridorNotEvaluable('CHELSA source is not EPSG:4326')
        w=native_bounding_window(first.transform,lon,lat,padding,shape=first.shape,maximum_cells=max_cells)
        reference=(first.width,first.height,tuple(first.transform),str(first.crs))
        raw={}
        for time in TIMES:
            cols=[]
            for v in VARS:
                ds=sources[(time,v)]
                if (ds.width,ds.height,tuple(ds.transform),str(ds.crs))!=reference:
                    raise ClimateCorridorNotEvaluable('CHELSA variable or time has different native grid')
                z=ds.read(1,window=w).astype(float)
                if ds.nodata is not None:
                    z[z==float(ds.nodata)]=np.nan
                cols.append(z)
            raw[time]=np.stack(cols,axis=-1)
        t=win_transform(w,first.transform)
        return raw,t,{'rows':int(w.height),'cols':int(w.width),
                      'row_offset':int(w.row_off),'col_offset':int(w.col_off),
                      'native_pixel_lon_degrees':float(first.transform.a),
                      'native_pixel_lat_degrees':float(abs(first.transform.e))}


def main():
    import rasterio
    from rasterio.transform import rowcol
    from ttf.historical_host_memory import great_circle_quadrature_latlon, historical_analog_memory
    # The original response-blind builder handles the exact source-native conventions.
    try:
        from scripts.build_historical_host_memory_predictor import sample_matrix, norm_id, load_assets
    except ModuleNotFoundError:
        from build_historical_host_memory_predictor import sample_matrix,norm_id,load_assets
    p=argparse.ArgumentParser()
    for name in ('contract','candidates','localities','edges','host_pairs','native_units','wgsrpd_support','frozen_predictors','manifest'):
        p.add_argument('--'+name.replace('_','-'),required=True,type=Path)
    p.add_argument('--historical-asset',type=Path,action='append',required=True)
    p.add_argument('--output-summary',type=Path,required=True)
    p.add_argument('--output-edges',type=Path,required=True)
    a=p.parse_args()
    rule=json.loads(a.contract.read_text())
    if rule['status']!='FIXED_PILOT_RULES_BEFORE_NEW_CORRIDOR_OUTCOME_NOT_GENETIC_PREREGISTRATION' or rule['one_fixed_pilot_species']!=SPECIES:
        raise RuntimeError('exploratory fixed pilot contract status/species drift')
    if sha256_file(a.contract)!='a46e8c61bbc07e3b093c4664b12825c7b193455be99d457d8d96343a358d5ef7':
        raise RuntimeError('frozen pilot cost/ROI/algorithm contract bytes changed')
    verify_pilot_source_bindings(rule)
    for name,expected in BASE_SOURCE_SHA.items():
        pth=getattr(a,name)
        if sha256_file(pth)!=expected:
            raise RuntimeError('original response-blind source SHA mismatch: '+name)
    if len(a.historical_asset)!=8:
        raise RuntimeError('exact eight CHELSA climate variables are required')
    assets=load_assets(a.historical_asset)
    original=load_csv(a.candidates)
    r=[x for x in original if x['species']==SPECIES]
    if len(r)!=1 or int(r[0]['n_localities'])!=26 or int(r[0]['edges'])!=64:
        raise RuntimeError('frozen selected species geo identity drift')
    candidate=r[0]
    local=sorted([x for x in load_csv(a.localities) if x['species']==SPECIES],key=lambda x:int(x['locality_index']))
    edges=sorted([x for x in load_csv(a.edges) if x['species']==SPECIES],key=lambda x:int(x['edge_index']))
    if len(local)!=26 or len(edges)!=64 or [int(x['locality_index']) for x in local]!=list(range(26)) or [int(x['edge_index']) for x in edges]!=list(range(64)):
        raise RuntimeError('frozen pilot locality/edge index drift')
    pair=np.asarray([[int(x['node_left']),int(x['node_right'])] for x in edges],dtype=int)
    if np.any(pair<0) or np.any(pair>=26) or np.any(pair[:,0]==pair[:,1]) or len(set(map(tuple,np.sort(pair,axis=1))))!=64:
        raise RuntimeError('invalid original fixed graph topology')
    lat=np.array([float(x['latitude']) for x in local]);lon=np.array([float(x['longitude']) for x in local]);xy=list(zip(lon,lat))
    # Crosscheck exact insect-accepted host associations from sidecar before using host ranges.
    selected_ids={norm_id(q) for q in candidate['accepted_host_ids'].split(';')}
    actual_ids={norm_id(x['accepted_plant_name_id']) for x in load_csv(a.host_pairs) if x['insect_species']==SPECIES}
    if selected_ids!=actual_ids or len(actual_ids)!=3:
        raise RuntimeError('selected host plant identity drift')
    native=defaultdict(set)
    for rec in load_csv(a.native_units):
        h=norm_id(rec['accepted_plant_name_id'])
        if h in selected_ids:
            native[h].add(rec['area_code_l3'].strip())
    units=set().union(*(native[h] for h in selected_ids))
    support=load_csv(a.wgsrpd_support)
    coords_support=[(float(z['longitude']),float(z['latitude'])) for z in support]
    present_paths=[assets[(v,20)] for v in VARS]
    lgm_paths=[assets[(v,-190)] for v in VARS]
    orig_support=sample_matrix(coords_support,present_paths)
    finite=np.all(np.isfinite(orig_support),axis=1)
    if np.count_nonzero(finite)<100:
        raise ClimateCorridorNotEvaluable('global frozen climate support too small')
    mean=orig_support[finite].mean(axis=0);sd=orig_support[finite].std(axis=0,ddof=0)
    if np.any(sd<=0) or not np.all(np.isfinite(sd)):
        raise RuntimeError('invalid global response-blind climate standardization')
    stand=lambda x:(x-mean)/sd
    host_cloud=np.stack([stand(orig_support[j]) for j,z in enumerate(support) if finite[j] and z['area_code_l3'].strip() in units])
    if len(host_cloud)<10:
        raise ClimateCorridorNotEvaluable('host cloud inadequate')
    loc_climate=sample_matrix(xy,present_paths)
    validloc=np.all(np.isfinite(loc_climate),axis=1)
    self_cloud=stand(loc_climate[validloc])
    if len(self_cloud)<10:
        raise ClimateCorridorNotEvaluable('self climate cloud inadequate')
    raw,transform,shape=check_assets_and_read_native_window(assets,lon,lat)
    scenario={}
    for time in TIMES:
        z=raw[time]
        values=np.full(z.shape,np.nan,dtype=float)
        idx=np.all(np.isfinite(z),axis=-1)
        values[idx]=stand(z[idx]);
        for label,cloud in (('host',host_cloud),('self',self_cloud)):
            scenario[(label,time)]=analogue_cell_cost(values,cloud)
    pixels=np.asarray([rowcol(transform,float(lon[i]),float(lat[i])) for i in range(26)],dtype=int)
    result={k:cost_distance_between_pixels(cost,transform,pixels,pair,forbid_self_edges=True) for k,cost in scenario.items()}
    # Independently reconstruct the original straight-edge 5-point climate predictor.
    quad=great_circle_quadrature_latlon(np.column_stack([lat[pair[:,0]],lon[pair[:,0]]]),
                                      np.column_stack([lat[pair[:,1]],lon[pair[:,1]]]),segment_points=5)
    sample_xy=[(float(east),float(north)) for north,east in quad.reshape(-1,2)]
    climate_cur=stand(sample_matrix(sample_xy,present_paths)).reshape(64,5,4)
    climate_lgm=stand(sample_matrix(sample_xy,lgm_paths)).reshape(64,5,4)
    mh,dhcur,dhlgm=historical_analog_memory(climate_cur,climate_lgm,host_cloud)
    ms,dscur,dslgm=historical_analog_memory(climate_cur,climate_lgm,self_cloud)
    frozen=sorted([z for z in load_csv(a.frozen_predictors) if z['species']==SPECIES],key=lambda x:int(x['edge_index']))
    if len(frozen)!=64 or [int(q['edge_index']) for q in frozen]!=list(range(64)):
        raise RuntimeError('frozen species original 5-point predictor edge identity drift')
    err={}
    for k,arr in [('M_host',mh),('D_host_0BP',dhcur),('D_host_LGM',dhlgm),('M_self',ms),('D_self_0BP',dscur),('D_self_LGM',dslgm)]:
        err[k]=float(np.max(np.abs(arr-np.array([float(z[k]) for z in frozen]))))
    if any(v>1e-8 for v in err.values()):
        raise ClimateCorridorNotEvaluable('frozen original five-point climate predictor mismatch')
    h=result[('host','LGM')]-result[('host','0BP')]
    s=result[('self','LGM')]-result[('self','0BP')]
    d=h-s
    def statistics(v):
        return {'min':float(np.min(v)),'median':float(np.median(v)),
                'max':float(np.max(v)),'mean':float(np.mean(v))}
    zeros=int(np.sum(np.all(pixels[pair[:,0]]==pixels[pair[:,1]],axis=1)))
    output={
        'schema':'ttf_bungalotis_host_climate_corridor_pilot_v0.1',
        'status':'PASS_CHELSA_CLIMATE_ONLY_PILOT_NO_LAND_ICE_MASK',
        'source_sha256':BASE_SOURCE_SHA,
        'contract_sha256':sha256_file(a.contract),
        'species':SPECIES,'localities':26,'original_edges':64,
        'host_accepted_ids':sorted(selected_ids),'host_0bp_climate_cloud_points':len(host_cloud),
        'self_0bp_climate_cloud_points':len(self_cloud),
        'reference_climate_points':int(np.count_nonzero(finite)),
        'native_raster':shape,
        'joint_valid_cells':int(sum(np.all([np.isfinite(v) for v in raw.values()],axis=0).all(axis=-1).ravel())),
        'same_cell_original_edges':zeros,
        'frozen_5point_predictor_max_abs_error':err,
        'cost_distance_km':{f'{lab}_{t}':statistics(v) for (lab,t),v in result.items()},
        'change_LGM_minus_present_cost_km':{'host':statistics(h),'insect_self':statistics(s),'host_minus_insect_self':statistics(d)},
        'fraction_host_change_positive':float(np.mean(h>0)),
        'fraction_unique_change_positive':float(np.mean(d>0)),
        'ecological_boundary':{
           'COI_genetic_response_used':False,
           'fivepoint_straight_climate_average_is_not_two_dimensional_connectivity':True,
           'modern_host_climate_analogue_is_not_ancient_host_occupancy':True,
           'no_paleo_land_or_ice_mask_accepted':True,
           'finite_climate_pixels_not_verified_land':True,
           'ROI_boundary_can_artificially_limit_detours':True,
           'pilot_not_confirmatory_or_source_disjoint':True,
           'original_259_empirical_result':'NOT_EVALUABLE_EMPIRICAL_RESPONSE',
           'no_genetic_effect_p_value':True
        }
    }
    columns=['species','edge_index','node_left','node_right','geographic_chord_km','R_host_0BP_km','R_host_LGM_km',
             'R_self_0BP_km','R_self_LGM_km','M_host_lcp_km','M_self_lcp_km','M_host_minus_self_lcp_km']
    a.output_edges.parent.mkdir(parents=True,exist_ok=True);a.output_summary.parent.mkdir(parents=True,exist_ok=True)
    with a.output_edges.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,lineterminator='\n');writer.writeheader()
        xyz=np.asarray([[float(z[x]) for x in ('x_km','y_km','z_km')] for z in local])
        for i,e in enumerate(edges):
            values=[SPECIES,i,int(pair[i,0]),int(pair[i,1]),np.linalg.norm(xyz[pair[i,0]]-xyz[pair[i,1]]),
                    result[('host','0BP')][i],result[('host','LGM')][i],result[('self','0BP')][i],result[('self','LGM')][i],h[i],s[i],d[i]]
            writer.writerow(dict(zip(columns,[v if isinstance(v,(str,int,np.integer)) else repr(float(v)) for v in values])))
    output['output_edges_sha256']=sha256_file(a.output_edges)
    a.output_summary.write_text(json.dumps(output,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':output['status'],'lcp_edge_n':64,'summary':output['change_LGM_minus_present_cost_km'],'replay_max_err':max(err.values())},sort_keys=True))

if __name__=='__main__':main()
