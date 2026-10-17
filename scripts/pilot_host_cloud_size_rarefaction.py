#!/usr/bin/env python3
"""Response-free CHELSA LCP sensitivity to 665 versus 26 host-reference points.

All 26 original Bungalotis sampling locations, 64 frozen edges, 101x57
native CHELSA ROI, modern global reference scale and positive four-neighbor
km distance cost are unchanged from the *already evaluated* exploratory pilot.
After binding original host/self LCP edge output, exactly 64 host-cloud
rarefactions of 26 out of the original 665 reference rows are evaluated.

This is NOT a historical host distribution reconstruction, new response test,
or equivalently sampled host/self ecological domain. Never reads genotypes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

SPECIES='Bungalotis erythus'
HOST_ORIGINAL_N=665
RAREFIED_N=26
REPLICATES=64
SEED_NAME=(
    'hhm-host-26-rarefaction-v0.1|Bungalotis erythus|'
    '7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b'
)
RNG_SEED=int(hashlib.sha256(SEED_NAME.encode()).hexdigest()[:16],16)
BASE={
    'candidates':'c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b',
    'localities':'037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9',
    'edges':'ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9',
    'host_pairs':'0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9',
    'native_units':'c731906315f7452f83ce302c1d36c75ad94afc24ce7e24d461360312244ac558',
    'wgsrpd_support':'d0fc12f635ec56a442dd37dec06cf3aa06684a00b14f1bb1f2db2f6fe3c4a21f',
    'frozen_predictors':'70b748d588d9ab10309ace080f85b8053c6e14aa96d7a84288724c1000272a5b',
    'manifest':'7d82833996006fc1e99b86acb4ea7e80cde29c7f500df739ef48cb9ce145d450',
    'original_pilot_edges':'7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b',
    'original_pilot_summary':'6969f482a1967db2db159cbc5fd62b722268a0fc45d7dcd1b979f1973b9a3493',
}


def sha256_path(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):
            h.update(chunk)
    return h.hexdigest()


def rows_from_csv(path:Path):
    with path.open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))


def fixed_subsample_indices(n_host:int=665,n_match:int=26,n_replicates:int=64)->np.ndarray:
    if n_host!=HOST_ORIGINAL_N:
        raise ValueError('expected the exact original 665-host-point frozen cloud')
    if n_match!=RAREFIED_N or n_host<n_match:
        raise ValueError('rarefied size must equal 26 fixed self reference points')
    if n_replicates!=REPLICATES:
        raise ValueError('frozen 64 source-blind rarefaction replicates required')
    rng=np.random.default_rng(RNG_SEED)
    return np.stack([rng.choice(n_host,size=n_match,replace=False) for _ in range(n_replicates)])


def verify_original_edges(full:dict,original:list[dict],atol:float=1e-8)->float:
    columns={
        ('host','0BP'):'R_host_0BP_km',
        ('host','LGM'):'R_host_LGM_km',
        ('self','0BP'):'R_self_0BP_km',
        ('self','LGM'):'R_self_LGM_km',
    }
    if not original or any(int(r['edge_index'])!=i for i,r in enumerate(original)):
        raise RuntimeError('original frozen baseline edge indices are missing or out of order')
    errors=[]
    for k,field in columns.items():
        computed=np.asarray(full[k],dtype=float)
        observed=np.asarray([float(r[field]) for r in original])
        if computed.shape!=observed.shape or not np.isfinite(computed).all():
            raise RuntimeError('original baseline dimensions or values are nonfinite')
        errors.append(float(np.max(np.abs(computed-observed))))
    maximum=max(errors)
    if maximum>atol:
        raise RuntimeError('baseline replay differs from frozen source 64-edge LCP costs')
    return maximum


def summarize_rarefied(host_delta_by_draw:np.ndarray,original_self_delta:np.ndarray)->dict:
    a=np.asarray(host_delta_by_draw,dtype=float)
    b=np.asarray(original_self_delta,dtype=float)
    if a.ndim!=2 or b.ndim!=1 or a.shape[1]!=len(b):
        raise ValueError('shape mismatch in rarefaction draw/edge matrix')
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('nonfinite rarefied cost distances')
    if len(b)==0 or len(a)==0:
        raise ValueError('empty rarefaction matrix')
    hmean=np.mean(a,axis=1)
    diff=hmean-float(b.mean())
    def describe(v,*,positive=False):
        result={
            'means_per_draw':[float(x) for x in v],
            'mean':float(np.mean(v)),
            'min':float(np.min(v)),
            'max':float(np.max(v)),
            'q05':float(np.quantile(v,.05)),
            'q25':float(np.quantile(v,.25)),
            'median':float(np.median(v)),
            'q75':float(np.quantile(v,.75)),
            'q95':float(np.quantile(v,.95)),
        }
        if positive:
            result['positive_draws']=int(np.sum(v>0))
            result['fraction_positive_draws']=float(np.mean(v>0))
        return result
    return {
        'draw_n':int(len(a)),
        'original_edge_n':int(len(b)),
        'fixed_self_mean_delta_cost_km':float(b.mean()),
        'rarefied_host_delta_cost_km':describe(hmean),
        'rarefied_host_minus_self_delta_cost_km':describe(diff,positive=True),
    }


def verify_fixed_sources(a,rule):
    if (rule.get('schema')!='ttf_hhm_bungalotis_climate_host_cloud_26point_rarefaction_v0.1'
            or rule.get('status')!='POST_OUTCOME_EXPLORATORY_TRANSPORT_AND_METHOD_FREEZE_BEFORE_RAREFACTION_RESULTS'
            or rule['method']['replicates']!=REPLICATES
            or rule['method']['host_reference_original_n']!=HOST_ORIGINAL_N
            or rule['method']['host_rarefied_n']!=RAREFIED_N
            or rule['study_unit']['original_edges']!=64
            or rule['study_unit']['original_localities']!=26
            or rule['study_unit']['species']!=SPECIES):
        raise RuntimeError('scientific exploratory rarefaction contract drift')
    mapping={
        'candidates':'candidate_sha256',
        'localities':'localities_sha256',
        'edges':'edges_sha256',
        'host_pairs':'host_pairs_sha256',
        'native_units':'native_units_sha256',
        'wgsrpd_support':'wgsrpd_support_sha256',
        'frozen_predictors':'FIXED_FROM_LEGACY_PREDICTOR',
        'manifest':'climate_manifest_sha256',
        'original_pilot_edges':'original_pilot_edges_sha256',
        'original_pilot_summary':'original_pilot_summary_sha256',
    }
    for key,sourcekey in mapping.items():
        if sha256_path(getattr(a,key))!=BASE[key]:
            raise RuntimeError('source SHA drift: '+key)
        if sourcekey!='FIXED_FROM_LEGACY_PREDICTOR' and rule['external_sources'].get(sourcekey)!=BASE[key]:
            raise RuntimeError('contract source SHA drift: '+key)
    if rule['external_sources']['chelsa_eight_tif_archive_sha256']!='0a797c609d26a185ba08d209a71203a584fa60f1045c4a8c335934372cee30b2':
        raise RuntimeError('CHELSA transport source drift')


def main()->int:
    import rasterio
    from rasterio.transform import rowcol
    try:
        from scripts.pilot_historical_host_climate_corridor import (
            check_assets_and_read_native_window,analogue_cell_cost,
            cost_distance_between_pixels,VARS, TIMES
        )
        from scripts.build_historical_host_memory_predictor import norm_id,load_assets,sample_matrix
    except ModuleNotFoundError:
        from pilot_historical_host_climate_corridor import (
            check_assets_and_read_native_window,analogue_cell_cost,
            cost_distance_between_pixels,VARS,TIMES
        )
        from build_historical_host_memory_predictor import norm_id,load_assets,sample_matrix
    ap=argparse.ArgumentParser()
    for name in (*BASE.keys(),'contract'):
        ap.add_argument('--'+name.replace('_','-'),required=True,type=Path)
    ap.add_argument('--historical-asset',type=Path,action='append',required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    rule=json.loads(args.contract.read_text())
    verify_fixed_sources(args,rule)
    if len(args.historical_asset)!=8:
        raise RuntimeError('frozen eight raster climate layers required')
    assets=load_assets(args.historical_asset)
    cand=[c for c in rows_from_csv(args.candidates) if c['species']==SPECIES]
    if len(cand)!=1 or int(cand[0]['n_localities'])!=26 or int(cand[0]['edges'])!=64:
        raise RuntimeError('one selected source-blind original species expected')
    loc=sorted((c for c in rows_from_csv(args.localities) if c['species']==SPECIES),key=lambda c:int(c['locality_index']))
    edges=sorted((c for c in rows_from_csv(args.edges) if c['species']==SPECIES),key=lambda c:int(c['edge_index']))
    if (len(loc),len(edges))!=(26,64) or any(int(x['locality_index'])!=i for i,x in enumerate(loc)) or any(int(x['edge_index'])!=i for i,x in enumerate(edges)):
        raise RuntimeError('original 26 endpoint or 64 edge order changed')
    nodes=np.array([(int(z['node_left']),int(z['node_right'])) for z in edges],dtype=int)
    if np.any(nodes<0) or np.any(nodes>=26) or np.any(nodes[:,0]==nodes[:,1]):
        raise RuntimeError('unqualified frozen 64-edge graph')
    lat=np.array([float(z['latitude']) for z in loc]);lon=np.array([float(z['longitude']) for z in loc]);
    selected_ids={norm_id(x) for x in cand[0]['accepted_host_ids'].split(';')}
    sidecar_ids={norm_id(x['accepted_plant_name_id']) for x in rows_from_csv(args.host_pairs) if x['insect_species']==SPECIES}
    if len(selected_ids)!=3 or selected_ids!=sidecar_ids:
        raise RuntimeError('frozen 3 accepted host plant IDs do not align')
    host_units=set()
    for x in rows_from_csv(args.native_units):
        if norm_id(x['accepted_plant_name_id']) in selected_ids:
            host_units.add(x['area_code_l3'].strip())
    support=rows_from_csv(args.wgsrpd_support)
    coords=[(float(x['longitude']),float(x['latitude'])) for x in support]
    now_assets=[assets[(v,20)] for v in VARS]
    climate_support=sample_matrix(coords,now_assets)
    finite=np.all(np.isfinite(climate_support),axis=1)
    if np.sum(finite)!=6343:
        raise RuntimeError('original 6343 global reference climatic observations drift')
    mean=climate_support[finite].mean(axis=0)
    sd=climate_support[finite].std(axis=0,ddof=0)
    if not np.isfinite(sd).all() or np.any(sd<=0):
        raise RuntimeError('frozen 6343 global climate reference has no finite positive spread')
    standardize=lambda x:(x-mean)/sd
    host_cloud=np.stack([standardize(climate_support[i]) for i,z in enumerate(support)
                         if finite[i] and z['area_code_l3'].strip() in host_units])
    if len(host_cloud)!=HOST_ORIGINAL_N:
        raise RuntimeError('frozen 665 host union climate vectors not recovered')
    loc_raw=sample_matrix(list(zip(lon,lat)),now_assets)
    if not np.isfinite(loc_raw).all():
        raise RuntimeError('frozen self-cloud localities have invalid 0BP climate')
    self_cloud=standardize(loc_raw)
    if len(self_cloud)!=RAREFIED_N:
        raise RuntimeError('frozen self reference had other than 26 points')
    raw,transform,shape=check_assets_and_read_native_window(assets,lon,lat)
    if (shape['rows'],shape['cols'])!=(101,57):
        raise RuntimeError('original 101x57 native CHELSA ROI changed')
    grid={}
    for time in TIMES:
        climate=np.asarray(raw[time],dtype=float)
        ok=np.isfinite(climate).all(axis=2)
        if int(ok.sum())!=5757:
            raise RuntimeError('original ROI finite pixel eligibility changed')
        standardized=np.full(climate.shape,np.nan,dtype=float)
        standardized[ok]=standardize(climate[ok])
        grid[time]=standardized
    pixels=np.array([rowcol(transform,float(lon[i]),float(lat[i])) for i in range(26)],dtype=int)
    frozen_costs={}
    for name,cloud in (('host',host_cloud),('self',self_cloud)):
        for time,climate in grid.items():
            raster=analogue_cell_cost(climate,cloud)
            frozen_costs[(name,time)]=cost_distance_between_pixels(raster,transform,pixels,nodes,forbid_self_edges=True)
    pilot_edges=rows_from_csv(args.original_pilot_edges)
    replay_max=verify_original_edges(frozen_costs,pilot_edges,atol=1e-8)
    pilot_summary=json.loads(args.original_pilot_summary.read_text())
    if (pilot_summary.get('species')!=SPECIES or pilot_summary.get('status')!='PASS_CHELSA_CLIMATE_ONLY_PILOT_NO_LAND_ICE_MASK'
        or pilot_summary.get('host_0bp_climate_cloud_points')!=665
        or pilot_summary.get('self_0bp_climate_cloud_points')!=26
        or pilot_summary.get('joint_valid_cells')!=5757):
        raise RuntimeError('frozen original pilot scientific provenance mismatch')
    self_delta=frozen_costs[('self','LGM')]-frozen_costs[('self','0BP')]
    host_full_delta=frozen_costs[('host','LGM')]-frozen_costs[('host','0BP')]
    indices=fixed_subsample_indices()
    sample_deltas=np.empty((REPLICATES,len(edges)),dtype=float)
    for rep,selected in enumerate(indices):
        cloud=host_cloud[selected]
        landscape_now=analogue_cell_cost(grid['0BP'],cloud)
        landscape_lgm=analogue_cell_cost(grid['LGM'],cloud)
        dc=cost_distance_between_pixels(landscape_now,transform,pixels,nodes,forbid_self_edges=True)
        dl=cost_distance_between_pixels(landscape_lgm,transform,pixels,nodes,forbid_self_edges=True)
        sample_deltas[rep]=dl-dc
        if (rep+1)%8==0:
            print(json.dumps({'rarefaction_draws_completed':rep+1,'total':REPLICATES}),flush=True)
    summary=summarize_rarefied(sample_deltas,self_delta)
    baseline_mean=float(host_full_delta.mean())
    result={
        'schema':'ttf_hhm_bungalotis_26point_host_reference_rarefaction_result_v0.1',
        'status':'EXPLORATORY_CLIMATE_ONLY_REFERENCE_SIZE_SENSITIVITY_NOT_GENETIC_TEST',
        'contract_sha256':sha256_path(args.contract),
        'input_sha256':{k:sha256_path(getattr(args,k)) for k in BASE},
        'study':{'species':SPECIES,'localities':26,'frozen_edges':64,'same_pixel_edges':1,
                 'original_host_reference_n':665,'fixed_insect_self_reference_n':26,
                 'rarefied_host_reference_n':26,'global_reference_n':6343,'native_cells':5757,
                 'draws':REPLICATES},
        'frozen_replay_max_absolute_lcp_error':replay_max,
        'source_row_order_bounded_sample_indices_sha256':hashlib.sha256(indices.astype('<i8').tobytes()).hexdigest(),
        'original_full_host_mean_lgm_minus_present_cost_km':baseline_mean,
        'original_fixed_self_mean_lgm_minus_present_cost_km':float(self_delta.mean()),
        'original_full_host_minus_self_mean_cost_km':float(baseline_mean-self_delta.mean()),
        'rarefaction':summary,
        'boundary':{
            'source_panel_already_genetically_opened':True,
            'previous_confirmatory_genetic_status':'NOT_EVALUABLE_EMPIRICAL_RESPONSE',
            'new_genetic_sequences_opened':False,
            'gene_distance_or_p_value_computed':False,
            'rarefaction_quantiles_not_confidence_intervals':True,
            'same_number_of_reference_points_not_same_biological_niche_domain':True,
            'host_climate_reference_is_modern_not_ancient_plant_occupancy':True,
            'no_verified_21ka_land_or_ice_mask':True,
            'ROI_and_native_chelsa_grid_unchanged':True,
            'site_graph_edges_not_independent_replicates':True,
            'not_a_new_confirmatory_hypothesis':True,
        }
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],
                      'rarefied_host_mean':summary['rarefied_host_delta_cost_km']['mean'],
                      'rarefied_minus_self_mean':summary['rarefied_host_minus_self_delta_cost_km']['mean'],
                      'baseline_replay_error':replay_max},sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
