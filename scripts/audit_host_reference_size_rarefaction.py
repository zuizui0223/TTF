#!/usr/bin/env python3
"""Exploratory, post-exposure original Bungalotis host-cloud-size rarefaction.

This uses only the same 26 response-blind sites, 64 original geographical
edges and exact CHELSA/HOSTS/WCVP/WGSRPD3 sources as the original 2D
climate-only corridor pilot. It never reads COI identities/genetic distances.
Source-cloud size equality is NOT spatial or ecological sampling equivalence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from rasterio.transform import rowcol

try:
    from scripts.pilot_historical_host_climate_corridor import (
        SPECIES, VARS, TIMES, BASE_SOURCE_SHA,
        sha256_file, load_csv, check_assets_and_read_native_window,
        analogue_cell_cost, cost_distance_between_pixels,verify_pilot_source_bindings
    )
    from scripts.build_historical_host_memory_predictor import sample_matrix, norm_id, load_assets
except ModuleNotFoundError:
    from pilot_historical_host_climate_corridor import (
        SPECIES, VARS, TIMES, BASE_SOURCE_SHA,
        sha256_file, load_csv, check_assets_and_read_native_window,
        analogue_cell_cost, cost_distance_between_pixels,verify_pilot_source_bindings
    )
    # Source-builder import is performed by main() only, never by the synthetic fixtures.
    sample_matrix=norm_id=load_assets=None

CONTRACT_SCHEMA='ttf_hhm_bungalotis_host_cloud_rarefaction_v0.1'
ORIGINAL_PILOT_SUMMARY_SHA='6969f482a1967db2db159cbc5fd62b722268a0fc45d7dcd1b979f1973b9a3493'
ORIGINAL_PILOT_EDGES_SHA='7d255b208a485e7bc6d02e4e5748781c5e10170ed8bf161542dce3135c897b3b'
SEED_NAMESPACE='ttf-hhm-host-reference-26-of-665-50draws-v0.1'


def sampled_reference_indices(n_total:int,sample_size:int,repetitions:int,*,integer_seed:int) -> np.ndarray:
    if any(isinstance(x,bool) or not isinstance(x,(int,np.integer)) for x in (n_total,sample_size,repetitions,integer_seed)):
        raise ValueError('reference sample counts and seed must be integers')
    if n_total<2 or sample_size<1 or sample_size>n_total or repetitions<1 or integer_seed<0:
        raise ValueError('cannot draw the fixed number of distinct reference climate support vectors')
    rng=np.random.default_rng(int(integer_seed))
    return np.stack([rng.choice(n_total,size=sample_size,replace=False) for _ in range(repetitions)])


def reference_corridor_delta(current_grid,lgm_grid,reference_cloud,transform,pixel_nodes,edges):
    cur=np.asarray(current_grid,float);lgm=np.asarray(lgm_grid,float)
    cloud=np.asarray(reference_cloud,float)
    if cur.ndim!=3 or lgm.shape!=cur.shape or cloud.ndim!=2 or cloud.shape[1]!=cur.shape[2] or not np.isfinite(cloud).all():
        raise ValueError('cloud and two climate grids must have matching finite reference dimensions')
    cost_cur=analogue_cell_cost(cur,cloud)
    cost_lgm=analogue_cell_cost(lgm,cloud)
    r_cur=cost_distance_between_pixels(cost_cur,transform,pixel_nodes,edges,forbid_self_edges=True)
    r_lgm=cost_distance_between_pixels(cost_lgm,transform,pixel_nodes,edges,forbid_self_edges=True)
    return r_cur,r_lgm,r_lgm-r_cur


def verify_baseline(actual,original,*,atol=1e-8)->float:
    if len(actual)!=len(original) or len(actual)!=3:
        raise RuntimeError('original 64-edge pilot replay mismatch')
    error=0.
    for first,second in zip(actual,original):
        a=np.asarray(first,float);b=np.asarray(second,float)
        if a.shape!=b.shape or a.ndim!=1 or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            raise RuntimeError('original 64-edge pilot replay mismatch')
        error=max(error,float(np.max(np.abs(a-b))))
    if error>atol:
        raise RuntimeError('original 64-edge pilot replay mismatch')
    return error


def run_rarefaction(current_grid,lgm_grid,host_cloud,transform,pixel_nodes,edges,
                    *,sample_size=26,repetitions=50,integer_seed:int)->dict:
    host=np.asarray(host_cloud,dtype=float)
    if host.ndim!=2 or len(host)<2 or not np.isfinite(host).all():
        raise ValueError('original host climate cloud invalid')
    inds=sampled_reference_indices(len(host),sample_size,repetitions,integer_seed=integer_seed)
    rows=[]
    for row in inds:
        current,lgm,delta=reference_corridor_delta(current_grid,lgm_grid,host[row],transform,pixel_nodes,edges)
        rows.append({'host_current_mean':float(np.mean(current)),
                     'host_lgm_mean':float(np.mean(lgm)),
                     'host_delta_mean':float(np.mean(delta)),
                     'host_delta_median_original_edges':float(np.median(delta))})
    v=np.asarray([r['host_delta_mean'] for r in rows],float)
    if np.any(~np.isfinite(v)):
        raise RuntimeError('nonfinite climate-only original-edge rarefaction values')
    return {
        'draw_count':int(repetitions),
        'sample_size':int(sample_size),
        'host_delta_mean_by_draw':[float(x) for x in v],
        'host_delta_mean_across_draws':float(np.mean(v)),
        'host_delta_median_across_draws':float(np.median(v)),
        'host_delta_min_across_draws':float(np.min(v)),
        'host_delta_max_across_draws':float(np.max(v)),
        'host_delta_quantile_10':float(np.quantile(v,0.10)),
        'host_delta_quantile_90':float(np.quantile(v,0.90)),
        'new_empirical_genetic_test':False,
        'new_p_values':0,
        'quantiles_are_descriptive_rarefaction_not_statistical_confidence_intervals':True,
    }


def main():
    from ttf.historical_host_memory import standardize_from_reference
    # Import original source-native raster helper in source-compatible direct CLI mode.
    global sample_matrix,norm_id,load_assets
    if sample_matrix is None:
        from build_historical_host_memory_predictor import sample_matrix, norm_id, load_assets
    a=argparse.ArgumentParser()
    for key in ('contract','original_pilot_contract','original_summary','original_edges','candidates','localities','edges',
                'host_pairs','native_units','wgsrpd_support','frozen_predictors','manifest'):
        a.add_argument('--'+key.replace('_','-'),type=Path,required=True)
    a.add_argument('--historical-asset',type=Path,required=True,action='append')
    a.add_argument('--output',type=Path,required=True)
    args=a.parse_args()
    contract=json.loads(args.contract.read_text())
    if contract.get('schema')!=CONTRACT_SCHEMA or contract.get('status')!='POST_LCP_OUTCOME_EXPLORATORY_FIXED_BEFORE_FIRST_HOST_RAREFACTION_RESULT':
        raise RuntimeError('the exploratory host rarefaction contract has not been frozen')
    rr=contract['rarefaction']
    if (rr['sample_size'],rr['draws'],rr['replace_within_draw'],rr['seed_namespace'])!=(26,50,False,SEED_NAMESPACE):
        raise RuntimeError('frozen rarefaction sampling/seed definition changed')
    if (contract['original_host_cloud_n'],contract['original_self_cloud_n'])!=(665,26):
        raise RuntimeError('original host/self sample sizes changed')
    parent=json.loads(args.original_pilot_contract.read_text())
    if sha256_file(args.original_pilot_contract)!=contract['parent_climate_lcp_contract_sha256']:
        raise RuntimeError('original climate LCP contract changed')
    verify_pilot_source_bindings(parent)
    for key in BASE_SOURCE_SHA:
        p=getattr(args,key)
        if sha256_file(p)!=BASE_SOURCE_SHA[key]:
            raise RuntimeError('frozen original geography/host/CHELSA source hash mismatch: '+key)
    if sha256_file(args.original_summary)!=ORIGINAL_PILOT_SUMMARY_SHA or sha256_file(args.original_edges)!=ORIGINAL_PILOT_EDGES_SHA:
        raise RuntimeError('original complete CHELSA pilot source hashes mismatch')
    if len(args.historical_asset)!=8:
        raise RuntimeError('eight frozen original CHELSA variables required')
    assets=load_assets(args.historical_asset)
    cand=[r for r in load_csv(args.candidates) if r['species']==SPECIES]
    if len(cand)!=1 or int(cand[0]['edges'])!=64 or int(cand[0]['n_localities'])!=26:
        raise RuntimeError('frozen Bungalotis candidate drift')
    local=sorted((r for r in load_csv(args.localities) if r['species']==SPECIES),key=lambda q:int(q['locality_index']))
    edge=sorted((r for r in load_csv(args.edges) if r['species']==SPECIES),key=lambda q:int(q['edge_index']))
    if len(local)!=26 or len(edge)!=64 or [int(r['locality_index']) for r in local]!=list(range(26)) or [int(r['edge_index']) for r in edge]!=list(range(64)):
        raise RuntimeError('original 26-node 64-edge graph drift')
    pairs=np.array([[int(r['node_left']),int(r['node_right'])] for r in edge],int)
    if np.any(pairs<0) or np.any(pairs>=26) or len(set(map(tuple,np.sort(pairs,axis=1))))!=64:
        raise RuntimeError('original 64-edge topology drift')
    lat=np.array([float(r['latitude']) for r in local]);lon=np.array([float(r['longitude']) for r in local])
    original_support=load_csv(args.wgsrpd_support)
    current_paths=[assets[(v,20)] for v in VARS]
    support_current=sample_matrix([(float(r['longitude']),float(r['latitude'])) for r in original_support],current_paths)
    finite=np.all(np.isfinite(support_current),axis=1)
    if int(np.sum(finite))!=6343:
        raise RuntimeError('frozen modern WGSRPD3 reference changed')
    reference=support_current[finite]
    mean=reference.mean(axis=0);sd=reference.std(axis=0,ddof=0)
    stand=lambda x:(x-mean)/sd
    original_host_ids={norm_id(x) for x in cand[0]['accepted_host_ids'].split(';')}
    host_pairs={norm_id(r['accepted_plant_name_id']) for r in load_csv(args.host_pairs) if r['insect_species']==SPECIES}
    if original_host_ids!=host_pairs or len(host_pairs)!=3:
        raise RuntimeError('original accepted host plant identity mismatch')
    native=defaultdict(set)
    for row in load_csv(args.native_units):
        ident=norm_id(row['accepted_plant_name_id'])
        if ident in host_pairs:native[ident].add(row['area_code_l3'].strip())
    unit_set=set().union(*(native[x] for x in host_pairs))
    host_cloud=np.stack([stand(support_current[i]) for i,s in enumerate(original_support)
                         if finite[i] and s['area_code_l3'].strip() in unit_set])
    if len(host_cloud)!=665:
        raise RuntimeError('original frozen accepted host climate support 665 drift')
    insect_current=sample_matrix(list(zip(lon,lat)),current_paths)
    if not np.all(np.isfinite(insect_current)):
        raise RuntimeError('original 26 locality modern climate invalid')
    self_cloud=stand(insect_current)
    raw,transform,shape=check_assets_and_read_native_window(assets,lon,lat)
    if (shape['rows'],shape['cols'])!=(101,57):
        raise RuntimeError('original native climate grid 101x57 changed')
    grid={}
    for time in ('0BP','LGM'):
        z=raw[time]
        good=np.all(np.isfinite(z),axis=2)
        normalized=np.full(z.shape,np.nan,dtype=float)
        normalized[good]=stand(z[good])
        grid[time]=normalized
    pixels=np.asarray([rowcol(transform,float(lon[i]),float(lat[i])) for i in range(26)],int)
    original=load_csv(args.original_edges)
    if len(original)!=64 or any(int(r['edge_index'])!=i for i,r in enumerate(original)):
        raise RuntimeError('original 64-edge results order/row identity drift')
    cur=reference_corridor_delta(grid['0BP'],grid['LGM'],host_cloud,transform,pixels,pairs)
    ec=np.array([float(r['R_host_0BP_km']) for r in original])
    el=np.array([float(r['R_host_LGM_km']) for r in original])
    eh=np.array([float(r['M_host_lcp_km']) for r in original])
    replay_err=verify_baseline(cur,(ec,el,eh))
    own=reference_corridor_delta(grid['0BP'],grid['LGM'],self_cloud,transform,pixels,pairs)
    replay_err=max(replay_err,verify_baseline(own,(
        np.array([float(r['R_self_0BP_km']) for r in original]),
        np.array([float(r['R_self_LGM_km']) for r in original]),
        np.array([float(r['M_self_lcp_km']) for r in original])
    )))
    expected_key=f'{SEED_NAMESPACE}|{SPECIES}|{ORIGINAL_PILOT_EDGES_SHA}'
    seed=int(hashlib.sha256(expected_key.encode()).hexdigest()[:16],16)
    rare=run_rarefaction(grid['0BP'],grid['LGM'],host_cloud,transform,pixels,pairs,
                         sample_size=26,repetitions=50,integer_seed=seed)
    self_delta=float(np.mean(own[2]));host_delta=float(np.mean(cur[2]))
    v=np.asarray(rare['host_delta_mean_by_draw'])
    diff=v-self_delta
    result={
      'schema':CONTRACT_SCHEMA+'_result',
      'status':'EXPLORATORY_HOST_REFERENCE_SIZE_SENSITIVITY_ONLY',
      'contract_sha256':sha256_file(args.contract),
      'original_pilot_sha256':{'edges':ORIGINAL_PILOT_EDGES_SHA,'summary':ORIGINAL_PILOT_SUMMARY_SHA},
      'sources_sha256':BASE_SOURCE_SHA,
      'source_manifest_sha256':BASE_SOURCE_SHA['manifest'],
      'species':SPECIES,'sites':26,'edges':64,'host_accepted_plant_species':3,
      'original_reference_cloud_sizes':{'host':665,'self':26,'global_climate_standardization':6343},
      'native_roi':{'rows':101,'cols':57},
      'same_pixel_edges_intact':int(np.sum(np.all(pixels[pairs[:,0]]==pixels[pairs[:,1]],axis=1))),
      'maximum_original_LCP_replay_abs_error':float(replay_err),
      'original_full_host_delta_cost_km_mean':host_delta,
      'original_self_delta_cost_km_mean':self_delta,
      'original_full_host_minus_self_delta_cost_km_mean':host_delta-self_delta,
      'rarefaction':rare,
      'rarefied_host_minus_fixed_self':{
        'mean_across_draws':float(np.mean(diff)),
        'median_across_draws':float(np.median(diff)),
        'min_across_draws':float(diff.min()),
        'max_across_draws':float(diff.max()),
        'fraction_positive_draws_descriptive_only':float(np.mean(diff>0)),
      },
      'source_seed_namespace':SEED_NAMESPACE,
      'source_rng_seed_64bit':int(seed),
      'scientific_boundary':{
        'after_prior_climate_LCP_outcome_exposure':True,
        'no_genetic_response_or_coi_distance_reopened':True,
        'original_259_result':'NOT_EVALUABLE_EMPIRICAL_RESPONSE',
        'matched_reference_n_not_matched_plant_insect_sampling_breadth':True,
        'rarefaction_quantiles_are_not_CI':True,
        'no_new_effect_p_value':True,
        'does_not_reconstruct_LGM_plant_presence_or_barriers':True,
        'fixed_climate_R0_and_5pt_pilot_unchanged':True
      }
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'original_host_delta':host_delta,
        'original_self_delta':self_delta,'rarefied_host_delta_mean':rare['host_delta_mean_across_draws'],
        'rarefied_minus_self_delta_mean':result['rarefied_host_minus_fixed_self']['mean_across_draws']},sort_keys=True))

if __name__=='__main__':main()
