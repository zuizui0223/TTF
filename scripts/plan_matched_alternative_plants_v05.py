#!/usr/bin/env python3
"""Stage-I modern-climate matching of unlisted alternatives, without LGM climate.

Exploratory ecological predictor information only; no genetic response, no
reconstructed occupied LGM host ranges or inference on species barriers.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, math
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from reference_helpers_min import (csv_rows,sha_file,normalize_id,
    filter_host_points,join_hosts,FROZEN_MANIFEST_SHA256,
    FROZEN_LOCALITY_SHA256,FROZEN_EDGE_SHA256)
from host_geometry_yearstrat_v03 import load_host_snapshot

VARIABLES=('bio01','bio07','bio12','bio15')
DRAWS=32

def match_sampling_indices(n_h,n_s,name,draw):
    n=min(n_h,n_s)
    if n<30:raise ValueError('at least 30 external observations per reference required')
    seed=int.from_bytes(hashlib.sha256(f'ttf-climate-match-v0.1|{name}|{draw}'.encode()).digest()[:8],'big')
    rng=np.random.default_rng(seed)
    h=np.sort(rng.choice(n_h,size=n,replace=False)) if n_h!=n else np.arange(n)
    s=np.sort(rng.choice(n_s,size=n,replace=False)) if n_s!=n else np.arange(n)
    return h,s

def analogue(points,cloud):
    if len(cloud)==0 or not np.isfinite(cloud).all() or not np.isfinite(points).all():
        raise ValueError('nonfinite or empty climate point cloud')
    return np.asarray(cKDTree(cloud).query(points,k=1,workers=1)[0],float)

def residual_fraction(y,controls):
    y=np.asarray(y,float).reshape(-1)
    c=np.asarray(controls,float)
    if c.ndim!=2 or len(c)!=len(y) or len(y)<=c.shape[1]+2 or not np.isfinite(c).all() or not np.isfinite(y).all():
        return None
    scale=float(y.std())
    if scale <=1e-12:return None
    design=np.column_stack([np.ones(len(y)),c])
    if np.linalg.matrix_rank(design)<design.shape[1]:
        return None
    r=y-design@np.linalg.lstsq(design,y,rcond=None)[0]
    return float(np.var(r)/np.var(y))

def matched_metrics(current,lgm,host,self_,log_geodist,species,draw):
    current=np.asarray(current,float);lgm=np.asarray(lgm,float)
    if (current.ndim!=3 or lgm.shape!=current.shape or current.shape[2]!=4 or
        current.shape[1]!=5 or len(log_geodist)!=len(current)):
        raise ValueError('expected edge/5-interior-points/4-climate arrays')
    h_idx,s_idx=match_sampling_indices(len(host),len(self_),species,draw)
    h=host[h_idx];s=self_[s_idx]
    n_edges,n_positions,n_dim=current.shape
    flat_cur=current.reshape(-1,n_dim);flat_lgm=lgm.reshape(-1,n_dim)
    h_cur=analogue(flat_cur,h).reshape(n_edges,n_positions).mean(axis=1)
    h_lgm=analogue(flat_lgm,h).reshape(n_edges,n_positions).mean(axis=1)
    s_cur=analogue(flat_cur,s).reshape(n_edges,n_positions).mean(axis=1)
    s_lgm=analogue(flat_lgm,s).reshape(n_edges,n_positions).mean(axis=1)
    hm=h_lgm-h_cur;sm=s_lgm-s_cur
    ctrl=np.column_stack([sm,h_cur,s_cur,log_geodist])
    before_geo=residual_fraction(hm,ctrl[:,:3])
    return {'host_delta_mean':float(hm.mean()),'self_delta_mean':float(sm.mean()),
            'relative_host_minus_self':float(np.mean(hm-sm)),
            'unique_history_fraction':residual_fraction(hm,ctrl),
            'unique_before_geography':before_geo,
            'edge_n':int(n_edges),'reference_n':len(h_idx)}

def read_climate_points(asset_paths,coords):
    import rasterio
    out=[]
    for path in asset_paths:
        with rasterio.open(path) as ds:
            if ds.crs is None or not ds.crs.is_geographic:
                raise ValueError('TIFF not geographically referenced')
            a=np.asarray([float(vals[0]) for vals in ds.sample(coords,masked=False)],float)
            if ds.nodata is not None:
                a[a==ds.nodata]=np.nan
            out.append(a)
    return np.column_stack(out)

def edge_dist_km(pair):
    a,b=pair
    la1,lo1=np.deg2rad(a);la2,lo2=np.deg2rad(b)
    h=np.sin((la1-la2)/2)**2+np.cos(la1)*np.cos(la2)*np.sin((lo1-lo2)/2)**2
    return float(2*6371.0088*np.arcsin(np.sqrt(np.clip(h,0,1))))

from alternative_plant_match_v05 import (genus_of,stable_sample_indices,climate_signature,
    rank_alternatives,group_alternatives,observation_signature,tag_donor_reference)

def main():
    from ttf.historical_host_memory import great_circle_quadrature_latlon
    parser=argparse.ArgumentParser()
    for field in ('host-manifest','host-snapshot','insect-summary','insect-points',
                  'original-v04-eligibility','localities','edges','global-support',
                  'asset-root','out'):
        parser.add_argument('--'+field,required=True,type=Path)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if sha_file(args.original_v04_eligibility)!='7cdaa90b869fd7fe03d14c2928e9eea2c2c3df9bcf0ae9b609c7f81d567eed06':
        raise ValueError('closed v04 reference-split-eligible species cohort drift')
    for p,ref in [(args.host_manifest,FROZEN_MANIFEST_SHA256),
                  (args.localities,FROZEN_LOCALITY_SHA256),(args.edges,FROZEN_EDGE_SHA256)]:
        if sha_file(p)!=ref:raise ValueError('source identity mismatch: '+str(p))
    completed,plants,review,failures,host_plan=load_host_snapshot(args.host_snapshot)
    assert (len(completed),len(plants),len(review),len(failures))==(518,506,1,11)
    taxa={t['wcvp_id']:t for t in host_plan['taxa']}
    alias=defaultdict(set)
    for r in csv_rows(args.host_manifest):
        alias[r['insect_species']].add(normalize_id(r['wcvp_id_as_frozen']))
    insect=defaultdict(list)
    with gzip.open(args.insect_points,'rt',encoding='utf8',newline='') as f:
        for r in csv.DictReader(f):insect[r['species']].append(r)
    local=defaultdict(dict);edges=defaultdict(list)
    for r in csv_rows(args.localities):
        local[r['species']][int(r['locality_index'])]=(float(r['latitude']),float(r['longitude']))
    for r in csv_rows(args.edges):
        edges[r['species']].append((int(r['node_left']),int(r['node_right'])))
    base_elig=csv_rows(args.original_v04_eligibility)
    if len(base_elig)!=641 or len({r['species'] for r in base_elig})!=641:
        raise ValueError('v04 eligibility full domain drift')
    candidate_names=sorted(r['species'] for r in base_elig if r['state']=='READY_FOR_SPLIT_NULL')
    if len(candidate_names)!=348:raise ValueError('v04 split cohort changed')
    active_ids={k for k,v in plants.items() if len(v)>0}
    expected=[n for n in candidate_names if alias[n]<=active_ids]
    if len(expected)!=336:raise ValueError('336 source-complete v04 eligibility drift')
    # Stage I, response blind: build only current (0BP) climate signatures.
    allcoords=[];seen={}
    def indices(points):
        ix=[]
        for lat,lon in points:
            p=(round(float(lat),8),round(float(lon),8))
            if p not in seen:seen[p]=len(allcoords);allcoords.append((p[1],p[0]))
            ix.append(seen[p])
        return np.array(ix,dtype=int)
    sample_index={}
    for key,obs in plants.items():
        if len(obs)<30:continue
        sample=[obs[i] for i in stable_sample_indices(obs,cap=64)]
        sample_index[key]=indices([(p['latitude'],p['longitude']) for p in sample])
    tasks={};eligibility=[]
    for record in sorted(base_elig,key=lambda r:r['species']):
        name=record['species'];known=alias[name]
        status=('SOURCE_V04_NOT_ELIGIBLE' if name not in candidate_names
                else 'MISSING_RECORDED_HOST_TAXON' if name not in expected else 'PENDING_CURRENT_SIGNATURE')
        eligibility.append({'species':name,'status':status,'host_count':len(known),
                           'insect_n':len(insect.get(name,[]))})
        if status!='PENDING_CURRENT_SIGNATURE':continue
        loc=local[name]; pairs=edges[name]
        orig=np.array([loc[i] for i in sorted(loc)],float)
        ed=np.array([[loc[i],loc[j]] for i,j in pairs],float)
        quadrature=great_circle_quadrature_latlon(ed[:,0,:],ed[:,1,:],segment_points=5)
        support=np.vstack([orig,quadrature.reshape(-1,2)])
        actual=filter_host_points(join_hosts(known,plants),support,buffer_km=50,thinning_km=10,cap=200)
        if len(actual)<30 or len(insect[name])<60:raise ValueError('original v04 qualification geometry drift')
        tasks[name]={'known_host_ids':sorted(known),'support':support,
                    'actual_point_n':len(actual),
                    'actual_point_digest':observation_signature(actual),
                    'actual_idx':indices([(p['latitude'],p['longitude']) for p in actual])}
    support_global=csv_rows(args.global_support)
    if len(support_global)<6343:raise ValueError('global support incomplete')
    global_idx=indices([(r['latitude'],r['longitude']) for r in support_global])
    climate_assets=[args.asset_root/f'CHELSA_TraCE21K_{v}_20_V1.0.tif'
                    for v in VARIABLES]
    if any(not p.is_file() for p in climate_assets):raise FileNotFoundError('current CHELSA asset missing')
    modern=read_climate_points(climate_assets,allcoords)
    reference=modern[global_idx]
    finite=reference[np.isfinite(reference).all(axis=1)]
    if len(finite)!=6343:raise ValueError('original modern global support 6343 drift')
    ref_mu=finite.mean(axis=0);ref_sd=finite.std(axis=0)
    if not np.all(ref_sd>0):raise ValueError('reference standard deviation invalid')
    modern=(modern-ref_mu)/ref_sd
    signatures={key:climate_signature(modern[idx]) for key,idx in sample_index.items()}
    signatures={key:s for key,s in signatures.items() if s is not None}
    assignments={};stat=Counter();n_assess=0
    for item in eligibility:
        if item['status']!='PENDING_CURRENT_SIGNATURE':stat[item['status']]+=1;continue
        name=item['species'];task=tasks[name];known=set(task['known_host_ids'])
        sig=climate_signature(modern[task['actual_idx']])
        if sig is None:
            item['status']='NONFINITE_MODERN_ACTUAL_HOST_REFERENCE'
            stat[item['status']]+=1;continue
        ranked=rank_alternatives(sig,known,taxa,signatures)
        alternatives=[];rank_scores={}
        for score,key in ranked:
            donor=filter_host_points(tag_donor_reference(key,plants),task['support'],buffer_km=50,thinning_km=10,cap=200)
            n_assess+=1
            if len(donor)<30:continue
            alternatives.append(key);rank_scores[key]=score
            if len(alternatives)==24:break
        n_host=min(len(known),3)
        groups=group_alternatives(alternatives,n_host,max_group=8,max_taxa=24)
        viable=[]
        for group in groups:
            pts=filter_host_points(join_hosts(group,plants),task['support'],
                                    buffer_km=50,thinning_km=10,cap=200)
            if len(pts)<30:continue
            viable.append({'group_index':len(viable),'donor_taxa':group,
                'n':len(pts),'point_digest':observation_signature(pts),
                'mean_plant_signature_gap':float(np.mean([rank_scores[x] for x in group]))})
        item.update(alternatives_checked=len(alternatives),group_n=len(viable),
                    actual_host_n=task['actual_point_n'])
        if len(viable)<4:
            item['status']='FEWER_THAN_FOUR_VALID_ALTERNATIVE_PLANT_GROUPS'
            stat[item['status']]+=1;continue
        item['status']='MATCHED_ALTERNATIVE_PLANT_CONTROL'
        stat[item['status']]+=1
        assignments[name]={'known_host_taxa':task['known_host_ids'],
            'actual_n':task['actual_point_n'],
            'actual_point_digest':task['actual_point_digest'],
            'alternative_groups':viable}
    for k in assignments:
        if k not in expected:raise ValueError('unexpected target made eligible')
    result={'schema':'ttf_lgm_matched_unlisted_plants_selection_v05',
         'stage':'CURRENT_0BP_SELECTION_ONLY_NO_LGM_CLIMATE_OPENED',
         'source_348_reference_split_cohort':len(candidate_names),
         'source_336_with_all_known_hosts_observed':len(expected),
         'eligible_alternative_control_n':len(assignments),
         'plant_taxa_considered':len(signatures),
         'distinct_alternative_taxa_assessed':n_assess,
         'status_counts':dict(stat),
         'matched_species':assignments,
         'source_sha256':{'manifest':sha_file(args.host_manifest),'original_v04_eligibility':sha_file(args.original_v04_eligibility),
                          'original_localities':sha_file(args.localities),'original_edges':sha_file(args.edges)},
         'modern_reference_mean':ref_mu.tolist(),'modern_reference_sd':ref_sd.tolist(),
         'no_lgm_raster_or_genetic_data_used':True,
         'caveat':'Unlisted modern plant taxa are not confirmed nonhosts. Moment matching is modern-only and may be imperfect.'}
    (args.out/'source_bound_0bp_alternatives.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    with (args.out/'all_641_current_matching_eligibility.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in eligibility))))
        w.writeheader();w.writerows(eligibility)
    print(json.dumps({k:v for k,v in result.items() if k!='matched_species'},sort_keys=True))
if __name__=='__main__':main()
