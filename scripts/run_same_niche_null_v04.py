#!/usr/bin/env python3
"""Exploratory real-host versus same-insect-niche disjoint reference-split null.

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
from same_niche_null_v04 import draw_independent_clouds,score_real_and_placebo

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

def main():
    from ttf.historical_host_memory import great_circle_quadrature_latlon
    parser=argparse.ArgumentParser()
    for name in ('host-manifest','host-snapshot','insect-summary','insect-points',
                 'localities','edges','global-support','asset-root','out'):
        parser.add_argument('--'+name,required=True,type=Path)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    for path,expected in ((args.host_manifest,FROZEN_MANIFEST_SHA256),
       (args.localities,FROZEN_LOCALITY_SHA256),(args.edges,FROZEN_EDGE_SHA256)):
        if sha_file(path)!=expected:raise ValueError('original response-free geometry/host identity drift')
    completed,plants,review,failures,plan=load_host_snapshot(args.host_snapshot)
    if (len(plants),len(completed),len(review),len(failures))!=(506,518,1,11):
        raise ValueError('frozen 518-taxon modern v0.2 host inventory drift')
    present_with_records={taxon for taxon,points in plants.items() if points}
    panel=csv_rows(args.insect_summary)
    if len(panel)!=641 or len({r['species'] for r in panel})!=641:
        raise ValueError('641-species domain drift')
    associated=defaultdict(set)
    for r in csv_rows(args.host_manifest):
        associated[r['insect_species']].add(normalize_id(r['wcvp_id_as_frozen']))
    geo=defaultdict(dict);edges=defaultdict(list)
    for r in csv_rows(args.localities):
        geo[r['species']][int(r['locality_index'])]=(float(r['latitude']),float(r['longitude']))
    for r in csv_rows(args.edges):
        edges[r['species']].append((int(r['node_left']),int(r['node_right'])))
    insect_points=defaultdict(list)
    with gzip.open(args.insect_points,'rt',encoding='utf-8',newline='') as f:
        for r in csv.DictReader(f):insect_points[r['species']].append(r)
    if set(insect_points)-{r['species'] for r in panel}:
        raise ValueError('external insect observations outside panel')
    todo={};allcoords=[];coordinate_index={}
    def indices(rows):
        out=[]
        for lat,lon in rows:
            lat=float(lat);lon=float(lon)
            pair=(round(lat,8),round(lon,8))
            if pair not in coordinate_index:
                coordinate_index[pair]=len(allcoords);allcoords.append((lon,lat))
            out.append(coordinate_index[pair])
        return np.array(out,dtype=int)
    rows_status=[]
    for row in sorted(panel,key=lambda r:r['species']):
        name=row['species'];ids=associated[name];present=ids & plants.keys()
        hrecords=join_hosts(present,plants)
        loc=geo[name];pairs=edges[name]
        orig=np.array([loc[i] for i in sorted(loc)],float)
        ed=np.array([[loc[i],loc[j]] for i,j in pairs],float)
        quad=great_circle_quadrature_latlon(ed[:,0,:],ed[:,1,:],segment_points=5)
        support=np.vstack([orig,quad.reshape(-1,2)])
        selected=filter_host_points(hrecords,support,buffer_km=50,thinning_km=10,cap=200)
        insects=insect_points.get(name,[])
        good=len(selected)>=30 and len(insects)>=30
        null_eligible=good and len(insects)>=60
        all_known_host_records=ids<=present_with_records
        all_known_host_taxonomy=ids<=plants.keys()
        meta={'species':name,'insect_n':len(insects),'host_n':len(selected),
              'total_host_taxa':len(ids),'observed_host_taxa':len(present),
              'host_taxa_with_observations':len(ids & present_with_records),
              'missing_host_taxa':len(ids-present),
              'all_known_hosts_with_observations':all_known_host_records,
              'all_known_host_taxonomy_resolved':all_known_host_taxonomy,
              'state':'READY_FOR_SPLIT_NULL' if null_eligible else ('PAIRED_BUT_INSUFFICIENT_INSECT_SPLIT' if good else 'NO_PAIRED_REFERENCE')}
        rows_status.append(meta)
        if not null_eligible:continue
        todo[name]={'host_idx':indices([(p['latitude'],p['longitude']) for p in selected]),
                    'insect_idx':indices([(p['latitude'],p['longitude']) for p in insects]),
                    'query_idx':indices(quad.reshape(-1,2)).reshape(len(pairs),5),
                    'distance':np.log1p(np.array([edge_dist_km(pair) for pair in ed])),
                    'all_known_host_records':all_known_host_records,
                    'all_known_host_taxonomy':all_known_host_taxonomy}
    if len(todo)!=348:raise ValueError('same-niche reference-split qualification of 348 species changed')
    global_rows=csv_rows(args.global_support)
    support_idx=indices([(r['latitude'],r['longitude']) for r in global_rows])
    curr_paths=[args.asset_root/f'CHELSA_TraCE21K_{v}_20_V1.0.tif' for v in VARIABLES]
    lgm_paths=[args.asset_root/f'CHELSA_TraCE21K_{v}_-190_V1.0.tif' for v in VARIABLES]
    for file in curr_paths+lgm_paths:
        if not file.is_file():raise FileNotFoundError(str(file))
    current=read_climate_points(curr_paths,allcoords)
    lgm=read_climate_points(lgm_paths,allcoords)
    ref=current[support_idx]
    finite_ref=np.isfinite(ref).all(axis=1)
    if int(finite_ref.sum())!=6343:raise ValueError('common 6343 reference points not recovered')
    mean=ref[finite_ref].mean(axis=0);sd=ref[finite_ref].std(axis=0)
    if not np.all(sd>0):raise ValueError('invalid climate standardization')
    current=(current-mean)/sd;lgm=(lgm-mean)/sd
    results=[]
    for name,data in todo.items():
        hi,si,qi=data['host_idx'],data['insect_idx'],data['query_idx']
        x=current[qi];y=lgm[qi];h=current[hi];s=current[si]
        base={'species':name,'host_reference_n':int(len(hi)),
              'insect_reference_n':int(len(si)), 'edges':int(len(qi)),
              'all_known_hosts_with_observations':data['all_known_host_records'],
              'all_known_host_taxonomy_resolved':data['all_known_host_taxonomy']}
        if not all(np.isfinite(a).all() for a in (x,y,h,s)):
            results.append({**base,'status':'NONFINITE_CLIMATE','draws_completed':0})
            continue
        draws=[]
        unidentifiable=False
        for d in range(DRAWS):
            actual_idx,self_idx,pseudo_idx=draw_independent_clouds(len(s),len(h),name,d)
            real,null=score_real_and_placebo(x,y,h[actual_idx],s[self_idx],s[pseudo_idx],data['distance'])
            if not real['identifiable'] or not null['identifiable']:
                unidentifiable=True
                break
            if len({len(actual_idx),len(self_idx),len(pseudo_idx)})!=1:
                raise RuntimeError('not matched source count')
            draws.append((real,null))
        if unidentifiable:
            results.append({**base,'status':'UNIDENTIFIABLE_SPLIT_REFERENCE','draws_completed':len(draws)})
            continue
        actual_fraction=np.asarray([a['unique_fraction'] for a,b in draws])
        null_fraction=np.asarray([b['unique_fraction'] for a,b in draws])
        actual_raw=np.asarray([a['residual_variance'] for a,b in draws])
        null_raw=np.asarray([b['residual_variance'] for a,b in draws])
        actual_current=np.asarray([a['host_current_mean'] for a,b in draws])
        null_current=np.asarray([b['host_current_mean'] for a,b in draws])
        actual_history=np.asarray([a['history_variance'] for a,b in draws])
        null_history=np.asarray([b['history_variance'] for a,b in draws])
        results.append({**base,'status':'COMPLETE_SAME_NICHE_NULL','draws_completed':DRAWS,
            'reference_n_per_cloud':len(draws) and len(draw_independent_clouds(len(s),len(h),name,0)[0]),
            'actual_unique_fraction_median':float(np.median(actual_fraction)),
            'pseudo_unique_fraction_median':float(np.median(null_fraction)),
            'paired_unique_fraction_diff_median':float(np.median(actual_fraction-null_fraction)),
            'diff_of_medians_unique_fraction':float(np.median(actual_fraction)-np.median(null_fraction)),
            'draws_real_unique_greater':int(np.sum(actual_fraction>null_fraction)),
            'actual_raw_residual_variance_median':float(np.median(actual_raw)),
            'pseudo_raw_residual_variance_median':float(np.median(null_raw)),
            'paired_raw_residual_variance_diff_median':float(np.median(actual_raw-null_raw)),
            'actual_total_history_variance_median':float(np.median(actual_history)),
            'pseudo_total_history_variance_median':float(np.median(null_history)),
            'actual_current_climate_mismatch_median':float(np.median(actual_current)),
            'pseudo_current_climate_mismatch_median':float(np.median(null_current)),
            'actual_host_delta_mean_median':float(np.median([a['host_delta_mean'] for a,b in draws])),
            'pseudo_host_delta_mean_median':float(np.median([b['host_delta_mean'] for a,b in draws]))})
        if len(results)%20==0:print(json.dumps({'processed':len(results),'qualified':len(todo)}),flush=True)
    fields=sorted(set().union(*(r.keys() for r in results)))
    with (args.out/'species_same_niche_null_v0.4.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(results)
    fields=sorted(set().union(*(r.keys() for r in rows_status)))
    with (args.out/'all_641_split_eligibility_v0.4.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows_status)
    ok=[r for r in results if r['status']=='COMPLETE_SAME_NICHE_NULL']
    complete_hosts=[r for r in ok if r['all_known_hosts_with_observations']]
    out={'schema':'ttf_lgm_plant_climate_same_niche_split_null_v0.4',
         'status':'EXPLORATORY_SAMPLING_NULL_NOT_HISTORICAL_PLANT_OCCUPANCY_OR_GENETIC_RESPONSE',
         'insect_panel':641,'previously_paired_reference_qualified':398,
         'split_eligible_species':348,'split_draws_per_species':DRAWS,
         'status_counts':dict(Counter(r['status'] for r in results)),
         'null_computable_species':len(ok),
         'with_all_known_hosts_observed':len(complete_hosts),
         'actual_unique_fraction_median':float(np.median([r['actual_unique_fraction_median'] for r in ok])) if ok else None,
         'pseudo_unique_fraction_median':float(np.median([r['pseudo_unique_fraction_median'] for r in ok])) if ok else None,
         'median_species_paired_unique_fraction_difference':float(np.median([r['diff_of_medians_unique_fraction'] for r in ok])) if ok else None,
         'species_actual_fraction_exceeds_pseudo':int(sum(r['diff_of_medians_unique_fraction']>0 for r in ok)),
         'median_species_paired_raw_residual_difference':float(np.median([r['paired_raw_residual_variance_diff_median'] for r in ok])) if ok else None,
         'species_actual_raw_variance_exceeds_pseudo':int(sum(r['paired_raw_residual_variance_diff_median']>0 for r in ok)),
         'interpretation_boundaries':['The surrogate is another independent sample of the same insect climate niche, not another known host plant',
            'The unique fractions have different denominators; raw residual variances must be compared alongside fractions',
            'Modern host and insect observations do not reconstruct LGM plant occupation or historical corridors',
            'The 32 draws are dependent reference splits, not biological replicates',
            'Plant nativity and sampling completeness remain unknown; observational reference data are bounded',
            'No genetic response or genetic barrier effect was opened, no confirmatory p-value']}
    (args.out/'summary.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
