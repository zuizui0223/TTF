#!/usr/bin/env python3
"""Score modern-only preselected alternative plant groups on 0BP/LGM CHELSA without genotypes.

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

from same_niche_null_v04 import score_real_and_placebo
from alternative_plant_match_v05 import equal_reference_indices,observation_signature,genus_of

def main():
    from ttf.historical_host_memory import great_circle_quadrature_latlon
    parser=argparse.ArgumentParser()
    for field in ('matched-plan','host-manifest','host-snapshot','insect-points',
                  'localities','edges','global-support','asset-root','out'):
        parser.add_argument('--'+field,required=True,type=Path)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    assignment=json.loads(args.matched_plan.read_text())
    if (assignment['schema']!='ttf_lgm_matched_unlisted_plants_selection_v05'
        or assignment['stage']!='CURRENT_0BP_SELECTION_ONLY_NO_LGM_CLIMATE_OPENED'
        or assignment['source_348_reference_split_cohort']!=348
        or assignment['source_336_with_all_known_hosts_observed']!=336
        or not assignment['no_lgm_raster_or_genetic_data_used']):
        raise ValueError('modern-only alternate plant matching plan drift')
    src=assignment['source_sha256']
    for path,key in ((args.host_manifest,'manifest'),(args.localities,'original_localities'),
                     (args.edges,'original_edges')):
        if sha_file(path)!=src[key]:raise ValueError('original 0BP-selection source binding mismatch')
    completed,plants,review,failures,host_plan=load_host_snapshot(args.host_snapshot)
    if (len(completed),len(plants),len(review),len(failures))!=(518,506,1,11):
        raise ValueError('source v0.2 host data drift')
    taxa={t['wcvp_id']:t for t in host_plan['taxa']}
    insects=defaultdict(list)
    with gzip.open(args.insect_points,'rt',newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):insects[r['species']].append(r)
    local=defaultdict(dict);edges=defaultdict(list)
    for r in csv_rows(args.localities):
        local[r['species']][int(r['locality_index'])]=(float(r['latitude']),float(r['longitude']))
    for r in csv_rows(args.edges):
        edges[r['species']].append((int(r['node_left']),int(r['node_right'])))
    todo={};allcoords=[];index={}
    def indices(x):
        out=[]
        for lat,lon in x:
            key=(round(float(lat),8),round(float(lon),8))
            if key not in index:index[key]=len(allcoords);allcoords.append((key[1],key[0]))
            out.append(index[key])
        return np.array(out,dtype=int)
    for name,item in sorted(assignment['matched_species'].items()):
        loc=local[name];pairs=edges[name]
        arr=np.array([[loc[i],loc[j]] for i,j in pairs],float)
        quad=great_circle_quadrature_latlon(arr[:,0,:],arr[:,1,:],segment_points=5)
        support=np.vstack((np.array([loc[i] for i in sorted(loc)],float),quad.reshape(-1,2)))
        real=filter_host_points(join_hosts(item['known_host_taxa'],plants),support,
                                buffer_km=50,thinning_km=10,cap=200)
        if len(real)!=item['actual_n'] or observation_signature(real)!=item['actual_point_digest']:
            raise ValueError('actual host points changed relative to 0BP-only matching stage: '+name)
        known=set(item['known_host_taxa']);known_genera={genus_of(taxa[x]['name']) for x in known}
        alternatives=[]
        for group in item['alternative_groups']:
            ids=group['donor_taxa']
            if known.intersection(ids) or any(genus_of(taxa[t]['name']) in known_genera for t in ids):
                raise ValueError('an alternative host ID or host genus was included')
            pts=filter_host_points(join_hosts(ids,plants),support,buffer_km=50,thinning_km=10,cap=200)
            if (len(pts)!=group['n'] or
                observation_signature(pts)!=group['point_digest']):
                raise ValueError('alternate plant points changed relative to 0BP-only matching')
            alternatives.append({'donor_taxa':ids,'group_index':group['group_index'],
               'modern_matching_score':group['mean_plant_signature_gap'],
               'index':indices([(p['latitude'],p['longitude']) for p in pts])})
        if not 4<=len(alternatives)<=8:raise ValueError('not enough stage1 matched controls')
        if len(insects[name])<60 or len(real)<30:raise ValueError('source v04 disjoint-cohort drift')
        todo[name]={'real_idx':indices([(p['latitude'],p['longitude']) for p in real]),
            'insect_idx':indices([(p['latitude'],p['longitude']) for p in insects[name]]),
            'quad_idx':indices(quad.reshape(-1,2)).reshape(len(pairs),5),
            'geo':np.log1p(np.array([edge_dist_km(p) for p in arr])),
            'groups':alternatives,
            'num_real_host_taxa':len(known),'edge_n':len(pairs)}
    if len(todo)!=assignment['eligible_alternative_control_n']:
        raise ValueError('source-bound alternate group plan species count changed')
    global_rows=csv_rows(args.global_support)
    global_idx=indices([(r['latitude'],r['longitude']) for r in global_rows])
    cp=[args.asset_root/f'CHELSA_TraCE21K_{v}_20_V1.0.tif' for v in VARIABLES]
    lp=[args.asset_root/f'CHELSA_TraCE21K_{v}_-190_V1.0.tif' for v in VARIABLES]
    if any(not p.is_file() for p in cp+lp):raise FileNotFoundError('CHELSA 0BP/LGM input missing')
    curr=read_climate_points(cp,allcoords)
    past=read_climate_points(lp,allcoords)
    ref=curr[global_idx];ref=ref[np.isfinite(ref).all(axis=1)]
    if len(ref)!=6343:raise ValueError('common 6343 present 0BP climate support changed')
    mu,sd=ref.mean(axis=0),ref.std(axis=0)
    if np.any(sd<=0):raise ValueError('degenerate climate standardization')
    curr=(curr-mu)/sd;past=(past-mu)/sd
    results=[];group_rows=[]
    for i,(name,data) in enumerate(sorted(todo.items())):
        q_idx=data['quad_idx'];h_idx=data['real_idx'];i_idx=data['insect_idx']
        x=curr[q_idx];y=past[q_idx];host=curr[h_idx];insect=curr[i_idx]
        present_groups=[curr[g['index']] for g in data['groups']]
        # Modern observation reference clouds are evaluated at 0BP only;
        # the LGM grid is sampled at the five path quadrature positions.
        arrays=[x,y,host,insect,*present_groups]
        row={'species':name,'edge_n':data['edge_n'],'real_host_taxa':data['num_real_host_taxa'],
            'actual_host_n':len(h_idx),'insect_n':len(i_idx),
            'alternate_groups':len(data['groups'])}
        if not all(np.isfinite(a).all() for a in arrays):
            row.update(status='NONFINITE_CLIMATE',draws_completed=0)
            results.append(row);continue
        outputs=[];failed=False
        for group,gcloud in zip(data['groups'],present_groups):
            pair=[]
            for draw in range(8):
                si,ri,ci=equal_reference_indices(len(host),len(gcloud),len(insect),name,
                                                 group['group_index'],draw)
                act,alt=score_real_and_placebo(x,y,host[ri],insect[si],gcloud[ci],data['geo'])
                if not act['identifiable'] or not alt['identifiable']:
                    failed=True;break
                pair.append((act,alt))
            if failed:break
            outputs.extend(pair)
            group_rows.append({'species':name,'group_index':group['group_index'],
                'donor_plant_ids':';'.join(group['donor_taxa']),
                'modern_plant_climate_signature_gap':group['modern_matching_score'],
                'actual_unique_fraction_median':float(np.median([z[0]['unique_fraction'] for z in pair])),
                'alternate_unique_fraction_median':float(np.median([z[1]['unique_fraction'] for z in pair])),
                'actual_current_mismatch_median':float(np.median([z[0]['host_current_mean'] for z in pair])),
                'alternate_current_mismatch_median':float(np.median([z[1]['host_current_mean'] for z in pair]))})
        if failed:
            row.update(status='UNIDENTIFIABLE_NULL_COMPARISON',draws_completed=len(outputs))
            results.append(row);continue
        a=np.array([v[0]['unique_fraction'] for v in outputs]);b=np.array([v[1]['unique_fraction'] for v in outputs])
        av=np.array([v[0]['residual_variance'] for v in outputs]);bv=np.array([v[1]['residual_variance'] for v in outputs])
        ac=np.array([v[0]['host_current_mean'] for v in outputs]);bc=np.array([v[1]['host_current_mean'] for v in outputs])
        row.update(status='COMPLETE_MODERN_MATCHED_PLANT_CONTROL',draws_completed=len(outputs),
            'actual_fraction_median':float(np.median(a)),
            'alternate_fraction_median':float(np.median(b)),
            'fraction_diff_of_medians':float(np.median(a)-np.median(b)),
            'fraction_paired_diff_median':float(np.median(a-b)),
            'actual_raw_residual_variance_median':float(np.median(av)),
            'alternate_raw_residual_variance_median':float(np.median(bv)),
            'raw_residual_paired_diff_median':float(np.median(av-bv)),
            'current_climate_mismatch_actual_median':float(np.median(ac)),
            'current_climate_mismatch_alternate_median':float(np.median(bc)),
            'current_mismatch_gap':float(np.median(ac)-np.median(bc)),
            'actual_exceeds_alternate_draws':int(np.sum(a>b)))
        results.append(row)
        if (i+1)%25==0:print(json.dumps({'species':i+1,'matched_species':len(todo)}),flush=True)
    fields=sorted(set().union(*(r.keys() for r in results)))
    with (args.out/'species_v05.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(results)
    if group_rows:
        with (args.out/'plant_group_v05.csv').open('w',newline='',encoding='utf8') as f:
            w=csv.DictWriter(f,fieldnames=sorted(group_rows[0]))
            w.writeheader();w.writerows(group_rows)
    ok=[r for r in results if r['status']=='COMPLETE_MODERN_MATCHED_PLANT_CONTROL']
    out={'schema':'ttf_lgm_modern_climate_matched_unlisted_plant_control_v0.5',
         'status':'EXPLORATORY_PLANT_NEGATIVE_CONTROL_NOT_GENETIC_OR_OCCUPANCY',
         'source_matched_species':len(todo),
         'climate_completed_species':len(ok),
         'statuses':dict(Counter(r['status'] for r in results)),
         'median_actual_unique_fraction':float(np.median([r['actual_fraction_median'] for r in ok])) if ok else None,
         'median_altplant_unique_fraction':float(np.median([r['alternate_fraction_median'] for r in ok])) if ok else None,
         'median_real_minus_alt_fraction':float(np.median([r['fraction_diff_of_medians'] for r in ok])) if ok else None,
         'species_actual_unique_greater':int(sum(r['fraction_diff_of_medians']>0 for r in ok)),
         'species_actual_raw_greater':int(sum(r['raw_residual_paired_diff_median']>0 for r in ok)),
         'median_modern_mismatch_gap':float(np.median([r['current_mismatch_gap'] for r in ok])) if ok else None,
         'median_raw_residual_diff':float(np.median([r['raw_residual_paired_diff_median'] for r in ok])) if ok else None,
         'source_plan_sha256':sha_file(args.matched_plan),
         'no_genetic_response_opened':True,
         'not_historical_plant_occupancy':True,
         'caveats':['Nonlisted plant taxa are not known true nonhosts',
                    'Current climate moment matching may leave modern fit or spatial/habitat discrepancies',
                    'Source host plant observations are 3 nonrandom shallow pages per year-bin and nativity is unknown',
                    'Repeated plant groups/draws are not independent organismal replicates',
                    'Original NOT_EVALUABLE_EMPIRICAL_RESPONSE unchanged']}
    (args.out/'summary.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
