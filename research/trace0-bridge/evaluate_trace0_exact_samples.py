#!/usr/bin/env python3
"""Consume exact COG samples from cloud workflow and recompute TTF-Q C2.

No genetic response, no re-tuning, and no simulation values substituted for
actual SHA-confirmed CHELSA-TraCE21k v1.0 rasters.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ttf_q_same_product_0bp import (
    read_archives, build_species_densities, relation_pairs,
    historical_unique_fraction, zscore, absorb_two_way
)

EXPECTED_COORDINATES_SHA256='ae5d5b1e192500e7b12b42789b4e8ee7931f439e8c46905b267b0554ef25cef6'
FROZEN_C2={'development':0.29893812458713254,'confirmatory':0.3163113642780213}
VARIABLES=['bio01','bio07','bio12','bio15']


def file_hash(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def evaluate(*,c1:Path,c2:Path,values_file:Path,receipt_file:Path,out:Path):
    receipt=json.loads(receipt_file.read_text())
    if receipt['status']!='SAMPLED_ALL_FOUR_VERIFIED_COG_ASSETS':
        raise ValueError('sample receipt does not attest four verified original COGs')
    if receipt['coordinate_sha256']!=EXPECTED_COORDINATES_SHA256 or receipt['variables']!=VARIABLES:
        raise ValueError('sample coordinate or variable order mismatch')
    if receipt['matrix_sha256'] != file_hash(values_file):
        raise ValueError('sampled climate matrix SHA256 mismatch')
    if len(receipt['assets'])!=4 or any(
      not x.get('sha256') or x.get('variable')!=v for x,v in zip(receipt['assets'],VARIABLES)
    ):
        raise ValueError('incomplete raster provenance')
    m,obs,d,g=read_archives(c1,c2)
    species=d['species_order'].astype(str)
    obs=obs.loc[obs.species.isin(species)].copy().reset_index(drop=True)
    actual_coords=obs[['longitude','latitude']].to_numpy(dtype='<f8')
    if hashlib.sha256(actual_coords.tobytes()).hexdigest()!=EXPECTED_COORDINATES_SHA256:
        raise ValueError('C1 occurrence coordinate bytes changed')
    if not np.array_equal(obs.groupby('species').size().reindex(species).to_numpy(),d['retained_environment_occurrences']):
        raise ValueError('C1 per-species support changed')
    values=np.load(values_file,allow_pickle=False)
    if values.shape!=(len(obs),4) or values.dtype.kind!='f' or not np.isfinite(values).all():
        raise ValueError(f'bad 0BP sampled values {values.shape}, {values.dtype}')
    if receipt.get('matrix_shape')!=[len(obs),4]:raise ValueError('shape mismatch')
    frozen_receipt={x['logical_basename'].lower():x for x in m['assets']}
    for a,v in zip(receipt['assets'],VARIABLES):
        f=f'CHELSA_TraCE21K_{v}_20_V1.0.tif'.lower()
        required=frozen_receipt[f]
        if a['sha256']!=required['sha256'] or a['size_bytes']!=required['size_bytes']:
            raise ValueError('raster SHA/size mismatch against original C1')
    den,details=build_species_densities(obs,species,values)
    outcome={'schema':'ttf_q_trace0_exact_product_comparison_v0.1',
      'status':'COMPLETE_RESPONSE_BLIND_SAME_PRODUCT_COMPARISON',
      'n_species':len(species),'n_occurrences':len(obs),
      'climate_model':'CHELSA-TraCE21k v1.0',
      'historical_relation':'21ka minus 0BP at current occurrence sites',
      'current_comparator':'CHELSA-TraCE21k v1.0 0BP same product',
      'observational_genetic_outcome_used':False,
      'samples_sha256':file_hash(values_file),'receipt_sha256':file_hash(receipt_file),
      'pca':details,'panels':{}}
    for p in FROZEN_C2:
        s=d[p+'_source_index'];t=d[p+'_target_index'];hist=d[p+'_R_hist'];curr=d[p+'_R_current'];geo=g[p+'_coverage']
        if not (np.array_equal(s,g[p+'_source_index']) and np.array_equal(t,g[p+'_target_index'])):
            raise ValueError('C2 exact dyad geometry changed')
        args=(s,t,d['class_name'],d['order_name'],d['family_name'],d['retained_environment_occurrences'])
        c2_base=historical_unique_fraction(hist,curr,geo,*args)
        if abs(c2_base-FROZEN_C2[p])>1e-9:
            raise ValueError(f'{p}: frozen C2 control not reproduced: {c2_base}')
        trace=relation_pairs(den,s,t)
        same_model=historical_unique_fraction(hist,trace,geo,*args)
        n=d['retained_environment_occurrences'];cls=d['class_name'];ordr=d['order_name'];fam=d['family_name']
        X=np.column_stack((zscore(curr),zscore(trace),zscore(geo),
           (cls[s]==cls[t]).astype(float),(ordr[s]==ordr[t]).astype(float),
           (fam[s]==fam[t]).astype(float),zscore(np.abs(np.log(n[s]/n[t])))))
        a=absorb_two_way(np.column_stack((zscore(hist),X)),s,t)
        y=a[:,0]; b=a[:,1:]; residual=y-b@np.linalg.lstsq(b,y,rcond=None)[0]
        dual=float(np.mean(residual**2))
        if dual>min(c2_base,same_model)+1e-9:raise AssertionError('nested control nonmonotonic')
        outcome['panels'][p]={
            'dyads':int(len(s)), 'frozen_original_0bp_v2p1_unique_fraction':c2_base,
            'trace0_same_product_unique_fraction':same_model,
            'both_present_comparators_unique_fraction':dual,
            'old_new_current_relation_pearson':float(np.corrcoef(curr,trace)[0,1]),
            'old_new_current_relation_mean_absolute_difference':float(np.mean(np.abs(curr-trace))),
            'trace0_current_relation_quantiles':np.quantile(trace,[0,.25,.5,.75,1]).tolist(),
        }
    out.mkdir(parents=True,exist_ok=True)
    (out/'ttf_q_trace0_same_product_comparison.json').write_text(json.dumps(outcome,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    print(json.dumps({'panels':outcome['panels'],'status':outcome['status']},indent=2))
    return outcome

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for k in ('c1','c2','samples','receipt','output'):ap.add_argument('--'+k,type=Path,required=True)
    args=ap.parse_args()
    evaluate(c1=args.c1,c2=args.c2,values_file=args.samples,receipt_file=args.receipt,out=args.output)
