#!/usr/bin/env python3
"""Independent data audit of the exploratory host-vs-unlisted-plant LGM null.

The only outcomes here are already-frozen modern/LGM *climatic analogue* scores.
No genetic responses, LGM plant occupancies, new source queries or model tuning.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
import numpy as np

EXPECTED={
  'prior_species':'678c513b015af3514d2d54b5bbeea78f02a73d7c00e70c25276ce27ea0cce0a2',
  'current_species':'591a1383b739d73c8ffc2cbae5c220dae4b7730629e6d447be80542b88124446',
  'current_groups':'236931f93d2efcfebb6988d87291ed6628e63a6ce0ae2cb227cc0825d637cf73',
  'current_summary':'27c432f54918589e6a4e1a373ce0b1c8271486a0ca94552659528a489e60a5f5'
}
def file_hash(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def read_csv(path):
    with Path(path).open(encoding='utf-8',newline='') as f:
        return list(csv.DictReader(f))

def classify_transition(prior,current):
    if set(prior)!=set(current):
        raise ValueError('prior and current comparisons must use the exact same species')
    counts=Counter()
    for species in prior:
        if not np.isfinite(float(prior[species])) or not np.isfinite(float(current[species])):
            raise ValueError('nonfinite species comparison')
        a='positive' if float(prior[species])>0 else 'negative'
        b='positive' if float(current[species])>0 else 'negative'
        counts[a+'_to_'+b]+=1
    return {k:int(counts.get(k,0)) for k in
       ('positive_to_positive','positive_to_negative',
        'negative_to_positive','negative_to_negative')}

def audit_pair_records(rows):
    if not rows or len({r['species'] for r in rows})!=len(rows):
        raise ValueError('species grain is missing/duplicated')
    def v(r,key):
        x=float(r[key])
        if not np.isfinite(x):raise ValueError('nonfinite '+key)
        return x
    return {'species':len(rows),
        'median_actual_fraction':float(median(v(r,'actual_fraction_median') for r in rows)),
        'median_alt_fraction':float(median(v(r,'alternate_fraction_median') for r in rows)),
        'species_median_fraction_difference':float(median(v(r,'fraction_diff_of_medians') for r in rows)),
        'positive_species_fraction_count':sum(v(r,'fraction_diff_of_medians')>0 for r in rows),
        'positive_species_raw_variance_count':sum(v(r,'raw_residual_paired_diff_median')>0 for r in rows),
        'median_raw_residual_difference':float(median(v(r,'raw_residual_paired_diff_median') for r in rows)),
        'median_modern_climate_gap_signed':float(median(v(r,'current_mismatch_gap') for r in rows)),
        'median_modern_climate_gap_abs':float(median(abs(v(r,'current_mismatch_gap')) for r in rows))}

def genus_cluster_interval(rows,*,draws=4000,seed=20261019):
    if draws<2 or not rows:raise ValueError('invalid descriptive bootstrap')
    by_genus=defaultdict(list)
    for r in rows:
        genus=r['species'].split()[0]
        by_genus[genus].append(float(r['fraction_diff_of_medians']))
    groups=[np.asarray(v,float) for _,v in sorted(by_genus.items())]
    count=len(groups);rng=np.random.default_rng(seed)
    medians=np.empty(draws);pos=np.empty(draws)
    for draw in range(draws):
        take=rng.integers(0,count,size=count)
        sample=np.concatenate([groups[k] for k in take])
        medians[draw]=np.median(sample)
        pos[draw]=np.mean(sample>0)
    return {'descriptive_not_confirmatory':True,'genus_clusters':count,
            'draws':draws,'seed':seed,
            'median_diff_limits':np.quantile(medians,[.025,.975]).tolist(),
            'positive_fraction_limits':np.quantile(pos,[.025,.975]).tolist(),
            'shared_donor_across_genera_not_modelled':True}

def audit_sources(v04_dir,v05_dir):
    old_file=Path(v04_dir)/'species_same_niche_null_v0.4.csv'
    new_file=Path(v05_dir)/'species_v05.csv'
    group_file=Path(v05_dir)/'plant_group_v05.csv'
    summary_file=Path(v05_dir)/'summary.json'
    for p,key in ((old_file,'prior_species'),(new_file,'current_species'),
                  (group_file,'current_groups'),(summary_file,'current_summary')):
        if file_hash(p)!=EXPECTED[key]:raise ValueError('frozen source SHA changed: '+str(p))
    old=read_csv(old_file);new=read_csv(new_file);groups=read_csv(group_file)
    summary=json.loads(summary_file.read_text())
    if summary['statuses']!={'COMPLETE_MODERN_MATCHED_PLANT_CONTROL':333,'NONFINITE_CLIMATE':2}:
        raise ValueError('wrong v0.5 status/denominator')
    if summary['source_matched_species']!=335 or summary['climate_completed_species']!=333:
        raise ValueError('wrong v0.5 target qualification')
    prior={r['species']:float(r['diff_of_medians_unique_fraction'])
           for r in old if r['status']=='COMPLETE_SAME_NICHE_NULL'}
    current={r['species']:r for r in new if r['status']=='COMPLETE_MODERN_MATCHED_PLANT_CONTROL'}
    if len(old)!=348 or len(prior)!=347 or len(new)!=335 or len(current)!=333:
        raise ValueError('v04/v05 saved source counts drift')
    matched_old={name:prior[name] for name in current if name in prior}
    if len(matched_old)!=333:raise ValueError('cannot produce unchanged paired333 cohort')
    transition=classify_transition(matched_old,
                                  {k:float(v['fraction_diff_of_medians']) for k,v in current.items()})
    values=list(current.values());metrics=audit_pair_records(values)
    if metrics['positive_species_fraction_count']!=147 or metrics['positive_species_raw_variance_count']!=142:
        raise ValueError('v0.5 source summary does not reconcile')
    if transition['positive_to_negative']!=132 or transition['negative_to_positive']!=10:
        raise ValueError('v04-to-v05 sign-transition counts changed')
    by_species=Counter(g['species'] for g in groups)
    if (len(groups)!=333*8 or set(by_species)!=set(current) or
        any(v!=8 for v in by_species.values())):
        raise ValueError('not eight source-fixed plant comparison groups per species')
    group_vals=defaultdict(list);best_group=[]
    for group in groups:
        current_species=current[group['species']]
        diff=float(group['actual_unique_fraction_median'])-float(group['alternate_unique_fraction_median'])
        if not np.isfinite(diff):raise ValueError('nonfinite plant group comparison')
        group_vals[group['species']].append(diff)
        if int(group['group_index'])==0:best_group.append(diff)
    sensitivity=[]
    for limit in (.005,.01,.02,.05,.1,.2):
        subset=[r for r in values if abs(float(r['current_mismatch_gap']))<=limit]
        sensitivity.append({'posthoc_abs_modern_mismatch_limit':limit,'species':len(subset),
            'actual_fraction_gt_other_plant':sum(float(r['fraction_diff_of_medians'])>0 for r in subset),
            'raw_residual_actual_gt_other_plant':sum(float(r['raw_residual_paired_diff_median'])>0 for r in subset),
            'species_paired_difference_median':float(median(float(r['fraction_diff_of_medians'])
                                                          for r in subset)) if subset else None})
    return {'schema':'ttf_v05_source_bound_independent_host_identity_audit_v0.1',
         'status':'EXPLORATORY_HOST_IDENTITY_ADVANTAGE_NOT_SUPPORTED_BY_UNLISTED_PLANTS',
         'scientific_boundary':'modern climate analogy; unlisted plant != certified nonhost; no LGM occupied plants or genetic response',
         'source_sha256':{k:EXPECTED[k] for k in sorted(EXPECTED)},
         'source_v04_total':348,'source_v04_scored':347,
         'source_v05_total':335,'species_compared':333,'group_rows':len(groups),
         'source_v05_missing_climate':[r['species'] for r in new if r['status']!='COMPLETE_MODERN_MATCHED_PLANT_CONTROL'],
         'source_v04_to_v05_paired':333,
         'v04_actual_positive_on_shared_species':sum(v>0 for v in matched_old.values()),
         'v05_actual_positive_on_shared_species':metrics['positive_species_fraction_count'],
         'v05_raw_positive_on_shared_species':metrics['positive_species_raw_variance_count'],
         'transition_counts':transition,'paired_metrics':metrics,
         'group_0_positive_species':sum(v>0 for v in best_group),
         'group_0_median_fraction_difference':float(median(best_group)),
         'groups_positive_majority_species':sum(sum(d>0 for d in ds)>=5 for ds in group_vals.values()),
         'groups_negative_majority_species':sum(sum(d>0 for d in ds)<=3 for ds in group_vals.values()),
         'posthoc_balance_sensitivity':sensitivity,
         'genus_cluster_descriptive_interval':genus_cluster_interval(values),
         'limitations':['Outcome is historical climate-analogue predictor variance, not genetic differentiation',
           'Alternative plants are unlisted for focal insect in the accepted HOSTS data, not biologically proven nonhosts',
           'Modern climate moment matching is imperfect; distributions, nativity and collection bias remain confounded',
           'Plant donors shared across insect species; even genus-block intervals do not capture all dependencies',
           'Balance windows and taxon-breadth splits are ex post and cannot be confirmatory',
           'Technical 8x8 repeated draws do not create independent species or biological replicates',
           'Original LGM genetic confirmation remains NOT_EVALUABLE_EMPIRICAL_RESPONSE']}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old',type=Path,default=Path('benchmarks/exploratory/lgm_same_niche_split_null_v04'))
    parser.add_argument('--new',type=Path,default=Path('benchmarks/exploratory/lgm_modern_matched_alternative_plant_v05'))
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    data=audit_sources(args.old,args.new)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,sort_keys=True,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:data[k] for k in ('status','species_compared',
         'v04_actual_positive_on_shared_species','v05_actual_positive_on_shared_species',
         'v05_raw_positive_on_shared_species','transition_counts','paired_metrics')},sort_keys=True))
if __name__=='__main__':main()
