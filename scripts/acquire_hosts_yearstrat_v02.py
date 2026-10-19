#!/usr/bin/env python3
"""Exploratory v0.2: zero-offset three-year-bin sampling of modern plant observations.

Never mix with original 3-deep-offset host v0.1 observations or treat modern
plant records as ancient host presence. No genetic or climate responses opened.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,os
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import requests
import acquire_lgm_external_references as insect
import acquire_lgm_host_references as host

BINS=((2010,2015),(2016,2020),(2021,2026))
SCHEME='THREE_NONOVERLAPPING_YEAR_BINS_SHALLOW_OFFSET_ZERO_V02'
EXPECTED_MANIFEST='a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()

def make_query(key,start,end):
    return {'taxonKey':int(key),'datasetKey':insect.DATASET,
        'basisOfRecord':'HUMAN_OBSERVATION','occurrenceStatus':'PRESENT',
        'hasCoordinate':'true','hasGeospatialIssue':'false',
        'year':f'{start},{end}',
        'eventDate':f'{start}-01-01,{end}-12-31' if end<2026 else f'{start}-01-01,2026-10-18',
        'coordinateUncertaintyInMeters':'0,1000','limit':300,'offset':0}

def acquire_binned(taxon,*,fetch=insect.get_json):
    result={'wcvp_id':taxon['wcvp_id'],'name':taxon['name'],
        'rank_hint':taxon['rank_hint'],'scheme':SCHEME,
        'status':'REQUEST_ERROR','quality_records':None,'points':[],
        'request_ledger':[],'native_range_filter_applied':False,
        'spatial_filter_applied':False,'climate_or_genetic_response_used':False}
    rank=taxon['rank_hint']
    if rank not in {'SPECIES','SUBSPECIES'}:
        result['status']='TAXONOMIC_RANK_REVIEW'
        return result
    with requests.Session() as session:
        session.headers.update({'User-Agent':'TTF-shallow-binned-host-v02 (https://github.com/zuizui0223/TTF)'})
        ledger=result['request_ledger']
        try:
            match=fetch(session,'https://api.gbif.org/v1/species/match',
                {'name':taxon['name'],'strict':'true','rank':rank,'kingdom':'Plantae'},ledger)
            result['taxonomy_match']=match
            key=host.plant_key(match,taxon['name'],rank)
            result['accepted_taxon_key']=key
            if key is None:
                result['status']='TAXONOMY_NOT_EXACT'
                return result
            found=[];seen=set();rejected=Counter();bins=[]
            for start,end in BINS:
                query=make_query(key,start,end)
                page=fetch(session,'https://api.gbif.org/v1/occurrence/search',query,ledger)
                n=int(page['count']);rows=page['results']
                if int(page.get('offset',0))!=0:
                    raise ValueError('nonzero returned offset')
                if len(rows)<min(300,n):
                    raise ValueError('short shallow search page')
                bins.append({'years':[start,end],'api_count':n,'retrieved':len(rows),'offset':0})
                for source in rows:
                    item,reason=host.clean_host_record(source,key,rank)
                    if item is None:
                        rejected[reason]+=1
                        continue
                    id_=str(item['gbif_id'])
                    if id_ in seen:
                        rejected['duplicate_gbif_id']+=1
                        continue
                    seen.add(id_)
                    found.append(item)
            result.update(status='OBSERVATIONS_ACQUIRED',
                quality_records=len(found),points=found,bin_census=bins,
                rejection_counts=dict(rejected))
        except (requests.RequestException,ValueError,KeyError,RuntimeError,OSError) as exc:
            result.update(status='REQUEST_ERROR',quality_records=None,
                points=[],error=type(exc).__name__+': '+str(exc)[:400])
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',required=True,type=Path)
    ap.add_argument('--original-plan',required=True,type=Path)
    ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--max-seconds',type=float,default=1600)
    ap.add_argument('--execute',action='store_true')
    a=ap.parse_args()
    if sha(a.manifest)!=EXPECTED_MANIFEST:
        raise ValueError('original accepted-host relationships changed')
    q=json.loads(a.original_plan.read_text())
    if q['binding']['manifest_sha256']!=EXPECTED_MANIFEST:
        raise ValueError('original host acquisition plan binding changed')
    plan=q['taxa']
    if len(plan)!=518 or sum(int(t['linked_insect_count']) for t in plan)!=1037:
        raise ValueError('full 518 taxonomy plan drift')
    a.out.mkdir(parents=True,exist_ok=True)
    binding={'schema':'ttf_external_host_yearstrat_v0.2',
        'new_code_sha256':sha(Path(__file__)),
        'original_host_plan_sha256':sha(a.original_plan),
        'original_host_script_sha256':sha(Path(host.__file__)),
        'original_insect_filter_sha256':sha(Path(insect.__file__)),
        'manifest_sha256':EXPECTED_MANIFEST,
        'dataset_key':insect.DATASET,
        'bins':BINS,'max_per_bin':300,'offset':0,
        'not_equivalent_to_v01':True}
    (a.out/'source_bound_plan.json').write_text(json.dumps(
        {'binding':binding,'taxa':plan},indent=2,sort_keys=True)+'\n')
    if not a.execute:
        print('PLAN_ONLY_518_YEAR_BIN_TAXA')
        return 0
    complete,pending=host.run_queue(plan,a.out/'checkpoints',binding,
        acquire=acquire_binned,max_seconds=a.max_seconds,workers=4)
    result_by_id={r['wcvp_id']:r for r in complete}
    allrows=[{'wcvp_id':t['wcvp_id'],'name':t['name'],
        'rank_hint':t['rank_hint'],'linked_insect_count':t['linked_insect_count'],
        'status':result_by_id.get(t['wcvp_id'],{}).get('status','NOT_ATTEMPTED'),
        'quality_records':result_by_id.get(t['wcvp_id'],{}).get('quality_records')}
        for t in plan]
    with (a.out/'all_518_host_taxa_v0.2.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]))
        w.writeheader();w.writerows(allrows)
    counts=dict(Counter(r['status'] for r in allrows))
    summary={'schema':'ttf_external_host_yearstrat_v0.2_result',
        'status':'PARTIAL' if pending else 'ALL_TAXA_ATTEMPTED',
        'attempted_taxa':len(complete),'pending_taxa':len(pending),
        'pending_ids':pending,'status_counts':counts,
        'quality_records':sum(len(r['points']) for r in complete),
        'binding':binding,'old_v01_result_mixed':False,
        'native_range_filter_applied':False,'spatial_filter_applied':False,
        'climate_or_genetic_response_used':False,
        'completed_at':datetime.now(timezone.utc).isoformat()}
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in {'pending_ids','binding'}},sort_keys=True))
    return 0
if __name__=='__main__':raise SystemExit(main())
