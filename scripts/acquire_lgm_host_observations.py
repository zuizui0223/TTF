#!/usr/bin/env python3
"""Acquire the fixed 518 accepted host plants from external observations.

This is pre-spatial-screening retrieval, not matched host/insect climate,
ancient host presence, a representative occurrence census, or genetic evidence.
The original insect query, quality filter and retry helpers are reused unchanged.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
import acquire_lgm_external_references as base

MANIFEST_SHA='a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044'


def plant_key(match: dict, name: str) -> int | None:
    if (match.get('matchType')!='EXACT' or match.get('rank')!='SPECIES'
        or match.get('kingdom')!='Plantae' or match.get('confidence',0)<95
        or match.get('canonicalName','').strip().casefold()!=name.strip().casefold()):
        return None
    key=match.get('acceptedUsageKey') or match.get('usageKey')
    return int(key) if key is not None else None


def host_map(rows: list[dict]) -> dict[str,str]:
    result={}
    for row in rows:
        try: n=Decimal(str(row['wcvp_id_as_frozen']).strip())
        except (InvalidOperation,KeyError) as e: raise ValueError('invalid WCVP ID') from e
        if not n.is_finite() or n<=0 or n!=n.to_integral_value() or n.adjusted()>12:
            raise ValueError('non-positive, non-integral or nonfinite WCVP ID')
        key=str(int(n));name=row['accepted_host_name'].strip()
        if not name or key in result and result[key]!=name: raise ValueError('conflicting host identity')
        result[key]=name
    return dict(sorted(result.items(),key=lambda x:(x[1],x[0])))


def observation_query(key: int) -> dict:
    return dict(taxonKey=key,datasetKey=base.DATASET,basisOfRecord='HUMAN_OBSERVATION',
                occurrenceStatus='PRESENT',hasCoordinate='true',hasGeospatialIssue='false',
                year='2010,2026',eventDate='2010-01-01,2026-10-18',coordinateUncertaintyInMeters='0,1000')


def acquire_host(host_id: str,name: str) -> dict:
    import requests
    session=requests.Session()
    session.headers.update({'User-Agent':'TTF-exploratory-host-reference/0.1 (https://github.com/zuizui0223/TTF)'})
    result=dict(host_id=host_id,host_name=name,status='REQUEST_ERROR',points=[],ledger=[])
    try:
        match=base.get_json(session,'https://api.gbif.org/v1/species/match',dict(name=name,strict='true',rank='SPECIES'),result['ledger'])
        key=plant_key(match,name)
        result.update(taxonomy=match,accepted_species_key=key)
        if key is None: result['status']='TAXONOMY_NOT_EXACT'; return result
        query=observation_query(key)
        first=base.get_json(session,'https://api.gbif.org/v1/occurrence/search',dict(query,limit=0),result['ledger'])
        count=int(first['count']);offsets=base.window_offsets(count)
        result.update(api_count_initial=count,page_offsets=offsets,api_count_changed=False,
                      retrieval_mode='complete_query' if count<=900 else 'bounded_3_windows')
        fetched=0;seen=set();rejections=Counter()
        for offset in offsets:
            page=base.get_json(session,'https://api.gbif.org/v1/occurrence/search',dict(query,offset=offset,limit=300),result['ledger'])
            if int(page.get('offset',offset))!=offset: raise RuntimeError('pagination offset mismatch')
            result['api_count_changed']|=int(page['count'])!=count
            for raw in page.get('results',[]):
                fetched+=1;row,reason=base.clean_record(raw,key)
                if row is None: rejections[reason]+=1;continue
                if row['gbif_id'] in seen: rejections['duplicate_gbif_id']+=1;continue
                seen.add(row['gbif_id'])
                row.update(host_id=host_id,host_name=name,accepted_species_key=key,
                           establishment_means=raw.get('establishmentMeans'),degree_of_establishment=raw.get('degreeOfEstablishment'))
                result['points'].append(row)
        result['points'].sort(key=lambda r:r['gbif_id'])
        result.update(status='RETRIEVAL_COMPLETE',fetched_records=fetched,quality_unique_records=len(result['points']),rejection_counts=dict(rejections))
    except Exception as exc:
        result.update(status='REQUEST_ERROR',error=type(exc).__name__+': '+str(exc)[:350])
        # Partial records are retained for provenance, not counted as completed references.
    finally: session.close()
    return result


def canonical(value) -> bytes:
    return json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def validate_result(host_id: str,result: dict) -> None:
    if result.get('host_id')!=host_id or result.get('status') not in {'REQUEST_ERROR','TAXONOMY_NOT_EXACT','RETRIEVAL_COMPLETE'}:
        raise ValueError('host checkpoint identity/status mismatch')
    if not isinstance(result.get('points'),list) or not isinstance(result.get('ledger'),list):
        raise ValueError('invalid host result fields')
    if any(str(p.get('host_id'))!=host_id for p in result['points']): raise ValueError('point host ID mismatch')


def save_result(directory: Path,host_id: str,binding: dict,result: dict) -> None:
    validate_result(host_id,result)
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    payload=dict(schema='ttf_host_observation_checkpoint_v0.1',binding=binding,result=result,
                 result_sha256=hashlib.sha256(canonical(result)).hexdigest())
    fd,temp=tempfile.mkstemp(prefix='.partial-',dir=directory)
    try:
        with os.fdopen(fd,'wb') as f:
            f.write(gzip.compress(canonical(payload)+b'\n',mtime=0));f.flush();os.fsync(f.fileno())
        os.replace(temp,directory/(host_id+'.json.gz'))
    finally:
        if os.path.exists(temp):os.unlink(temp)


def load_result(directory: Path,host_id: str,binding: dict) -> dict | None:
    path=Path(directory)/(host_id+'.json.gz')
    if not path.exists():return None
    payload=json.loads(gzip.decompress(path.read_bytes()))
    if payload.get('schema')!='ttf_host_observation_checkpoint_v0.1' or payload.get('binding')!=binding:
        raise ValueError('host checkpoint source mismatch')
    result=payload['result'];validate_result(host_id,result)
    if hashlib.sha256(canonical(result)).hexdigest()!=payload.get('result_sha256'):
        raise ValueError('host checkpoint digest mismatch')
    return result


def collect(mapping: dict,directory: Path,binding: dict,acquire,*,workers: int=2,max_seconds: float=900) -> dict:
    if not 1<=workers<=2 or max_seconds<=0:raise ValueError('invalid bounded execution settings')
    done={}
    for key in mapping:
        cached=load_result(directory,key,binding)
        if cached is not None:done[key]=cached
    pending=iter(k for k in mapping if k not in done);start=time.monotonic()
    def worker(k):
        result=acquire(k,mapping[k]);save_result(directory,k,binding,result);return result
    with ThreadPoolExecutor(max_workers=workers) as pool:
        active={}
        def submit():
            if time.monotonic()-start>=max_seconds:return
            key=next(pending,None)
            if key is not None:active[pool.submit(worker,key)]=key
        for _ in range(workers):submit()
        while active:
            finished,_=wait(active,return_when=FIRST_COMPLETED)
            for future in finished:done[active.pop(future)]=future.result()
            if len(done)%20==0:print(json.dumps({'host_attempted':len(done),'host_total':len(mapping)}),flush=True)
            for _ in finished:submit()
    return {k:done[k] for k in mapping if k in done}


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--max-seconds',type=float,default=900)
    args=ap.parse_args()
    if base.digest(args.manifest)!=MANIFEST_SHA:raise RuntimeError('fixed 641-insect host manifest changed')
    with args.manifest.open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
    mapping=host_map(rows)
    if len(rows)!=1037 or len(mapping)!=518 or len({r['insect_species'] for r in rows})!=641:
        raise RuntimeError('fixed host universe changed')
    args.output.mkdir(parents=True,exist_ok=True)
    binding=dict(manifest_sha256=MANIFEST_SHA,host_script_sha256=base.digest(__file__),
                 insect_helper_sha256=base.digest(base.__file__),dataset=base.DATASET)
    done=collect(mapping,args.output/'checkpoints',binding,acquire_host,max_seconds=args.max_seconds)
    details=[];records=[];keys=defaultdict(list)
    for k,name in mapping.items():
        r=done.get(k,dict(host_id=k,host_name=name,status='NOT_ATTEMPTED_IN_CURRENT_SNAPSHOT',points=[],ledger=[]))
        details.append({key:value for key,value in r.items() if key not in {'points','ledger','taxonomy'}})
        if r['status']=='RETRIEVAL_COMPLETE':records.extend(r['points'])
        if r.get('accepted_species_key'):keys[str(r['accepted_species_key'])].append(k)
    (args.output/'host_details.json').write_bytes(canonical(details)+b'\n')
    (args.output/'host_observations_pre_spatial_screen.json.gz').write_bytes(gzip.compress(canonical(records)+b'\n',mtime=0))
    summary=dict(schema='ttf_host_observation_retrieval_v0.1',binding=binding,total_host_species=518,
                 attempted_host_species=len(done),pending_host_species=518-len(done),
                 status='HOST_RETRIEVAL_COMPLETE_WITH_EXPLICIT_FAILURE_STATES' if len(done)==518 else 'HOST_RETRIEVAL_PARTIAL',
                 status_counts=dict(Counter(r['status'] for r in details)),
                 fetched_records=sum(r.get('fetched_records',0) for r in done.values()),
                 usable_records_before_spatial_screen=len(records),
                 host_species_with_records=sum(r.get('quality_unique_records',0)>0 and r['status']=='RETRIEVAL_COMPLETE' for r in done.values()),
                 shared_accepted_taxon_keys={k:v for k,v in keys.items() if len(v)>1},
                 api_count_changed_species=sum(r.get('api_count_changed',False) for r in done.values()),
                 finished_at=datetime.now(timezone.utc).isoformat(),
                 spatial_separation_executed=False,climate_recalculation_executed=False,genetic_response_used=False,
                 boundary='Observed host plants only. Bounded nonrandom retrieval; no native/introduced or cultivation filtering, no proof of local insect use, no time-specific ancient occupancy. Per-insect separation, pooling, thinning and climate reconstruction remain undone.')
    (args.output/'summary.json').write_bytes(canonical(summary)+b'\n')
    print(json.dumps(summary,sort_keys=True),flush=True)
    return 0

if __name__=='__main__':raise SystemExit(main())
