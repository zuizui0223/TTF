#!/usr/bin/env python3
"""External iNaturalist/GBIF reference acquisition; no climate or genetic response.

All 641 previously climate-complete species are attempted under one bounded
retrieval rule. A readiness count is not a biological absence or an inference
opening gate. Spatial separation concerns frozen localities and the original
five internal query points, NOT every possible two-dimensional corridor.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
import inspect
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import date, datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
import numpy as np
from scipy.spatial import cKDTree

DATASET='50c9509d-22c7-4a22-a47d-8c48425ef4a7'
EARTH_KM=6371.0088
BUFFERS=(0,25,50,100)
SOURCE_HASHES={
 'candidates':'c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b',
 'localities':'037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9',
 'edges':'ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9'}
FATAL_ISSUES={'ZERO_COORDINATE','COORDINATE_OUT_OF_RANGE','COORDINATE_INVALID','COUNTRY_COORDINATE_MISMATCH'}
POINT_FIELDS=['species','gbif_id','occurrence_id','dataset_key','license','event_date','latitude','longitude','uncertainty_m','nearest_evaluation_km']


def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,obj): Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
def window_offsets(count):
    if count<0: raise ValueError('negative occurrence count')
    if count==0: return []
    if count<=900: return list(range(0,count,300))
    last=min(count,100000)-300
    return sorted(set((0,last//2,last)))

def unit_vectors(coords):
    a=np.asarray(coords,dtype=float)
    if a.ndim!=2 or a.shape[1]!=2 or not np.isfinite(a).all():
        raise ValueError('finite latitude/longitude matrix required')
    if np.any(np.abs(a[:,0])>90) or np.any(np.abs(a[:,1])>180):
        raise ValueError('coordinate outside bounds')
    lat,lon=np.radians(a).T
    return np.column_stack((np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)))

def nearest_km(coords,support):
    q=unit_vectors(coords); s=unit_vectors(support)
    if not len(s): raise ValueError('empty evaluation support')
    chord=cKDTree(s).query(q,k=1,workers=1)[0]
    return 2*EARTH_KM*np.arcsin(np.clip(chord/2,0,1))

def accepted_key(match,original_name):
    if (match.get('matchType')!='EXACT' or match.get('rank')!='SPECIES'
        or match.get('class')!='Insecta' or match.get('confidence',0)<95
        or match.get('canonicalName','').strip().casefold()!=original_name.strip().casefold()):
        return None
    key=match.get('acceptedUsageKey') or match.get('usageKey')
    return int(key) if key is not None else None

def clean_record(record,key):
    r=record
    if (r.get('datasetKey')!=DATASET or r.get('basisOfRecord')!='HUMAN_OBSERVATION'
        or r.get('occurrenceStatus')!='PRESENT' or r.get('speciesKey')!=key):
        return None,'source_taxon_basis_status'
    if FATAL_ISSUES.intersection(r.get('issues',[])): return None,'coordinate_issue'
    try:
        lat=float(r['decimalLatitude']); lon=float(r['decimalLongitude'])
        unc=float(r['coordinateUncertaintyInMeters']); day=date.fromisoformat(r['eventDate'][:10])
        oid=str(r['occurrenceID']); parsed=urlparse(oid)
        if (not all(math.isfinite(x) for x in (lat,lon,unc)) or abs(lat)>90 or abs(lon)>180
            or (lat==0 and lon==0) or not 0<=unc<=1000): return None,'coordinate_quality'
        if not date(2010,1,1)<=day<=date(2026,10,18): return None,'date'
        if parsed.hostname not in ('www.inaturalist.org','inaturalist.org') or not parsed.path.startswith('/observations/'):
            return None,'observation_identity'
        gbif_id=int(r['key'])
    except (KeyError,TypeError,ValueError,OverflowError): return None,'missing_or_invalid_field'
    return dict(gbif_id=gbif_id,occurrence_id=oid,dataset_key=DATASET,
                license=r.get('license',''),event_date=str(day),latitude=lat,longitude=lon,
                uncertainty_m=unc),'retained'

def thin_records(records,clearance_km,cap=200):
    if clearance_km<0 or cap<1: raise ValueError('invalid separation or cap')
    ordered=sorted(records,key=lambda r:hashlib.sha256(('ttf-external-reference-v01|'+str(r['gbif_id'])).encode()).digest())
    kept=[]; coords=[]; seen_ids=set(); seen_coords=set()
    threshold=2*math.sin(10/(2*EARTH_KM))
    for r in ordered:
        if r['nearest_evaluation_km']-r['uncertainty_m']/1000<clearance_km: continue
        xy=(r['latitude'],r['longitude'])
        if r['occurrence_id'] in seen_ids or xy in seen_coords: continue
        seen_ids.add(r['occurrence_id']); seen_coords.add(xy)
        v=unit_vectors(np.array([xy]))[0]
        if coords and np.min(np.linalg.norm(np.asarray(coords)-v,axis=1))<threshold: continue
        kept.append(r); coords.append(v)
        if len(kept)==cap: break
    return kept

def get_json(session,url,params,ledger):
    import requests
    for attempt in range(3):
        try:
            response=session.get(url,params=params,timeout=(10,35))
            ledger.append(dict(url=response.url,status=response.status_code,attempt=attempt+1,
                               retrieved_at=datetime.now(timezone.utc).isoformat(),
                               response_sha256=hashlib.sha256(response.content).hexdigest()))
            if response.status_code==429 or response.status_code>=500:
                delay=float(response.headers.get('Retry-After',2**(attempt+1)))
                if delay>120: raise RuntimeError('server retry delay exceeds bounded acquisition')
                if attempt<2: time.sleep(max(1,delay)); continue
            response.raise_for_status()
            return response.json()
        except (requests.RequestException,ValueError):
            if attempt==2: raise
            time.sleep(2**(attempt+1))
    raise RuntimeError('request attempts exhausted')

def acquire_one(name,support):
    import requests
    session=requests.Session(); session.headers.update({'User-Agent':'TTF-exploratory-reference/0.1 (https://github.com/zuizui0223/TTF)'})
    ledger=[]; summary=dict(species=name,status='REQUEST_ERROR',external_reference_points_50km=0)
    points=[]; ledger_name=name
    try:
        match=get_json(session,'https://api.gbif.org/v1/species/match',{'name':name,'strict':'true','rank':'SPECIES'},ledger)
        key=accepted_key(match,name)
        summary.update(taxonomy_match=match.get('matchType'),canonical_name=match.get('canonicalName'),
                       accepted_species_key=key,match_confidence=match.get('confidence'))
        if key is None:
            summary['status']='TAXONOMY_NOT_EXACT'; return summary,points,ledger
        query=dict(taxonKey=key,datasetKey=DATASET,basisOfRecord='HUMAN_OBSERVATION',
                   occurrenceStatus='PRESENT',hasCoordinate='true',hasGeospatialIssue='false',
                   year='2010,2026',eventDate='2010-01-01,2026-10-18',coordinateUncertaintyInMeters='0,1000')
        count_response=get_json(session,'https://api.gbif.org/v1/occurrence/search',dict(query,limit=0),ledger)
        count=int(count_response['count']); offsets=window_offsets(count)
        summary.update(api_count_initial=count,retrieval_mode='complete_query' if count<=900 else 'bounded_3_windows',
                       page_offsets=';'.join(map(str,offsets)),api_count_changed=False)
        clean=[]; seen=set(); rejects=Counter(); fetched=0
        for off in offsets:
            page=get_json(session,'https://api.gbif.org/v1/occurrence/search',dict(query,offset=off,limit=300),ledger)
            if int(page.get('offset',off))!=off: raise RuntimeError('pagination offset mismatch')
            summary['api_count_changed']|=int(page['count'])!=count
            for record in page.get('results',[]):
                fetched+=1
                item,reason=clean_record(record,key)
                if item is None: rejects[reason]+=1; continue
                if item['gbif_id'] in seen: rejects['duplicate_gbif_id']+=1; continue
                seen.add(item['gbif_id']); clean.append(item)
        if clean:
            distances=nearest_km(np.array([[r['latitude'],r['longitude']] for r in clean]),support)
            for r,dist in zip(clean,distances): r['nearest_evaluation_km']=float(dist); r['species']=name
        summary.update(fetched_records=fetched,quality_unique_records=len(clean),rejection_counts=dict(rejects))
        for width in BUFFERS:
            eligible=sum(r['nearest_evaluation_km']-r['uncertainty_m']/1000>=width for r in clean)
            selected=thin_records(clean,width)
            summary[f'eligible_before_thinning_{width}km']=eligible
            summary[f'thinned_capped_{width}km']=len(selected)
            if width==50: points=selected
        summary['external_reference_points_50km']=len(points)
        summary['status']='GEOMETRY_REFERENCE_GE30' if len(points)>=30 else 'INSUFFICIENT_IN_RETRIEVED_REFERENCE'
    except Exception as e:
        summary['status']='REQUEST_ERROR';summary['error']=type(e).__name__+': '+str(e)[:350]
        points=[]; summary['external_reference_points_50km']=0
    finally: session.close()
    return summary,points,ledger



def _canonical_json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')


def _validate_saved_result(name, result):
    if not isinstance(result, (tuple, list)) or len(result) != 3:
        raise ValueError('checkpoint result must contain summary, points, ledger')
    summary, points, ledger = result
    statuses = {'REQUEST_ERROR','TAXONOMY_NOT_EXACT','GEOMETRY_REFERENCE_GE30','INSUFFICIENT_IN_RETRIEVED_REFERENCE'}
    if summary.get('species') != name or summary.get('status') not in statuses:
        raise ValueError('checkpoint species or status mismatch')
    if not isinstance(points, list) or not isinstance(ledger, list):
        raise ValueError('invalid checkpoint points or request ledger')
    if any(point.get('species') != name for point in points):
        raise ValueError('checkpoint point species mismatch')
    if summary.get('external_reference_points_50km') != len(points):
        raise ValueError('checkpoint point count mismatch')
    if summary['status'] == 'GEOMETRY_REFERENCE_GE30' and len(points) < 30:
        raise ValueError('checkpoint support status contradicts count')
    return list(result)


def _checkpoint_path(directory, name):
    return Path(directory) / (hashlib.sha256(name.encode('utf-8')).hexdigest()+'.json')


def _save_checkpoint(directory, name, binding, result):
    result = _validate_saved_result(name, result)
    payload = {'schema':'ttf_external_reference_species_checkpoint_v0.1',
               'species':name, 'binding':binding, 'result':result,
               'result_sha256':hashlib.sha256(_canonical_json(result)).hexdigest()}
    target = _checkpoint_path(directory, name)
    fd, temporary = tempfile.mkstemp(prefix='.partial-', dir=directory)
    try:
        with os.fdopen(fd, 'wb') as handle:
            handle.write(_canonical_json(payload)+b'\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)
    return result


def _load_checkpoint(directory, name, binding):
    path = _checkpoint_path(directory, name)
    if not path.exists(): return None
    payload = json.loads(path.read_text(encoding='utf-8'))
    if (payload.get('schema') != 'ttf_external_reference_species_checkpoint_v0.1'
            or payload.get('binding') != binding):
        raise ValueError('checkpoint input binding mismatch')
    if payload.get('species') != name:
        raise ValueError('checkpoint species mismatch')
    if hashlib.sha256(_canonical_json(payload['result'])).hexdigest() != payload.get('result_sha256'):
        raise ValueError('checkpoint result digest mismatch')
    return _validate_saved_result(name, payload['result'])


def checkpointed_acquisition(names, support, directory, binding, acquire, *,
                            workers=4, max_new_species=None, max_seconds=None):
    """Resume only unattempted taxa; successful AND failed attempts are immutable.

    This controls transport, not ecological eligibility. A soft deadline stops
    submitting new taxa and lets the at-most-four active taxa finish their
    original bounded retry rule. Finished species are persisted inside workers,
    even if another worker raises. A corrupted/source-mismatched cache is fatal.
    """
    names = list(names)
    if len(names) != len(set(names)): raise ValueError('duplicate species')
    if not 1 <= workers <= 4: raise ValueError('workers must remain between 1 and 4')
    if max_new_species is not None and max_new_species < 0:
        raise ValueError('negative batch size')
    if max_seconds is not None and max_seconds <= 0: raise ValueError('invalid soft deadline')
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    completed = {}
    for name in names:
        cached = _load_checkpoint(directory, name, binding)
        if cached is not None: completed[name] = cached
    todo = [name for name in names if name not in completed]
    if max_new_species is not None: todo = todo[:max_new_species]
    remaining = iter(todo)
    start = time.monotonic()
    def worker(name):
        result = acquire(name, support[name])
        return _save_checkpoint(directory, name, binding, result)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        active = {}
        def submit_next():
            if max_seconds is not None and time.monotonic()-start >= max_seconds: return False
            name = next(remaining, None)
            if name is None: return False
            active[pool.submit(worker, name)] = name
            return True
        for _ in range(workers): submit_next()
        while active:
            done, _ = wait(active, return_when=FIRST_COMPLETED)
            for future in done:
                name = active.pop(future)
                completed[name] = future.result()
            print(json.dumps({'checkpointed_species':len(completed),'panel_species':len(names),
                              'status_counts':dict(Counter(r[0]['status'] for r in completed.values()))}), flush=True)
            for _ in done: submit_next()
    return [completed[n] for n in names if n in completed], [n for n in names if n not in completed]


def panel_summary_rows(names, completed):
    mapping = {row[0]['species']:row[0] for row in completed}
    return [mapping.get(n, {'species':n,'status':'NOT_ATTEMPTED_IN_CURRENT_SNAPSHOT'}) for n in names]


def acquisition_code_fingerprint():
    # Binding is invariant to checkpoint/export-only changes, not to query,
    # filtering, thinning, retry, or taxonomy changes.
    functions = (window_offsets,unit_vectors,nearest_km,accepted_key,clean_record,thin_records,get_json,acquire_one)
    return hashlib.sha256(('\n'.join(inspect.getsource(f) for f in functions)+
        json.dumps({'dataset':DATASET,'earth_km':EARTH_KM,'buffers':BUFFERS,
                    'fatal_issues':sorted(FATAL_ISSUES)},sort_keys=True)).encode()).hexdigest()


def main():
    ap=argparse.ArgumentParser()
    for name in ('candidates','localities','edges','contract','output'):
        ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--max-new-species',type=int,default=None)
    ap.add_argument('--max-seconds',type=float,default=None)
    args=ap.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    contract=json.loads(args.contract.read_text())
    if (contract['schema']!='ttf_external_insect_reference_v0.1' or contract['buffers_km']!=list(BUFFERS)
        or contract['dataset_key']!=DATASET or contract['minimum_descriptive_reference_n']!=30):
        raise RuntimeError('acquisition contract drift')
    for k,h in SOURCE_HASHES.items():
        if digest(getattr(args,k))!=h: raise RuntimeError('source hash mismatch: '+k)
    def read(p):
        with p.open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))
    candidates=read(args.candidates)
    if len(candidates)!=642: raise RuntimeError('expected original 642 candidates')
    panel=[r for r in candidates if r['species'].strip()!='Polyommatus coridon']
    names=sorted(r['species'].strip() for r in panel)
    names_sha=hashlib.sha256(('\n'.join(names)+'\n').encode()).hexdigest()
    if len(names)!=641 or len(set(names))!=641 or names_sha!=contract['641_species_sha256']:
        raise RuntimeError('641-species ecological domain differs')
    # Reuse the original, frozen five INTERNAL great-circle quadrature points.
    module_path=Path(__file__).parents[1]/'src/ttf/historical_host_memory.py'
    spec=importlib.util.spec_from_file_location('original_hhm',module_path)
    original=importlib.util.module_from_spec(spec); spec.loader.exec_module(original)
    local=defaultdict(dict); edges=defaultdict(list)
    for r in read(args.localities): local[r['species']][int(r['locality_index'])]=(float(r['latitude']),float(r['longitude']))
    for r in read(args.edges): edges[r['species']].append((int(r['node_left']),int(r['node_right'])))
    support={}
    for name in names:
        loc=local[name]; pairs=edges[name]
        quad=original.great_circle_quadrature_latlon(np.array([loc[i] for i,j in pairs]),np.array([loc[j] for i,j in pairs]),segment_points=5)
        support[name]=np.vstack((np.array(list(loc.values())),quad.reshape(-1,2)))
    host_manifest=[]
    for r in panel:
        hs=r['accepted_host_names'].split(';'); ids=r['accepted_host_ids'].split(';')
        if len(hs)!=len(ids): raise RuntimeError('unaligned accepted-host manifest')
        for name,host_id in zip(hs,ids): host_manifest.append({'insect_species':r['species'].strip(),'accepted_host_name':name,'wcvp_id_as_frozen':host_id})
    with (args.output/'host_reference_request_manifest.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(host_manifest[0]));w.writeheader();w.writerows(host_manifest)
    binding=dict(contract_sha256=digest(args.contract), source_hashes=SOURCE_HASHES,
                 species_sha256=names_sha,original_quadrature_module_sha256=digest(module_path),
                 acquisition_code_sha256=acquisition_code_fingerprint())
    completed,pending=checkpointed_acquisition(names,support,args.output/'checkpoints',binding,acquire_one,
                         workers=4,max_new_species=args.max_new_species,max_seconds=args.max_seconds)
    results=panel_summary_rows(names,completed); all_points=[];all_ledger=[]
    for item,points,ledger in completed:
        all_points.extend(points)
        all_ledger.extend({'species':item['species'],**r} for r in ledger)
    results.sort(key=lambda r:r['species']);all_points.sort(key=lambda r:(r['species'],r['gbif_id']))
    fields=sorted(set().union(*(r.keys() for r in results)))
    with (args.output/'all_641_species_external_reference.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(results)
    with gzip.open(args.output/'external_insect_points_50km.csv.gz','wt',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=POINT_FIELDS);w.writeheader();w.writerows(all_points)
    dump(args.output/'species_details.json',results)
    dump(args.output/'request_ledger.json',all_ledger)
    key_groups=defaultdict(list)
    for r in results:
        if r.get('accepted_species_key'): key_groups[str(r['accepted_species_key'])].append(r['species'])
    summary=dict(schema='ttf_external_insect_reference_result_v0.1',status='ACQUISITION_PARTIAL_WITH_CHECKPOINTS' if pending else 'ACQUISITION_COMPLETE_WITH_EXPLICIT_FAILURE_STATES',
        pending_species=len(pending),pending_species_names=pending,
        checkpoint_binding=binding,
        attempted_species=len(completed),original_candidate_species=642,previous_nonfinite_species_excluded='Polyommatus coridon',
        species_sha256=names_sha,contract_sha256=digest(args.contract),source_hashes=SOURCE_HASHES,
        original_quadrature_module_sha256=digest(module_path),dataset_key=DATASET,
        completed_at=datetime.now(timezone.utc).isoformat(),status_counts=dict(Counter(r['status'] for r in results)),
        retained_reference_points_50km=len(all_points),quality_unique_records=sum(r.get('quality_unique_records',0) for r in results),
        fetched_records=sum(r.get('fetched_records',0) for r in results),
        species_ge30_by_buffer={str(b):sum(r.get(f'thinned_capped_{b}km',0)>=30 for r in results) for b in BUFFERS},
        shared_accepted_taxon_keys={k:v for k,v in key_groups.items() if len(v)>1},
        retrieval_modes=dict(Counter(r.get('retrieval_mode','unresolved') for r in results)),
        api_index_changed_species=sum(r.get('api_count_changed',False) for r in results),
        host_reference_manifest_rows=len(host_manifest),host_reference_acquisition_executed=False,
        real_climate_recalculation_executed=False,genetic_response_used=False,
        genetic_source_disjointness_proven=False,
        interpretation='New external observation references, spatially separated from original localities and five evaluated internal positions. Not full-route spatial separation, source-disjoint genetic validation, matched host niche coverage, historical occupancy, or a genetic effect.')
    dump(args.output/'summary.json',summary)
    print(json.dumps(summary,ensure_ascii=False,sort_keys=True),flush=True)
    return 0

if __name__=='__main__': raise SystemExit(main())
