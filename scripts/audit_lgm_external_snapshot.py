#!/usr/bin/env python3
"""Independent audit of saved external-reference geometry; never reads genomes."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path
import numpy as np

RADIUS_KM = 6371.0088

def nearest_haversine(query, reference):
    q=np.asarray(query,float); r=np.asarray(reference,float)
    if q.ndim!=2 or r.ndim!=2 or q.shape[1]!=2 or r.shape[1]!=2 or len(r)==0:
        raise ValueError('nonempty lat/lon reference required')
    if not np.isfinite(q).all() or not np.isfinite(r).all():raise ValueError('nonfinite coordinates')
    r=np.radians(r);result=[]
    for start in range(0,len(q),64):
        b=np.radians(q[start:start+64]);dlat=b[:,0,None]-r[None,:,0];dlon=b[:,1,None]-r[None,:,1]
        h=np.sin(dlat/2)**2+np.cos(b[:,0,None])*np.cos(r[None,:,0])*np.sin(dlon/2)**2
        distance=2*RADIUS_KM*np.arcsin(np.sqrt(np.clip(h,0,1)))
        result.extend(np.min(distance,axis=1).tolist())
    return np.asarray(result)

def audit_points(records,support):
    groups=defaultdict(list)
    for r in records:groups[r['species']].append(r)
    error=0.;clearances=[];min_spacing=[]
    for name,rows in groups.items():
        if name not in support:raise ValueError('point species outside panel')
        if len(rows)>200:raise ValueError('reference cap violated')
        if len({str(r['gbif_id']) for r in rows})!=len(rows):raise ValueError('duplicate observation')
        xy=np.array([[float(r['latitude']),float(r['longitude'])] for r in rows])
        if np.any(np.abs(xy[:,0])>90) or np.any(np.abs(xy[:,1])>180):raise ValueError('coordinate range')
        distances=nearest_haversine(xy,support[name])
        uncertainty=np.array([float(r['uncertainty_m']) for r in rows])
        if not np.isfinite(uncertainty).all() or np.any((uncertainty<0)|(uncertainty>1000)):raise ValueError('uncertainty')
        clearance=distances-uncertainty/1000
        if np.min(clearance)<50-1e-6:raise ValueError('50-km separation violation')
        stored=np.array([float(r['nearest_evaluation_km']) for r in rows])
        if not np.isfinite(stored).all():raise ValueError('nonfinite stored distance')
        e=float(np.max(np.abs(stored-distances)))
        if e>1e-6:raise ValueError('stored nearest distance mismatch')
        error=max(error,e);clearances.extend(clearance.tolist())
        if len(xy)>1:
            spacing=min(float(nearest_haversine(xy[i:i+1],np.delete(xy,i,axis=0))[0]) for i in range(len(xy)))
            if spacing<10-1e-6:raise ValueError('10-km thinning violation')
            min_spacing.append(spacing)
    return dict(retained_points=len(records),species_with_retained_points=len(groups),
        maximum_independent_distance_error_km=error,
        minimum_verified_clearance_km=min(clearances) if clearances else None,
        minimum_verified_within_species_spacing_km=min(min_spacing) if min_spacing else None,
        counts_by_species={k:len(v) for k,v in sorted(groups.items())})

def read_csv(p):
    with Path(p).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))

def main():
    ap=argparse.ArgumentParser()
    for name in ('snapshot','localities','edges','quadrature-module','output'):
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args();summary=json.loads((args.snapshot/'summary.json').read_text())
    for key in ('localities','edges'):
        path=getattr(args,key)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=summary['source_hashes'][key]:raise ValueError('source geometry hash drift')
    if hashlib.sha256(args.quadrature_module.read_bytes()).hexdigest()!=summary['original_quadrature_module_sha256']:raise ValueError('quadrature hash drift')
    spec=importlib.util.spec_from_file_location('audit_original_quadrature',args.quadrature_module)
    original=importlib.util.module_from_spec(spec);spec.loader.exec_module(original)
    local=defaultdict(dict);edges=defaultdict(list)
    for r in read_csv(args.localities):local[r['species']][int(r['locality_index'])]=(float(r['latitude']),float(r['longitude']))
    for r in read_csv(args.edges):edges[r['species']].append((int(r['node_left']),int(r['node_right'])))
    rows=read_csv(args.snapshot/'all_641_species_external_reference.csv')
    names=[r['species'] for r in rows]
    if len(names)!=641 or len(set(names))!=641:raise ValueError('all-641 snapshot identity invalid')
    if hashlib.sha256(('\n'.join(sorted(names))+'\n').encode()).hexdigest()!=summary['species_sha256']:raise ValueError('panel name hash drift')
    support={}
    for name in names:
        loc=local[name];pairs=edges[name]
        quad=original.great_circle_quadrature_latlon(np.array([loc[i] for i,j in pairs]),np.array([loc[j] for i,j in pairs]),segment_points=5)
        support[name]=np.vstack((list(loc.values()),quad.reshape(-1,2)))
    with gzip.open(args.snapshot/'external_insect_points_50km.csv.gz','rt',newline='',encoding='utf-8') as f:points=list(csv.DictReader(f))
    audit=audit_points(points,support)
    if len(points)!=summary['retained_reference_points_50km']:raise ValueError('summary point count mismatch')
    for row in rows:
        n=audit['counts_by_species'].get(row['species'],0)
        if row['status']=='NOT_ATTEMPTED_IN_CURRENT_SNAPSHOT':
            if row.get('external_reference_points_50km') or n:raise ValueError('pending species has observations')
        elif n!=int(row['external_reference_points_50km']):raise ValueError('species point count mismatch')
        if (n>=30)!=(row['status']=='GEOMETRY_REFERENCE_GE30'):raise ValueError('support label mismatch')
    counts=dict(Counter(r['status'] for r in rows))
    if counts!=summary['status_counts']:raise ValueError('status count mismatch')
    hosts=read_csv(args.snapshot/'host_reference_request_manifest.csv')
    host_ids=set();host_names=set();insects=defaultdict(set)
    for r in hosts:
        host_id=Decimal(r['wcvp_id_as_frozen'])
        if not host_id.is_finite() or host_id!=host_id.to_integral_value():raise ValueError('host ID not integral')
        host_ids.add(str(int(host_id)));host_names.add(r['accepted_host_name']);insects[r['insect_species']].add(str(int(host_id)))
    audit.update(schema='ttf_external_reference_independent_snapshot_audit_v0.1',
        status='PASS_SAVED_SNAPSHOT_GEOMETRY_ONLY',panel_species=641,
        attempted_species=summary['attempted_species'],pending_species=summary.get('pending_species',0),
        status_counts=counts,host_association_rows=len(hosts),unique_accepted_host_ids=len(host_ids),
        unique_accepted_host_names=len(host_names),host_manifest_insect_species=len(insects),
        snapshot_summary_sha256=hashlib.sha256((args.snapshot/'summary.json').read_bytes()).hexdigest(),
        genetic_response_used=False,climate_recalculation_executed=False,
        boundary='Partial acquisition is not a representative species sample. Haversine distances independently verify separation from discrete evaluated points, not every continuous route or genetic-source disjointness. Host taxa are a request manifest only, not acquired observations.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(audit,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='counts_by_species'},sort_keys=True))

if __name__=='__main__':main()
