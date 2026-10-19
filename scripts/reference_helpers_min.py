"""Source-audited modern reference helper for exploratory LGM climate comparison.
No genetic identities, past host occupancy or biological effect estimation.
"""
import csv,hashlib,json
from pathlib import Path
from decimal import Decimal,InvalidOperation
import numpy as np
from scipy.spatial import cKDTree
FROZEN_MANIFEST_SHA256='a28fba9efa8a9dc79c162336f62c883059c5dddbfabec2861f554d038dd98044'
FROZEN_LOCALITY_SHA256='037cd8fa1f059fb67c349a465540d3d5fac469b5d14a2a4d2658d74c036c0ae9'
FROZEN_EDGE_SHA256='ab10a876895cf00817e8ce555665ebb78a0f2ac64323d2e93c213e1857738ba9'
R=6371.0088

def normalize_id(raw):
    try:d=Decimal(str(raw).strip())
    except InvalidOperation as e:raise ValueError('bad taxon id') from e
    if not d.is_finite() or d<=0 or d!=d.to_integral_value() or d.adjusted()>12:
        raise ValueError('bad WCVP taxon ID')
    return str(int(d))

def csv_rows(p):
    with Path(p).open(encoding='utf8',newline='') as f:return list(csv.DictReader(f))

def sha_file(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
    return h.hexdigest()

def unit_sphere(xy):
    xy=np.asarray(xy,float)
    if xy.ndim!=2 or xy.shape[1]!=2 or not np.isfinite(xy).all() or np.any(abs(xy[:,0])>90) or np.any(abs(xy[:,1])>180):raise ValueError('invalid coords')
    lat=np.deg2rad(xy[:,0]);lon=np.deg2rad(xy[:,1]);c=np.cos(lat)
    return np.column_stack([c*np.cos(lon),c*np.sin(lon),np.sin(lat)])

def minimum_greatcircle_km(xy,support):
    if len(support)==0:raise ValueError('missing eval geography')
    d=cKDTree(unit_sphere(support)).query(unit_sphere(xy),k=1,workers=1)[0]
    return 2*R*np.arcsin(np.clip(d/2,0,1))

def join_hosts(ids,data):
    seen={}
    for host in sorted(set(ids)):
        for raw in data.get(host,[]):
            p=dict(raw);key=str(p['gbif_id'])
            if key in seen:
                old=seen[key]
                if (abs(float(old['latitude'])-float(p['latitude']))>1e-8 or abs(float(old['longitude'])-float(p['longitude']))>1e-8):raise ValueError('conflicting duplicate observation')
                continue
            p['wcvp_id']=host;seen[key]=p
    return list(seen.values())

def filter_host_points(points,support,*,buffer_km=50,thinning_km=10,cap=200):
    if not points:return []
    xy=np.array([[float(p['latitude']),float(p['longitude'])] for p in points]);d=minimum_greatcircle_km(xy,support)
    rows=[]
    for p,dist in zip(points,d):
        u=float(p['uncertainty_m'])
        if not np.isfinite(u) or not 0<=u<=1000:raise ValueError('bad uncertainty')
        if dist-u/1000>=buffer_km-1e-9:
            v=dict(p);v['nearest_evaluation_km']=float(dist);rows.append(v)
    rows.sort(key=lambda p:(hashlib.sha256(f"ttf-host-paired-v0.1|{p['wcvp_id']}|{p['gbif_id']}".encode()).hexdigest(),str(p['gbif_id'])))
    if thinning_km<=0:return rows[:cap]
    selected=[];kept=[];threshold=2*np.sin(thinning_km/(2*R))
    for p in rows:
        v=unit_sphere([[float(p['latitude']),float(p['longitude'])]])[0]
        if kept and np.min(np.linalg.norm(np.asarray(kept)-v,axis=1))<threshold-1e-12:continue
        selected.append(p);kept.append(v)
        if len(selected)>=cap:break
    return selected

def load_host_snapshot(directory):
    directory=Path(directory);plan=json.loads((directory/'host_acquisition_plan.json').read_text());binding=plan['binding']
    completed={};points={};review=set();failed=set()
    for f in sorted((directory/'checkpoints').glob('*.json')):
        q=json.loads(f.read_text());r=q['result'];ident=normalize_id(r['wcvp_id'])
        if q.get('schema')!='ttf_host_observation_checkpoint_v0.1' or q['binding']!=binding or ident!=f.stem:raise ValueError('wrong checkpoint identity')
        sha=hashlib.sha256(json.dumps(r,ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
        if sha!=q['result_sha256'] or ident in completed:raise ValueError('checkpoint hash/duplicate')
        completed[ident]=r
        if r['status']=='OBSERVATIONS_ACQUIRED':
            if r['quality_records']!=len(r['points']):raise ValueError('checkpoint count mismatch')
            points[ident]=r['points']
        elif r['status']=='TAXONOMIC_RANK_REVIEW':review.add(ident)
        else:failed.add(ident)
    return completed,points,review,failed,plan
