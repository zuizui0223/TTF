"""Pure functions for a modern-climate-matched *unlisted plant* negative control.

All matching is computed from current 0BP plant climate only, never LGM.
An unlisted plant is not necessarily biologically unused by the insect.
"""
import hashlib
import numpy as np

def genus_of(name: str)->str:
    pieces=name.strip().split()
    if len(pieces)<2:return ""
    return pieces[0].casefold()

def stable_sample_indices(observations,cap=64):
    if cap<1:raise ValueError('nonpositive cap')
    hashes=[]
    for i,p in enumerate(observations):
        if 'gbif_id' not in p:raise ValueError('missing GBIF id')
        h=hashlib.sha256(f"ttf-v05-plant-signature|{p['gbif_id']}".encode()).hexdigest()
        hashes.append((h,str(p['gbif_id']),i))
    return [r[2] for r in sorted(hashes)[:cap]]

def climate_signature(x):
    a=np.asarray(x,dtype=float)
    if a.ndim!=2 or a.shape[1]!=4 or len(a)<30 or not np.isfinite(a).all():
        return None
    return {'center':a.mean(axis=0).tolist(),
            'spread':a.std(axis=0).tolist(),
            'n':len(a)}

def distance_to_signature(actual, candidate):
    if actual is None or candidate is None: return float('inf')
    a=np.array(actual['center']);b=np.array(candidate['center'])
    sa=np.array(actual['spread']);sb=np.array(candidate['spread'])
    if len(a)!=4 or len(b)!=4 or not np.isfinite([*a,*b,*sa,*sb]).all():
        return float('inf')
    return float(np.linalg.norm(a-b)+.5*np.linalg.norm(sa-sb))

def rank_alternatives(actual_signature,known_host_ids,all_taxa,signatures):
    host_ids=set(known_host_ids)
    host_genera={genus_of(all_taxa[x]['name']) for x in host_ids}
    result=[]
    for key,taxon in all_taxa.items():
        if key in host_ids or genus_of(taxon['name']) in host_genera:
            continue
        if key not in signatures:continue
        score=distance_to_signature(actual_signature,signatures[key])
        if np.isfinite(score):result.append((score,key))
    return sorted(result,key=lambda x:(x[0],x[1]))

def group_alternatives(eligible,actual_host_count,max_group=8,max_taxa=24):
    if not 1<=actual_host_count<=3:raise ValueError('unsupported source host richness')
    groups=[];distinct=[]
    for item in eligible:
        if item not in distinct:distinct.append(item)
        if len(distinct)>=max_taxa:break
    for start in range(0,len(distinct),actual_host_count):
        group=distinct[start:start+actual_host_count]
        if len(group)==actual_host_count:groups.append(group)
        if len(groups)>=max_group:break
    return groups

def observation_signature(points):
    ids=[str(v['gbif_id']) for v in points]
    if len(ids)!=len(set(ids)):raise ValueError('duplicate reference observation')
    return hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()

def equal_reference_indices(n_actual,n_alternative,n_insect,name,group,draw):
    if min(n_actual,n_alternative)<30 or n_insect<60:
        raise ValueError('insufficient independent sampling')
    k=min(100,n_actual,n_alternative,n_insect//2)
    seed=int.from_bytes(hashlib.sha256(
        f'ttf-v05-real-vs-unlistedplant|{name}|{group}|{draw}'.encode()).digest()[:8],'big')
    rng=np.random.default_rng(seed)
    i=np.sort(rng.choice(n_insect,size=k,replace=False))
    r=np.sort(rng.choice(n_actual,size=k,replace=False))
    c=np.sort(rng.choice(n_alternative,size=k,replace=False))
    return i,r,c
