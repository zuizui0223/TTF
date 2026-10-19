"""Exploratory response-blind same-insect-niche sampling-noise reference null.

NOT a null model for host association, LGM plant occupation, dispersal or genetics.
All actual-host, self and disjoint pseudo-host clouds are equal-size.
"""
from __future__ import annotations
import hashlib
import numpy as np
from scipy.spatial import cKDTree

def draw_independent_clouds(insect_n: int, host_n: int, species: str, draw: int):
    if insect_n < 60 or host_n < 30 or draw < 0 or not species:
        raise ValueError('requires >=60 external insect and >=30 host points')
    k=min(100,insect_n//2,host_n)
    seed=int.from_bytes(hashlib.sha256(
        f'ttf-v04-same-niche-split|{species}|{draw}'.encode('utf-8')).digest()[:8],'big')
    rng=np.random.default_rng(seed)
    perm=rng.permutation(insect_n)
    self_idx=np.sort(perm[:k])
    pseudo_idx=np.sort(perm[k:2*k])
    host_idx=np.sort(rng.choice(host_n,size=k,replace=False))
    if set(self_idx).intersection(pseudo_idx):
        raise RuntimeError('insect references overlapped')
    return host_idx,self_idx,pseudo_idx

def _distance_score(query: np.ndarray,cloud: np.ndarray)->np.ndarray:
    dims=query.shape[-1]
    if query.ndim!=3 or cloud.ndim!=2 or cloud.shape[1]!=dims:
        raise ValueError('invalid climate cube/cloud dimensions')
    if not np.isfinite(query).all() or not np.isfinite(cloud).all() or len(cloud)<30:
        raise ValueError('nonfinite or inadequate climate reference')
    n,p,_=query.shape
    result,_=cKDTree(cloud).query(query.reshape(-1,dims),k=1,workers=1)
    return result.reshape(n,p).mean(axis=1)

def score_analog_history(current:np.ndarray,lgm:np.ndarray,host:np.ndarray,
                         insect_self:np.ndarray,geo:np.ndarray)->dict:
    current=np.asarray(current,float);lgm=np.asarray(lgm,float)
    host=np.asarray(host,float);insect_self=np.asarray(insect_self,float)
    geo=np.asarray(geo,float)
    if current.ndim!=3 or current.shape!=lgm.shape or current.shape[1:]!=(5,4):
        raise ValueError('expected edge x 5 point x 4 climate cube')
    if geo.shape!=(len(current),) or not np.isfinite(geo).all():
        raise ValueError('invalid geographic controls')
    hc=_distance_score(current,host);hl=_distance_score(lgm,host)
    sc=_distance_score(current,insect_self);sl=_distance_score(lgm,insect_self)
    hm=hl-hc;sm=sl-sc
    design=np.column_stack([np.ones(len(current)),sm,hc,sc,geo])
    output={'history_variance':float(np.var(hm)),
            'host_delta_mean':float(hm.mean()),
            'host_current_mean':float(hc.mean()),
            'self_delta_mean':float(sm.mean()),
            'unique_fraction':None,'residual_variance':None,'identifiable':False}
    if (len(current)<=design.shape[1]+1 or not np.isfinite(design).all() or
        output['history_variance']<=1e-12 or
        np.linalg.matrix_rank(design)<design.shape[1]):
        return output
    beta=np.linalg.lstsq(design,hm,rcond=None)[0]
    resid=float(np.var(hm-design@beta))
    output.update(unique_fraction=resid/output['history_variance'],
                  residual_variance=resid,identifiable=True)
    return output

def score_real_and_placebo(current:np.ndarray,lgm:np.ndarray,actual_host:np.ndarray,
                           self_cloud:np.ndarray,placebo_host:np.ndarray,geo:np.ndarray):
    if len(actual_host)!=len(self_cloud) or len(self_cloud)!=len(placebo_host):
        raise ValueError('all three clouds must have equal reference count')
    return (score_analog_history(current,lgm,actual_host,self_cloud,geo),
            score_analog_history(current,lgm,placebo_host,self_cloud,geo))
