"""Synthetic spatial-null calibration, species-heldout AUC. Never fit real lineage breaks."""
from __future__ import annotations
import numpy as np
from sklearn.linear_model import LogisticRegression
from scipy.stats import rankdata
BASELINE_COEFF=np.array([.30,-.14,.12,-.09,-.22,.18,-.13,.12,.12,.22,-.18])

def _moments(values,species,role):
    ids=np.unique(species[role==0])
    if len(ids)<2:raise ValueError("at least two training species required")
    first=np.asarray([values[species==s].mean(axis=0) for s in ids])
    second=np.asarray([(values[species==s]**2).mean(axis=0) for s in ids])
    mu=first.mean(axis=0)
    sd=np.sqrt(np.maximum(second.mean(axis=0)-mu*mu,0))
    return mu,np.where(sd>1e-10,sd,1.)

def standardize_predictors(climate_barrier,lonlat,species,role):
    x=np.asarray(climate_barrier,float);xy=np.asarray(lonlat,float)
    s=np.asarray(species);role=np.asarray(role)
    if x.ndim!=2 or x.shape[1]!=13 or xy.shape!=(len(x),2) or len(s)!=len(x) or len(role)!=len(x):
        raise ValueError("aligned 13 climate/barrier fields + geographic sites required")
    if not np.isfinite(x).all() or not np.isfinite(xy).all() or (x[:,8:13]<0).any():
        raise ValueError("nonfinite climate/barrier or negative geographic distance")
    _,climate_sd=_moments(x[:,4:8],s,role)
    hist=np.linalg.norm((x[:,:4]-x[:,4:8])/climate_sd,axis=1)
    raw=np.column_stack([x[:,4:8],np.sqrt(x[:,8:12]),np.log1p(x[:,12]),
                         xy[:,1],xy[:,0],hist])
    mu,scale=_moments(raw,s,role)
    return (raw-mu)/scale,mu,scale,climate_sd

def training_history_unique_fraction(p,species,role):
    p=np.asarray(p,float);s=np.asarray(species);r=np.asarray(role)
    X=[];Y=[]
    for sid in np.unique(s[r==0]):
        a=p[s==sid];w=1/np.sqrt(len(a))
        X.append(a[:,:-1]*w);Y.append(a[:,-1]*w)
    X=np.vstack(X);Y=np.concatenate(Y)
    norm=float(Y@Y)
    if norm<=1e-12:raise ValueError("historical feature has zero training variance")
    resid=Y-X@np.linalg.lstsq(X,Y,rcond=1e-12)[0]
    return float((resid@resid)/norm)

def _spatial_effect(lonlat,rng):
    xy=np.asarray(lonlat,float);ctr=xy.mean(axis=0)
    km=np.column_stack([(xy[:,0]-ctr[0])*111.2*np.cos(np.deg2rad(ctr[1])),
                        (xy[:,1]-ctr[1])*111.2])
    field=np.zeros(len(xy),float)
    for length in (100.,250.,550.):
        angle=rng.uniform(-np.pi,np.pi);phase=rng.uniform(0,2*np.pi)
        direction=np.array([np.cos(angle),np.sin(angle)])
        field+=rng.normal()*np.cos(2*np.pi*(km@direction)/length+phase)
    std=float(field.std())
    return (field-field.mean())/std if std>1e-10 else np.zeros(len(xy))

def draw_synthetic_labels(p,lonlat,species,case_counts,rng,history_log_odds=0.,spatial_sigma=.75):
    p=np.asarray(p,float);sp=np.asarray(species);xy=np.asarray(lonlat,float)
    if p.ndim!=2 or p.shape[1]!=12 or len(sp)!=len(p) or len(xy)!=len(p):
        raise ValueError("12 locked standardized predictors required")
    y=np.zeros(len(p),dtype='u1')
    for sid in np.unique(sp):
        indices=np.flatnonzero(sp==sid);k=case_counts.get(int(sid))
        if k is None or not (5<=k<len(indices)):raise ValueError("bad frozen case count")
        odds=p[indices,:11]@BASELINE_COEFF+history_log_odds*p[indices,11]
        if spatial_sigma:odds+=spatial_sigma*_spatial_effect(xy[indices],rng)
        odds+=rng.gumbel(size=len(indices))
        selected=np.argpartition(odds,-k)[-k:]
        y[indices[selected]]=1
    return y

def species_auc(y,score):
    y=np.asarray(y,dtype='u1');v=np.asarray(score,float)
    n1=int(y.sum());n0=len(y)-n1
    if n0<=0 or n1<=0 or len(y)!=len(v) or not np.isfinite(v).all():
        raise ValueError("AUC needs finite scores and both classes")
    return float((rankdata(v,method='average')[y==1].sum()-n1*(n1+1)/2)/(n1*n0))

def evaluate_species_holdout(p,y,species,role,bootstrap_reps=199,seed=20261029):
    p=np.asarray(p,float);y=np.asarray(y,dtype='u1')
    sp=np.asarray(species);role=np.asarray(role)
    train_ids=np.unique(sp[role==0]);test_ids=np.unique(sp[role==1])
    if not len(train_ids) or not len(test_ids) or len(np.intersect1d(train_ids,test_ids)):
        raise ValueError("species-disjoint train/eval groups required")
    w=np.zeros(len(y),float)
    for sid in train_ids:
        ii=np.flatnonzero(sp==sid);a=ii[y[ii]==1];b=ii[y[ii]==0]
        if not len(a) or not len(b):raise ValueError("a species lacks a class")
        w[a]=.5/len(a);w[b]=.5/len(b)
    models=[]
    for n in (11,12):
        model=LogisticRegression(C=2.,solver="lbfgs",max_iter=150,tol=1e-5,random_state=0)
        model.fit(p[role==0,:n],y[role==0],sample_weight=w[role==0])
        models.append(model)
    auc=np.empty((len(test_ids),2),float)
    for i,sid in enumerate(test_ids):
        ii=sp==sid
        for j,(n,model) in enumerate(zip((11,12),models)):
            auc[i,j]=species_auc(y[ii],model.decision_function(p[ii,:n]))
    diff=auc[:,1]-auc[:,0]
    rng=np.random.default_rng(seed)
    draw=rng.integers(len(test_ids),size=(bootstrap_reps,len(test_ids)))
    limits=np.quantile(diff[draw].mean(axis=1),[.008333333333,.991666666667])
    return {"train_species":len(train_ids),"heldout_species":len(test_ids),
            "mean_baseline_auc":float(auc[:,0].mean()),
            "mean_augmented_auc":float(auc[:,1].mean()),
            "mean_delta_auc":float(diff.mean()),"ci_9833":[float(v) for v in limits],
            "positive_decision":bool(limits[0]>0)}
