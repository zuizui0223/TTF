"""Locked ridge/MSE and synthetic-world operations. No original trait response imported."""
from __future__ import annotations
import numpy as np
from collections import Counter
ALPHA=0.1
SEED=20261114
CI_PROBS=(.00833333333333333,.991666666666667)

def equal_species_weights(species,role,wanted="train"):
    sp=np.asarray(species);role=np.asarray(role)
    ids=np.unique(sp[role==wanted])
    if not len(ids):raise ValueError("No eligible species for "+wanted)
    w=np.zeros(len(sp))
    for sid in ids:
        at=(sp==sid)&(role==wanted)
        w[at]=1/(len(ids)*int(at.sum()))
    return w

def scale_train(X,w):
    X=np.asarray(X,float);w=np.asarray(w,float)
    if X.ndim!=2 or len(X)!=len(w) or not np.isfinite(X).all() or np.any(w<0) or not np.isclose(w.sum(),1,atol=1e-9):
        raise ValueError("Finite columns and species-equal normalized source train weights required")
    mu=w@X
    sd=np.sqrt(np.maximum(w@((X-mu)**2),0))
    return mu,np.where(sd>1e-10,sd,1.)

def ridge_predict(X,y,species,role,history=None,alpha=ALPHA):
    X=np.asarray(X,float);y=np.asarray(y,float)
    sp=np.asarray(species);role=np.asarray(role)
    if X.ndim!=2 or X.shape!=(len(y),21) or len(sp)!=len(y) or len(role)!=len(y):
        raise ValueError("One model requires the 21 frozen baseline columns")
    if history is not None:
        h=np.asarray(history,float)
        if h.shape!=(len(y),):raise ValueError("Exactly one historical increment")
        X=np.column_stack([X,h])
    if not np.isfinite(X).all() or not np.isfinite(y).all() or alpha!=ALPHA:
        raise ValueError("Nonfinite predictors/outcome or respecified alpha")
    w=equal_species_weights(sp,role)
    train=role=="train"
    mu,sd=scale_train(X[train],w[train])
    Z=(X-mu)/sd
    mean_y=w[train]@y[train]
    A=Z[train].T@(Z[train]*w[train,None])+alpha*np.eye(Z.shape[1])
    coef=np.linalg.solve(A,Z[train].T@(w[train]*(y[train]-mean_y)))
    return mean_y+Z@coef

def species_mse_comparison(y,pred_base,pred_hist,species,role,reps=1999,seed=SEED):
    y=np.asarray(y,float);b=np.asarray(pred_base,float);a=np.asarray(pred_hist,float)
    sp=np.asarray(species);r=np.asarray(role)
    if not (len(y)==len(b)==len(a)==len(sp)==len(r)) or not np.isfinite(y).all() or not np.isfinite(b).all() or not np.isfinite(a).all():
        raise ValueError("Unequal or nonfinite source response/predictions")
    ids=np.unique(sp[r=="heldout"])
    if len(ids)<1:raise ValueError("No species disjoint holdout")
    mb=np.array([np.mean((y[(sp==s)&(r=="heldout")]-b[(sp==s)&(r=="heldout")])**2) for s in ids])
    ma=np.array([np.mean((y[(sp==s)&(r=="heldout")]-a[(sp==s)&(r=="heldout")])**2) for s in ids])
    M=float(mb.mean());A=float(ma.mean())
    if M<=1e-12:raise ValueError("Baseline MSE not defined")
    draw=np.random.default_rng(seed).integers(0,len(ids),size=(reps,len(ids)))
    B=mb[draw].mean(axis=1);C=ma[draw].mean(axis=1)
    if (B<=1e-12).any():raise ValueError("Bootstrap MSE denominator zero")
    ci=np.quantile((B-C)/B,CI_PROBS)
    return {"heldout_species":len(ids),"baseline_species_equal_MSE":M,
            "plus_history_species_equal_MSE":A,
            "relative_MSE_gain":(M-A)/M,
            "species_bootstrap_9833":[float(x) for x in ci],
            "positive":bool(ci[0]>0)}

def known_truth_labels(X,H,species,studies,role,rng,coefficient):
    X=np.asarray(X,float);H=np.asarray(H,float)
    sp=np.asarray(species);st=np.asarray(studies);r=np.asarray(role)
    if X.shape!=(len(H),21) or len(sp)!=len(H) or len(st)!=len(H):
        raise ValueError("Frozen synthetic predictor geometry required")
    w=equal_species_weights(sp,r);train=r=="train"
    mu,sd=scale_train(X[train],w[train])
    Z=(X-mu)/sd
    muH,sdH=scale_train(H[train,None],w[train])
    h=(H-muH[0])/sdH[0]
    design=np.column_stack([np.ones(len(H)),Z])
    beta=np.linalg.lstsq(design[train]*np.sqrt(w[train])[:,None],
                         h[train]*np.sqrt(w[train]),rcond=1e-12)[0]
    residual=h-design@beta
    sigma=np.sqrt(w[train]@(residual[train]**2))
    if sigma<1e-10:raise ValueError("No source historical predictor residual")
    residual/=sigma
    base_coef=np.array([(0.12 if i%2==0 else -0.12) if i<11 else
                        (0.06 if i%2==0 else -0.06) for i in range(21)])
    source_species={x:rng.normal(0,.5) for x in np.unique(sp)}
    source_study={x:rng.normal(0,.5) for x in np.unique(st)}
    noise=np.array([source_species[x]+source_study[z] for x,z in zip(sp,st)])
    return Z@base_coef+coefficient*residual+noise+rng.normal(0,.5,len(H))

def wilson95(k,n):
    if n<1 or not 0<=k<=n:raise ValueError("No proper source simulations")
    z=1.959963984540054;p=k/n;d=1+z*z/n
    center=(p+z*z/(2*n))/d
    spread=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [float(max(0,center-spread)),float(min(1,center+spread))]
