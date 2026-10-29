"""Pre-registered one-shot He contrast estimator: fixed linear nested models, species holdout."""
import collections
import numpy as np

def evaluate(groups, draws=1999, seed=20261029):
    train=[g for g in groups.values() if g["role"]==0]
    test=[g for g in groups.values() if g["role"]==1]
    if min(len(train),len(test))<20:raise ValueError("too few study groups")
    now=np.vstack([g["x"][:,4:8] for g in train])
    clim_sd=now.std(axis=0)
    if np.any(clim_sd<=1e-10):raise ValueError("degenerate training climate")
    nn=np.concatenate([g["x"][:,10] for g in train])
    known=nn[np.isfinite(nn)&(nn>0)]
    logn_median=float(np.median(np.log(known))) if len(known) else 0.
    def transform(g):
        a=g["x"];n=a[:,10]
        historical=np.linalg.norm((a[:,:4]-a[:,4:8])/clim_sd,axis=1)
        ln=np.full(len(n),logn_median)
        good=np.isfinite(n)&(n>0)
        ln[good]=np.log(n[good])
        x=np.column_stack([a[:,4:8],a[:,8:10],ln,~good,historical]).astype(float)
        x-=x.mean(axis=0)
        y=g["y"].copy();y-=y.mean()
        return x,y
    trainxy=[transform(g) for g in train]
    scale=np.sqrt(np.sum([np.sum(x*x,axis=0)/len(x) for x,y in trainxy],axis=0)/len(trainxy))
    scale[scale<1e-10]=1.
    X=np.vstack([x/scale/np.sqrt(len(x)) for x,y in trainxy])
    Y=np.concatenate([y/np.sqrt(len(y)) for x,y in trainxy])
    baseline=np.linalg.lstsq(X[:,:-1],Y,rcond=1e-10)[0]
    enriched=np.linalg.lstsq(X,Y,rcond=1e-10)[0]
    d=collections.defaultdict(list)
    for g in test:
        x,y=transform(g);x=x/scale;n=len(y)
        e0=y-x[:,:-1]@baseline
        e1=y-x@enriched
        d[g["species"]].append((2*n*np.mean(e0**2)/(n-1),
                                 2*n*np.mean(e1**2)/(n-1)))
    all_scores=np.asarray([p for z in d.values() for p in z])
    loss0=float(all_scores[:,0].mean());loss1=float(all_scores[:,1].mean())
    gain=(loss0-loss1)/loss0
    rng=np.random.default_rng(seed);sp=list(d);vals=[]
    for k in range(draws):
        sampled=rng.integers(len(sp),size=len(sp))
        a=np.asarray([p for i in sampled for p in d[sp[i]]])
        vals.append((a[:,0].mean()-a[:,1].mean())/a[:,0].mean())
    bounds=np.quantile(vals,[.025,.975])
    return {"status":"OPENED_NEW_INDEPENDENT_GENDIVRANGE_OUTCOME",
      "n_train_species":len({g["species"] for g in train}),
      "n_eval_species":len(d),
      "n_train_groups":len(train),"n_eval_groups":len(test),
      "n_train_sites":sum(len(g["y"]) for g in train),
      "n_eval_sites":sum(len(g["y"]) for g in test),
      "baseline_pair_MSE":loss0,"history_pair_MSE":loss1,
      "relative_heldout_MSE_improvement":gain,
      "species_cluster_bootstrap95":[float(x) for x in bounds],
      "bootstrap_draws":draws,"bootstrap_seed":seed,
      "historical_predictor_normalized_slope":float(enriched[-1]),
      "positive_prediction_gain_95lower_gt_zero":bool(gain>0 and bounds[0]>0),
      "prior_TTF_genetic_response_reopened":False}
