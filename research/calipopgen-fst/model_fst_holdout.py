"""Frozen species-heldout study-centered Fst contrast estimator. No source IO."""
from collections import defaultdict
import numpy as np

def evaluate(groups, bootstrap_draws=1999, seed=20261029):
    train=[g for g in groups if g["role"]==0]
    test=[g for g in groups if g["role"]==1]
    if not train or not test:raise ValueError("empty species role")
    if {g["species"] for g in train}&{g["species"] for g in test}:
        raise ValueError("species leakage")
    for g in groups:
        X=np.asarray(g["x"],float);y=np.asarray(g["y"],float)
        if len(y)<5 or X.shape!=(len(y),5) or not np.isfinite(X).all() or not np.isfinite(y).all():
            raise ValueError("invalid pairwise group")
        if abs(float(y.mean()))>1e-8 or np.max(np.abs(X.mean(axis=0)))>1e-8:
            raise ValueError("group must be centered consistently")
    x0=np.vstack([g["x"][:,:4]/np.sqrt(len(g["y"])) for g in train])
    x1=np.vstack([g["x"]/np.sqrt(len(g["y"])) for g in train])
    ytrain=np.concatenate([g["y"]/np.sqrt(len(g["y"])) for g in train])
    b0=np.linalg.lstsq(x0,ytrain,rcond=1e-12)[0]
    b1=np.linalg.lstsq(x1,ytrain,rcond=1e-12)[0]
    scores=defaultdict(list)
    for g in test:
        X=g["x"];y=g["y"]
        loss0=float(np.mean((y-X[:,:4]@b0)**2))
        loss1=float(np.mean((y-X@b1)**2))
        scores[g["species"]].append((loss0,loss1))
    pair_by_group=np.asarray([p for v in scores.values() for p in v],float)
    baseline=float(pair_by_group[:,0].mean())
    newer=float(pair_by_group[:,1].mean())
    gain=(baseline-newer)/baseline if baseline>0 else float("nan")
    rng=np.random.default_rng(seed);species=list(scores)
    bootstrap=np.empty(bootstrap_draws,float)
    for k in range(bootstrap_draws):
        selected=rng.integers(len(species),size=len(species))
        drawn=np.asarray([p for i in selected for p in scores[species[i]]])
        bootstrap[k]=(drawn[:,0].mean()-drawn[:,1].mean())/drawn[:,0].mean()
    bounds=np.quantile(bootstrap,[.025,.975])
    return {"status":"INDEPENDENT_CALIPOPGEN_FST_FIRST_SHOT",
       "train_species":len({g["species"] for g in train}),
       "heldout_species":len(scores),
       "train_study_marker_groups":len(train),
       "heldout_study_marker_groups":len(test),
       "train_pairwise_comparisons":sum(len(g["y"]) for g in train),
       "heldout_pairwise_comparisons":sum(len(g["y"]) for g in test),
       "baseline_group_equal_MSE":baseline,
       "plus_history_group_equal_MSE":newer,
       "relative_heldout_MSE_gain":float(gain),
       "history_coefficient":float(b1[-1]),
       "species_cluster_bootstrap95":bounds.tolist(),
       "bootstrap_draws":bootstrap_draws,"bootstrap_seed":seed,
       "positive_gain_95_lower_gt_zero":bool(gain>0 and bounds[0]>0),
       "within_study_only_contrast":True}
