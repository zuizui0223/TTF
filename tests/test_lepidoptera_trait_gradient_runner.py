from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np

SCRIPT=Path("scripts/run_lepidoptera_trait_gradient_qualification_batch.py")
spec=importlib.util.spec_from_file_location("lep_runner",SCRIPT)
mod=importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(mod)

def test_two_sided_target_bootstrap_rejects_clear_positive_mean():
    slopes=np.asarray([0.7,0.9,1.0,1.1,1.2,0.8,1.3,0.95,1.05,1.15])
    mean,p=mod.centered_two_sided_target_bootstrap(slopes,seed=3,n_bootstrap=999)
    assert mean>0.0
    assert p<=0.05

def test_two_sided_target_bootstrap_keeps_exact_zero_mean():
    slopes=np.asarray([-2.,-1.5,-1.,-.5,0.,.5,1.,1.5,2.])
    mean,p=mod.centered_two_sided_target_bootstrap(slopes,seed=4,n_bootstrap=999)
    assert abs(mean)<1e-12
    assert p>0.5

def test_pair_surface_residual_is_geometry_orthogonal_within_target():
    pairs=[]
    for target,shift in (("a",0.0),("b",0.3)):
        for i in range(10):
            pairs.append({
                "target":target,
                "trait_similarity":0.2+0.03*i+shift+0.01*(i%3),
                "coverage":0.5+0.02*i,
                "centroid_distance":100+20*i,
                "edge_count_ratio":0.8+0.05*(i%4),
                "locality_count_ratio":0.9+0.04*(i%5),
            })
    target,residual=mod.prepare_pair_surface(pairs)
    geometry=np.asarray([[p[k] for k in ("coverage","centroid_distance","edge_count_ratio","locality_count_ratio")] for p in pairs],float)
    for label in ("a","b"):
        idx=np.flatnonzero(target==label)
        X=np.column_stack([np.ones(len(idx)),geometry[idx]])
        assert np.max(np.abs(X.T@residual[idx]))<1e-8

def test_wilson_bounds_contain_rate():
    lo,hi=mod.wilson(25,500)
    assert lo<0.05<hi


def test_score_pair_contributions_uses_implicit_canonical_target_edges():
    class World:
        edge_response={
            "t":np.asarray([0.,1.,3.,7.]),
            "s1":np.asarray([7.,3.,1.,0.]),
            "s2":np.asarray([0.,3.,1.,7.]),
            "s3":np.asarray([0.,1.,3.,7.]),
        }
    pairs=[
        {"target":"t","source":"s1","source_edge_index":[0,1,2,3]},
        {"target":"t","source":"s2","source_edge_index":[0,1,2,3]},
        {"target":"t","source":"s3","source_edge_index":[0,1,2,3]},
    ]
    target=np.asarray(["t","t","t"])
    residual=np.asarray([-1.,0.5,1.5])
    stat,slopes=mod.score_pair_contributions(World(),pairs,target,residual)
    assert np.isfinite(stat)
    assert slopes.shape==(1,)
