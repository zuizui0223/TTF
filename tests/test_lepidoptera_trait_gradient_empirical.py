from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.lepidoptera_trait_gradient_empirical import (
    pair_transfer_congruence,
    post_ibd_edge_response,
)


def test_post_ibd_edge_response_is_finite_rank01():
    theta=np.linspace(0,2*np.pi,16,endpoint=False)
    xyz=np.column_stack([np.cos(theta),np.sin(theta),np.zeros(len(theta))])
    g=prepare_density_scaled_genetic_geometry(xyz,neighbor_fraction=.15)
    nodes=np.asarray(g.edge_nodes,int)
    geographic=np.linalg.norm(g.coordinates[nodes[:,0]]-g.coordinates[nodes[:,1]],axis=1)
    genetic=0.3*geographic+np.linspace(0,1,len(geographic))
    r=post_ibd_edge_response(genetic,g,min_training_edges=5)
    assert r.shape==(g.n_edges,)
    assert np.isfinite(r).all()
    assert np.all((r>=0)&(r<=1))


def test_pair_transfer_congruence_recovers_identical_aligned_response():
    names=("a","b")
    responses={"a":np.array([0.,1.,2.,3.]),"b":np.array([0.,1.,2.,3.])}
    pairs=np.array([[1,0,0.5,1,1,1,1,0,4]],float)
    pair_id=np.zeros(4,dtype=int)
    a=np.arange(4,dtype=int)
    score=pair_transfer_congruence(responses,names,pairs,pair_id,a,a)
    assert np.allclose(score,[1.0])


def test_empirical_rule_is_frozen_before_identity_opening():
    rule=json.loads(Path("docs/supporting/lepidoptera_trait_gradient_empirical_rule_v0.1.json").read_text())
    assert rule["status"]=="FROZEN_BEFORE_NEW_LEPIDOPTERA_NUCLEOTIDE_IDENTITY_OPENING"
    assert rule["survivor_design"]["species"]==240
    assert rule["survivor_design"]["geographically_supported_eval_species"]==110
    assert rule["survivor_references"]["worlds_each"]==499
    assert all(value is False for value in rule["outcome_firewall"].values())
