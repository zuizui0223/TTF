import json
from pathlib import Path

import numpy as np

from scripts.run_codistributed_recurrence_fixed_tst_formal_shard import (
    expected_ranges,
    seed_range,
)
from ttf.codistributed_geometry_null import (
    frozen_uint64_seed,
    simulate_fixed_dyad_tst_batch,
)
from ttf.genetic_geometry import prepare_genetic_sampling_geometry


ROOT=Path(__file__).resolve().parents[1]


def _complete_geometry(dx: float, dy: float):
    base=np.asarray([
        [0.0,0.0],
        [1.0,0.1],
        [0.2,1.1],
        [1.3,1.4],
        [2.0,0.7],
        [0.7,2.2],
    ],dtype=float)
    return prepare_genetic_sampling_geometry(
        base+np.asarray([dx,dy]),
        k=5,
    )


def test_formal_seed_slice_preserves_original_replicate_indices():
    master=20261003
    namespace="formal-test"
    cell="A2"
    full=[
        frozen_uint64_seed(master,namespace,cell,i)
        for i in range(9)
    ]
    assert seed_range(master,namespace,cell,2,6)==full[2:6]
    assert seed_range(master,namespace,cell,6,9)==full[6:9]


def test_frozen_formal_ranges_match_contract():
    p=json.loads(
        (ROOT/"docs/supporting/genetic_codistributed_recurrence_formal_sharding_v0.1.json").read_text()
    )
    assert expected_ranges(p,"center")==[(0,199)]
    assert expected_ranges(p,"reference")==[(0,200),(200,400),(400,600),(600,800),(800,999)]
    assert expected_ranges(p,"evaluation")==[(0,200),(200,400),(400,500)]
    assert expected_ranges(p,"positive")==[(0,200),(200,400),(400,500)]


def test_reduced_fixed_tst_shards_match_monolithic_world_order():
    geometries={
        "s0":_complete_geometry(0.0,0.0),
        "s1":_complete_geometry(0.3,0.2),
        "t0":_complete_geometry(-0.2,0.4),
        "t1":_complete_geometry(0.5,-0.3),
    }
    pairs=[("s0","t0"),("s0","t1"),("s1","t0"),("s1","t1")]
    seeds=seed_range(20261003,"reduced-shard-equivalence","A2",0,6)
    kwargs=dict(
        shared_fraction=0.0,
        residual_amplitude=2.0,
        ibd_strength=1.0,
        noise_sd=0.10,
        transition_width=0.20,
        noise_dimensions=2,
        min_training_edges=5,
        bandwidth_km=2.0,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )
    monolithic=simulate_fixed_dyad_tst_batch(
        geometries,pairs,seeds,**kwargs
    )
    sharded=np.concatenate([
        simulate_fixed_dyad_tst_batch(geometries,pairs,seeds[0:2],**kwargs),
        simulate_fixed_dyad_tst_batch(geometries,pairs,seeds[2:4],**kwargs),
        simulate_fixed_dyad_tst_batch(geometries,pairs,seeds[4:6],**kwargs),
    ],axis=1)
    np.testing.assert_allclose(sharded,monolithic,rtol=0,atol=1e-12)
    np.testing.assert_allclose(
        np.sum(sharded,axis=1)/sharded.shape[1],
        np.sum(monolithic,axis=1)/monolithic.shape[1],
        rtol=0,atol=1e-14,
    )
