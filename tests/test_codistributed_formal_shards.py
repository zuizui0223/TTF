import itertools

import numpy as np

from ttf.codistributed_formal_shards import frozen_seed_range, validate_exact_ranges
from ttf.codistributed_geometry_null import (
    GeometryNullCenter,
    center_dyad_scores,
    prepare_fixed_dyad_transfer_cache,
    simulate_fixed_dyad_tst_batch,
    template_edges_from_geometries,
)
from ttf.genetic_geometry import GeneticSamplingGeometry, endpoint_disjoint_training_counts
from ttf.relational_dyadic import batch_primary_test, prepare_dyadic_regression


def _geometry(offset):
    base = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ]
    )
    xyz = base + np.asarray(offset, dtype=float)
    nodes = np.asarray(list(itertools.combinations(range(6), 2)), dtype=np.int64)
    return GeneticSamplingGeometry(
        coordinates=xyz,
        record_to_locality=np.arange(6, dtype=np.int64),
        records_per_locality=np.ones(6, dtype=np.int64),
        edge_nodes=nodes,
        endpoint_disjoint_training_edges=endpoint_disjoint_training_counts(nodes),
        graph_k=5,
    )


def test_declared_shard_ranges_reject_gap_or_overlap():
    want = [(0, 2), (2, 5), (5, 7)]
    assert validate_exact_ranges(want, want) == tuple(want)
    for bad in (
        [(0, 2), (3, 5), (5, 7)],
        [(0, 3), (2, 5), (5, 7)],
        [(0, 2), (2, 5)],
    ):
        try:
            validate_exact_ranges(bad, want)
        except RuntimeError:
            pass
        else:
            raise AssertionError("invalid shard layout was accepted")


def test_reduced_world_shards_reproduce_monolithic_fixed_tst_and_beta():
    geometries = {
        "s1": _geometry([0.0, 0.0, 0.0]),
        "s2": _geometry([0.2, 0.1, 0.0]),
        "s3": _geometry([-0.2, 0.2, 0.0]),
        "t1": _geometry([0.1, 0.3, 0.0]),
        "t2": _geometry([0.4, -0.1, 0.0]),
        "t3": _geometry([-0.3, -0.2, 0.0]),
    }
    sources = ("s1", "s2", "s3")
    targets = ("t1", "t2", "t3")
    pairs = tuple((s, t) for s in sources for t in targets)
    cache = prepare_fixed_dyad_transfer_cache(
        template_edges_from_geometries(geometries),
        pairs,
        bandwidth_km=2.0,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
        min_training_edges=5,
    )

    master = 20261003
    namespace = "codistributed-recurrence-shard-reduced-fixture"
    cell = "A2"
    all_seeds = frozen_seed_range(master, namespace, cell, 0, 7)
    monolithic = simulate_fixed_dyad_tst_batch(
        geometries,
        pairs,
        all_seeds,
        shared_fraction=0.4,
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
        prepared_cache=cache,
    )

    ranges = ((0, 2), (2, 5), (5, 7))
    sharded = np.concatenate(
        [
            simulate_fixed_dyad_tst_batch(
                geometries,
                pairs,
                frozen_seed_range(master, namespace, cell, start, stop),
                shared_fraction=0.4,
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
                prepared_cache=cache,
            )
            for start, stop in ranges
        ],
        axis=1,
    )
    np.testing.assert_allclose(sharded, monolithic, rtol=0, atol=0)

    source_index = np.repeat(np.arange(3, dtype=np.int64), 3)
    target_index = np.tile(np.arange(3, dtype=np.int64), 3)
    predictor = np.asarray([0.1, -0.3, 0.5, -0.2, 0.4, 0.8, -0.7, 0.6, -0.1])
    prepared = prepare_dyadic_regression(
        source_index,
        target_index,
        predictor[:, None],
        primary_index=0,
    )
    center = GeometryNullCenter(
        mu0=np.linspace(-0.05, 0.05, len(pairs)),
        amplitudes=(0.5, 1.0, 2.0, 3.0),
        worlds_per_amplitude=(199, 199, 199, 199),
    )
    mono_beta = batch_primary_test(
        prepared, center_dyad_scores(monolithic, center)
    ).coefficient
    shard_beta = np.concatenate(
        [
            batch_primary_test(
                prepared, center_dyad_scores(sharded[:, start:stop], center)
            ).coefficient
            for start, stop in ranges
        ]
    )
    np.testing.assert_allclose(shard_beta, mono_beta, rtol=0, atol=1e-14)
