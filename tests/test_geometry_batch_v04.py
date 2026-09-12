import numpy as np

from ttf.core import SpeciesSample, split_species
from ttf.geometry import SpeciesGeometry
from ttf.geometry_batch_v04 import (
    V04_LOCKED_CANDIDATE,
    V04_PROPENSITY_DIRECTIONS,
    V04_PROPENSITY_SEED,
    V04_PROPENSITY_TRANSITION_WIDTH,
    run_geometry_calibration_v04_locked_batched,
    score_v04_prepared_batch,
)
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.private_geometry_control import partial_spearman_controls
from ttf.transfer import prepare_transfer


def make_geometries(n_species=12, n_records=18):
    rng = np.random.default_rng(9917)
    return tuple(
        SpeciesGeometry(
            species=f"sp_{i:02d}",
            coordinates=rng.normal(size=(n_records + i % 3, 2)) + np.array([i * 0.03, -i * 0.02]),
        )
        for i in range(n_species)
    )


def test_v04_lock_constants_are_explicit():
    assert V04_LOCKED_CANDIDATE == "eval_geometry_partial"
    assert V04_PROPENSITY_DIRECTIONS == 2048
    assert V04_PROPENSITY_SEED == 20260908
    assert V04_PROPENSITY_TRANSITION_WIDTH == 0.2


def test_v04_batch_score_matches_manual_partial_spearman():
    geometries = make_geometries()
    labels = [g.species for g in geometries]
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=44)
    gmap = {g.species: g for g in geometries}
    samples = [
        SpeciesSample(name, gmap[name].coordinates, np.arange(len(gmap[name].coordinates)))
        for name in train + evaluation
    ]
    graphs = fixed_graphs(samples, k=3)
    edges = edges_on_fixed_graphs(samples, graphs)
    prepared = prepare_transfer(
        [edges[name] for name in train],
        [edges[name] for name in evaluation],
        bandwidth=1.1,
        prior_mean=0.0,
    )
    rng = np.random.default_rng(7)
    width = 3
    train_values = {
        name: rng.normal(size=(edges[name].n_edges, width)) for name in train
    }
    eval_values = {
        name: rng.normal(size=(edges[name].n_edges, width)) for name in evaluation
    }
    eval_length = {name: edges[name].length for name in evaluation}
    eval_propensity = {
        name: rng.normal(size=edges[name].n_edges) for name in evaluation
    }
    got = score_v04_prepared_batch(
        prepared,
        train_values,
        eval_values,
        eval_length=eval_length,
        eval_propensity=eval_propensity,
    )

    n_train = max(sl.stop for sl in prepared.train_slices.values())
    flat = np.empty((n_train, width))
    for name in train:
        flat[prepared.train_slices[name], :] = train_values[name]
    manual_stats = []
    for column in range(width):
        species_scores = []
        for name in evaluation:
            pred = (
                prepared.eval_projection[name] @ flat[:, column]
                + prepared.eval_prior_offset[name]
            )
            score = partial_spearman_controls(
                pred,
                eval_values[name][:, column],
                [eval_length[name], eval_propensity[name]],
            )
            assert np.isclose(got.species_scores[name][column], score, atol=1e-12)
            species_scores.append(score)
        manual_stats.append(np.mean(species_scores))
    assert np.allclose(got.statistics, manual_stats, atol=1e-12, rtol=0)


def test_v04_small_calibration_smoke():
    cells = run_geometry_calibration_v04_locked_batched(
        make_geometries(),
        shared_fractions=(0.0, 1.0),
        amplitudes=(2.0,),
        n_replicates=3,
        n_bootstrap=99,
        eval_fraction=0.5,
        k=3,
        bandwidth=1.0,
        seed=1234,
        world_batch_size=2,
        propensity_directions=64,
        propensity_seed=20260908,
    )
    assert len(cells) == 2
    assert all(np.isfinite(cell.mean_statistic) for cell in cells)
    assert all(0 <= cell.rejection_rate <= 1 for cell in cells)
