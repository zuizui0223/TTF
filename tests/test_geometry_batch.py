import numpy as np

from ttf.batch import score_prepared_batch
from ttf.calibration import run_geometry_calibration
from ttf.core import SpeciesSample, edge_turnover, split_species
from ttf.geometry import SpeciesGeometry, simulate_fixed_geometry_boundary_world
from ttf.geometry_batch import run_geometry_calibration_batched
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.transfer import prepare_transfer


def make_geometries(n_species: int = 12, n_records: int = 16):
    rng = np.random.default_rng(20260908)
    return tuple(
        SpeciesGeometry(
            species=f"sp_{i:02d}",
            coordinates=(
                np.array([0.3 * np.cos(i), 0.3 * np.sin(i)])
                + rng.normal(0.0, 0.6, size=(n_records + i % 2, 2))
            ),
        )
        for i in range(n_species)
    )


def prepared_fixture():
    geometries = make_geometries()
    labels = [item.species for item in geometries]
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=77)
    geometry_map = {item.species: item for item in geometries}
    names = train + evaluation
    templates = [
        SpeciesSample(
            species=name,
            coordinates=geometry_map[name].coordinates,
            trait=np.arange(len(geometry_map[name].coordinates), dtype=float),
        )
        for name in names
    ]
    graphs = fixed_graphs(templates, k=3)
    edges = edges_on_fixed_graphs(templates, graphs)
    prepared = prepare_transfer(
        [edges[name] for name in train],
        [edges[name] for name in evaluation],
        bandwidth=0.5,
        segment_points=5,
    )
    return geometries, train, evaluation, graphs, prepared


def test_batched_prepared_scores_match_scalar_worlds():
    geometries, train, evaluation, graphs, prepared = prepared_fixture()
    names = train + evaluation
    train_columns = {name: [] for name in train}
    eval_columns = {name: [] for name in evaluation}
    scalar = []

    for seed in (10, 11, 12, 13):
        world = simulate_fixed_geometry_boundary_world(
            geometries,
            shared_fraction=0.5,
            amplitude=2.0,
            noise_sd=0.6,
            seed=seed,
        )
        sample_map = {sample.species: sample for sample in world.samples}
        turnover = {
            name: edge_turnover(sample_map[name], graphs[name]) for name in names
        }
        scalar.append(
            prepared.score(
                {name: turnover[name] for name in train},
                {name: turnover[name] for name in evaluation},
            )
        )
        for name in train:
            train_columns[name].append(turnover[name])
        for name in evaluation:
            eval_columns[name].append(turnover[name])

    batch = score_prepared_batch(
        prepared,
        {name: np.column_stack(train_columns[name]) for name in train},
        {name: np.column_stack(eval_columns[name]) for name in evaluation},
    )
    for column, result in enumerate(scalar):
        assert np.isclose(batch.statistics[column], result.statistic, atol=1e-12, rtol=0.0)
        assert batch.n_eval_species[column] == result.n_eval_species
        for name in evaluation:
            assert np.isclose(
                batch.species_scores[name][column],
                result.species_scores[name],
                atol=1e-12,
                rtol=0.0,
            )


def test_batched_geometry_calibration_matches_scalar_fixture():
    geometries = make_geometries()
    labels = [item.species for item in geometries]
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=123)
    common = dict(
        geometries=geometries,
        shared_fractions=(0.0, 1.0),
        amplitudes=(2.0,),
        n_replicates=4,
        n_bootstrap=99,
        train_species=train,
        eval_species=evaluation,
        k=3,
        bandwidth=0.5,
        noise_sd=0.6,
        seed=31415,
    )
    scalar = run_geometry_calibration(**common)
    batched = run_geometry_calibration_batched(**common, world_batch_size=2)
    assert len(scalar) == len(batched)
    for left, right in zip(scalar, batched):
        assert left.shared_fraction == right.shared_fraction
        assert left.amplitude == right.amplitude
        assert left.n_replicates == right.n_replicates
        assert left.rejection_rate == right.rejection_rate
        assert left.median_p_value == right.median_p_value
        assert np.isclose(left.mean_statistic, right.mean_statistic, atol=1e-12, rtol=0.0)
        assert np.isclose(
            left.mean_null_statistic,
            right.mean_null_statistic,
            atol=1e-12,
            rtol=0.0,
        )
