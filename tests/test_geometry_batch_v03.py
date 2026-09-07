import numpy as np

from ttf.core import split_species
from ttf.geometry import SpeciesGeometry
from ttf.geometry_batch_v03 import (
    V03_LOCKED_CANDIDATE,
    run_geometry_calibration_v03_locked_batched,
)


def make_geometries(n_species: int = 12, n_records: int = 18):
    rng = np.random.default_rng(20260908)
    return tuple(
        SpeciesGeometry(
            species=f"fresh_{i:02d}",
            coordinates=(
                np.array([0.2 * np.cos(i), 0.2 * np.sin(i)])
                + rng.normal(0.0, 0.7, size=(n_records, 2))
            ),
        )
        for i in range(n_species)
    )


def test_v03_candidate_identity_is_locked():
    assert V03_LOCKED_CANDIDATE == "train_orthogonalized_only"


def test_v03_world_batch_size_does_not_change_results():
    geometries = make_geometries()
    labels = [item.species for item in geometries]
    train, evaluation = split_species(labels, eval_fraction=0.5, seed=77)
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
        seed=314159,
    )
    a = run_geometry_calibration_v03_locked_batched(**common, world_batch_size=2)
    b = run_geometry_calibration_v03_locked_batched(**common, world_batch_size=4)
    assert len(a) == len(b)
    for left, right in zip(a, b):
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
