import numpy as np

from ttf.q5_factor_diagnostic import simulate_factor_world


def _geom(profile: str):
    world = simulate_factor_world(
        profile=profile,
        response_mode="null_private_relation",
        amplitude=2.0,
        noise_sd=0.0,
        transition_width=0.2,
        seed=11,
    )
    return world.geometry


def test_balanced_profile_is_fully_uniform_geometry():
    g = _geom("balanced")
    assert {x.nominal_n for x in g} == {60}
    assert {x.clustering for x in g} == {"uniform"}
    assert {x.truncation for x in g} == {"full"}
    assert {x.missingness for x in g} == {0.0}


def test_n_only_changes_sample_size_only():
    g = _geom("n_only")
    assert {x.nominal_n for x in g} == {30, 45, 60, 90, 120}
    assert {x.clustering for x in g} == {"uniform"}
    assert {x.truncation for x in g} == {"full"}
    assert {x.missingness for x in g} == {0.0}


def test_clustering_only_changes_clustering_only():
    g = _geom("clustering_only")
    assert {x.nominal_n for x in g} == {60}
    assert {x.clustering for x in g} == {"uniform", "moderate", "strong"}
    assert {x.truncation for x in g} == {"full"}
    assert {x.missingness for x in g} == {0.0}


def test_truncation_only_changes_support_only():
    g = _geom("truncation_only")
    assert {x.nominal_n for x in g} == {60}
    assert {x.clustering for x in g} == {"uniform"}
    assert {x.truncation for x in g} == {"full", "one_sided", "interior_gap"}
    assert {x.missingness for x in g} == {0.0}


def test_missingness_only_changes_row_loss_only():
    g = _geom("missingness_only")
    assert {x.nominal_n for x in g} == {60}
    assert {x.clustering for x in g} == {"uniform"}
    assert {x.truncation for x in g} == {"full"}
    assert {x.missingness for x in g} == {0.0, 0.1, 0.25}
    assert min(x.effective_n for x in g) >= 8


def test_n_shift_only_reverses_density_weight_between_splits():
    g = _geom("n_shift_only")
    train = [x.nominal_n for x in g if x.split == "train"]
    evaluation = [x.nominal_n for x in g if x.split == "eval"]
    assert np.mean(train) > np.mean(evaluation)
    assert {x.clustering for x in g} == {"uniform"}


def test_clustering_shift_only_makes_evaluation_more_clustered():
    g = _geom("clustering_shift_only")
    weight = {"uniform": 0, "moderate": 1, "strong": 2}
    train = [weight[x.clustering] for x in g if x.split == "train"]
    evaluation = [weight[x.clustering] for x in g if x.split == "eval"]
    assert np.mean(evaluation) > np.mean(train)
    assert {x.nominal_n for x in g} == {60}
