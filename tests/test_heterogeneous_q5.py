import numpy as np

from ttf.heterogeneous_inference import (
    density_scaled_k,
    heterogeneous_paired_heldout_species_bootstrap_test,
)
from ttf.heterogeneous_simulate import simulate_q5_world
from ttf.mismatch import pointwise_mismatch


def test_density_scaled_k_uses_effective_n():
    assert density_scaled_k(30) == 5
    assert density_scaled_k(45) == 7
    assert density_scaled_k(60) == 9
    assert density_scaled_k(90) == 14
    assert density_scaled_k(120) == 18
    assert density_scaled_k(8) == 2


def test_q5_geometry_is_response_blind_for_same_seed():
    a = simulate_q5_world(
        response_mode="null_component",
        geometry_profile="matched",
        seed=123,
    )
    b = simulate_q5_world(
        response_mode="power_shared_mismatch",
        geometry_profile="matched",
        seed=123,
    )
    assert a.geometry == b.geometry
    for sa, sb in zip(a.samples, b.samples):
        assert sa.species == sb.species
        assert np.array_equal(sa.coordinates, sb.coordinates)


def test_matched_geometry_represents_all_nominal_n_levels_in_each_split():
    world = simulate_q5_world(
        response_mode="null_component",
        geometry_profile="matched",
        seed=77,
    )
    expected = {30, 45, 60, 90, 120}
    train = {g.nominal_n for g in world.geometry if g.split == "train"}
    evaluation = {g.nominal_n for g in world.geometry if g.split == "eval"}
    assert train == expected
    assert evaluation == expected


def test_shifted_geometry_is_lower_density_and_more_clustered_in_evaluation():
    world = simulate_q5_world(
        response_mode="null_component",
        geometry_profile="shifted",
        seed=91,
    )
    train = [g for g in world.geometry if g.split == "train"]
    evaluation = [g for g in world.geometry if g.split == "eval"]
    assert np.mean([g.nominal_n for g in train]) > np.mean([g.nominal_n for g in evaluation])
    assert sum(g.clustering == "strong" for g in train) == 2
    assert sum(g.clustering == "strong" for g in evaluation) == 10
    assert {g.nominal_n for g in train} == {30, 45, 60, 90, 120}
    assert {g.nominal_n for g in evaluation} == {30, 45, 60, 90, 120}


def test_shared_relation_world_keeps_mismatch_magnitude_exactly_constant():
    world = simulate_q5_world(
        response_mode="power_shared_relation",
        geometry_profile="matched",
        amplitude=2.0,
        noise_sd=0.0,
        seed=222,
    )
    for sample in world.samples:
        mismatch = pointwise_mismatch(sample)
        assert np.all(mismatch == 2.0)


def test_variable_k_inference_smoke_and_reports_realized_graph_rule():
    world = simulate_q5_world(
        response_mode="power_shared_mismatch",
        geometry_profile="matched",
        amplitude=2.0,
        noise_sd=0.4,
        seed=404,
    )
    train = [f"sp_{i:03d}" for i in range(20)]
    evaluation = [f"sp_{i:03d}" for i in range(20, 40)]
    result = heterogeneous_paired_heldout_species_bootstrap_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        graph_fraction=0.15,
        bandwidth=0.2,
        n_bootstrap=99,
        seed=505,
        edge_chunk_size=32,
        train_chunk_size=2048,
    )
    assert np.isfinite(result.mismatch_statistic)
    assert np.isfinite(result.coupling_statistic)
    assert set(result.graph_k) == set(train + evaluation)
    for name, n in result.effective_n.items():
        assert result.graph_k[name] == density_scaled_k(n, 0.15)
