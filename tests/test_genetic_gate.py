import numpy as np

from ttf.genetic_gate import prepare_genetic_ttf_design, score_genetic_world
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_simulate import simulate_genetic_distance_world


def _panel(n_species=12, n_localities=24):
    theta = np.linspace(0.0, 2.0 * np.pi, n_localities, endpoint=False)
    base = np.column_stack([np.cos(theta), np.sin(theta)])
    out = {}
    for i in range(n_species):
        shift = np.array([0.03 * (i % 3), 0.03 * (i // 3)])
        out[f"sp{i:02d}"] = prepare_density_scaled_genetic_geometry(base + shift)
    return out


def test_genetic_design_is_species_disjoint_and_endpoint_safe():
    geometries = _panel()
    design = prepare_genetic_ttf_design(
        geometries,
        eval_fraction=0.5,
        split_seed=17,
        bandwidth=1.5,
        min_training_edges=5,
    )
    assert len(design.train_species) == 6
    assert len(design.eval_species) == 6
    assert not (set(design.train_species) & set(design.eval_species))
    assert all(
        geometries[name].min_endpoint_disjoint_training_edges >= 5
        for name in design.train_species + design.eval_species
    )


def test_genetic_world_scores_through_ibd_and_v011_geometry_layers():
    geometries = _panel()
    design = prepare_genetic_ttf_design(
        geometries,
        eval_fraction=0.5,
        split_seed=19,
        bandwidth=1.5,
        min_training_edges=5,
    )
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=1.0,
        residual_amplitude=1.5,
        ibd_strength=1.0,
        noise_sd=0.05,
        seed=23,
    )
    score = score_genetic_world(design, world)
    assert np.isfinite(score.statistic)
    assert np.isfinite(score.training_strength)
    assert set(score.species_scores) == set(design.eval_species)
    assert all(np.isfinite(list(score.species_scores.values())))
    assert set(score.ibd) == set(design.train_species + design.eval_species)
    for name, result in score.ibd.items():
        assert len(result.residual_turnover) == geometries[name].n_edges
        assert np.all(np.isfinite(result.residual_turnover))
