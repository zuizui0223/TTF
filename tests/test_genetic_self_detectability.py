import numpy as np

from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_self_detectability import (
    prepare_genetic_self_detectability,
    score_genetic_self_world_batch,
)
from ttf.genetic_simulate import simulate_genetic_distance_world


def _panel(n_species=12, n_localities=24):
    theta = np.linspace(0.0, 2.0 * np.pi, n_localities, endpoint=False)
    base = np.column_stack([np.cos(theta), np.sin(theta)])
    out = {}
    for i in range(n_species):
        shift = np.array([0.03 * (i % 3), 0.03 * (i // 3)])
        out[f"sp{i:02d}"] = prepare_density_scaled_genetic_geometry(base + shift)
    return out


def _design():
    geometries = _panel()
    design = prepare_genetic_ttf_design(
        geometries,
        eval_fraction=0.5,
        split_seed=19,
        bandwidth=1.5,
        min_training_edges=5,
    )
    self_design = prepare_genetic_self_detectability(
        design,
        bandwidth=1.5,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    return geometries, design, self_design


def test_self_field_projection_excludes_both_target_endpoints():
    _, design, self_design = _design()
    for name in self_design.species:
        nodes = design.template_edges[name].nodes
        projection = self_design.by_species[name].projection
        for target, (left, right) in enumerate(nodes):
            incident = (
                (nodes[:, 0] == left)
                | (nodes[:, 1] == left)
                | (nodes[:, 0] == right)
                | (nodes[:, 1] == right)
            )
            assert np.array_equal(projection[target, incident], np.zeros(np.count_nonzero(incident)))
            assert self_design.by_species[name].eligible_counts[target] >= 5


def test_self_field_preserves_constant_half_response_with_half_prior():
    _, design, self_design = _design()
    for name in self_design.species:
        prepared = self_design.by_species[name]
        response = np.full(design.template_edges[name].n_edges, 0.5)
        predicted = prepared.prior_offset + prepared.projection @ response
        assert np.allclose(predicted, 0.5, atol=1e-12, rtol=0.0)


def test_self_detectability_batch_returns_finite_species_level_scores():
    geometries, design, self_design = _design()
    worlds = [
        simulate_genetic_distance_world(
            geometries,
            shared_fraction=0.0,
            residual_amplitude=2.0,
            ibd_strength=1.0,
            noise_sd=0.10,
            seed=seed,
        )
        for seed in (101, 102, 103)
    ]
    scored = score_genetic_self_world_batch(design, self_design, worlds)
    assert scored.statistics.shape == (3,)
    assert np.all(np.isfinite(scored.statistics))
    assert np.array_equal(scored.n_species, np.full(3, len(design.eval_species)))
    assert set(scored.species_scores) == set(design.eval_species)
    assert all(np.all(np.isfinite(v)) for v in scored.species_scores.values())


def test_self_field_prediction_for_target_does_not_use_incident_responses():
    _, design, self_design = _design()
    name = self_design.species[0]
    edges = design.template_edges[name]
    prepared = self_design.by_species[name]
    target = 0
    left, right = edges.nodes[target]
    incident = (
        (edges.nodes[:, 0] == left)
        | (edges.nodes[:, 1] == left)
        | (edges.nodes[:, 0] == right)
        | (edges.nodes[:, 1] == right)
    )
    base = np.linspace(0.0, 1.0, edges.n_edges)
    altered = base.copy()
    altered[incident] += 1000.0
    pred0 = float(prepared.prior_offset[target] + prepared.projection[target] @ base)
    pred1 = float(prepared.prior_offset[target] + prepared.projection[target] @ altered)
    assert pred0 == pred1
