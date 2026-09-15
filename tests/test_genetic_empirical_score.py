from __future__ import annotations

import json

import numpy as np
import pytest

from ttf.genetic_empirical_score import (
    score_genetic_distance_mapping,
    score_self_detectability_mapping,
    score_total_genetic_distance_mapping,
)
from ttf.genetic_gate import prepare_genetic_ttf_design, score_genetic_world
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_self_detectability import (
    prepare_genetic_self_detectability,
    score_genetic_self_world_batch,
)
from ttf.genetic_simulate import simulate_genetic_distance_world


def _geometries() -> dict:
    out = {}
    for species_index in range(12):
        coordinates = np.column_stack(
            [
                np.arange(12, dtype=float) * 25.0 + species_index * 3.0,
                np.sin(np.arange(12, dtype=float) / 2.0) * 10.0,
                np.full(12, species_index * 2.0, dtype=float),
            ]
        )
        out[f"species_{species_index:02d}"] = prepare_density_scaled_genetic_geometry(
            coordinates, neighbor_fraction=0.15
        )
    return out


def _design(geometries: dict):
    names = tuple(sorted(geometries))
    return names, prepare_genetic_ttf_design(
        geometries,
        train_species=names[:6],
        eval_species=names[6:],
        bandwidth=500.0,
        prior_strength=0.25,
        segment_points=5,
        min_training_edges=5,
        strength_neighbours=4,
    )


def test_empirical_cross_adapter_matches_synthetic_scorer_exactly() -> None:
    geometries = _geometries()
    names, design = _design(geometries)
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.5,
        residual_amplitude=2.0,
        ibd_strength=1.0,
        noise_sd=0.1,
        transition_width=0.2,
        noise_dimensions=2,
        seed=12345,
    )
    expected = score_genetic_world(design, world)
    observed = score_genetic_distance_mapping(design, world.genetic_distance)
    assert observed.statistic == expected.statistic
    assert observed.training_strength == expected.training_strength
    assert observed.species_scores == expected.species_scores
    for name in names:
        np.testing.assert_array_equal(
            observed.ibd[name].residual_turnover,
            expected.ibd[name].residual_turnover,
        )


def test_empirical_self_adapter_matches_synthetic_self_scorer_exactly() -> None:
    geometries = _geometries()
    _, design = _design(geometries)
    self_design = prepare_genetic_self_detectability(
        design,
        bandwidth=500.0,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.0,
        residual_amplitude=2.0,
        ibd_strength=1.0,
        noise_sd=0.1,
        transition_width=0.2,
        noise_dimensions=2,
        seed=67890,
    )
    expected = score_genetic_self_world_batch(design, self_design, [world])
    observed = score_self_detectability_mapping(
        design, self_design, world.genetic_distance
    )
    assert observed.statistic == float(expected.statistics[0])
    assert observed.species_scores == {
        name: float(expected.species_scores[name][0]) for name in self_design.species
    }


def _independent_rank(values: np.ndarray) -> np.ndarray:
    # Pairwise counting, independent of the production sorting/ranking helpers.
    return np.array([np.sum(values < value) + (np.sum(values == value) + 1) / 2 for value in values])


def _dense_total_oracle(design, distances: dict) -> dict[str, float]:
    positions, weights, responses = [], [], []
    for name in design.train_species:
        edges = design.template_edges[name]
        y = _independent_rank(distances[name])
        y -= y.mean()
        z = _independent_rank(edges.length)
        z -= z.mean()
        residual = y if z @ z <= np.finfo(float).eps else y - z * ((z @ y) / (z @ z))
        sd = residual.std()
        target_sd = np.sqrt((len(y)**2 - 1) / (12 * len(y)**2))
        residual = np.zeros_like(y) if sd <= np.sqrt(np.finfo(float).eps) * len(y) else residual * target_sd / sd
        positions.extend(edges.midpoint)
        weights.extend([1 / edges.n_edges] * edges.n_edges)
        responses.extend(residual)
    positions, weights, responses = map(np.array, (positions, weights, responses))
    config = design.prepared
    scores = {}
    for name in design.eval_species:
        edges = design.template_edges[name]
        fractions = (np.arange(config.segment_points) + 0.5) / config.segment_points
        points = edges.start[:, None] + fractions[None, :, None] * (edges.end - edges.start)[:, None]
        distance2 = np.sum((points[:, :, None] - positions[None, None])**2, axis=-1)
        kernel = np.exp(-distance2 / (2 * config.bandwidth**2)) * weights
        prediction = ((kernel @ responses + config.prior_strength * config.prior_mean) /
                      (kernel.sum(axis=-1) + config.prior_strength)).mean(axis=1)
        x, y = _independent_rank(prediction), _independent_rank(distances[name])
        x, y = x - x.mean(), y - y.mean()
        denom = np.sqrt((x @ x) * (y @ y))
        scores[name] = 0.0 if denom <= np.finfo(float).eps else float((x @ y) / denom)
    return scores


def _total_fixture():
    geometries = _geometries()
    # Unequal edge counts distinguish species-equal weighting from edge pooling.
    for index, name in enumerate(sorted(geometries)):
        n = 12 + 4 * (index % 3)
        coordinates = np.column_stack([
            np.arange(n) * 17.0 + index * 3.0,
            np.sin(np.arange(n) / 2.0) * 10.0,
            np.full(n, index * 2.0),
        ])
        geometries[name] = prepare_density_scaled_genetic_geometry(coordinates, neighbor_fraction=0.15)
    _, design = _design(geometries)
    rng = np.random.default_rng(7139)
    # Ties are deliberate; no empirical sequence or archive is read.
    distances = {name: rng.integers(0, 7, len(geometry.edge_nodes)).astype(float)
                 for name, geometry in geometries.items()}
    return design, distances


@pytest.mark.parametrize("edge_chunk,train_chunk", [(1, 7), (32, 4096)])
def test_total_matches_independent_dense_oracle_with_ties(edge_chunk, train_chunk):
    design, distances = _total_fixture()
    expected = _dense_total_oracle(design, distances)
    actual = score_total_genetic_distance_mapping(
        design, distances, edge_chunk_size=edge_chunk, train_chunk_size=train_chunk)
    for name, value in expected.items():
        assert actual.species_scores[name] == pytest.approx(value, abs=1e-12)
    assert actual.statistic == pytest.approx(np.mean(list(expected.values())), abs=1e-12)
    payload = actual.as_dict()
    assert payload["inferentially_qualified"] is False
    assert payload["used_for_primary_decision"] is False
    assert "p_value" not in payload and "positive" not in payload
    json.dumps(payload, allow_nan=False)


def test_total_preserves_primary_scores_and_inputs():
    design, distances = _total_fixture()
    copied = {name: values.copy() for name, values in distances.items()}
    before = score_genetic_distance_mapping(design, distances)
    total = score_total_genetic_distance_mapping(design, distances)
    after = score_genetic_distance_mapping(design, distances)
    assert before.statistic == after.statistic
    assert before.training_strength == after.training_strength
    assert before.species_scores == after.species_scores
    assert total.statistic != before.statistic  # distinct responses, not relabelled primary
    for name in distances:
        np.testing.assert_array_equal(distances[name], copied[name])
        np.testing.assert_array_equal(before.ibd[name].residual, after.ibd[name].residual)


def test_total_constant_vectors_keep_zero_scores_and_all_species():
    design, distances = _total_fixture()
    for values in distances.values():
        values.fill(1.0)
    result = score_total_genetic_distance_mapping(design, distances)
    assert result.statistic == 0.0
    assert result.species_scores == {name: 0.0 for name in design.eval_species}
    assert result.as_dict()["finite_evaluation_species_count"] == len(design.eval_species)


@pytest.mark.parametrize("invalid", ["missing", "extra", "shape", "nan", "negative"])
def test_total_rejects_invalid_response_mapping(invalid):
    design, distances = _total_fixture()
    name = design.eval_species[0]
    if invalid == "missing":
        del distances[name]
    elif invalid == "extra":
        distances["unexpected"] = distances[name].copy()
    elif invalid == "shape":
        distances[name] = distances[name][:-1]
    elif invalid == "nan":
        distances[name][0] = np.nan
    else:
        distances[name][0] = -1
    with pytest.raises(ValueError):
        score_total_genetic_distance_mapping(design, distances)


def test_total_nonfinite_species_does_not_shrink_aggregate(monkeypatch):
    import ttf.genetic_empirical_score as module
    from types import SimpleNamespace
    design, distances = _total_fixture()
    scores = {name: np.array([0.5]) for name in design.eval_species}
    scores[design.eval_species[0]][0] = np.nan
    monkeypatch.setattr(module, "score_chunked_batch", lambda *args, **kwargs:
                        SimpleNamespace(statistics=np.array([0.5]), species_scores=scores))
    result = score_total_genetic_distance_mapping(design, distances)
    assert result.statistic is None
    assert len(result.species_scores) == len(design.eval_species)
    assert result.as_dict()["status"] == "NOT_EVALUABLE_DESCRIPTIVE"
    json.dumps(result.as_dict(), allow_nan=False)
