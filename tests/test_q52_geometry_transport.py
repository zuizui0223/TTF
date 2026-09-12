import numpy as np

from ttf.core import SpeciesEdges
from ttf.geometry_transport import geometry_transport_train_weights
from ttf.mismatch import PairedSpeciesSample
from ttf.q52_transport import DEVELOPMENT_MODES, q52_geometry_transport_development_test


def _toy_edges(name: str, shift: float = 0.0, turnover_offset: float = 0.0) -> SpeciesEdges:
    theta = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False) + float(shift)
    coordinates = np.column_stack((np.cos(theta), np.sin(theta)))
    nodes = np.asarray([[i, (i + 1) % 8] for i in range(8)], dtype=np.int64)
    nodes = np.sort(nodes, axis=1)
    nodes = np.unique(nodes, axis=0)
    start = coordinates[nodes[:, 0]]
    end = coordinates[nodes[:, 1]]
    midpoint = 0.5 * (start + end)
    length = np.linalg.norm(end - start, axis=1)
    turnover = np.linspace(0.1, 0.9, len(nodes)) + float(turnover_offset)
    return SpeciesEdges(
        species=name,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=midpoint,
        length=length,
        turnover=turnover,
    )


def _species_slices(edge_sets):
    out = {}
    cursor = 0
    for edges in edge_sets:
        out[edges.species] = slice(cursor, cursor + edges.n_edges)
        cursor += edges.n_edges
    return out


def test_transport_weights_are_response_blind_and_preserve_equal_system_mass():
    train_a = [_toy_edges(f"tr_{i}", shift=0.07 * i) for i in range(3)]
    train_b = [
        _toy_edges(f"tr_{i}", shift=0.07 * i, turnover_offset=10.0 + i)
        for i in range(3)
    ]
    evaluation = [_toy_edges(f"ev_{i}", shift=0.25 + 0.11 * i) for i in range(3)]
    slices = _species_slices(train_a)

    for mode in ("pooled_inverse_density", "self_inverse_density", "target_density_ratio"):
        w_a = geometry_transport_train_weights(
            train_a, evaluation, mode=mode, bandwidth=0.2
        )
        w_b = geometry_transport_train_weights(
            train_b, evaluation, mode=mode, bandwidth=0.2
        )
        np.testing.assert_allclose(w_a, w_b, atol=0.0, rtol=0.0)
        assert np.isfinite(w_a).all()
        assert np.all(w_a > 0)
        for sl in slices.values():
            assert np.isclose(w_a[sl].sum(), 1.0, atol=1e-12, rtol=0.0)


def _toy_paired_samples():
    samples = []
    theta = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    for i in range(12):
        phase = 0.15 if i < 6 else 0.18
        local = theta + 0.025 * i
        coords = np.column_stack((np.cos(local), np.sin(local)))
        h = 0.5 * (1.0 + np.tanh(np.sin(local - phase) / 0.2))
        a = np.zeros(len(local), dtype=float)
        b = 2.0 * h + 0.05 * np.sin(3.0 * local + i)
        samples.append(
            PairedSpeciesSample(
                species=f"sp_{i:03d}",
                coordinates=coords,
                state_a=a,
                state_b=b,
            )
        )
    return samples


def test_q52_development_scorer_reports_frozen_baseline_and_three_candidates():
    samples = _toy_paired_samples()
    result = q52_geometry_transport_development_test(
        samples,
        train_species=[f"sp_{i:03d}" for i in range(6)],
        eval_species=[f"sp_{i:03d}" for i in range(6, 12)],
        graph_fraction=0.15,
        bandwidth=0.2,
        n_bootstrap=99,
        seed=909,
        edge_chunk_size=8,
        train_chunk_size=128,
        density_chunk_size=64,
    )
    assert tuple(result.modes) == DEVELOPMENT_MODES
    for mode_result in result.modes.values():
        assert np.isfinite(mode_result.statistic)
        assert 0.0 <= mode_result.bootstrap.p_value <= 1.0
        assert len(mode_result.species_scores) == 6
        assert mode_result.mean_train_effective_edge_count > 0
        assert mode_result.min_train_effective_edge_count > 0
