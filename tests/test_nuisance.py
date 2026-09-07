import numpy as np

from ttf.core import SpeciesEdges, rank01, spearman_rho
from ttf.nuisance import residualize_edge_length


def edge_fixture():
    length = np.linspace(0.1, 2.0, 20)
    nodes = np.column_stack([np.arange(20), np.arange(1, 21)])
    start = np.column_stack([np.arange(20, dtype=float), np.zeros(20)])
    end = start + np.column_stack([length, np.zeros(20)])
    # Strong deterministic IBD trend plus a small non-distance pattern.
    raw = 3.0 * rank01(length) + 0.05 * np.sin(np.arange(20))
    return SpeciesEdges(
        species="sp",
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=length,
        turnover=rank01(raw),
    )


def test_residualization_preserves_geometry_and_reranks_turnover():
    edges = edge_fixture()
    out = residualize_edge_length(edges, degree=1)
    assert np.array_equal(out.nodes, edges.nodes)
    assert np.allclose(out.start, edges.start)
    assert np.allclose(out.end, edges.end)
    assert np.allclose(out.length, edges.length)
    assert np.all((out.turnover > 0) & (out.turnover < 1))
    # Average ranks deliberately retain exact residual ties rather than breaking
    # them by edge order.  Rank-standardized turnover remains centered at 0.5.
    assert np.isclose(out.turnover.mean(), 0.5, atol=1e-12)


def test_residualization_removes_dominant_monotone_edge_length_signal():
    edges = edge_fixture()
    before = abs(spearman_rho(edges.length, edges.turnover))
    after = abs(spearman_rho(edges.length, residualize_edge_length(edges).turnover))
    assert before > 0.95
    assert after < 0.35
