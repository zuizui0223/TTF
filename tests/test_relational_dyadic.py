import numpy as np

from ttf.relational_dyadic import batch_primary_test, prepare_dyadic_regression, prepare_two_way_absorber


def complete_bipartite(ns=8, nt=7):
    source = np.repeat(np.arange(ns), nt)
    target = np.tile(np.arange(nt), ns)
    return source, target


def test_exact_two_way_absorption_removes_cluster_means():
    source, target = complete_bipartite()
    rng = np.random.default_rng(4)
    values = rng.normal(size=(len(source), 3))
    absorber = prepare_two_way_absorber(source, target)
    residual = absorber.residualize(values)
    for s in np.unique(source):
        assert np.max(np.abs(residual[source == s].mean(axis=0))) < 1e-10
    for t in np.unique(target):
        assert np.max(np.abs(residual[target == t].mean(axis=0))) < 1e-10


def test_primary_cluster_test_recovers_positive_relation():
    source, target = complete_bipartite(12, 10)
    rng = np.random.default_rng(8)
    x1 = rng.normal(size=len(source))
    x2 = rng.normal(size=len(source))
    x = np.column_stack((x1, x2))
    prepared = prepare_dyadic_regression(source, target, x)
    source_effect = rng.normal(scale=0.3, size=12)
    target_effect = rng.normal(scale=0.3, size=10)
    y = source_effect[source] + target_effect[target] + 0.6 * x1 + 0.1 * x2 + rng.normal(scale=0.2, size=len(source))
    result = batch_primary_test(prepared, y)
    assert result.coefficient[0] > 0.4
    assert result.p_value_one_sided[0] < 0.05
