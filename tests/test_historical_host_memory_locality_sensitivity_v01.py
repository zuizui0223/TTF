"""Artificial-data checks only; never opens an empirical sequence or distance file."""
from pathlib import Path
import importlib.util
import itertools

import numpy as np
import pytest


def api():
    path = Path(__file__).resolve().parents[1] / "scripts/historical_host_memory_locality_sensitivity.py"
    assert path.is_file(), "locality sensitivity implementation is not yet present"
    spec = importlib.util.spec_from_file_location("hhm_locality_audit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def toy(n=8):
    rng = np.random.default_rng(91204)
    xyz = rng.normal(size=(n, 3)) * 120
    nodes = np.array(list(itertools.combinations(range(n), 2)), dtype=int)
    distance = np.linalg.norm(xyz[nodes[:, 0]] - xyz[nodes[:, 1]], axis=1)
    predictors = rng.normal(size=(len(nodes), 3))
    y = rng.uniform(0.001, 0.09, size=len(nodes))
    return nodes, distance, predictors, y, n


def explicit_residual(nodes, distance, y):
    x = np.log1p(distance)
    result = np.empty(len(y))
    for i, (a, b) in enumerate(nodes):
        train = np.all((nodes != a) & (nodes != b), axis=1)
        assert train.sum() >= 5
        design = np.column_stack((np.ones(train.sum()), x[train]))
        coef = np.linalg.lstsq(design, y[train], rcond=None)[0]
        result[i] = y[i] - np.array([1.0, x[i]]) @ coef
    return result


def expected_slope(nodes, distance, standardized_predictors, y):
    residual = explicit_residual(nodes, distance, y)
    design = np.column_stack((np.ones(len(y)), standardized_predictors))
    return np.linalg.lstsq(design, residual, rcond=None)[0][1]


def test_endpoint_safe_operator_matches_explicit_fits():
    m = api()
    nodes, distance, _, y, _ = toy()
    np.testing.assert_allclose(m.endpoint_disjoint_residual(nodes, distance, y),
                               explicit_residual(nodes, distance, y), rtol=1e-9, atol=1e-12)


def test_every_locality_recomputes_ibd_and_preserves_full_predictor_units():
    m = api()
    nodes, distance, x, y, n = toy()
    z = (x - x.mean(axis=0)) / x.std(axis=0)
    report = m.diagnose_locality_deletions(nodes, distance, x, y, n)
    assert report["baseline"]["unscaled_slope"] == pytest.approx(expected_slope(nodes, distance, z, y), abs=1e-11)
    assert report["n_deletions"] == n
    assert report["n_estimable_deletions"] == n
    for row in report["deletions"]:
        keep = np.all(nodes != row["removed_locality"], axis=1)
        expected = expected_slope(nodes[keep], distance[keep], z[keep], y[keep])
        assert row["status"] == "ESTIMABLE"
        assert row["unscaled_slope"] == pytest.approx(expected, abs=1e-11)
        assert row["retained_edges"] == int(keep.sum())
    assert report["primary_result_replaced"] is False
    assert report["p_value"] is None


def test_removed_locality_cannot_leak_through_full_sample_ibd_residuals():
    m = api()
    nodes, distance, x, y, n = toy()
    changed = y.copy()
    changed[np.any(nodes == 0, axis=1)] += 0.2
    a = m.diagnose_locality_deletions(nodes, distance, x, y, n)
    b = m.diagnose_locality_deletions(nodes, distance, x, changed, n)
    assert a["deletions"][0]["unscaled_slope"] == pytest.approx(b["deletions"][0]["unscaled_slope"], abs=1e-12)
    assert abs(a["baseline"]["unscaled_slope"] - b["baseline"]["unscaled_slope"]) > 1e-8


def test_zero_response_remains_zero_unscaled_and_undefined_standardized():
    m = api()
    nodes, distance, x, y, n = toy()
    report = m.diagnose_locality_deletions(nodes, distance, x, y * 0, n)
    for row in [report["baseline"], *report["deletions"]]:
        assert row["unscaled_slope"] == 0.0
        assert row["standardized_slope"] is None
    assert report["n_sign_reversals"] is None


def test_all_nonestimable_deletions_reported_not_hidden():
    m = api()
    nodes, distance, x, y, n = toy(n=6)
    report = m.diagnose_locality_deletions(nodes, distance, x, y, n)
    assert report["baseline"]["status"] == "ESTIMABLE"
    assert report["n_estimable_deletions"] == 0
    assert report["n_nonestimable_deletions"] == 6
    assert all(r["status"] == "NOT_ESTIMABLE" for r in report["deletions"])
    assert all(r["reason"] == "INSUFFICIENT_ENDPOINT_DISJOINT_EDGES" for r in report["deletions"])
    assert report["min_unscaled_slope"] is None
    assert report["all_deletions_same_nonzero_sign"] is None


def test_predictor_units_do_not_change_when_one_predictor_is_rescaled():
    m = api()
    nodes, distance, x, y, n = toy()
    a = m.diagnose_locality_deletions(nodes, distance, x, y, n)
    xx = x.copy()
    xx[:, 0] = 100 * xx[:, 0] + 37
    b = m.diagnose_locality_deletions(nodes, distance, xx, y, n)
    np.testing.assert_allclose([r["unscaled_slope"] for r in a["deletions"]],
                               [r["unscaled_slope"] for r in b["deletions"]], atol=1e-11)


def test_edge_order_does_not_change_result():
    m = api()
    nodes, distance, x, y, n = toy()
    a = m.diagnose_locality_deletions(nodes, distance, x, y, n)
    b = m.diagnose_locality_deletions(nodes[::-1], distance[::-1], x[::-1], y[::-1], n)
    np.testing.assert_allclose([r["unscaled_slope"] for r in a["deletions"]],
                               [r["unscaled_slope"] for r in b["deletions"]], atol=1e-11)


@pytest.mark.parametrize("problem", ["duplicate", "self", "bad_id", "nan_response", "constant_predictor"])
def test_bad_inputs_fail_instead_of_repairing_the_graph(problem):
    m = api()
    nodes, distance, x, y, n = toy()
    if problem == "duplicate":
        nodes[1] = nodes[0][::-1]
    elif problem == "self":
        nodes[1] = [2, 2]
    elif problem == "bad_id":
        nodes[1, 0] = n
    elif problem == "nan_response":
        y[0] = np.nan
    else:
        x[:, 1] = 1.0
    with pytest.raises((ValueError, RuntimeError)):
        m.diagnose_locality_deletions(nodes, distance, x, y, n)
