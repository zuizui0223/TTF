import numpy as np

from ttf.mismatch import (
    PairedSpeciesSample,
    build_coupling_edges,
    build_mismatch_edges,
    default_mismatch,
    pointwise_mismatch,
    pointwise_relative_state,
)


def _coords(n: int) -> np.ndarray:
    return np.column_stack([np.arange(n, dtype=float), np.zeros(n, dtype=float)])


def test_default_mismatch_scalar_and_vector() -> None:
    assert default_mismatch(np.array(2.0), np.array(5.0)) == 3.0
    assert np.isclose(default_mismatch(np.array([0.0, 0.0]), np.array([3.0, 4.0])), 5.0)


def test_common_shift_has_no_coupling_change() -> None:
    a = np.array([0.0, 1.0, 2.0, 3.0])
    b = np.array([0.0, 1.0, 2.0, 3.0])
    sample = PairedSpeciesSample("s", _coords(4), a, b)

    mismatch = pointwise_mismatch(sample)
    relation = pointwise_relative_state(sample)
    assert np.allclose(mismatch, 0.0)
    assert np.allclose(relation, 0.0)

    nodes = np.array([[0, 1], [1, 2], [2, 3]])
    mismatch_edges = build_mismatch_edges(sample, edge_nodes=nodes)
    coupling_edges = build_coupling_edges(sample, edge_nodes=nodes)

    # Constant raw turnover is represented by a constant mid-rank field.
    assert np.allclose(mismatch_edges.turnover, 0.5)
    assert np.allclose(coupling_edges.turnover, 0.5)


def test_asymmetric_shift_changes_mismatch_and_coupling() -> None:
    a = np.array([0.0, 1.0, 2.0, 4.0])
    b = np.array([0.0, 1.0, 2.0, 3.0])
    sample = PairedSpeciesSample("s", _coords(4), a, b)
    nodes = np.array([[0, 1], [1, 2], [2, 3]])

    mismatch_edges = build_mismatch_edges(sample, edge_nodes=nodes)
    coupling_edges = build_coupling_edges(sample, edge_nodes=nodes)

    assert mismatch_edges.turnover[-1] > mismatch_edges.turnover[0]
    assert coupling_edges.turnover[-1] > coupling_edges.turnover[0]


def test_relative_coupling_detects_rotation_at_constant_mismatch_magnitude() -> None:
    # ||A-B|| is one at both sites, but the relative vector rotates 90 degrees.
    a = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 1.0]])
    b = np.zeros_like(a)
    sample = PairedSpeciesSample("s", _coords(3), a, b)
    nodes = np.array([[0, 1], [1, 2]])

    mismatch_edges = build_mismatch_edges(sample, edge_nodes=nodes)
    coupling_edges = build_coupling_edges(sample, edge_nodes=nodes)

    assert np.allclose(pointwise_mismatch(sample), 1.0)
    assert np.allclose(mismatch_edges.turnover, 0.5)
    assert coupling_edges.turnover[0] > coupling_edges.turnover[1]


def test_custom_mismatch_supports_noncommensurable_state_shapes() -> None:
    a = np.array([[1.0, 2.0], [2.0, 2.0], [3.0, 1.0]])
    b = np.array([1.0, 2.0, 4.0])
    sample = PairedSpeciesSample("s", _coords(3), a, b)

    values = pointwise_mismatch(
        sample,
        mismatch=lambda aa, bb: abs(float(np.sum(aa)) - float(bb)),
    )
    assert np.allclose(values, [2.0, 2.0, 0.0])


def test_pairing_shape_drift_is_rejected() -> None:
    try:
        PairedSpeciesSample(
            "s",
            _coords(3),
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0, 2.0]),
        )
    except ValueError as exc:
        assert "share the first dimension" in str(exc)
    else:
        raise AssertionError("expected paired-support validation failure")
