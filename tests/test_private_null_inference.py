import numpy as np

from ttf.private_null_inference import (
    envelope_upper_pvalue,
    envelope_upper_pvalues,
    upper_monte_carlo_pvalue,
)


def test_upper_monte_carlo_pvalue_uses_plus_one_correction():
    reference = np.arange(99, dtype=float)
    assert upper_monte_carlo_pvalue(1000.0, reference) == 0.01
    assert upper_monte_carlo_pvalue(-1.0, reference) == 1.0


def test_envelope_takes_least_favourable_private_configuration():
    refs = {
        "low": np.linspace(-1.0, 1.0, 199),
        "high": np.linspace(0.0, 2.0, 199),
    }
    p, label, component = envelope_upper_pvalue(0.75, refs)
    assert label == "high"
    assert p == max(component.values())
    assert component["high"] > component["low"]


def test_vectorized_envelope_matches_scalar_calls():
    refs = {
        "a": np.linspace(-2.0, 2.0, 199),
        "b": np.linspace(-1.0, 3.0, 199),
    }
    observed = np.array([-0.2, 0.5, 1.7])
    vector, labels = envelope_upper_pvalues(observed, refs)
    for i, value in enumerate(observed):
        scalar, label, _ = envelope_upper_pvalue(float(value), refs)
        assert vector[i] == scalar
        assert labels[i] == label
