from __future__ import annotations

import numpy as np
import pytest

from ttf.mismatch import PairedSpeciesSample
from ttf.q52_formal import q52_selected_transport_test


def _sample(name: str, phase: float) -> PairedSpeciesSample:
    theta = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    coordinates = np.column_stack((np.cos(theta), np.sin(theta)))
    state_a = np.sin(theta - phase)
    state_b = np.zeros_like(state_a)
    return PairedSpeciesSample(
        species=name,
        coordinates=coordinates,
        state_a=state_a,
        state_b=state_b,
    )


def test_selected_formal_scorer_returns_only_requested_candidate() -> None:
    samples = [
        _sample("sp_000", 0.0),
        _sample("sp_001", 0.1),
        _sample("sp_002", 0.0),
        _sample("sp_003", 0.1),
    ]
    result = q52_selected_transport_test(
        samples,
        train_species=["sp_000", "sp_001"],
        eval_species=["sp_002", "sp_003"],
        candidate="pooled_inverse_density",
        graph_fraction=0.25,
        bandwidth=0.5,
        n_bootstrap=49,
        seed=7,
        density_chunk_size=16,
    )
    assert result.candidate == "pooled_inverse_density"
    assert np.isfinite(result.statistic)
    assert set(result.species_scores) == {"sp_002", "sp_003"}
    assert result.mean_train_effective_edge_count > 0
    assert result.min_train_effective_edge_count > 0


def test_selected_formal_scorer_rejects_unfrozen_candidate() -> None:
    samples = [
        _sample("sp_000", 0.0),
        _sample("sp_001", 0.1),
        _sample("sp_002", 0.0),
        _sample("sp_003", 0.1),
    ]
    with pytest.raises(ValueError, match="not an allowed frozen Q5.2 mode"):
        q52_selected_transport_test(
            samples,
            train_species=["sp_000", "sp_001"],
            eval_species=["sp_002", "sp_003"],
            candidate="q5_1_uniform",  # type: ignore[arg-type]
            graph_fraction=0.25,
            bandwidth=0.5,
            n_bootstrap=49,
            seed=7,
        )
