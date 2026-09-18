from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pytest

from ttf.genetic_empirical_score import (
    score_genetic_distance_mapping,
    score_self_detectability_mapping,
    score_total_genetic_distance_mapping,
)
from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_self_detectability import prepare_genetic_self_detectability
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_compact_execution import (
    prepare_phylogatr_compact_cached_transfer,
    prepare_phylogatr_compact_ttf_design,
)


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


def _designs(geometries: dict):
    names = tuple(sorted(geometries))
    kwargs = dict(
        train_species=names[:6],
        eval_species=names[6:],
        bandwidth=500.0,
        prior_strength=0.25,
        segment_points=5,
        min_training_edges=5,
        strength_neighbours=4,
    )
    return (
        prepare_genetic_ttf_design(geometries, **kwargs),
        prepare_phylogatr_compact_ttf_design(geometries, **kwargs),
    )


def test_compact_empirical_primary_and_self_match_historical_reference() -> None:
    from ttf.phylogatr_compact_empirical import (
        score_phylogatr_compact_distance_mapping,
        score_phylogatr_compact_self_mapping,
    )

    geometries = _geometries()
    historical, compact = _designs(geometries)
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.5,
        residual_amplitude=2.0,
        ibd_strength=1.0,
        noise_sd=0.1,
        transition_width=0.2,
        noise_dimensions=2,
        seed=82471,
    )

    expected_primary = score_genetic_distance_mapping(
        historical, world.genetic_distance
    )
    cached = prepare_phylogatr_compact_cached_transfer(
        compact, edge_chunk_size=3, train_chunk_size=17
    )
    observed_primary = score_phylogatr_compact_distance_mapping(
        compact, world.genetic_distance, cached
    )
    assert observed_primary.statistic == pytest.approx(expected_primary.statistic, abs=1e-12)
    assert observed_primary.training_strength == pytest.approx(
        expected_primary.training_strength, abs=1e-12
    )
    assert observed_primary.species_scores == pytest.approx(
        expected_primary.species_scores, abs=1e-12
    )

    historical_self_design = prepare_genetic_self_detectability(
        historical,
        bandwidth=500.0,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    compact_self_design = prepare_genetic_self_detectability(
        compact,
        bandwidth=500.0,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    expected_self = score_self_detectability_mapping(
        historical, historical_self_design, world.genetic_distance
    )
    observed_self = score_phylogatr_compact_self_mapping(
        compact, compact_self_design, world.genetic_distance
    )
    assert observed_self.statistic == pytest.approx(expected_self.statistic, abs=1e-12)
    assert observed_self.species_scores == pytest.approx(
        expected_self.species_scores, abs=1e-12
    )

    expected_total = score_total_genetic_distance_mapping(
        historical, world.genetic_distance
    )
    compact_total = score_total_genetic_distance_mapping(
        compact, world.genetic_distance
    )
    assert compact_total.statistic == pytest.approx(expected_total.statistic, abs=1e-12)
    assert compact_total.species_scores == pytest.approx(
        expected_total.species_scores, abs=1e-12
    )


def test_phase4_empirical_cli_and_authorizer_freeze_compact_execution_surface() -> None:
    empirical = Path("scripts/run_phylogatr_phase4_empirical_test.py").read_text()
    authorizer = Path("scripts/authorize_phylogatr_phase4_identity_opening.py").read_text()

    assert "prepare_phylogatr_compact_ttf_design" in empirical
    assert "prepare_phylogatr_compact_cached_transfer" in empirical
    assert "score_phylogatr_compact_distance_mapping" in empirical
    assert "score_phylogatr_compact_self_mapping" in empirical
    assert "prepare_genetic_ttf_design(" not in empirical
    assert "score_genetic_distance_mapping(" not in empirical
    assert "score_self_detectability_mapping(" not in empirical

    for required in (
        "src/ttf/phylogatr_compact_empirical.py",
        "src/ttf/phylogatr_compact_execution.py",
        "src/ttf/phylogatr_compact_ibd.py",
        "src/ttf/phylogatr_compact_self_detectability.py",
        "scripts/run_phylogatr_phase3_self_reference_shard.py",
        "scripts/run_phylogatr_phase3_self_evaluation_shard.py",
        "docs/supporting/genetic_phylogatr_phase4_compact_execution_v0.1.json",
    ):
        assert required in authorizer

    receipt = json.loads(
        Path("docs/supporting/genetic_phylogatr_phase4_compact_execution_v0.1.json").read_text()
    )
    assert receipt["schema"] == "ttf_genetic_phylogatr_phase4_compact_execution_v0.1"
    assert receipt["status"] == "FROZEN_BEFORE_PHASE4_IDENTITY_OPENING_OR_EMPIRICAL_RESULT"
    assert receipt["trigger"]["fresh_sequence_identity_opened"] is False
    assert receipt["trigger"]["fresh_empirical_ttf_statistic_opened"] is False
    assert all(value is False for value in receipt["firewall"].values())
