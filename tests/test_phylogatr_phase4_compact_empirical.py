from __future__ import annotations

import importlib
import importlib.util
import json
from pathlib import Path

import numpy as np

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


def _panel() -> dict[str, object]:
    out = {}
    for i in range(12):
        x = np.arange(14, dtype=float)
        coordinates = np.column_stack(
            [x * 21.0 + i * 2.0, np.sin(x / 2.0) * 9.0, np.full(14, i)]
        )
        out[f"sp{i:02d}"] = prepare_density_scaled_genetic_geometry(
            coordinates, neighbor_fraction=0.15
        )
    return out


def test_compact_empirical_primary_self_and_total_match_existing_execution() -> None:
    spec = importlib.util.find_spec("ttf.phylogatr_compact_empirical")
    assert spec is not None, "compact empirical adapter must exist before identity opening"
    module = importlib.import_module("ttf.phylogatr_compact_empirical")

    geometries = _panel()
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
    generic = prepare_genetic_ttf_design(geometries, **kwargs)
    compact = prepare_phylogatr_compact_ttf_design(geometries, **kwargs)
    cached = prepare_phylogatr_compact_cached_transfer(
        compact, edge_chunk_size=4, train_chunk_size=32
    )
    world = simulate_genetic_distance_world(
        geometries,
        shared_fraction=0.5,
        residual_amplitude=2.0,
        ibd_strength=1.0,
        noise_sd=0.1,
        transition_width=0.2,
        noise_dimensions=2,
        seed=98231,
    )

    expected_primary = score_genetic_distance_mapping(generic, world.genetic_distance)
    compact_primary = module.score_phylogatr_compact_genetic_distance_mapping(
        compact, world.genetic_distance, cached
    )
    assert np.isclose(compact_primary.statistic, expected_primary.statistic, rtol=0.0, atol=2e-13)
    assert np.isclose(
        compact_primary.training_strength,
        expected_primary.training_strength,
        rtol=0.0,
        atol=2e-13,
    )
    for name in compact.eval_species:
        assert np.isclose(
            compact_primary.species_scores[name],
            expected_primary.species_scores[name],
            rtol=0.0,
            atol=2e-13,
        )
    for name in names:
        assert np.array_equal(
            compact_primary.ibd[name].residual_turnover,
            expected_primary.ibd[name].residual_turnover,
        )

    generic_self = prepare_genetic_self_detectability(
        generic,
        bandwidth=500.0,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    compact_self = prepare_genetic_self_detectability(
        compact,
        bandwidth=500.0,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    expected_self = score_self_detectability_mapping(
        generic, generic_self, world.genetic_distance
    )
    compact_self_score = module.score_phylogatr_compact_self_detectability_mapping(
        compact, compact_self, world.genetic_distance
    )
    assert np.isclose(compact_self_score.statistic, expected_self.statistic, rtol=0.0, atol=2e-13)
    for name in compact.eval_species:
        assert np.isclose(
            compact_self_score.species_scores[name],
            expected_self.species_scores[name],
            rtol=0.0,
            atol=2e-13,
        )

    generic_total = score_total_genetic_distance_mapping(generic, world.genetic_distance)
    compact_total = score_total_genetic_distance_mapping(compact, world.genetic_distance)
    assert compact_total.as_dict() == generic_total.as_dict()


def test_phase4_empirical_cli_is_wired_to_compact_execution() -> None:
    source = Path("scripts/run_phylogatr_phase4_empirical_test.py").read_text()
    assert "prepare_phylogatr_compact_ttf_design" in source
    assert "prepare_phylogatr_compact_cached_transfer" in source
    assert "score_phylogatr_compact_genetic_distance_mapping" in source
    assert "score_phylogatr_compact_self_detectability_mapping" in source
    assert "prepare_genetic_ttf_design(" not in source
    assert "score_genetic_distance_mapping(" not in source
    assert "score_self_detectability_mapping(" not in source


def test_phase4_compact_execution_is_frozen_and_authorized_code_is_pinned() -> None:
    rule = json.loads(
        Path("docs/supporting/genetic_phylogatr_phase4_compact_execution_v0.1.json").read_text()
    )
    assert rule["schema"] == "ttf_genetic_phylogatr_phase4_compact_execution_v0.1"
    assert rule["status"] == (
        "FROZEN_AFTER_GATE_D_AND_SELF_PASS_BEFORE_PHASE4_IDENTITY_OPENING_OR_EMPIRICAL_RESULT"
    )
    assert rule["trigger"]["fresh_sequence_identity_opened"] is False
    assert rule["trigger"]["fresh_empirical_ttf_result_opened"] is False
    assert all(value is False for value in rule["firewall"].values())

    authorizer = Path("scripts/authorize_phylogatr_phase4_identity_opening.py").read_text()
    for required in (
        "src/ttf/phylogatr_compact_ibd.py",
        "src/ttf/phylogatr_compact_execution.py",
        "src/ttf/phylogatr_compact_empirical.py",
        "docs/supporting/genetic_phylogatr_phase4_compact_execution_v0.1.json",
    ):
        assert required in authorizer
