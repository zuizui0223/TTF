from __future__ import annotations

import numpy as np

from ttf.genetic_gate import prepare_genetic_ttf_design
from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry
from ttf.genetic_self_detectability import (
    prepare_genetic_self_detectability,
    score_genetic_self_world_batch,
)
from ttf.genetic_simulate import simulate_genetic_distance_world
from ttf.phylogatr_compact_execution import prepare_phylogatr_compact_ttf_design
from ttf.phylogatr_compact_self_detectability import (
    score_phylogatr_compact_self_world_batch,
)


def _panel() -> dict[str, object]:
    theta = np.linspace(0.0, 2.0 * np.pi, 18, endpoint=False)
    base = np.column_stack([np.cos(theta), np.sin(theta)])
    return {
        f"sp{i:02d}": prepare_density_scaled_genetic_geometry(
            base + np.array([0.03 * (i % 3), 0.03 * (i // 3)])
        )
        for i in range(12)
    }


def test_compact_self_scoring_matches_existing_reference_exactly() -> None:
    geometries = _panel()
    names = tuple(sorted(geometries))
    train = names[:6]
    evaluation = names[6:]
    kwargs = dict(
        train_species=train,
        eval_species=evaluation,
        bandwidth=1.5,
        prior_strength=0.25,
        segment_points=5,
        min_training_edges=5,
        strength_neighbours=4,
    )
    reference_design = prepare_genetic_ttf_design(geometries, **kwargs)
    compact_design = prepare_phylogatr_compact_ttf_design(geometries, **kwargs)

    reference_self = prepare_genetic_self_detectability(
        reference_design,
        bandwidth=1.5,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )
    compact_self = prepare_genetic_self_detectability(
        compact_design,
        bandwidth=1.5,
        prior_strength=0.25,
        prior_mean=0.5,
        segment_points=5,
    )

    worlds = [
        simulate_genetic_distance_world(
            geometries,
            shared_fraction=0.0,
            residual_amplitude=amplitude,
            ibd_strength=1.0,
            noise_sd=0.10,
            transition_width=0.20,
            noise_dimensions=2,
            seed=seed,
        )
        for amplitude, seed in ((0.0, 7001), (2.0, 7002), (2.0, 7003))
    ]

    reference = score_genetic_self_world_batch(
        reference_design,
        reference_self,
        worlds,
    )
    compact = score_phylogatr_compact_self_world_batch(
        compact_design,
        compact_self,
        worlds,
    )

    assert np.allclose(compact.statistics, reference.statistics, rtol=0.0, atol=2e-13)
    assert np.array_equal(compact.n_species, reference.n_species)
    for name in evaluation:
        assert np.allclose(
            compact.species_scores[name],
            reference.species_scores[name],
            rtol=0.0,
            atol=2e-13,
        )


def test_self_compact_execution_rule_is_frozen_before_self_results() -> None:
    import json
    from pathlib import Path

    payload = json.loads(
        Path(
            "docs/supporting/genetic_phylogatr_phase3_self_compact_execution_v0.1.json"
        ).read_text()
    )
    assert payload["schema"] == "ttf_genetic_phylogatr_phase3_self_compact_execution_v0.1"
    assert payload["status"] == (
        "FROZEN_AFTER_GATE_D_PASS_BEFORE_ANY_SELF_DETECTABILITY_RESULT_OR_NUCLEOTIDE_IDENTITY_OPENING"
    )
    assert payload["trigger"]["formal_gate_d_status"] == "PASS"
    assert payload["trigger"]["self_detectability_result_seen_before_amendment"] is False
    assert payload["trigger"]["fresh_sequence_identity_opened"] is False
    assert all(value is False for value in payload["firewall"].values())
