import numpy as np

from ttf.core import SpeciesEdges
from ttf.mesoscopic import cut_evidence
from ttf.mesoscopic_consensus import (
    ConsensusBoundaryField,
    CutPosterior,
    consensus_species_score_from_posterior,
    fit_consensus_boundary_from_posteriors,
    species_cut_posterior,
)


def _posterior(
    species: str,
    observable_ids: list[int],
    probabilities: list[float],
    *,
    n_cuts: int = 4,
) -> CutPosterior:
    cuts = np.arange(n_cuts, dtype=float) + 0.5
    observable = np.zeros(n_cuts, dtype=bool)
    observable[observable_ids] = True
    q = np.zeros(n_cuts, dtype=float)
    q[observable_ids] = np.asarray(probabilities, dtype=float)
    private = np.zeros(n_cuts, dtype=float)
    private[observable_ids] = 1.0 / len(observable_ids)
    return CutPosterior(
        species=species,
        cuts=cuts,
        observable=observable,
        probability=q,
        private_probability=private,
        total_log_evidence=np.zeros(n_cuts),
        n_edges=10,
    )


def _edges(species: str) -> SpeciesEdges:
    coords = np.asarray(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
        dtype=float,
    )
    nodes = np.asarray(
        [[0, 1], [1, 2], [2, 3], [0, 2], [1, 3]],
        dtype=np.int64,
    )
    start = coords[nodes[:, 0]]
    end = coords[nodes[:, 1]]
    return SpeciesEdges(
        species=species,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=np.asarray([0.10, 0.90, 0.20, 0.75, 0.85]),
    )


def test_uniform_private_soft_labels_fit_uniform_field_despite_missing_cuts() -> None:
    posteriors = [
        _posterior("a", [0, 1, 2], [1 / 3, 1 / 3, 1 / 3]),
        _posterior("b", [1, 2, 3], [1 / 3, 1 / 3, 1 / 3]),
        _posterior("c", [0, 2, 3], [1 / 3, 1 / 3, 1 / 3]),
        _posterior("d", [0, 1, 3], [1 / 3, 1 / 3, 1 / 3]),
    ]
    field = fit_consensus_boundary_from_posteriors(posteriors, prior_strength=1.0)
    assert np.allclose(field.probability, np.full(4, 0.25), atol=1e-12, rtol=0.0)
    assert field.iterations == 1


def test_shared_soft_labels_concentrate_consensus_without_species_count_weighting() -> None:
    posteriors = [
        _posterior(f"s{i}", [0, 1, 2, 3], [0.02, 0.94, 0.02, 0.02])
        for i in range(8)
    ]
    field = fit_consensus_boundary_from_posteriors(posteriors, prior_strength=1.0)
    assert int(np.argmax(field.probability)) == 1
    assert field.probability[1] > 0.80
    assert field.n_training_species == 8
    assert field.n_informative_species == 8


def test_geometry_private_soft_label_has_nonpositive_expected_log_gain() -> None:
    cuts = np.arange(4, dtype=float) + 0.5
    field = ConsensusBoundaryField(
        cuts=cuts,
        probability=np.asarray([0.70, 0.10, 0.10, 0.10]),
        logits=np.log(np.asarray([0.70, 0.10, 0.10, 0.10])),
        prior_strength=1.0,
        n_training_species=8,
        n_informative_species=8,
        iterations=10,
    )
    private = _posterior("heldout", [0, 1, 2, 3], [0.25, 0.25, 0.25, 0.25])
    score = consensus_species_score_from_posterior(field, private)
    assert score < 0.0
    expected = -np.sum(0.25 * np.log(0.25 / field.probability))
    assert np.isclose(score, expected, atol=1e-12, rtol=0.0)


def test_uniform_consensus_scores_zero_for_any_heldout_soft_label() -> None:
    cuts = np.arange(4, dtype=float) + 0.5
    field = ConsensusBoundaryField(
        cuts=cuts,
        probability=np.full(4, 0.25),
        logits=np.zeros(4),
        prior_strength=1.0,
        n_training_species=8,
        n_informative_species=8,
        iterations=1,
    )
    heldout = _posterior("heldout", [0, 1, 2], [0.05, 0.90, 0.05])
    assert abs(consensus_species_score_from_posterior(field, heldout)) < 1e-12


def test_species_posterior_restores_total_profile_evidence_before_normalization() -> None:
    edges = _edges("s")
    cuts = np.asarray([0.5, 1.5, 2.5])
    per_edge = cut_evidence(edges, cuts)
    posterior = species_cut_posterior(edges, cuts)
    mask = posterior.observable
    assert np.allclose(
        posterior.total_log_evidence[mask],
        edges.n_edges * per_edge.log_evidence[mask],
        atol=1e-12,
        rtol=0.0,
    )
    assert np.isclose(posterior.probability.sum(), 1.0, atol=1e-12, rtol=0.0)
    assert int(np.argmax(posterior.probability)) == 1
