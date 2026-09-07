import numpy as np

from ttf.core import SpeciesEdges
from ttf.mesoscopic import (
    cut_evidence,
    fit_mesoscopic_boundary,
    mesoscopic_species_score,
)
from ttf.mesoscopic_randomization import (
    fit_mesoscopic_boundary_from_evidence,
    mesoscopic_alignment_randomization_test,
    mesoscopic_species_score_from_evidence,
    permute_cut_evidence_labels,
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


def test_alignment_permutation_preserves_observability_and_evidence_multiset() -> None:
    cuts = np.asarray([0.5, 1.5, 2.5])
    evidence = cut_evidence(_edges("s"), cuts)
    randomized = permute_cut_evidence_labels(evidence, np.random.default_rng(4))
    assert np.array_equal(randomized.observable, evidence.observable)
    mask = evidence.observable
    assert np.allclose(
        np.sort(randomized.log_evidence[mask]),
        np.sort(evidence.log_evidence[mask]),
        atol=0.0,
        rtol=0.0,
    )
    assert np.all(randomized.log_evidence[~mask] == 0.0)


def test_evidence_fast_path_matches_direct_mesoscopic_field_and_score() -> None:
    cuts = np.asarray([0.5, 1.5, 2.5])
    train = [_edges(f"train{i}") for i in range(4)]
    heldout = _edges("heldout")
    direct = fit_mesoscopic_boundary(train, cuts, prior_strength=1.0)
    evidence = [cut_evidence(edges, cuts) for edges in train]
    fast = fit_mesoscopic_boundary_from_evidence(evidence, prior_strength=1.0)
    assert np.allclose(direct.probability, fast.probability, atol=1e-12, rtol=0.0)
    assert np.allclose(direct.mean_log_evidence, fast.mean_log_evidence, atol=1e-12, rtol=0.0)
    heldout_evidence = cut_evidence(heldout, cuts)
    assert np.isclose(
        mesoscopic_species_score(direct, heldout),
        mesoscopic_species_score_from_evidence(fast, heldout_evidence),
        atol=1e-12,
        rtol=0.0,
    )


def test_shared_cut_alignment_beats_relabelled_training_null() -> None:
    cuts = np.asarray([0.5, 1.5, 2.5])
    train = [_edges(f"train{i}") for i in range(8)]
    evaluation = [_edges(f"eval{i}") for i in range(8)]
    result = mesoscopic_alignment_randomization_test(
        train,
        evaluation,
        cuts,
        prior_strength=0.0,
        n_randomizations=199,
        seed=20260912,
    )
    assert result.observed_statistic > result.null_mean
    assert result.p_value <= 0.05
    assert result.n_eval_species == 8
