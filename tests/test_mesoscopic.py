import numpy as np

from ttf.core import SpeciesEdges
from ttf.mesoscopic import (
    MesoscopicBoundaryField,
    cut_evidence,
    fit_mesoscopic_boundary,
    mesoscopic_bootstrap_test,
    mesoscopic_species_score,
    mesoscopic_transfer,
)


def _edges(species: str, *, signal_cut: float = 1.5) -> SpeciesEdges:
    coords = np.asarray([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    nodes = np.asarray([[0, 1], [1, 2], [2, 3], [0, 2], [1, 3]], dtype=np.int64)
    start = coords[nodes[:, 0]]
    end = coords[nodes[:, 1]]
    midpoint = 0.5 * (start + end)
    length = np.linalg.norm(end - start, axis=1)
    crosses = ((start[:, 0] < signal_cut) & (end[:, 0] > signal_cut)) | (
        (end[:, 0] < signal_cut) & (start[:, 0] > signal_cut)
    )
    # Strictly ordered rank-like scores with a large gap between boundary and
    # non-boundary edges.  The exact values are immaterial to the contract.
    turnover = np.asarray([0.10, 0.90, 0.20, 0.75, 0.85])
    if signal_cut != 1.5:
        base = np.linspace(0.1, 0.5, len(nodes))
        turnover = base + 0.45 * crosses.astype(float)
    return SpeciesEdges(
        species=species,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=midpoint,
        length=length,
        turnover=turnover,
    )


def test_cut_evidence_peaks_at_supported_transition() -> None:
    cuts = np.asarray([0.5, 1.5, 2.5])
    ev = cut_evidence(_edges("s"), cuts)
    assert ev.observable.tolist() == [True, True, True]
    assert int(np.argmax(ev.log_evidence)) == 1
    assert ev.log_evidence[1] > ev.log_evidence[0]
    assert ev.log_evidence[1] > ev.log_evidence[2]


def test_unobservable_cuts_are_not_negative_evidence() -> None:
    cuts = np.asarray([-0.5, 0.5, 1.5, 2.5, 3.5])
    ev = cut_evidence(_edges("s"), cuts)
    assert ev.observable.tolist() == [False, True, True, True, False]
    assert ev.log_evidence[0] == 0.0
    assert ev.log_evidence[-1] == 0.0


def test_uniform_private_baseline_scores_exactly_zero() -> None:
    cuts = np.asarray([0.5, 1.5, 2.5])
    field = MesoscopicBoundaryField(
        cuts=cuts,
        probability=np.full(3, 1.0 / 3.0),
        mean_log_evidence=np.zeros(3),
        opportunity_species=np.full(3, 10),
        prior_strength=1.0,
    )
    assert abs(mesoscopic_species_score(field, _edges("heldout"))) < 1e-12


def test_shared_training_species_concentrate_boundary_probability() -> None:
    cuts = np.asarray([0.5, 1.5, 2.5])
    train = [_edges(f"train-{i}") for i in range(8)]
    field = fit_mesoscopic_boundary(train, cuts, prior_strength=1.0)
    assert int(np.argmax(field.probability)) == 1
    assert field.probability[1] > 1.0 / 3.0
    assert mesoscopic_species_score(field, _edges("eval")) > 0.0


def test_mesoscopic_transfer_and_bootstrap_use_species_macro_average() -> None:
    cuts = np.asarray([0.5, 1.5, 2.5])
    train = [_edges(f"train-{i}") for i in range(8)]
    evaluation = [_edges(f"eval-{i}") for i in range(6)]
    transfer = mesoscopic_transfer(train, evaluation, cuts)
    scores = np.asarray(list(transfer.species_scores.values()), dtype=float)
    assert transfer.n_eval_species == 6
    assert np.isclose(transfer.statistic, scores.mean(), atol=1e-12, rtol=0.0)
    result = mesoscopic_bootstrap_test(
        train,
        evaluation,
        cuts,
        n_bootstrap=99,
        seed=17,
    )
    assert np.isclose(result.bootstrap.observed_mean, transfer.statistic, atol=1e-12, rtol=0.0)
    assert result.p_value <= 0.05
