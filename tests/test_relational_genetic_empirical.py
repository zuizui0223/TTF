import numpy as np

from ttf.core import SpeciesEdges
from ttf.relational_genetic_empirical import source_only_transfer_scores


def edges(name, shift, turnover):
    start = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]) + shift
    end = start + np.asarray([[0.4, 0.0, 0.0], [0.4, 0.0, 0.0], [0.4, 0.0, 0.0]])
    nodes = np.asarray([[0, 1], [1, 2], [2, 3]])
    return SpeciesEdges(
        species=name,
        nodes=nodes,
        start=start,
        end=end,
        midpoint=0.5 * (start + end),
        length=np.linalg.norm(end - start, axis=1),
        turnover=np.asarray(turnover, dtype=float),
    )


class Response:
    def __init__(self, species, train, evaluation):
        self.species = species
        self.train_turnover = np.asarray(train, dtype=float)
        self.eval_turnover = np.asarray(evaluation, dtype=float)


def test_source_only_transfer_score_is_directed_and_finite():
    source = edges("source", 0.0, [0.0, 0.0, 0.0])
    target = edges("target", 0.1, [0.0, 0.0, 0.0])
    score = source_only_transfer_scores(
        {"source": source, "target": target},
        {
            "source": Response("source", [-0.2, 0.0, 0.2], [-0.2, 0.0, 0.2]),
            "target": Response("target", [-0.2, 0.0, 0.2], [-0.2, 0.0, 0.2]),
        },
        [("source", "target")],
        bandwidth_km=1.0,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )
    assert score.shape == (1,)
    assert np.isfinite(score[0])
    assert score[0] > 0.9
