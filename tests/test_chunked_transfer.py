from __future__ import annotations

import numpy as np

from ttf.batch import score_prepared_batch
from ttf.chunked_transfer import prepare_chunked_transfer, score_chunked_batch
from ttf.core import SpeciesSample, build_species_edges
from ttf.transfer import prepare_transfer


def _edge_set(name: str, seed: int):
    rng = np.random.default_rng(seed)
    coords = rng.normal(size=(9, 3)) * 300.0
    sample = SpeciesSample(
        species=name,
        coordinates=coords,
        trait=np.arange(len(coords), dtype=float),
    )
    return build_species_edges(sample, k=3)


def test_chunked_batch_matches_dense_prepared_transfer():
    train = [_edge_set("train_a", 1), _edge_set("train_b", 2), _edge_set("train_c", 3)]
    evaluation = [_edge_set("eval_a", 4), _edge_set("eval_b", 5)]
    dense = prepare_transfer(
        train,
        evaluation,
        bandwidth=500.0,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )
    chunked = prepare_chunked_transfer(
        train,
        evaluation,
        bandwidth=500.0,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )

    rng = np.random.default_rng(99)
    worlds = 7
    train_turnover = {
        edges.species: rng.normal(size=(edges.n_edges, worlds)) for edges in train
    }
    eval_turnover = {
        edges.species: rng.normal(size=(edges.n_edges, worlds)) for edges in evaluation
    }

    expected = score_prepared_batch(dense, train_turnover, eval_turnover)
    actual = score_chunked_batch(
        chunked,
        train_turnover,
        eval_turnover,
        query_chunk_size=3,
        train_chunk_size=5,
    )

    np.testing.assert_allclose(actual.statistics, expected.statistics, atol=1e-12, rtol=0.0)
    np.testing.assert_array_equal(actual.n_eval_species, expected.n_eval_species)
    assert set(actual.species_scores) == set(expected.species_scores)
    for species in expected.species_scores:
        np.testing.assert_allclose(
            actual.species_scores[species],
            expected.species_scores[species],
            atol=1e-12,
            rtol=0.0,
        )


def test_chunked_batch_rejects_invalid_chunk_sizes():
    train = [_edge_set("train_a", 11)]
    evaluation = [_edge_set("eval_a", 12)]
    prepared = prepare_chunked_transfer(train, evaluation, bandwidth=500.0)
    train_turnover = {train[0].species: np.ones((train[0].n_edges, 2))}
    eval_turnover = {evaluation[0].species: np.ones((evaluation[0].n_edges, 2))}

    for q, p in ((0, 1), (1, 0)):
        try:
            score_chunked_batch(
                prepared,
                train_turnover,
                eval_turnover,
                query_chunk_size=q,
                train_chunk_size=p,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("invalid chunk size was accepted")
