import numpy as np

from ttf import (
    SpeciesSample,
    balanced_schedule,
    build_species_edges,
    inclusion_counts,
    rank01,
    recurrence_probability,
    split_species,
)


def test_graphs_are_species_local_by_construction():
    a = SpeciesSample(
        "a",
        np.array([[0.0, 0.0], [0.1, 0.0], [0.2, 0.0], [0.3, 0.0]]),
        np.array([0.0, 1.0, 0.0, 1.0]),
    )
    b = SpeciesSample(
        "b",
        np.array([[10.0, 0.0], [10.1, 0.0], [10.2, 0.0], [10.3, 0.0]]),
        np.array([1.0, 0.0, 1.0, 0.0]),
    )
    ea = build_species_edges(a, k=2)
    eb = build_species_edges(b, k=2)
    assert ea.species == "a"
    assert eb.species == "b"
    assert np.all(ea.nodes < len(a.coordinates))
    assert np.all(eb.nodes < len(b.coordinates))


def test_within_species_rank_is_monotone_scale_invariant():
    x = np.array([0.1, 0.4, 0.2, 0.9, 0.7])
    assert np.array_equal(rank01(x), rank01(np.exp(x)))


def test_balanced_schedule_hits_exact_long_run_balance():
    labels = [f"sp{i}" for i in range(13)]
    schedule = balanced_schedule(
        labels,
        n_realizations=17,
        per_realization=7,
        seed=42,
    )
    counts = inclusion_counts(schedule, labels)
    assert len(schedule) == 17
    assert all(len(row) == 7 and len(set(row)) == 7 for row in schedule)
    assert max(counts.values()) - min(counts.values()) <= 1


def test_species_split_is_deterministic_and_disjoint():
    labels = [f"sp{i}" for i in range(12)]
    a = split_species(labels, eval_fraction=0.5, seed=8)
    b = split_species(labels, eval_fraction=0.5, seed=8)
    assert a == b
    assert set(a[0]).isdisjoint(a[1])
    assert set(a[0]) | set(a[1]) == set(labels)


def test_recurrence_is_descriptive_probability_not_single_map():
    fields = np.array(
        [
            [0.9, 0.1, 0.2, 0.3],
            [0.8, 0.7, 0.1, 0.2],
            [0.1, 0.9, 0.8, 0.2],
        ]
    )
    result = recurrence_probability(fields, top_fraction=0.25)
    assert np.allclose(result.probability, [2 / 3, 1 / 3, 0.0, 0.0])
    assert np.array_equal(result.valid_realizations, np.array([3, 3, 3, 3]))
