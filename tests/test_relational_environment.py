import numpy as np

from ttf.relational_environment import (
    canonical_pca,
    deterministic_page_offsets,
    filter_and_thin_occurrences,
    frozen_grid_bounds,
    grid_centres,
    normalized_kde_grid,
    project_whiten,
    relation_pairs,
    schoener_d,
)


def test_page_offsets_small_are_sequential_and_large_are_bounded():
    assert deterministic_page_offsets(0) == ()
    assert deterministic_page_offsets(650) == (0, 300, 600)
    large = deterministic_page_offsets(50_000)
    assert len(large) == 20
    assert large == tuple(sorted(set(large)))
    assert large[0] == 0
    assert large[-1] <= 49_700


def test_occurrence_thinning_is_hash_deterministic_and_caps():
    rows = [
        {"source_key": i, "latitude": 10 + i * 0.2, "longitude": 20.0}
        for i in range(1, 20)
    ]
    a = filter_and_thin_occurrences("Example species", rows, minimum_distance_km=10, maximum_retained=5)
    b = filter_and_thin_occurrences("Example species", list(reversed(rows)), minimum_distance_km=10, maximum_retained=5)
    assert [(x["source_key"], x["priority_rank"]) for x in a] == [(x["source_key"], x["priority_rank"]) for x in b]
    assert len(a) == 5


def test_pca_sign_and_schoener_d_contracts():
    rng = np.random.default_rng(3)
    x = rng.normal(size=(100, 4)) @ np.asarray([
        [1.0, .2, 0, 0],
        [.1, 1.0, .2, 0],
        [0, .2, 1.0, .1],
        [0, 0, .1, 1.0],
    ])
    mean, sd, eigval, eigvec = canonical_pca(x)
    for axis in range(eigvec.shape[1]):
        col = eigvec[:, axis]
        assert col[np.argmax(np.abs(col))] >= 0
    white = project_whiten(x, mean, sd, eigval, eigvec)
    low, high = frozen_grid_bounds(white)
    centres = grid_centres(low, high)
    p = normalized_kde_grid(white[:50], centres)
    q = normalized_kde_grid(white[50:], centres)
    assert np.isclose(schoener_d(p, p), 1.0)
    assert 0 <= schoener_d(p, q) <= 1
    density = np.vstack((p, q))
    s, t, r = relation_pairs(density, [0], [0, 1])
    assert list(s) == [0, 0]
    assert list(t) == [0, 1]
    assert np.isclose(r[0], 1.0)
