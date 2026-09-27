import numpy as np

from ttf.historical_geographic_information import (
    nested_historical_geographic_information,
)


def make_panel(ns=12, nt=10, seed=61):
    rng = np.random.default_rng(seed)
    source = np.repeat(np.arange(ns), nt)
    target = np.tile(np.arange(ns, ns + nt), ns)
    n_species = ns + nt
    same_class = np.ones(len(source))
    same_order = (
        np.asarray([i % 5 for i in range(n_species)])[source]
        == np.asarray([i % 5 for i in range(n_species)])[target]
    ).astype(float)
    same_family = (
        np.asarray([i % 8 for i in range(n_species)])[source]
        == np.asarray([i % 8 for i in range(n_species)])[target]
    ).astype(float)
    counts = np.asarray([35 + i % 20 for i in range(n_species)], float)
    ratio = counts[source] / counts[target]
    return rng, source, target, same_class, same_order, same_family, ratio


def test_c2_localizes_geography_when_history_tracks_geographic_opportunity():
    rng, s, t, sc, so, sf, ratio = make_panel(seed=62)
    current = rng.normal(size=len(s))
    geography = rng.normal(size=len(s))
    history = (
        0.2 * current
        + 1.5 * geography
        + 0.15 * rng.normal(size=len(s))
    )
    out = nested_historical_geographic_information(
        source_index=s,
        target_index=t,
        r_hist=history,
        r_current=current,
        geographic_coverage=geography,
        same_class=sc,
        same_order=so,
        same_family=sf,
        occurrence_count_ratio=ratio,
    )
    assert out["nested_information"]["geography_retained_fraction"] < 0.2
    assert (
        out["nested_information"][
            "geography_attributable_fraction_of_C1_unique_information"
        ]
        > 0.8
    )


def test_c2_retains_history_when_geography_is_independent():
    rng, s, t, sc, so, sf, ratio = make_panel(seed=63)
    current = rng.normal(size=len(s))
    geography = rng.normal(size=len(s))
    history = rng.normal(size=len(s))
    out = nested_historical_geographic_information(
        source_index=s,
        target_index=t,
        r_hist=history,
        r_current=current,
        geographic_coverage=geography,
        same_class=sc,
        same_order=so,
        same_family=sf,
        occurrence_count_ratio=ratio,
    )
    assert out["nested_information"]["geography_retained_fraction"] > 0.8
    c1 = out["C1_recomputed"]["total_unique_variance_fraction"]
    c2 = out["C2_geography_adjusted"]["total_unique_variance_fraction"]
    assert c2 <= c1 + 1e-12
    assert out["interpretation"]["genetic_response_used"] is False
