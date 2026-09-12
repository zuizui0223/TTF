import numpy as np

from ttf.profiled_private_null import profiled_private_pvalue, robust_profile_distance


def test_robust_profile_distance_centres_on_median():
    x = np.linspace(-1.0, 1.0, 101)
    d0, med, scale = robust_profile_distance(0.0, x)
    assert med == 0.0
    assert scale > 0.0
    assert d0 == 0.0
    d1, _, _ = robust_profile_distance(1.0, x)
    assert d1 > d0


def test_profiled_pvalue_uses_strength_not_observed_statistic_to_choose_configs():
    n = 1999
    refs = {
        "weak": (np.linspace(-0.2, 0.2, n), np.linspace(-1.0, 1.0, n)),
        "medium": (np.linspace(0.8, 1.2, n), np.linspace(-0.5, 1.5, n)),
        "strong": (np.linspace(1.8, 2.2, n), np.linspace(0.0, 2.0, n)),
    }
    first = profiled_private_pvalue(
        0.0, 1.0, refs, profile_draws=999, selected_configs=2
    )
    second = profiled_private_pvalue(
        100.0, 1.0, refs, profile_draws=999, selected_configs=2
    )
    assert first[1] == second[1]
    assert first[1][0] == "medium"
    assert set(first[1]) == {"medium", "strong"} or set(first[1]) == {"medium", "weak"}
    assert second[0] <= first[0]


def test_profiled_pvalue_is_max_over_selected_component_pvalues():
    n = 1999
    refs = {
        "a": (np.zeros(n), np.linspace(-1.0, 1.0, n)),
        "b": (np.ones(n), np.linspace(0.0, 2.0, n)),
    }
    p, selected, component, _ = profiled_private_pvalue(
        0.75, 0.1, refs, profile_draws=999, selected_configs=2
    )
    assert set(selected) == {"a", "b"}
    assert p == max(component.values())
