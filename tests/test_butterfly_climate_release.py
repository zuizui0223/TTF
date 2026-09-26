from __future__ import annotations

from ttf.butterfly_climate_release import (
    ExpansionDescriptor,
    freedman_lane_partial_spearman_permutation,
    host_breadth_stratum,
    one_sided_partial_spearman_permutation,
    partial_spearman,
    select_independent_panel,
)


def _d(
    species: str,
    host_families: float,
    host_units: int,
    resolved: int,
    wing: float,
):
    return ExpansionDescriptor(
        species=species,
        host_family_count=host_families,
        host_wgsrpd3_unit_count=host_units,
        resolved_host_species=resolved,
        wing_size_proxy=wing,
        voltinism="M",
    )


def test_host_breadth_strata():
    assert host_breadth_stratum(1) == "1_family"
    assert host_breadth_stratum(2) == "2_families"
    assert host_breadth_stratum(3) == "3_to_5_families"
    assert host_breadth_stratum(5) == "3_to_5_families"
    assert host_breadth_stratum(6) == "6plus_families"


def test_independent_panel_is_balanced_order_invariant_and_excludes_pilot():
    rows = []
    for label, host_families in [
        ("a", 1),
        ("b", 2),
        ("c", 3),
        ("d", 7),
    ]:
        for i in range(6):
            rows.append(
                _d(
                    f"{label}{i} species",
                    host_families,
                    10 + i * 20,
                    max(int(host_families), 1) + i,
                    2.0 + i,
                )
            )
    excluded = {"a0 species", "d0 species"}
    selected = select_independent_panel(
        rows,
        excluded_species=excluded,
        per_stratum=3,
        minimum_native_units=10,
    )
    selected_rev = select_independent_panel(
        list(reversed(rows)),
        excluded_species=excluded,
        per_stratum=3,
        minimum_native_units=10,
    )
    assert selected == selected_rev
    assert not (set(selected) & excluded)
    assert len(selected) == 12
    counts = {name: 0 for name in (
        "1_family",
        "2_families",
        "3_to_5_families",
        "6plus_families",
    )}
    by_species = {row.species: row for row in rows}
    for name in selected:
        counts[host_breadth_stratum(by_species[name].host_family_count)] += 1
    assert set(counts.values()) == {3}


def test_panel_requires_host_taxonomy_lower_bound():
    rows = []
    for host_families in (1, 2, 3, 6):
        for i in range(3):
            rows.append(
                _d(
                    f"{host_families}-{i}",
                    host_families,
                    50 + i,
                    int(host_families),
                    3 + i,
                )
            )
    rows.append(_d("bad", 6, 200, 2, 10))
    selected = select_independent_panel(
        rows,
        excluded_species=(),
        per_stratum=3,
        minimum_native_units=10,
    )
    assert "bad" not in selected


def test_partial_spearman_detects_negative_relation_after_control():
    response = [0.95, 0.85, 0.72, 0.60, 0.45, 0.30]
    predictor = [1, 1, 2, 3, 6, 12]
    control = [20, 100, 40, 80, 30, 120]
    value = partial_spearman(response, predictor, control)
    assert value < -0.8


def test_permutation_test_is_deterministic():
    response = [0.95, 0.85, 0.72, 0.60, 0.45, 0.30]
    predictor = [1, 1, 2, 3, 6, 12]
    control = [20, 100, 40, 80, 30, 120]
    a = one_sided_partial_spearman_permutation(
        response,
        predictor,
        control,
        iterations=199,
    )
    b = one_sided_partial_spearman_permutation(
        response,
        predictor,
        control,
        iterations=199,
    )
    assert a == b
    assert a["observed_partial_spearman"] < 0
    assert 0 < a["one_sided_p_value"] <= 1


def test_freedman_lane_observed_statistic_matches_partial_spearman():
    species = ["a", "b", "c", "d", "e", "f", "g", "h"]
    response = [0.91, 0.82, 0.74, 0.61, 0.55, 0.44, 0.31, 0.20]
    predictor = [1, 1, 2, 2, 3, 5, 7, 10]
    control = [20, 120, 55, 200, 80, 140, 95, 260]
    expected = partial_spearman(response, predictor, control)
    result = freedman_lane_partial_spearman_permutation(
        response,
        predictor,
        control,
        species,
        iterations=199,
    )
    assert abs(result["observed_partial_spearman"] - expected) < 1e-12
    assert result["method"].startswith("Freedman-Lane")


def test_freedman_lane_permutation_is_row_order_invariant():
    species = ["a", "b", "c", "d", "e", "f", "g", "h"]
    response = [0.91, 0.82, 0.74, 0.61, 0.55, 0.44, 0.31, 0.20]
    predictor = [1, 1, 2, 2, 3, 5, 7, 10]
    control = [20, 120, 55, 200, 80, 140, 95, 260]
    a = freedman_lane_partial_spearman_permutation(
        response,
        predictor,
        control,
        species,
        iterations=199,
    )
    order = [6, 2, 7, 0, 5, 3, 1, 4]
    b = freedman_lane_partial_spearman_permutation(
        [response[i] for i in order],
        [predictor[i] for i in order],
        [control[i] for i in order],
        [species[i] for i in order],
        iterations=199,
    )
    assert a == b


def test_freedman_lane_keeps_predictor_control_structure_fixed():
    species = [f"sp{i}" for i in range(10)]
    response = [0.95, 0.84, 0.81, 0.70, 0.62, 0.54, 0.43, 0.34, 0.25, 0.12]
    predictor = [1, 1, 2, 2, 3, 4, 6, 7, 9, 12]
    control = [10, 40, 25, 70, 55, 120, 90, 180, 160, 250]
    result = freedman_lane_partial_spearman_permutation(
        response,
        predictor,
        control,
        species,
        iterations=99,
    )
    assert result["observed_partial_spearman"] < 0
    assert 0 < result["one_sided_p_value"] <= 1
