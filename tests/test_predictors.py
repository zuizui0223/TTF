import numpy as np

from ttf import (
    CalibrationCell,
    build_species_edges,
    compete_predictor_spaces,
    qualify_calibration,
    simulate_circular_boundary_world,
    split_species,
)


def test_predictor_space_increments_are_paired_on_same_species_split():
    world = simulate_circular_boundary_world(
        n_species=12,
        records_per_species=30,
        shared_fraction=1.0,
        amplitude=2.0,
        noise_sd=0.5,
        transition_width=0.15,
        seed=88,
    )
    train, evaluation = split_species(
        [s.species for s in world.samples], eval_fraction=0.5, seed=9
    )
    edge_map = {s.species: build_species_edges(s, k=4) for s in world.samples}
    features = {}
    for name, edges in edge_map.items():
        nuisance = edges.length
        informative = edges.turnover + 0.01 * np.arange(edges.n_edges)
        features[name] = np.column_stack((nuisance, informative))
    result = compete_predictor_spaces(
        [edge_map[s] for s in train],
        [edge_map[s] for s in evaluation],
        features,
        spaces={"G": [0], "GE": [0, 1]},
        increments={"E|G": ("GE", "G")},
        ridge=0.1,
    )
    assert np.isclose(
        result.increments["E|G"],
        result.scores["GE"] - result.scores["G"],
    )
    assert result.scores["GE"] > result.scores["G"]


def test_qualification_requires_zero_shared_adversarial_arm_and_power():
    cells = (
        CalibrationCell(0.0, 1.0, 100, 0.05, 0.04, 0.0, 0.0, 0.5),
        CalibrationCell(0.0, 2.0, 100, 0.05, 0.06, 0.0, 0.0, 0.5),
        CalibrationCell(1.0, 2.0, 100, 0.05, 0.85, 0.2, 0.0, 0.01),
    )
    report = qualify_calibration(
        cells,
        moderate_amplitude=2.0,
        type1_ceiling=0.10,
        power_floor=0.80,
    )
    assert report.passed
    assert report.max_zero_shared_rejection == 0.06
    assert report.full_shared_moderate_power == 0.85
