from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .calibration import CalibrationCell, seed_for
from .core import SpeciesSample, edge_turnover, split_species
from .geometry import SpeciesGeometry, simulate_fixed_geometry_boundary_world
from .geometry_control import length_orthogonalized_turnover
from .inference import centered_species_bootstrap_mean_test
from .nulls import edges_on_fixed_graphs, fixed_graphs
from .private_geometry_control import (
    partial_spearman_controls,
    pooled_affine_frame,
    random_private_transition_propensity,
)
from .transfer import PreparedTransfer, prepare_transfer

V04_LOCKED_CANDIDATE = "eval_geometry_partial"
V04_PROPENSITY_DIRECTIONS = 2048
V04_PROPENSITY_SEED = 20260908
V04_PROPENSITY_TRANSITION_WIDTH = 0.2


@dataclass(frozen=True)
class V04BatchScore:
    statistics: np.ndarray
    species_scores: Mapping[str, np.ndarray]


def score_v04_prepared_batch(
    prepared: PreparedTransfer,
    train_turnover: Mapping[str, np.ndarray],
    eval_turnover: Mapping[str, np.ndarray],
    *,
    eval_length: Mapping[str, np.ndarray],
    eval_propensity: Mapping[str, np.ndarray],
) -> V04BatchScore:
    """Score worlds using the locked v0.4 evaluation-only geometry partial."""
    widths = set()
    for name in prepared.train_species:
        values = np.asarray(train_turnover[name], dtype=float)
        if values.ndim != 2:
            raise ValueError("training turnover must be edges x worlds")
        widths.add(values.shape[1])
    for name in prepared.eval_species:
        values = np.asarray(eval_turnover[name], dtype=float)
        if values.ndim != 2:
            raise ValueError("evaluation turnover must be edges x worlds")
        widths.add(values.shape[1])
    if len(widths) != 1:
        raise ValueError("batch widths disagree")
    width = int(widths.pop())
    if width < 1:
        raise ValueError("empty batch")

    n_train_edges = max(sl.stop for sl in prepared.train_slices.values())
    train_values = np.empty((n_train_edges, width), dtype=float)
    for name in prepared.train_species:
        sl = prepared.train_slices[name]
        values = np.asarray(train_turnover[name], dtype=float)
        if values.shape != (sl.stop - sl.start, width):
            raise ValueError(f"training turnover shape drift for {name}")
        train_values[sl, :] = values

    score_map: dict[str, np.ndarray] = {}
    sums = np.zeros(width, dtype=float)
    counts = np.zeros(width, dtype=np.int64)
    for name in prepared.eval_species:
        projection = prepared.eval_projection[name]
        predicted = projection @ train_values + prepared.eval_prior_offset[name][:, None]
        target = np.asarray(eval_turnover[name], dtype=float)
        if target.shape != predicted.shape:
            raise ValueError(f"evaluation turnover shape drift for {name}")
        length = np.asarray(eval_length[name], dtype=float)
        propensity = np.asarray(eval_propensity[name], dtype=float)
        if length.shape != (predicted.shape[0],) or propensity.shape != length.shape:
            raise ValueError(f"evaluation geometry-control shape drift for {name}")
        scores = np.asarray(
            [
                partial_spearman_controls(
                    predicted[:, column],
                    target[:, column],
                    [length, propensity],
                )
                for column in range(width)
            ],
            dtype=float,
        )
        score_map[name] = scores
        finite = np.isfinite(scores)
        sums[finite] += scores[finite]
        counts[finite] += 1
    if np.any(counts == 0):
        raise ValueError("one or more worlds had no finite held-out species scores")
    return V04BatchScore(statistics=sums / counts, species_scores=score_map)


def run_geometry_calibration_v04_locked_batched(
    geometries: Sequence[SpeciesGeometry],
    *,
    shared_fractions: Sequence[float],
    amplitudes: Sequence[float],
    n_replicates: int,
    n_bootstrap: int,
    train_species: Sequence[str] | None = None,
    eval_species: Sequence[str] | None = None,
    eval_fraction: float = 0.5,
    split_seed: int | None = None,
    k: int = 4,
    max_distance: float | None = None,
    bandwidth: float = 0.2,
    prior_strength: float = 0.25,
    segment_points: int = 5,
    noise_sd: float = 0.8,
    transition_width: float = 0.2,
    alpha: float = 0.05,
    seed: int = 20260909,
    world_batch_size: int = 25,
    propensity_directions: int = V04_PROPENSITY_DIRECTIONS,
    propensity_seed: int = V04_PROPENSITY_SEED,
    propensity_transition_width: float = V04_PROPENSITY_TRANSITION_WIDTH,
) -> tuple[CalibrationCell, ...]:
    """Run the locked v0.4 evaluation-only private-geometry correction.

    Training and field projection are exactly v0.3: training turnover rank is
    orthogonalized only against edge-length rank and the opportunity-corrected
    dense kernel field is unchanged. The sole v0.4 change is the held-out score:
    prediction and raw held-out turnover are rank-residualized against two
    outcome-free geometry controls, edge length and the expected edge contrast
    under random species-private median hyperplanes, before correlation.
    """
    if n_replicates < 1 or n_bootstrap < 99 or world_batch_size < 2:
        raise ValueError("invalid replicate/bootstrap/batch configuration")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0,1)")
    if propensity_directions < 32 or propensity_transition_width <= 0:
        raise ValueError("invalid propensity configuration")

    data = tuple(geometries)
    labels = tuple(sorted(item.species for item in data))
    if len(set(labels)) != len(labels):
        raise ValueError("species labels must be unique")
    if (train_species is None) != (eval_species is None):
        raise ValueError("train_species and eval_species must be supplied together")
    if train_species is None:
        use_split_seed = (
            seed_for(seed, "fixed_geometry_split") if split_seed is None else int(split_seed)
        )
        train, evaluation = split_species(
            labels, eval_fraction=eval_fraction, seed=use_split_seed
        )
    else:
        train = tuple(map(str, train_species))
        evaluation = tuple(map(str, eval_species or ()))
        if not train or not evaluation or set(train) & set(evaluation):
            raise ValueError("non-empty species-disjoint train/evaluation sets required")
        missing = (set(train) | set(evaluation)) - set(labels)
        if missing:
            raise ValueError(f"unknown species in split: {sorted(missing)}")
    if len(evaluation) < 6:
        raise ValueError("v0.4 calibration requires at least six held-out species")

    geometry_map = {item.species: item for item in data}
    used_names = train + evaluation
    templates = [
        SpeciesSample(
            species=name,
            coordinates=geometry_map[name].coordinates,
            trait=np.arange(len(geometry_map[name].coordinates), dtype=float),
            blocks=geometry_map[name].blocks,
        )
        for name in used_names
    ]
    graphs = fixed_graphs(templates, k=k, max_distance=max_distance)
    template_edges = edges_on_fixed_graphs(templates, graphs)
    train_length = {name: template_edges[name].length for name in train}
    eval_length = {name: template_edges[name].length for name in evaluation}

    center, scale, active = pooled_affine_frame(data)
    eval_propensity = {
        name: random_private_transition_propensity(
            geometry_map[name].coordinates,
            graphs[name],
            center=center,
            scale=scale,
            active=active,
            transition_width=propensity_transition_width,
            n_directions=propensity_directions,
            seed=propensity_seed,
        )
        for name in evaluation
    }
    prepared = prepare_transfer(
        [template_edges[name] for name in train],
        [template_edges[name] for name in evaluation],
        bandwidth=bandwidth,
        prior_strength=prior_strength,
        prior_mean=0.0,
        segment_points=segment_points,
    )

    cells: list[CalibrationCell] = []
    for shared in map(float, shared_fractions):
        for amplitude in map(float, amplitudes):
            p_values: list[float] = []
            observed: list[float] = []
            null_means: list[float] = []
            for batch_start in range(0, int(n_replicates), int(world_batch_size)):
                batch_stop = min(batch_start + int(world_batch_size), int(n_replicates))
                train_columns = {name: [] for name in train}
                eval_columns = {name: [] for name in evaluation}
                null_seeds: list[int] = []
                for replicate in range(batch_start, batch_stop):
                    world = simulate_fixed_geometry_boundary_world(
                        data,
                        shared_fraction=shared,
                        amplitude=amplitude,
                        noise_sd=noise_sd,
                        transition_width=transition_width,
                        seed=seed_for(
                            seed,
                            "fixed_geometry_world",
                            shared,
                            amplitude,
                            replicate,
                        ),
                    )
                    sample_map = {sample.species: sample for sample in world.samples}
                    turnover = {
                        name: edge_turnover(sample_map[name], graphs[name])
                        for name in used_names
                    }
                    for name in train:
                        train_columns[name].append(
                            length_orthogonalized_turnover(
                                turnover[name], train_length[name]
                            )
                        )
                    for name in evaluation:
                        eval_columns[name].append(turnover[name])
                    null_seeds.append(
                        seed_for(
                            seed,
                            "fixed_geometry_bootstrap",
                            shared,
                            amplitude,
                            replicate,
                        )
                    )

                scored = score_v04_prepared_batch(
                    prepared,
                    {name: np.column_stack(train_columns[name]) for name in train},
                    {name: np.column_stack(eval_columns[name]) for name in evaluation},
                    eval_length=eval_length,
                    eval_propensity=eval_propensity,
                )
                for column, null_seed in enumerate(null_seeds):
                    species_scores = np.asarray(
                        [
                            scored.species_scores[name][column]
                            for name in evaluation
                            if np.isfinite(scored.species_scores[name][column])
                        ],
                        dtype=float,
                    )
                    if len(species_scores) < 6:
                        raise ValueError("fewer than six finite held-out species scores")
                    bootstrap = centered_species_bootstrap_mean_test(
                        species_scores,
                        n_bootstrap=n_bootstrap,
                        seed=null_seed,
                    )
                    statistic = float(scored.statistics[column])
                    if not np.isclose(
                        bootstrap.observed_mean, statistic, atol=1e-12, rtol=0.0
                    ):
                        raise RuntimeError("v0.4 score mean drifted from species bootstrap")
                    p_values.append(bootstrap.p_value)
                    observed.append(statistic)
                    null_means.append(float(bootstrap.null_means.mean()))

            p = np.asarray(p_values, dtype=float)
            cells.append(
                CalibrationCell(
                    shared_fraction=shared,
                    amplitude=amplitude,
                    n_replicates=int(n_replicates),
                    alpha=float(alpha),
                    rejection_rate=float(np.mean(p <= alpha)),
                    mean_statistic=float(np.mean(observed)),
                    mean_null_statistic=float(np.mean(null_means)),
                    median_p_value=float(np.median(p)),
                )
            )
    return tuple(cells)
