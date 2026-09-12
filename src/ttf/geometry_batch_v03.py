from __future__ import annotations

from typing import Sequence

import numpy as np

from .batch import score_prepared_batch
from .calibration import CalibrationCell, seed_for
from .core import SpeciesSample, edge_turnover, split_species
from .geometry import SpeciesGeometry, simulate_fixed_geometry_boundary_world
from .geometry_control import length_orthogonalized_turnover
from .inference import centered_species_bootstrap_mean_test
from .nulls import edges_on_fixed_graphs, fixed_graphs
from .transfer import prepare_transfer

V03_LOCKED_CANDIDATE = "train_orthogonalized_only"


def run_geometry_calibration_v03_locked_batched(
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
    seed: int = 20260908,
    world_batch_size: int = 25,
) -> tuple[CalibrationCell, ...]:
    """Run the prospectively locked v0.3 train-only edge-length correction.

    The only change from v0.2 is on the training response side: within each
    training species, turnover rank is linearly orthogonalized against kNN
    edge-length rank before fitting the shared boundary field. The held-out
    species target and held-out Spearman score are unchanged. The correction is
    centered, so the frozen field prior mean is zero. No coefficient is estimated
    from the fresh validation panel and no evaluation-side residualization is
    performed.
    """
    if n_replicates < 1:
        raise ValueError("n_replicates must be >= 1")
    if n_bootstrap < 99:
        raise ValueError("heldout species bootstrap requires at least 99 resamples")
    if world_batch_size < 2:
        raise ValueError("world_batch_size must be >= 2")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1)")

    data = tuple(geometries)
    labels = tuple(sorted(item.species for item in data))
    if len(set(labels)) != len(labels):
        raise ValueError("species labels must be unique")

    if (train_species is None) != (eval_species is None):
        raise ValueError("train_species and eval_species must be supplied together")
    if train_species is None:
        use_split_seed = (
            seed_for(seed, "fixed_geometry_split")
            if split_seed is None
            else int(split_seed)
        )
        train, evaluation = split_species(
            labels,
            eval_fraction=eval_fraction,
            seed=use_split_seed,
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
        raise ValueError("v0.3 calibration requires at least six held-out species")

    geometry_map = {item.species: item for item in data}
    used_names = train + evaluation
    template_samples = [
        SpeciesSample(
            species=name,
            coordinates=geometry_map[name].coordinates,
            trait=np.arange(len(geometry_map[name].coordinates), dtype=float),
            blocks=geometry_map[name].blocks,
        )
        for name in used_names
    ]
    graphs = fixed_graphs(template_samples, k=k, max_distance=max_distance)
    template_edges = edges_on_fixed_graphs(template_samples, graphs)
    train_length = {name: template_edges[name].length for name in train}
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
                                turnover[name],
                                train_length[name],
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

                scored = score_prepared_batch(
                    prepared,
                    {name: np.column_stack(train_columns[name]) for name in train},
                    {name: np.column_stack(eval_columns[name]) for name in evaluation},
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
                        bootstrap.observed_mean,
                        statistic,
                        atol=1e-12,
                        rtol=0.0,
                    ):
                        raise RuntimeError(
                            "v0.3 batched species-score mean drifted from transfer statistic"
                        )
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
