"""Transferable Turnover Fields (TTF).

A cross-taxon method for testing whether within-species trait-transition
structure learned from some species predicts transition structure in unseen
species.
"""

from .atlas import RecurrenceEstimate, recurrence_probability
from .calibration import (
    CalibrationCell,
    QualificationReport,
    qualify_calibration,
    run_calibration,
)
from .core import (
    SpeciesEdges,
    SpeciesSample,
    balanced_schedule,
    build_species_edges,
    inclusion_counts,
    knn_edges,
    rank01,
    spearman_rho,
    split_species,
)
from .nulls import (
    PermutationResult,
    fixed_graphs,
    permutation_test,
    permute_trait_within_species,
)
from .predictors import (
    PredictorCompetitionResult,
    PredictorTransferResult,
    compete_predictor_spaces,
    ridge_predictor_transfer,
)
from .simulate import SyntheticWorld, simulate_circular_boundary_world
from .transfer import (
    FieldEstimate,
    KernelBoundaryModel,
    PreparedTransfer,
    TransferResult,
    fit_boundary_model,
    prepare_transfer,
    transfer_statistic,
)

__all__ = [
    "CalibrationCell",
    "FieldEstimate",
    "KernelBoundaryModel",
    "PermutationResult",
    "PredictorCompetitionResult",
    "PredictorTransferResult",
    "PreparedTransfer",
    "QualificationReport",
    "RecurrenceEstimate",
    "SpeciesEdges",
    "SpeciesSample",
    "SyntheticWorld",
    "TransferResult",
    "balanced_schedule",
    "build_species_edges",
    "compete_predictor_spaces",
    "fit_boundary_model",
    "fixed_graphs",
    "inclusion_counts",
    "knn_edges",
    "permutation_test",
    "permute_trait_within_species",
    "prepare_transfer",
    "qualify_calibration",
    "rank01",
    "recurrence_probability",
    "ridge_predictor_transfer",
    "run_calibration",
    "simulate_circular_boundary_world",
    "spearman_rho",
    "split_species",
    "transfer_statistic",
]

__version__ = "0.1.0"
