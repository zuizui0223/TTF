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
    edge_turnover,
    inclusion_counts,
    knn_edges,
    rank01,
    spearman_rho,
    split_species,
)
from .genetics import (
    PairwiseDistanceSample,
    build_pairwise_distance_edges,
    pairwise_sequence_distance_matrix,
    sequence_p_distance,
)
from .inference import (
    MeanBootstrapResult,
    SpeciesBootstrapResult,
    centered_species_bootstrap_mean_test,
    heldout_species_bootstrap_test,
)
from .mesoscopic import (
    CutEvidence,
    MesoscopicBootstrapResult,
    MesoscopicBoundaryField,
    MesoscopicTransferResult,
    cut_evidence,
    fit_mesoscopic_boundary,
    mesoscopic_bootstrap_test,
    mesoscopic_species_score,
    mesoscopic_transfer,
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
    "CutEvidence",
    "FieldEstimate",
    "KernelBoundaryModel",
    "MeanBootstrapResult",
    "MesoscopicBootstrapResult",
    "MesoscopicBoundaryField",
    "MesoscopicTransferResult",
    "PairwiseDistanceSample",
    "PermutationResult",
    "PredictorCompetitionResult",
    "PredictorTransferResult",
    "PreparedTransfer",
    "QualificationReport",
    "RecurrenceEstimate",
    "SpeciesBootstrapResult",
    "SpeciesEdges",
    "SpeciesSample",
    "SyntheticWorld",
    "TransferResult",
    "balanced_schedule",
    "build_pairwise_distance_edges",
    "build_species_edges",
    "centered_species_bootstrap_mean_test",
    "compete_predictor_spaces",
    "cut_evidence",
    "edge_turnover",
    "fit_boundary_model",
    "fit_mesoscopic_boundary",
    "fixed_graphs",
    "heldout_species_bootstrap_test",
    "inclusion_counts",
    "knn_edges",
    "mesoscopic_bootstrap_test",
    "mesoscopic_species_score",
    "mesoscopic_transfer",
    "pairwise_sequence_distance_matrix",
    "permutation_test",
    "permute_trait_within_species",
    "prepare_transfer",
    "qualify_calibration",
    "rank01",
    "recurrence_probability",
    "ridge_predictor_transfer",
    "run_calibration",
    "sequence_p_distance",
    "simulate_circular_boundary_world",
    "spearman_rho",
    "split_species",
    "transfer_statistic",
]

__version__ = "0.2.0"
