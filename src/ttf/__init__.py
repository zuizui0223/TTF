"""Transferable Turnover Fields (TTF).

A cross-system method for testing whether within-system transition structure
learned from some systems predicts transition structure in unseen systems.
"""

from .atlas import RecurrenceEstimate, recurrence_probability
from .calibration import (
    CalibrationCell,
    QualificationReport,
    qualify_calibration,
    run_calibration,
    run_geometry_calibration,
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
from .geometry import (
    GeometrySyntheticWorld,
    SpeciesGeometry,
    geometry_fingerprint,
    simulate_fixed_geometry_boundary_world,
)
from .inference import (
    MeanBootstrapResult,
    SpeciesBootstrapResult,
    centered_species_bootstrap_mean_test,
    heldout_species_bootstrap_test,
)
from .mismatch import (
    PairedSpeciesSample,
    build_coupling_edges,
    build_mismatch_edges,
    default_mismatch,
    default_relative_state,
    mismatch_species_sample,
    pointwise_mismatch,
    pointwise_relative_state,
    relative_species_sample,
)
from .mismatch_inference import (
    PairedHeldoutInferenceResult,
    paired_heldout_species_bootstrap_test,
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
    "GeometrySyntheticWorld",
    "KernelBoundaryModel",
    "MeanBootstrapResult",
    "PairedHeldoutInferenceResult",
    "PairedSpeciesSample",
    "PermutationResult",
    "PredictorCompetitionResult",
    "PredictorTransferResult",
    "PreparedTransfer",
    "QualificationReport",
    "RecurrenceEstimate",
    "SpeciesBootstrapResult",
    "SpeciesEdges",
    "SpeciesGeometry",
    "SpeciesSample",
    "SyntheticWorld",
    "TransferResult",
    "balanced_schedule",
    "build_coupling_edges",
    "build_mismatch_edges",
    "build_species_edges",
    "centered_species_bootstrap_mean_test",
    "edge_turnover",
    "compete_predictor_spaces",
    "default_mismatch",
    "default_relative_state",
    "fit_boundary_model",
    "fixed_graphs",
    "geometry_fingerprint",
    "heldout_species_bootstrap_test",
    "inclusion_counts",
    "knn_edges",
    "mismatch_species_sample",
    "paired_heldout_species_bootstrap_test",
    "permutation_test",
    "permute_trait_within_species",
    "pointwise_mismatch",
    "pointwise_relative_state",
    "prepare_transfer",
    "qualify_calibration",
    "rank01",
    "recurrence_probability",
    "relative_species_sample",
    "ridge_predictor_transfer",
    "run_calibration",
    "run_geometry_calibration",
    "simulate_circular_boundary_world",
    "simulate_fixed_geometry_boundary_world",
    "spearman_rho",
    "split_species",
    "transfer_statistic",
]

__version__ = "0.3.0.dev0"
