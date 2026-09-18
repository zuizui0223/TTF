#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from ttf.genetic_empirical_score import score_total_genetic_distance_mapping
from ttf.genetic_self_detectability import prepare_genetic_self_detectability
from ttf.phylogatr_character_mask import canonical_mask_sha256, read_canonical_mask_alignment
from ttf.phylogatr_compact_empirical import (
    score_phylogatr_compact_distance_mapping,
    score_phylogatr_compact_self_mapping,
)
from ttf.phylogatr_compact_execution import (
    prepare_phylogatr_compact_cached_transfer,
    prepare_phylogatr_compact_ttf_design,
)
from ttf.phylogatr_confirmatory import fasta_header_sha256, fasta_headers_only, sha256_path
from ttf.phylogatr_empirical import extract_species_frozen_edge_distances
from ttf.phylogatr_phase4 import load_phylogatr_phase4_context
from ttf.private_null_inference import upper_monte_carlo_pvalue
from ttf.profiled_private_null import profiled_private_pvalue


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run the one authorized fresh phylogatR empirical genetic TTF test."
    )
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--geometry", type=Path, required=True)
    ap.add_argument("--phase1-manifest", type=Path, required=True)
    ap.add_argument("--phase2-manifest", type=Path, required=True)
    ap.add_argument("--phase3-rule", type=Path, required=True)
    ap.add_argument("--phase3-authorization", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--qualification", type=Path, required=True)
    ap.add_argument("--self-rule", type=Path, required=True)
    ap.add_argument("--self-references", type=Path, required=True)
    ap.add_argument("--self-qualification", type=Path, required=True)
    ap.add_argument("--phase4-rule", type=Path, required=True)
    ap.add_argument("--phase4-authorization", type=Path, required=True)
    ap.add_argument("--opening-state", type=Path, required=True)
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    context = load_phylogatr_phase4_context(
        args.geometry,
        args.phase1_manifest,
        args.phase2_manifest,
        args.phase3_rule,
        args.phase3_authorization,
        args.references,
        args.qualification,
        args.self_rule,
        args.self_references,
        args.self_qualification,
        args.phase4_rule,
        args.phase4_authorization,
        args.opening_state,
        verify_code=True,
        repo_root=args.repo_root,
    )

    root = args.root.resolve()
    if sha256_path(root / "cite.txt") != context.phase1["provenance"]["cite_sha256"]:
        raise RuntimeError("fresh cite.txt drift before empirical opening")
    if sha256_path(root / "genes.txt") != context.phase1["provenance"]["genes_sha256"]:
        raise RuntimeError("fresh genes.txt drift before empirical opening")

    used = context.train_species + context.eval_species
    selected = context.phase1["selected_panels"]
    ledger = context.phase2["species_ledger"]
    minimum_fraction = float(
        context.phase4_rule["sequence_contract"][
            "minimum_comparable_fraction_of_frozen_alignment_length"
        ]
    )
    genetic_distance: dict[str, np.ndarray] = {}
    source_checks: dict[str, dict[str, int | bool]] = {}

    for species in used:
        meta = selected[species]
        if species not in ledger or ledger[species].get("survives") is not True:
            raise RuntimeError(f"authorized survivor lacks Phase-2 admissibility ledger: {species}")
        fasta_path = root / str(meta["aligned_fasta_relative"])
        occurrence_path = root / str(meta["occurrence_relative"])
        if not fasta_path.is_file() or not occurrence_path.is_file():
            raise RuntimeError(f"fresh empirical source file missing for {species}")
        if sha256_path(occurrence_path) != str(meta["occurrence_sha256"]):
            raise RuntimeError(f"fresh occurrence source drift for {species}")
        headers = fasta_headers_only(fasta_path)
        if fasta_header_sha256(headers) != str(meta["aligned_header_sha256"]):
            raise RuntimeError(f"fresh aligned FASTA header drift for {species}")

        mask_alignment = read_canonical_mask_alignment(fasta_path)
        mask_sha = canonical_mask_sha256(mask_alignment)
        if mask_sha != str(ledger[species]["canonical_mask_sha256"]):
            raise RuntimeError(f"fresh canonical-mask drift since Phase 2 for {species}")

        distances = extract_species_frozen_edge_distances(
            fasta_path,
            occurrence_path,
            context.frozen_latlon[species],
            context.geometries[species],
            minimum_comparable_fraction=minimum_fraction,
        )
        genetic_distance[species] = np.asarray(distances.genetic_distance, dtype=float)
        source_checks[species] = {
            "header_sha_verified": True,
            "occurrence_sha_verified": True,
            "phase2_mask_sha_verified": True,
            "frozen_edges": int(len(distances.genetic_distance)),
            "minimum_valid_sequence_pairs_per_edge": int(np.min(distances.valid_pair_counts)),
        }

    core = context.phase3_rule["core_method"]
    design = prepare_phylogatr_compact_ttf_design(
        context.geometries,
        train_species=context.train_species,
        eval_species=context.eval_species,
        bandwidth=float(core["bandwidth_km"]),
        prior_strength=float(core["prior_strength"]),
        segment_points=int(core["segment_points"]),
        min_training_edges=int(
            context.phase3_rule["geometry_contract"][
                "minimum_endpoint_disjoint_ibd_training_edges"
            ]
        ),
        strength_neighbours=int(core["strength_neighbours"]),
    )
    execution = context.phase3_authorization["execution"]
    cached_transfer = prepare_phylogatr_compact_cached_transfer(
        design,
        edge_chunk_size=int(execution["edge_chunk_size"]),
        train_chunk_size=int(execution["train_chunk_size"]),
    )
    primary = score_phylogatr_compact_distance_mapping(
        design,
        genetic_distance,
        cached_transfer,
    )
    reference_family = {
        label: (
            np.asarray(payload["training_strength"], dtype=float),
            np.asarray(payload["statistic"], dtype=float),
        )
        for label, payload in context.references["references"].items()
    }
    primary_p, selected_configs, component_p, profile_distances = profiled_private_pvalue(
        primary.statistic,
        primary.training_strength,
        reference_family,
        profile_draws=int(context.phase3_rule["qualification"]["profile_strength_draws"]),
        selected_configs=int(core["profiled_private_selected_configurations"]),
    )
    alpha = float(context.phase4_rule["primary_estimand"]["alpha"])
    primary_positive = bool(float(primary_p) <= alpha)

    self_geometry = context.self_rule["geometry"]
    self_design = prepare_genetic_self_detectability(
        design,
        bandwidth=float(self_geometry["bandwidth_km"]),
        prior_strength=float(self_geometry["prior_strength"]),
        prior_mean=float(self_geometry["prior_mean"]),
        segment_points=int(self_geometry["segment_points"]),
    )
    empirical_self = score_phylogatr_compact_self_mapping(
        design, self_design, genetic_distance
    )
    self_reference = np.asarray(context.self_references["statistics"], dtype=float)
    self_p = float(upper_monte_carlo_pvalue(empirical_self.statistic, self_reference))
    self_alpha = float(context.self_rule["inference"]["alpha"])
    self_method_qualified = bool(context.self_qualification.get("passed") is True)
    empirical_self_positive = bool(self_p <= self_alpha)

    if primary_positive:
        decision = "TRANSFERABLE_PLACE_COMPONENT"
    elif self_method_qualified and empirical_self_positive:
        decision = "LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN"
    else:
        decision = "NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING"

    # The frozen primary decision is complete before this unqualified descriptor.
    total = score_total_genetic_distance_mapping(
        design,
        genetic_distance,
        edge_chunk_size=int(execution["edge_chunk_size"]),
        train_chunk_size=int(execution["train_chunk_size"]),
    )

    out = {
        "schema": "ttf_genetic_phylogatr_phase4_empirical_result_v0.1",
        "status": "EMPIRICAL_RESULT_OPENED_UNDER_FROZEN_PHASE4_AUTHORIZATION",
        "phase4_authorization_sha256": sha256_path(args.phase4_authorization),
        "geometry_fingerprint_sha256": context.phase3_authorization[
            "geometry_fingerprint_sha256"
        ],
        "source_integrity": {
            "cite_sha256_verified": True,
            "genes_sha256_verified": True,
            "all_survivor_headers_occurrences_masks_verified": True,
            "species": source_checks,
        },
        "primary_place_beyond_ibd": {
            "statistic": float(primary.statistic),
            "training_strength": float(primary.training_strength),
            "profiled_private_p_value": float(primary_p),
            "alpha": alpha,
            "positive": primary_positive,
            "selected_configurations": list(selected_configs),
            "component_p_values": {key: float(value) for key, value in component_p.items()},
            "profile_distances": {key: float(value) for key, value in profile_distances.items()},
            "heldout_species_scores": {
                key: float(value) for key, value in sorted(primary.species_scores.items())
            },
        },
        "within_species_self_diagnostic": {
            "synthetic_method_qualified": self_method_qualified,
            "synthetic_qualification_status": context.self_qualification["status"],
            "statistic": float(empirical_self.statistic),
            "p_value": self_p,
            "alpha": self_alpha,
            "positive": empirical_self_positive,
            "species_scores": {
                key: float(value)
                for key, value in sorted(empirical_self.species_scores.items())
            },
        },
        "secondary_total_genetic_transfer": total.as_dict(),
        "decision": decision,
        "interpretation": {
            "TRANSFERABLE_PLACE_COMPONENT": "Evidence that geographic location predicts post-IBD mitochondrial intraspecific differentiation in unseen species within the fresh confirmatory domain.",
            "LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN": "Cross-species transfer is not detected despite qualified and positive within-species spatial detectability; this is consistent with lineage-conditioned spatial structure within the tested domain, not proof of zero transfer or a lineage-specific historical cause.",
            "NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING": "Cross-species transfer is not detected, but within-species evidence is insufficient for a lineage-conditioned interpretation.",
        }[decision],
        "confirmatory_sequence_identity_opened": True,
        "confirmatory_pairwise_genetic_distances_opened": True,
        "confirmatory_ttf_statistic_opened": True,
        "decker_empirical_genetic_outcomes_opened": False,
        "serialized_sequence_identity": False,
        "serialized_edge_genetic_distance_vectors": False,
        "claim_boundary": "This result applies only to the exact frozen fresh phylogatR panel, marker, graph, split, response definition, and prequalified reference family.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "status": out["status"],
                "decision": decision,
                "primary_statistic": primary.statistic,
                "primary_p_value": primary_p,
                "self_method_qualified": self_method_qualified,
                "self_statistic": empirical_self.statistic,
                "self_p_value": self_p,
                "decker_empirical_genetic_outcomes_opened": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
