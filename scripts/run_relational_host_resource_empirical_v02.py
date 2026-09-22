#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from ttf.core import SpeciesEdges
from ttf.phylogatr_confirmatory import (
    choose_one_panel_per_species,
    read_genes_rows,
    scan_phylogatr_phase1,
)
from ttf.phylogatr_empirical import (
    canonical_mask_sha256_from_identity,
    frozen_edge_mean_p_distances,
    read_aligned_nucleotide_identity,
    sequences_by_frozen_locality,
)
from ttf.relational_dyadic import batch_primary_test, prepare_dyadic_regression
from ttf.relational_genetic_empirical import post_ibd_responses, source_only_transfer_scores


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def species_digest(names: list[str]) -> str:
    return hashlib.sha256(("\n".join(sorted(names)) + "\n").encode()).hexdigest()


def zscore(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    sd = float(x.std())
    if not np.isfinite(sd) or sd <= np.finfo(float).eps:
        raise ValueError("cannot z-score constant or non-finite predictor")
    return (x - float(x.mean())) / sd


def read_survivors(path: Path) -> set[str]:
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    if not rows or "species" not in rows[0]:
        raise ValueError("survivor CSV must contain species column")
    names = [str(row["species"]).strip() for row in rows]
    if any(not x for x in names) or len(names) != len(set(names)):
        raise ValueError("survivor species must be unique and non-empty")
    return set(names)


def remap(values: np.ndarray) -> np.ndarray:
    unique = sorted(map(int, np.unique(values)))
    mapping = {value: index for index, value in enumerate(unique)}
    return np.asarray([mapping[int(value)] for value in values], dtype=np.int64)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--design-npz", type=Path, required=True)
    ap.add_argument("--mask-result", type=Path, required=True)
    ap.add_argument("--survivors", type=Path, required=True)
    ap.add_argument("--repair", type=Path, required=True)
    ap.add_argument("--requalification", type=Path, required=True)
    ap.add_argument("--authorization", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    authorization = json.loads(args.authorization.read_text())
    if authorization.get("schema") != "ttf_relational_host_resource_identity_opening_authorization_v0.2":
        raise RuntimeError("unexpected repaired relational identity-opening authorization")
    if authorization.get("status") != "AUTHORIZE_ONE_CONTINUATION_AFTER_RESPONSE_INDEPENDENT_MASK_REPAIR":
        raise RuntimeError("repaired relational empirical continuation is not authorized")

    repair = json.loads(args.repair.read_text())
    if repair.get("schema") != "ttf_relational_host_resource_mask_parser_repair_v0.1":
        raise RuntimeError("unexpected mask-parser repair contract")
    if repair.get("status") != "FROZEN_RESPONSE_INDEPENDENT_MASK_ADMISSIBILITY_REPAIR_BEFORE_ANY_T_ST_OR_BETA_R":
        raise RuntimeError("mask-parser repair is not frozen at the pre-relational-response boundary")
    if repair["response_state_at_repair_freeze"] != {
        "T_st_opened": False,
        "beta_R_opened": False,
        "relational_decision_opened": False,
    }:
        raise RuntimeError("repair was frozen after relational response opening")

    requalification = json.loads(args.requalification.read_text())
    if requalification.get("schema") != "ttf_relational_host_resource_repaired_survivor_requalification_v0.1":
        raise RuntimeError("unexpected repaired survivor requalification")
    if requalification.get("status") != "PASS_TO_REPAIRED_EMPIRICAL_AUTHORIZATION":
        raise RuntimeError("repaired survivor geometry did not pass requalification")
    if requalification.get("gates", {}).get("overall_pass") is not True:
        raise RuntimeError("repaired survivor requalification gate is not PASS")
    if any(bool(v) for v in requalification.get("response_state", {}).values()):
        raise RuntimeError("requalification artifact indicates relational response was opened")
    if requalification["geometry"] != authorization["repaired_geometry"]:
        raise RuntimeError("authorized repaired survivor geometry drift")
    if repair["repaired_survivor_geometry_contract"] != authorization["repair_geometry_contract"]:
        raise RuntimeError("authorized repair geometry contract drift")

    expected = authorization["input_sha256"]
    actual = {
        "design_npz": sha256_path(args.design_npz),
        "mask_result": sha256_path(args.mask_result),
        "survivor_species_csv": sha256_path(args.survivors),
        "genes_txt": sha256_path(args.root / "genes.txt"),
        "cite_txt": sha256_path(args.root / "cite.txt"),
    }
    for key, value in actual.items():
        if value != expected[key]:
            raise RuntimeError(f"authorized input hash drift for {key}")

    mask_result = json.loads(args.mask_result.read_text())
    if mask_result.get("status") != "PASS_TO_SURVIVOR_REQUALIFICATION":
        raise RuntimeError("character-mask result did not pass")
    if any(bool(v) for v in mask_result["response_firewall"].values()):
        raise RuntimeError("character-mask response firewall was already open")
    mask_by_species = {str(row["species"]): row for row in mask_result["species"]}

    original_survivors = read_survivors(args.survivors)
    excluded_by_repair = {
        str(row["species"])
        for row in repair["response_blind_header_audit"]["affected_species"]
    }
    survivors = original_survivors - excluded_by_repair
    if len(survivors) != int(authorization["species"]["mask_survivors"]):
        raise RuntimeError("repaired survivor species count drift")
    if species_digest(list(survivors)) != authorization["species"]["survivor_digest_sha256"]:
        raise RuntimeError("repaired survivor species digest drift")

    data = np.load(args.design_npz, allow_pickle=False)
    species_order = np.asarray(data["species_order"]).astype(str)
    family = np.asarray(data["family"]).astype(str)
    source_index = np.asarray(data["source_index"], dtype=np.int64)
    target_index = np.asarray(data["target_index"], dtype=np.int64)

    pair_mask = np.asarray([
        species_order[s] in survivors and species_order[t] in survivors
        for s, t in zip(source_index, target_index)
    ], dtype=bool)
    target_counts = Counter(map(int, target_index[pair_mask]))
    eligible_targets = {target for target, count in target_counts.items() if count >= 5}
    pair_mask &= np.asarray([int(t) in eligible_targets for t in target_index], dtype=bool)

    source_index = source_index[pair_mask]
    target_index = target_index[pair_mask]
    source_names = species_order[source_index]
    target_names = species_order[target_index]
    pairs = list(zip(map(str, source_names), map(str, target_names)))
    empirical_species = sorted(set(map(str, source_names)) | set(map(str, target_names)))
    sources = sorted(set(map(str, source_names)))
    targets = sorted(set(map(str, target_names)))

    frozen_species = authorization["species"]
    if len(pairs) != int(frozen_species["directed_pairs"]):
        raise RuntimeError("authorized dyad count drift")
    if species_digest(sources) != frozen_species["source_digest_sha256"]:
        raise RuntimeError("authorized source species digest drift")
    if species_digest(targets) != frozen_species["target_digest_sha256"]:
        raise RuntimeError("authorized target species digest drift")
    if species_digest(empirical_species) != frozen_species["empirical_digest_sha256"]:
        raise RuntimeError("authorized empirical species digest drift")

    # Reconstruct only the authorized species panels from response-blind metadata.
    gene_rows = read_genes_rows(args.root / "genes.txt")
    all_species = {str(row.get("species", "")).strip() for row in gene_rows}
    excluded = all_species - set(empirical_species)
    aliases = tuple(authorization["response_contract"]["marker_aliases"])
    scan = scan_phylogatr_phase1(
        args.root,
        aliases=aliases,
        excluded_species=excluded,
        min_localities=12,
        min_endpoint_training_edges=5,
        neighbor_fraction=0.15,
    )
    selected = choose_one_panel_per_species(scan.candidates)
    panel_map = {panel.species: panel for panel in selected}
    if set(panel_map) != set(empirical_species):
        raise RuntimeError("authorized empirical species panels could not be reconstructed exactly")

    # Verify geometry against the response-blind host-resource design before identity use.
    offsets = np.asarray(data["coordinate_offsets"], dtype=np.int64)
    coordinates = np.asarray(data["coordinates"], dtype=float)
    species_position = {name: i for i, name in enumerate(species_order)}
    edge_map: dict[str, SpeciesEdges] = {}
    for name in empirical_species:
        panel = panel_map[name]
        pos = species_position[name]
        expected_coord = coordinates[offsets[pos] : offsets[pos + 1]]
        actual_coord = np.asarray(panel.geometry.coordinates, dtype=float)
        if expected_coord.shape != actual_coord.shape or not np.allclose(
            expected_coord, actual_coord, rtol=0.0, atol=1e-9
        ):
            raise RuntimeError(f"authorized geometry coordinate drift for {name}")
        nodes = np.asarray(panel.geometry.edge_nodes, dtype=np.int64)
        start = actual_coord[nodes[:, 0]]
        end = actual_coord[nodes[:, 1]]
        edge_map[name] = SpeciesEdges(
            species=name,
            nodes=nodes,
            start=start,
            end=end,
            midpoint=0.5 * (start + end),
            length=np.linalg.norm(end - start, axis=1),
            turnover=np.zeros(len(nodes), dtype=float),
        )

    # This is the single authorized continuation after the response-independent
    # duplicate-header mask repair. No T_st or beta_R was produced by the stopped
    # v0.1 attempt, and the repaired species set is already fixed above.
    response_map = {}
    pair_count_summary = {}
    for name in empirical_species:
        panel = panel_map[name]
        mask_row = mask_by_species[name]
        if mask_row.get("status") != "PASS_MASK":
            raise RuntimeError(f"authorized empirical species lacks PASS mask: {name}")
        alignment = read_aligned_nucleotide_identity(panel.fasta_path)
        if canonical_mask_sha256_from_identity(alignment) != str(mask_row["mask_sha256"]):
            raise RuntimeError(f"canonical mask digest drift at identity opening for {name}")
        occurrences = __import__(
            "ttf.phylogatr_confirmatory", fromlist=["read_occurrence_rows"]
        ).read_occurrence_rows(panel.occurrence_path)
        grouped = sequences_by_frozen_locality(
            alignment,
            occurrences,
            panel.canonical_latlon,
        )
        distance = frozen_edge_mean_p_distances(
            grouped,
            panel.geometry.edge_nodes,
            alignment_length=alignment.alignment_length,
            minimum_comparable_fraction=0.50,
        )
        response_map[name] = post_ibd_responses(
            edge_map[name],
            distance.genetic_distance,
            min_training_edges=5,
        )
        pair_count_summary[name] = {
            "edges": int(len(distance.genetic_distance)),
            "valid_sequence_pairs_total": int(np.sum(distance.valid_pair_counts)),
            "valid_sequence_pairs_min_per_edge": int(np.min(distance.valid_pair_counts)),
        }
        # No sequence identity or edge-level genetic distance is serialized.

    t_st = source_only_transfer_scores(
        edge_map,
        response_map,
        pairs,
        bandwidth_km=500.0,
        prior_strength=0.25,
        prior_mean=0.0,
        segment_points=5,
    )

    host = zscore(np.asarray(data["host_resource_jaccard"], dtype=float)[pair_mask])
    coverage = zscore(np.asarray(data["coverage"], dtype=float)[pair_mask])
    centroid = zscore(np.log1p(np.asarray(data["centroid_distance_km"], dtype=float)[pair_mask]))
    locality = zscore(np.abs(np.log(np.asarray(data["locality_count_ratio"], dtype=float)[pair_mask])))
    same_family = (family[source_index] == family[target_index]).astype(float)
    same_family -= float(same_family.mean())
    predictors = np.column_stack((host, coverage, centroid, locality, same_family))

    source_cluster = remap(source_index)
    target_cluster = remap(target_index)
    prepared = prepare_dyadic_regression(
        source_cluster,
        target_cluster,
        predictors,
        primary_index=0,
    )
    expected_condition = float(authorization["relational_model"]["predictor_condition_number"])
    if not np.isclose(prepared.condition_number, expected_condition, rtol=0.0, atol=1e-9):
        raise RuntimeError("authorized survivor predictor condition number drift")
    tested = batch_primary_test(prepared, t_st[:, None])
    beta = float(tested.coefficient[0])
    se = float(tested.standard_error[0])
    z = float(tested.z_score[0])
    p = float(tested.p_value_one_sided[0])
    alpha = float(authorization["relational_model"]["alpha"])
    positive = bool(p <= alpha and beta > 0.0)
    decision = (
        "RELATIONAL_HOST_RESOURCE_POSITIVE"
        if positive
        else "RELATIONAL_HOST_RESOURCE_NULL_WITH_QUALIFIED_POWER"
    )

    quantiles = np.quantile(t_st, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    dyads = [
        {
            "source": source,
            "target": target,
            "T_st": float(score),
            "host_resource_jaccard": float(raw_host),
            "geographic_coverage": float(raw_coverage),
            "centroid_distance_km": float(raw_centroid),
            "locality_count_ratio": float(raw_locality),
            "same_family": bool(sf),
        }
        for source, target, score, raw_host, raw_coverage, raw_centroid, raw_locality, sf in zip(
            source_names,
            target_names,
            t_st,
            np.asarray(data["host_resource_jaccard"], dtype=float)[pair_mask],
            np.asarray(data["coverage"], dtype=float)[pair_mask],
            np.asarray(data["centroid_distance_km"], dtype=float)[pair_mask],
            np.asarray(data["locality_count_ratio"], dtype=float)[pair_mask],
            family[source_index] == family[target_index],
        )
    ]

    payload = {
        "schema": "ttf_relational_host_resource_empirical_result_v0.2",
        "status": "EMPIRICAL_RELATIONAL_RESULT_OPENED_UNDER_REPAIRED_FROZEN_AUTHORIZATION",
        "authorization_sha256": sha256_path(args.authorization),
        "input_sha256": actual,
        "species": {
            "sources": len(sources),
            "targets": len(targets),
            "empirical_union": len(empirical_species),
            "directed_pairs": len(pairs),
        },
        "primary": {
            "estimand": "beta_R for z(host-resource Jaccard)",
            "coefficient": beta,
            "standard_error_two_way_cluster": se,
            "z_score": z,
            "p_value_one_sided": p,
            "alpha": alpha,
            "positive": positive,
        },
        "T_st_distribution": {
            "mean": float(np.mean(t_st)),
            "sd": float(np.std(t_st)),
            "min": float(quantiles[0]),
            "q10": float(quantiles[1]),
            "q25": float(quantiles[2]),
            "median": float(quantiles[3]),
            "q75": float(quantiles[4]),
            "q90": float(quantiles[5]),
            "max": float(quantiles[6]),
        },
        "decision": decision,
        "pair_count_summary_by_species": pair_count_summary,
        "dyads": dyads,
        "outcome_state": {
            "sequence_identity_opened": True,
            "pairwise_genetic_distances_computed_in_memory": True,
            "serialized_sequence_identity": False,
            "serialized_edge_genetic_distance_vectors": False,
            "T_st_computed": True,
            "beta_R_computed": True,
        },
        "repair_provenance": {
            "parent_authorization_v0.1_stopped_before_T_st": True,
            "excluded_duplicate_header_species": sorted(excluded_by_repair),
            "repaired_survivor_requalification_passed": True,
        },
        "one_shot": {
            "continuation_after_response_independent_structural_repair": True,
            "post_result_retuning_allowed": False,
            "result_selection_rerun_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "decision": decision,
        "beta_R": beta,
        "se": se,
        "z": z,
        "p": p,
        "dyads": len(pairs),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
