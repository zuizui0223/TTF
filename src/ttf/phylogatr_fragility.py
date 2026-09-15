from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median
from typing import Mapping, Sequence

from .genetic_geometry import GeneticSamplingGeometry, prepare_genetic_sampling_geometry
from .genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from .geometry import SpeciesGeometry, geometry_fingerprint
from .heterogeneous_inference import density_scaled_k


_PHASE2_SCHEMA = "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1"
_PHASE3_RULE_SCHEMA = "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1"
_FRAGILITY_RULE_SCHEMA = "ttf_genetic_phylogatr_gate_d_fragility_diagnostic_rule_v0.1"


@dataclass(frozen=True)
class FragilityLevel:
    retention_fraction: float
    geometries: dict[str, GeneticSamplingGeometry]
    source_indices: dict[str, tuple[int, ...]]
    geometry_fingerprint_sha256: str
    metrics: dict[str, float | int]


@dataclass(frozen=True)
class FragilityPlan:
    phase2: dict
    phase3_rule: dict
    fragility_rule: dict
    phase2_manifest_sha256: str
    phase3_rule_sha256: str
    fragility_rule_sha256: str
    train_species: tuple[str, ...]
    eval_species: tuple[str, ...]
    levels: tuple[FragilityLevel, ...]


def _load_json(path: Path, schema: str) -> dict:
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _git_blob_sha1(path: Path) -> str:
    data = Path(path).read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _assert_phase2_blindness(phase2: Mapping[str, object]) -> None:
    if phase2.get("character_mask_opened") is not True:
        raise RuntimeError("fragility diagnostic requires frozen Phase-2 character mask")
    if phase2.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("fragility diagnostic forbidden after fresh nucleotide identity opening")
    if phase2.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("fragility diagnostic forbidden after fresh genetic-distance opening")
    if phase2.get("confirmatory_ttf_statistic_opened") is not False:
        raise RuntimeError("fragility diagnostic forbidden after fresh empirical TTF opening")
    mask = phase2.get("mask_contract")
    if not isinstance(mask, Mapping):
        raise RuntimeError("Phase-2 mask contract missing")
    if mask.get("nucleotide_identity_persisted") is not False:
        raise RuntimeError("Phase-2 mask persisted nucleotide identity")
    if mask.get("pairwise_nucleotide_differences_computed") is not False:
        raise RuntimeError("Phase-2 mask computed pairwise nucleotide differences")


def _rank_localities(
    dataset_digest_sha256: str,
    species: str,
    n_localities: int,
) -> tuple[int, ...]:
    scored: list[tuple[bytes, int]] = []
    for index in range(n_localities):
        payload = (
            f"{dataset_digest_sha256}|fragility|v0.1|{species}|{index}"
        ).encode("utf-8")
        scored.append((hashlib.sha256(payload).digest(), index))
    return tuple(index for _, index in sorted(scored))


def retained_locality_indices(
    dataset_digest_sha256: str,
    species: str,
    n_localities: int,
    retention_fraction: float,
    *,
    minimum_localities: int = 12,
) -> tuple[int, ...]:
    if n_localities < minimum_localities:
        raise RuntimeError(f"{species} has fewer than {minimum_localities} localities")
    if not 0.0 < retention_fraction <= 1.0:
        raise ValueError("retention_fraction must be in (0, 1]")
    target = min(
        n_localities,
        max(minimum_localities, int(math.ceil(retention_fraction * n_localities))),
    )
    ranking = _rank_localities(dataset_digest_sha256, species, n_localities)
    return tuple(sorted(ranking[:target]))


def _fingerprint(geometries: Mapping[str, GeneticSamplingGeometry]) -> str:
    return geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=geometries[name].coordinates)
            for name in sorted(geometries)
        ]
    )


def _metrics(geometries: Mapping[str, GeneticSamplingGeometry]) -> dict[str, float | int]:
    locality_counts = [int(len(geometry.coordinates)) for geometry in geometries.values()]
    graph_k = [int(geometry.graph_k) for geometry in geometries.values()]
    endpoint_support = [
        int(geometry.min_endpoint_disjoint_training_edges)
        for geometry in geometries.values()
    ]
    return {
        "species_count": int(len(geometries)),
        "minimum_localities": int(min(locality_counts)),
        "median_localities": float(median(locality_counts)),
        "maximum_localities": int(max(locality_counts)),
        "minimum_graph_k": int(min(graph_k)),
        "median_graph_k": float(median(graph_k)),
        "maximum_graph_k": int(max(graph_k)),
        "minimum_endpoint_disjoint_ibd_training_edges": int(min(endpoint_support)),
        "median_endpoint_disjoint_ibd_training_edges": float(median(endpoint_support)),
        "maximum_endpoint_disjoint_ibd_training_edges": int(max(endpoint_support)),
    }


def _build_level(
    original: Mapping[str, GeneticSamplingGeometry],
    *,
    dataset_digest_sha256: str,
    retention_fraction: float,
    minimum_localities: int,
    neighbor_fraction: float,
) -> FragilityLevel:
    geometries: dict[str, GeneticSamplingGeometry] = {}
    source_indices: dict[str, tuple[int, ...]] = {}
    for species in sorted(original):
        source = original[species]
        keep = retained_locality_indices(
            dataset_digest_sha256,
            species,
            len(source.coordinates),
            retention_fraction,
            minimum_localities=minimum_localities,
        )
        coordinates = source.coordinates[list(keep)]
        graph_k = density_scaled_k(len(coordinates), fraction=neighbor_fraction)
        geometries[species] = prepare_genetic_sampling_geometry(coordinates, k=graph_k)
        source_indices[species] = keep
    return FragilityLevel(
        retention_fraction=float(retention_fraction),
        geometries=geometries,
        source_indices=source_indices,
        geometry_fingerprint_sha256=_fingerprint(geometries),
        metrics=_metrics(geometries),
    )


def _validate_nested(levels: Sequence[FragilityLevel]) -> None:
    ordered = sorted(levels, key=lambda level: level.retention_fraction, reverse=True)
    for higher, lower in zip(ordered, ordered[1:]):
        if set(higher.source_indices) != set(lower.source_indices):
            raise RuntimeError("fragility stress level species set drift")
        for species in sorted(higher.source_indices):
            if not set(lower.source_indices[species]).issubset(higher.source_indices[species]):
                raise RuntimeError(f"non-nested locality thinning for {species}")


def build_fragility_plan(
    geometry_csv: Path,
    phase2_manifest_path: Path,
    phase3_rule_path: Path,
    fragility_rule_path: Path,
) -> FragilityPlan:
    phase2 = _load_json(phase2_manifest_path, _PHASE2_SCHEMA)
    phase3_rule = _load_json(phase3_rule_path, _PHASE3_RULE_SCHEMA)
    fragility_rule = _load_json(fragility_rule_path, _FRAGILITY_RULE_SCHEMA)

    expected_phase3_blob = str(fragility_rule["formal_phase3_rule_git_blob_sha"])
    observed_phase3_blob = _git_blob_sha1(phase3_rule_path)
    if observed_phase3_blob != expected_phase3_blob:
        raise RuntimeError(
            "formal Phase-3 rule Git blob SHA drift before fragility diagnostic"
        )

    if phase2.get("status") != fragility_rule["input_contract"]["required_phase2_status"]:
        raise RuntimeError("fragility diagnostic requires PASS_TO_SYNTHETIC_GATE Phase 2")
    _assert_phase2_blindness(phase2)
    if sha256_path(geometry_csv) != phase2.get("geometry_csv_sha256"):
        raise RuntimeError("Phase-2 geometry CSV SHA256 drift before fragility diagnostic")

    neighbor_fraction = float(phase3_rule["geometry_contract"]["neighbor_fraction"])
    table = load_frozen_genetic_geometry_csv(
        geometry_csv,
        expected_sha256=str(phase2["geometry_csv_sha256"]),
        expected_neighbor_fraction=neighbor_fraction,
    )
    if set(table.species) != set(map(str, phase2["split"]["train_species"])) | set(
        map(str, phase2["split"]["eval_species"])
    ):
        raise RuntimeError("Phase-2 split does not exactly cover fragility geometry")
    if set(map(str, phase2["split"]["train_species"])) & set(
        map(str, phase2["split"]["eval_species"])
    ):
        raise RuntimeError("Phase-2 train/eval split overlaps")
    if int(phase2["species"]["survivors"]) != len(table.species):
        raise RuntimeError("Phase-2 survivor count disagrees with fragility geometry")
    if len(table.species) < int(phase3_rule["geometry_contract"]["minimum_surviving_species"]):
        raise RuntimeError("Phase-2 survivor panel is below the formal Phase-3 minimum")

    full_fingerprint = _fingerprint(table.geometries)
    if full_fingerprint != phase2.get("geometry_fingerprint_sha256"):
        raise RuntimeError("Phase-2 geometry fingerprint drift before fragility diagnostic")

    minimum_phase3_support = int(
        phase3_rule["geometry_contract"]["minimum_endpoint_disjoint_ibd_training_edges"]
    )
    if min(
        geometry.min_endpoint_disjoint_training_edges for geometry in table.geometries.values()
    ) < minimum_phase3_support:
        raise RuntimeError("full Phase-2 geometry lacks the formal Phase-3 IBD support minimum")

    stress = fragility_rule["stress_design"]
    minimum_localities = int(stress["minimum_localities_per_species"])
    fractions = tuple(float(value) for value in stress["retention_fractions"])
    if fractions != tuple(sorted(set(fractions), reverse=True)):
        raise RuntimeError("fragility retention fractions must be unique and descending")
    if fractions[0] != 1.0:
        raise RuntimeError("fragility design must begin at retention_fraction=1.0")

    dataset_digest = str(phase2["phase1"]["dataset_digest_sha256"])
    levels = tuple(
        _build_level(
            table.geometries,
            dataset_digest_sha256=dataset_digest,
            retention_fraction=fraction,
            minimum_localities=minimum_localities,
            neighbor_fraction=neighbor_fraction,
        )
        for fraction in fractions
    )
    _validate_nested(levels)
    if levels[0].geometry_fingerprint_sha256 != full_fingerprint:
        raise RuntimeError("retention=1.0 fragility geometry is not the exact Phase-2 geometry")

    return FragilityPlan(
        phase2=phase2,
        phase3_rule=phase3_rule,
        fragility_rule=fragility_rule,
        phase2_manifest_sha256=sha256_path(phase2_manifest_path),
        phase3_rule_sha256=sha256_path(phase3_rule_path),
        fragility_rule_sha256=sha256_path(fragility_rule_path),
        train_species=tuple(map(str, phase2["split"]["train_species"])),
        eval_species=tuple(map(str, phase2["split"]["eval_species"])),
        levels=levels,
    )


def write_fragility_plan(plan: FragilityPlan, output_dir: Path) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    level_receipts: list[dict] = []

    for level in plan.levels:
        label = f"{level.retention_fraction:.3f}".rstrip("0").rstrip(".").replace(".", "p")
        csv_path = output_dir / f"geometry_retention_{label}.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "species",
                    "locality_index",
                    "source_locality_index",
                    "x_km",
                    "y_km",
                    "z_km",
                    "graph_k",
                ],
            )
            writer.writeheader()
            for species in sorted(level.geometries):
                geometry = level.geometries[species]
                source = level.source_indices[species]
                for local_index, original_index in enumerate(source):
                    x_km, y_km, z_km = map(float, geometry.coordinates[local_index])
                    writer.writerow(
                        {
                            "species": species,
                            "locality_index": local_index,
                            "source_locality_index": original_index,
                            "x_km": repr(x_km),
                            "y_km": repr(y_km),
                            "z_km": repr(z_km),
                            "graph_k": int(geometry.graph_k),
                        }
                    )
        level_receipts.append(
            {
                "retention_fraction": level.retention_fraction,
                "geometry_csv": csv_path.name,
                "geometry_csv_sha256": sha256_path(csv_path),
                "geometry_fingerprint_sha256": level.geometry_fingerprint_sha256,
                "metrics": level.metrics,
                "source_locality_indices": {
                    species: list(indices)
                    for species, indices in sorted(level.source_indices.items())
                },
            }
        )

    receipt = {
        "schema": "ttf_genetic_phylogatr_gate_d_fragility_geometry_plan_v0.1",
        "status": "response_blind_diagnostic_geometry_plan_only",
        "phase2_manifest_sha256": plan.phase2_manifest_sha256,
        "phase3_rule_sha256": plan.phase3_rule_sha256,
        "fragility_rule_sha256": plan.fragility_rule_sha256,
        "phase2_geometry_csv_sha256": plan.phase2["geometry_csv_sha256"],
        "phase2_geometry_fingerprint_sha256": plan.phase2["geometry_fingerprint_sha256"],
        "dataset_digest_sha256": plan.phase2["phase1"]["dataset_digest_sha256"],
        "species": {
            "survivors": int(plan.phase2["species"]["survivors"]),
            "train": len(plan.train_species),
            "eval": len(plan.eval_species),
            "train_species": list(plan.train_species),
            "eval_species": list(plan.eval_species),
        },
        "levels": level_receipts,
        "authority_firewall": dict(plan.fragility_rule["authority_firewall"]),
        "formal_decision_rule": plan.fragility_rule["formal_decision_rule"],
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "formal_gate_d_decision_made_by_this_receipt": False,
        "phase4_identity_opening_authorized_by_this_receipt": False,
    }
    receipt_path = output_dir / "fragility_geometry_plan_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt_path


__all__ = [
    "FragilityLevel",
    "FragilityPlan",
    "build_fragility_plan",
    "retained_locality_indices",
    "write_fragility_plan",
]
