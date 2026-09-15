from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Mapping

from .genetic_geometry import GeneticSamplingGeometry
from .genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path
from .geometry import SpeciesGeometry, geometry_fingerprint


_PHASE2_SCHEMA = "ttf_genetic_phylogatr_confirmatory_phase2_mask_v0.1"
_PHASE3_RULE_SCHEMA = "ttf_genetic_phylogatr_phase3_gate_d_rule_v0.1"
_GEOMETRY_PLAN_SCHEMA = "ttf_genetic_phylogatr_gate_d_fragility_geometry_plan_v0.1"
_EXECUTION_RULE_SCHEMA = "ttf_genetic_phylogatr_gate_d_fragility_execution_rule_v0.1"
_FORMAL_QUALIFICATION_SCHEMA = "ttf_genetic_phylogatr_phase3_qualification_v0.1"
_EXECUTION_PLAN_SCHEMA = "ttf_genetic_phylogatr_gate_d_fragility_execution_plan_v0.1"


@dataclass(frozen=True)
class FragilityExecutionContext:
    geometries: dict[str, GeneticSamplingGeometry]
    train_species: tuple[str, ...]
    eval_species: tuple[str, ...]
    retention_fraction: float
    retention_label: str
    master_seed: int
    level: dict
    plan: dict
    phase3_rule: dict


def _load_json(path: Path, schema: str) -> dict:
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") != schema:
        raise RuntimeError(f"unexpected schema for {path}: {payload.get('schema')!r}")
    return payload


def _assert_phase2_blindness(phase2: Mapping[str, object]) -> None:
    if phase2.get("character_mask_opened") is not True:
        raise RuntimeError("fragility execution requires frozen Phase-2 character mask")
    if phase2.get("confirmatory_sequence_identity_opened") is not False:
        raise RuntimeError("fragility execution forbidden after fresh nucleotide identity opening")
    if phase2.get("confirmatory_pairwise_genetic_distances_opened") is not False:
        raise RuntimeError("fragility execution forbidden after fresh genetic-distance opening")
    if phase2.get("confirmatory_ttf_statistic_opened") is not False:
        raise RuntimeError("fragility execution forbidden after fresh empirical TTF opening")
    mask = phase2.get("mask_contract")
    if not isinstance(mask, Mapping):
        raise RuntimeError("Phase-2 mask contract missing")
    if mask.get("nucleotide_identity_persisted") is not False:
        raise RuntimeError("Phase-2 mask persisted nucleotide identity")
    if mask.get("pairwise_nucleotide_differences_computed") is not False:
        raise RuntimeError("Phase-2 mask computed pairwise nucleotide differences")


def canonical_retention_label(retention_fraction: float) -> str:
    value = float(retention_fraction)
    if not 0.0 < value <= 1.0:
        raise ValueError("retention fraction must be in (0, 1]")
    return f"{value:.6f}".rstrip("0").rstrip(".")


def derive_fragility_level_master_seed(
    dataset_digest_sha256: str,
    full_geometry_fingerprint_sha256: str,
    retention_fraction: float,
    thinned_geometry_fingerprint_sha256: str,
) -> int:
    label = canonical_retention_label(retention_fraction)
    payload = (
        "ttf|phylogatr|fragility|v0.1|"
        + str(dataset_digest_sha256)
        + "|"
        + str(full_geometry_fingerprint_sha256)
        + "|"
        + label
        + "|"
        + str(thinned_geometry_fingerprint_sha256)
    ).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _ranges(total: int, size: int) -> list[tuple[int, int]]:
    if total < 1 or size < 1:
        raise ValueError("total and shard size must be positive")
    return [(start, min(start + size, total)) for start in range(0, total, size)]


def _formal_anchor(qualification: Mapping[str, object]) -> dict:
    cells = qualification.get("cells")
    if not isinstance(cells, list):
        raise RuntimeError("formal Phase-3 qualification lacks cells")
    a3 = None
    shared_a2 = None
    for cell in cells:
        if not isinstance(cell, Mapping):
            continue
        shared = float(cell["shared_fraction"])
        amplitude = float(cell["residual_amplitude"])
        if shared == 0.0 and amplitude == 3.0:
            a3 = cell
        if shared == 1.0 and amplitude == 2.0:
            shared_a2 = cell
    if a3 is None or shared_a2 is None:
        raise RuntimeError("formal qualification lacks A3 private-null or shared-A2 cell")
    type1 = qualification.get("type1_gate")
    power = qualification.get("power_gate")
    if not isinstance(type1, Mapping) or not isinstance(power, Mapping):
        raise RuntimeError("formal qualification lacks frozen gate summaries")
    return {
        "retention_fraction": 1.0,
        "source": "formal_phase3_qualification_receipt",
        "formal_status": str(qualification.get("status")),
        "A3_private_null_rejection_rate": float(a3["rejection_rate"]),
        "A3_private_null_Wilson95_upper": float(a3["wilson95_high"]),
        "max_private_null_Wilson95_upper": float(type1["max_observed_wilson95_upper"]),
        "shared_A2_rejection_rate": float(shared_a2["rejection_rate"]),
        "shared_A2_Wilson95_lower": float(shared_a2["wilson95_low"]),
        "distance_of_max_private_null_upper_to_formal_0.10_ceiling": (
            float(type1["max_observed_wilson95_upper"])
            - float(type1["wilson95_upper_ceiling"])
        ),
        "distance_of_shared_A2_lower_to_formal_0.80_floor": (
            float(power["observed_wilson95_lower"])
            - float(power["wilson95_lower_floor"])
        ),
        "formal_gate_d_decision_authority": True,
        "diagnostic_decision_authority": False,
    }


def build_fragility_execution_plan(
    geometry_plan_receipt_path: Path,
    phase2_manifest_path: Path,
    phase3_rule_path: Path,
    formal_qualification_path: Path,
    execution_rule_path: Path,
) -> dict:
    geometry_plan_path = Path(geometry_plan_receipt_path)
    phase2_path = Path(phase2_manifest_path)
    phase3_path = Path(phase3_rule_path)
    formal_path = Path(formal_qualification_path)
    execution_path = Path(execution_rule_path)

    geometry_plan = _load_json(geometry_plan_path, _GEOMETRY_PLAN_SCHEMA)
    phase2 = _load_json(phase2_path, _PHASE2_SCHEMA)
    phase3_rule = _load_json(phase3_path, _PHASE3_RULE_SCHEMA)
    formal = _load_json(formal_path, _FORMAL_QUALIFICATION_SCHEMA)
    execution_rule = _load_json(execution_path, _EXECUTION_RULE_SCHEMA)

    if geometry_plan.get("status") != "response_blind_diagnostic_geometry_plan_only":
        raise RuntimeError("fragility geometry plan is not frozen diagnostic-only input")
    _assert_phase2_blindness(phase2)
    if geometry_plan.get("phase2_geometry_csv_sha256") != phase2.get("geometry_csv_sha256"):
        raise RuntimeError("fragility plan/Phase-2 geometry SHA drift")
    if geometry_plan.get("phase2_geometry_fingerprint_sha256") != phase2.get(
        "geometry_fingerprint_sha256"
    ):
        raise RuntimeError("fragility plan/Phase-2 geometry fingerprint drift")
    if geometry_plan.get("phase2_manifest_sha256") != sha256_path(phase2_path):
        raise RuntimeError("fragility geometry plan Phase-2 manifest provenance drift")
    if geometry_plan.get("phase3_rule_sha256") != sha256_path(phase3_path):
        raise RuntimeError("fragility geometry plan formal Phase-3 rule provenance drift")

    if formal.get("status") not in {"PASS", "NOT_EVALUABLE"}:
        raise RuntimeError("formal full-geometry Phase-3 qualification is incomplete")
    if formal.get("geometry_fingerprint_sha256") != phase2.get("geometry_fingerprint_sha256"):
        raise RuntimeError("formal qualification is not for exact full Phase-2 geometry")
    for key in (
        "confirmatory_sequence_identity_opened",
        "confirmatory_pairwise_genetic_distances_opened",
        "confirmatory_ttf_statistic_opened",
    ):
        if formal.get(key) is not False:
            raise RuntimeError(f"formal qualification firewall is open for {key}")

    qualification = phase3_rule["qualification"]
    execution = phase3_rule["execution_defaults"]
    reference_n = int(qualification["reference_worlds_per_configuration"])
    observed_n = int(qualification["observed_worlds_per_cell"])
    reference_ranges = _ranges(reference_n, int(execution["reference_shard_size"]))
    observed_ranges = _ranges(observed_n, int(execution["observed_shard_size"]))
    expected_retentions = [
        float(value)
        for value in execution_rule["thinned_geometry_policy"]["run_retention_fractions"]
    ]

    by_retention = {
        float(level["retention_fraction"]): level for level in geometry_plan["levels"]
    }
    if 1.0 not in by_retention:
        raise RuntimeError("fragility geometry plan lacks retention=1.0 anchor")
    minimum_support = int(
        phase3_rule["geometry_contract"]["minimum_endpoint_disjoint_ibd_training_edges"]
    )
    dataset_digest = str(geometry_plan["dataset_digest_sha256"])
    full_fingerprint = str(geometry_plan["phase2_geometry_fingerprint_sha256"])

    levels: list[dict] = []
    for retention in expected_retentions:
        if retention not in by_retention:
            raise RuntimeError(f"fragility geometry plan lacks frozen retention={retention}")
        geometry = dict(by_retention[retention])
        observed_support = int(
            geometry["metrics"]["minimum_endpoint_disjoint_ibd_training_edges"]
        )
        runnable = observed_support >= minimum_support
        label = canonical_retention_label(retention)
        master_seed = derive_fragility_level_master_seed(
            dataset_digest,
            full_fingerprint,
            retention,
            str(geometry["geometry_fingerprint_sha256"]),
        )
        reference_jobs = []
        observed_jobs = []
        if runnable:
            reference_jobs = [
                {"configuration": config, "start": start, "stop": stop}
                for config in phase3_rule["synthetic_world"][
                    "private_reference_configurations"
                ]
                for start, stop in reference_ranges
            ]
            observed_jobs = [
                {
                    "shared_fraction": float(cell[0]),
                    "residual_amplitude": float(cell[1]),
                    "start": start,
                    "stop": stop,
                }
                for cell in qualification["mandatory_primary_cells"]
                for start, stop in observed_ranges
            ]
        levels.append(
            {
                "retention_fraction": retention,
                "retention_label": label,
                "status": (
                    "SYNTHETIC_DIAGNOSTIC_RUNNABLE"
                    if runnable
                    else "STRUCTURAL_SUPPORT_BELOW_FORMAL_MINIMUM"
                ),
                "geometry_csv": geometry["geometry_csv"],
                "geometry_csv_sha256": geometry["geometry_csv_sha256"],
                "geometry_fingerprint_sha256": geometry[
                    "geometry_fingerprint_sha256"
                ],
                "metrics": geometry["metrics"],
                "minimum_formal_endpoint_disjoint_support": minimum_support,
                "master_seed": master_seed,
                "reference_jobs": reference_jobs,
                "observed_jobs": observed_jobs,
                "reference_job_count": len(reference_jobs),
                "observed_job_count": len(observed_jobs),
                "formal_gate_d_decision_authority": False,
                "phase4_identity_opening_authority": False,
            }
        )

    return {
        "schema": _EXECUTION_PLAN_SCHEMA,
        "status": "DIAGNOSTIC_SYNTHETIC_EXECUTION_PLAN_ONLY",
        "geometry_plan_receipt_sha256": sha256_path(geometry_plan_path),
        "phase2_manifest_sha256": sha256_path(phase2_path),
        "phase3_rule_sha256": sha256_path(phase3_path),
        "formal_qualification_sha256": sha256_path(formal_path),
        "execution_rule_sha256": sha256_path(execution_path),
        "dataset_digest_sha256": dataset_digest,
        "full_geometry_fingerprint_sha256": full_fingerprint,
        "species": geometry_plan["species"],
        "formal_anchor": _formal_anchor(formal),
        "levels": levels,
        "authority_firewall": execution_rule["authority_firewall"],
        "confirmatory_sequence_identity_opened": False,
        "confirmatory_pairwise_genetic_distances_opened": False,
        "confirmatory_ttf_statistic_opened": False,
        "formal_gate_d_decision_made_by_this_plan": False,
        "phase4_identity_opening_authorized_by_this_plan": False,
    }


def load_fragility_execution_context(
    execution_plan_path: Path,
    geometry_dir: Path,
    phase2_manifest_path: Path,
    phase3_rule_path: Path,
    retention_fraction: float,
) -> FragilityExecutionContext:
    plan = _load_json(execution_plan_path, _EXECUTION_PLAN_SCHEMA)
    phase2 = _load_json(phase2_manifest_path, _PHASE2_SCHEMA)
    phase3_rule = _load_json(phase3_rule_path, _PHASE3_RULE_SCHEMA)
    _assert_phase2_blindness(phase2)
    if sha256_path(phase2_manifest_path) != plan["phase2_manifest_sha256"]:
        raise RuntimeError("fragility execution Phase-2 manifest SHA drift")
    if sha256_path(phase3_rule_path) != plan["phase3_rule_sha256"]:
        raise RuntimeError("fragility execution formal Phase-3 rule SHA drift")
    retention = float(retention_fraction)
    candidates = [
        level for level in plan["levels"] if float(level["retention_fraction"]) == retention
    ]
    if len(candidates) != 1:
        raise RuntimeError(f"fragility execution plan lacks unique retention={retention}")
    level = candidates[0]
    if level["status"] != "SYNTHETIC_DIAGNOSTIC_RUNNABLE":
        raise RuntimeError(f"fragility retention={retention} is not synthetic-runnable")

    geometry_path = Path(geometry_dir) / str(level["geometry_csv"])
    if sha256_path(geometry_path) != level["geometry_csv_sha256"]:
        raise RuntimeError("fragility thinned geometry CSV SHA drift")
    neighbor_fraction = float(phase3_rule["geometry_contract"]["neighbor_fraction"])
    table = load_frozen_genetic_geometry_csv(
        geometry_path,
        expected_sha256=str(level["geometry_csv_sha256"]),
        expected_neighbor_fraction=neighbor_fraction,
    )
    fingerprint = geometry_fingerprint(
        [
            SpeciesGeometry(species=name, coordinates=table.geometries[name].coordinates)
            for name in sorted(table.geometries)
        ]
    )
    if fingerprint != level["geometry_fingerprint_sha256"]:
        raise RuntimeError("fragility thinned geometry fingerprint drift")

    train = tuple(map(str, plan["species"]["train_species"]))
    evaluation = tuple(map(str, plan["species"]["eval_species"]))
    if not train or not evaluation or set(train) & set(evaluation):
        raise RuntimeError("fragility inherited split is empty or overlapping")
    if set(train) | set(evaluation) != set(table.species):
        raise RuntimeError("fragility inherited split does not cover thinned geometry")

    return FragilityExecutionContext(
        geometries=table.geometries,
        train_species=train,
        eval_species=evaluation,
        retention_fraction=retention,
        retention_label=str(level["retention_label"]),
        master_seed=int(level["master_seed"]),
        level=level,
        plan=plan,
        phase3_rule=phase3_rule,
    )


__all__ = [
    "FragilityExecutionContext",
    "build_fragility_execution_plan",
    "canonical_retention_label",
    "derive_fragility_level_master_seed",
    "load_fragility_execution_context",
]
