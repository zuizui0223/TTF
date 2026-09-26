#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from ttf.butterfly_climate_release import (
    average_ranks,
    freedman_lane_partial_spearman_permutation,
    one_sided_partial_spearman_permutation,
)


PROTOCOL_SCHEMAS = {
    "ttf_butterfly_climate_release_independent_test_v0.1": (
        "FROZEN_PILOT_DERIVED_HYPOTHESIS_BEFORE_INDEPENDENT_GBIF_OR_CLIMATE"
    ),
    "ttf_butterfly_climate_release_independent_test_v0.2": (
        "FROZEN_RESPONSE_BLIND_PRIMARY_INFERENCE_CORRECTION_BEFORE_INDEPENDENT_PRECLIMATE_OR_CLIMATE_RESULT"
    ),
    "ttf_butterfly_climate_release_independent_test_v0.2.1": (
        "FROZEN_RESPONSE_BLIND_ROW_ORDER_INVARIANT_PERMUTATION_FIX_BEFORE_INDEPENDENT_PRECLIMATE_OR_CLIMATE_RESULT"
    ),
}
PANEL_SCHEMA = "ttf_butterfly_climate_release_independent_panel_v0.1"
CLIMATE_SCHEMA = "ttf_butterfly_resource_envelope_climate_crossfit_v0.1"


def spearman(a, b) -> float:
    x = average_ranks(a)
    y = average_ranks(b)
    if np.std(x) <= np.sqrt(np.finfo(float).eps):
        raise ValueError("first Spearman variable has zero variance")
    if np.std(y) <= np.sqrt(np.finfo(float).eps):
        raise ValueError("second Spearman variable has zero variance")
    return float(np.corrcoef(x, y)[0, 1])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--protocol-json", type=Path, required=True)
    ap.add_argument("--panel-json", type=Path, required=True)
    ap.add_argument("--climate-summary-json", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    args = ap.parse_args()

    protocol = json.loads(args.protocol_json.read_text(encoding="utf-8"))
    protocol_schema = str(protocol.get("schema") or "")
    expected_status = PROTOCOL_SCHEMAS.get(protocol_schema)
    if expected_status is None:
        raise RuntimeError("unexpected independent-test protocol schema")
    if protocol.get("status") != expected_status:
        raise RuntimeError("independent-test protocol is not frozen")

    panel = json.loads(args.panel_json.read_text(encoding="utf-8"))
    if panel.get("schema") != PANEL_SCHEMA:
        raise RuntimeError("unexpected independent-panel schema")
    panel_species = tuple(map(str, panel.get("species", [])))
    excluded = set(
        map(str, protocol["independent_panel"]["excluded_species"])
    )
    if set(panel_species) & excluded:
        raise RuntimeError("pilot species leaked into independent test panel")

    panel_descriptors = {
        str(row["species"]): row
        for row in panel.get("descriptors", [])
    }
    if set(panel_descriptors) != set(panel_species):
        raise RuntimeError("panel descriptor coverage drift")

    climate = json.loads(
        args.climate_summary_json.read_text(encoding="utf-8")
    )
    if climate.get("schema") != CLIMATE_SCHEMA:
        raise RuntimeError("unexpected climate cross-fit schema")
    if climate.get("analysis_program_schema") != protocol_schema:
        raise RuntimeError("climate summary is not from the selected independent protocol")
    if climate.get("analysis_scope") != "independent_pilot_derived_test":
        raise RuntimeError("climate summary scope drift")

    climate_rows = {
        str(row["species"]): row
        for row in climate.get("species", [])
    }
    if not set(climate_rows).issubset(set(panel_species)):
        raise RuntimeError("climate result contains species outside independent panel")

    informative = []
    for name in panel_species:
        row = climate_rows.get(name)
        if not row or not bool(row.get("climate_crossfit_informative")):
            continue
        score = row.get("probability_never_observed_mismatch_greater")
        host_units = row.get("host_units")
        descriptor = panel_descriptors[name]
        if score is None or host_units is None:
            raise RuntimeError(f"informative species missing primary fields: {name}")
        informative.append(
            {
                "species": name,
                "climate_filtering_score": float(score),
                "host_family_count": float(descriptor["host_family_count"]),
                "contemporary_host_resource_units": int(host_units),
                "wing_size_proxy": float(descriptor["wing_size_proxy"]),
                "voltinism": str(descriptor["voltinism"]),
                "host_breadth_stratum": str(
                    descriptor["host_breadth_stratum"]
                ),
                "mismatch_delta_never_minus_eval_observed": (
                    None
                    if row.get("mismatch_delta_never_minus_eval_observed") is None
                    else float(row["mismatch_delta_never_minus_eval_observed"])
                ),
                "effort_supported_eval_observed_units": int(
                    row["effort_supported_eval_observed_units"]
                ),
                "effort_supported_never_observed_units": int(
                    row["effort_supported_never_observed_units"]
                ),
            }
        )

    minimum = int(
        protocol["climate_crossfit"]["minimum_climate_informative_species"]
    )
    if len(informative) < minimum:
        payload = {
            "schema": (
                (
                    "ttf_butterfly_climate_release_independent_result_v0.2.1"
                    if protocol_schema.endswith("_v0.2.1")
                    else "ttf_butterfly_climate_release_independent_result_v0.2"
                )
                if protocol_schema in {
                    "ttf_butterfly_climate_release_independent_test_v0.2",
                    "ttf_butterfly_climate_release_independent_test_v0.2.1",
                }
                else "ttf_butterfly_climate_release_independent_result_v0.1"
            ),
            "status": "NOT_EVALUABLE_INDEPENDENT_CLIMATE_RELEASE_TEST",
            "reason": "CLIMATE_INFORMATIVE_SPECIES_BELOW_FROZEN_MINIMUM",
            "independent_panel_species": len(panel_species),
            "climate_informative_species": len(informative),
            "minimum_required": minimum,
            "species": informative,
            "primary_test_run": False,
            "pilot_species_used": False,
            "genetic_response_used": False,
        }
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    response = [row["climate_filtering_score"] for row in informative]
    predictor = [row["host_family_count"] for row in informative]
    control = [
        np.log1p(row["contemporary_host_resource_units"])
        for row in informative
    ]
    primary = protocol["primary_test"]
    iterations = int(primary["permutation"]["iterations"])
    if protocol_schema in {
                    "ttf_butterfly_climate_release_independent_test_v0.2",
                    "ttf_butterfly_climate_release_independent_test_v0.2.1",
                }:
        test = freedman_lane_partial_spearman_permutation(
            response,
            predictor,
            control,
            [row["species"] for row in informative],
            iterations=iterations,
            tag=str(primary["permutation"]["tag"]),
        )
    else:
        test = one_sided_partial_spearman_permutation(
            response,
            predictor,
            control,
            iterations=iterations,
            tag=str(primary["permutation"]["tag"]),
        )
    alpha = float(primary["alpha"])
    supported = (
        float(test["observed_partial_spearman"]) < 0
        and float(test["one_sided_p_value"]) < alpha
    )

    by_voltinism = defaultdict(list)
    for row in informative:
        by_voltinism[row["voltinism"]].append(
            row["climate_filtering_score"]
        )

    payload = {
        "schema": (
                "ttf_butterfly_climate_release_independent_result_v0.2"
                if protocol_schema in {
                    "ttf_butterfly_climate_release_independent_test_v0.2",
                    "ttf_butterfly_climate_release_independent_test_v0.2.1",
                }
                else "ttf_butterfly_climate_release_independent_result_v0.1"
            ),
        "status": "INDEPENDENT_TEST_COMPLETE",
        "protocol_schema": protocol_schema,
        "hypothesis": protocol["ecological_hypothesis"]["statement"],
        "independent_panel_species": len(panel_species),
        "climate_informative_species": len(informative),
        "minimum_required": minimum,
        "primary_test": {
            "statistic": primary["statistic"],
            "alternative": primary["alternative"],
            "alpha": alpha,
            **test,
            "decision": (
                "SUPPORT_PILOT_DERIVED_CLIMATE_RELEASE_HYPOTHESIS"
                if supported
                else "PILOT_DERIVED_CLIMATE_RELEASE_HYPOTHESIS_NOT_SUPPORTED"
            ),
        },
        "secondary_descriptive": {
            "spearman_filtering_vs_host_family_count": spearman(
                response, predictor
            ),
            "spearman_filtering_vs_log_contemporary_resource_units": spearman(
                response, control
            ),
            "spearman_filtering_vs_wing_size_proxy": spearman(
                response,
                [row["wing_size_proxy"] for row in informative],
            ),
            "median_filtering_score_by_voltinism": {
                name: float(np.median(values))
                for name, values in sorted(by_voltinism.items())
            },
        },
        "species": informative,
        "claim_boundary": {
            "pilot_species_used_in_primary_test": False,
            "panel_species_frozen_before_independent_gbif": True,
            "primary_direction_frozen_before_independent_gbif": True,
            "quality_thresholds_frozen_before_independent_gbif": True,
            "primary_inference_method": (
                (
                    "response_blind_row_order_invariant_residual_permutation_fix_after_transport_started_before_ecological_response"
                    if protocol_schema.endswith("_v0.2.1")
                    else "response_blind_residual_permutation_correction_after_transport_started_before_ecological_response"
                )
                if protocol_schema in {
                    "ttf_butterfly_climate_release_independent_test_v0.2",
                    "ttf_butterfly_climate_release_independent_test_v0.2.1",
                }
                else "raw_predictor_label_permutation_v0.1"
            ),
            "independent_ecological_response_seen_before_v0_2_correction": False,
            "never_observed_resource_unit_is_not_true_absence": True,
            "sampling_effort_proxy_is_other_independent_panel_butterflies": True,
            "genetic_response_used": False,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
