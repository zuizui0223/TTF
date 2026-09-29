#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from ttf.palearctic_insect_screen import (
    apply_gbif_admissibility,
    deterministic_role_split,
    load_geojson,
    select_realm_geometries,
    summarize_species_realm,
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_gbif(path: Path) -> dict[str, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"species", "retained_occurrences"}
        if not required.issubset(reader.fieldnames or ()):
            raise RuntimeError("GBIF ledger must contain species,retained_occurrences")
        return {str(row["species"]): int(row["retained_occurrences"]) for row in reader}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--localities", type=Path, required=True)
    ap.add_argument("--realm-geojson", type=Path, required=True)
    ap.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/palearctic_insect_lgm_eligibility_v0.2.json"),
    )
    ap.add_argument("--gbif-ledger", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    contract = json.loads(args.contract.read_text())
    if contract.get("schema") != "ttf_palearctic_insect_lgm_eligibility_v0.2":
        raise RuntimeError("unexpected eligibility contract")
    if contract.get("status") != "FROZEN_BEFORE_AUTHORITATIVE_WWF_REALM_ASSIGNMENT_AND_BEFORE_GBIF_CENSUS":
        raise RuntimeError("eligibility contract was not frozen before realm/GBIF results")

    rule = contract["scientific_filter_unchanged"]
    geojson = load_geojson(args.realm_geojson)
    geoms = select_realm_geometries(
        geojson,
        realm=str(rule["required_realm"]),
        realmcode="PA",
    )
    summaries = summarize_species_realm(
        read_rows(args.localities),
        geometries=geoms,
        allowed_orders=rule["terrestrial_allowed_orders"],
        excluded_coleoptera_families=rule["excluded_coleoptera_families"],
        minimum_fraction=float(rule["minimum_fraction_of_genetic_panel_localities_in_realm"]),
        minimum_localities=int(rule["minimum_number_of_genetic_panel_localities_in_realm"]),
    )

    eligible = tuple(sorted(row.species for row in summaries if row.biologically_eligible))
    payload = {
        "schema": "ttf_palearctic_insect_lgm_screen_result_v0.1",
        "status": "REALM_SCREEN_COMPLETE_GBIF_NOT_YET_APPLIED",
        "realm_source_contract": contract["stage_A_biological_eligibility"]["realm_source"],
        "realm_features": len(geoms),
        "species_screened": len(summaries),
        "biologically_eligible_species": len(eligible),
        "biologically_eligible_names": list(eligible),
        "species": [
            {
                "species": row.species,
                "class": row.class_name,
                "order": row.order,
                "family": row.family,
                "total_localities": row.total_localities,
                "palearctic_localities": row.palearctic_localities,
                "palearctic_fraction": row.palearctic_fraction,
                "terrestrial_core": row.terrestrial_core,
                "realm_eligible": row.realm_eligible,
                "biologically_eligible": row.biologically_eligible,
            }
            for row in summaries
        ],
        "response_firewall": contract["response_firewall"],
    }

    if args.gbif_ledger is not None:
        retained = read_gbif(args.gbif_ledger)
        admissible = apply_gbif_admissibility(
            summaries,
            retained,
            minimum_occurrences=int(contract["stage_B_sdm_occurrence_admissibility"]["minimum_after_thinning"]),
        )
        namespace = str(contract["stage_C_role_split"]["namespace"])
        source, target = deterministic_role_split(admissible, namespace=namespace)
        min_source = int(contract["stage_C_role_split"]["minimum_sources"])
        min_target = int(contract["stage_C_role_split"]["minimum_targets"])
        passed = len(source) >= min_source and len(target) >= min_target
        payload.update(
            {
                "status": "PASS_TO_SDM_FIT" if passed else "NOT_EVALUABLE_SUBPANEL_GEOMETRY",
                "gbif_admissible_species": len(admissible),
                "gbif_admissible_names": list(admissible),
                "source_species": list(source),
                "target_species": list(target),
                "source_count": len(source),
                "target_count": len(target),
                "role_split_pass": bool(passed),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        key: payload[key]
        for key in (
            "status",
            "species_screened",
            "biologically_eligible_species",
        )
    } | {
        key: payload[key]
        for key in ("gbif_admissible_species", "source_count", "target_count")
        if key in payload
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
