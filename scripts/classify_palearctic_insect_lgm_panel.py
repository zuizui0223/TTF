#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

ROLE_NAMESPACE = "palearctic-insect-lgm-role-v0.1"
ARCHIVE_SHA = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"


def role_order(species: list[str]) -> list[str]:
    return sorted(
        species,
        key=lambda name: (
            hashlib.sha256(f"{ROLE_NAMESPACE}|{ARCHIVE_SHA}|{name}".encode()).hexdigest(),
            name,
        ),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--occurrences", type=Path, required=True)
    ap.add_argument("--gbif-ledger", type=Path, required=True)
    ap.add_argument("--realm-geojson", type=Path, required=True)
    ap.add_argument("--output-panel", type=Path, required=True)
    args = ap.parse_args()

    from shapely.geometry import Point, shape
    from shapely.ops import unary_union

    contract = json.loads(args.contract.read_text())
    candidates = json.loads(args.candidates.read_text())
    gbif = json.loads(args.gbif_ledger.read_text())
    realm = json.loads(args.realm_geojson.read_text())

    if contract.get("schema") != "ttf_palearctic_insect_lgm_refugia_v0.1":
        raise RuntimeError("unexpected contract")
    if candidates.get("schema") != "ttf_palearctic_insect_lgm_taxonomic_candidates_v0.1":
        raise RuntimeError("unexpected candidates")
    if gbif.get("schema") != "ttf_palearctic_insect_lgm_gbif_aggregate_v0.1":
        raise RuntimeError("unexpected GBIF ledger")

    realm_field = contract["response_blind_panel_eligibility"]["palearctic_realm"]["realm_field"]
    realm_value = contract["response_blind_panel_eligibility"]["palearctic_realm"]["required_value"]
    features = [
        feature
        for feature in realm.get("features", [])
        if str(feature.get("properties", {}).get(realm_field)) == str(realm_value)
    ]
    if not features:
        raise RuntimeError("Palearctic realm GeoJSON contains no matching features")
    polygon = unary_union([shape(feature["geometry"]) for feature in features])

    by_species: dict[str, list[tuple[float, float]]] = {}
    with args.occurrences.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            by_species.setdefault(str(row["species"]), []).append(
                (float(row["longitude"]), float(row["latitude"]))
            )

    ledger = {str(row["species"]): row for row in gbif["species"]}
    threshold = float(
        contract["response_blind_panel_eligibility"]["palearctic_realm"][
            "minimum_fraction_retained_occurrences_in_realm"
        ]
    )
    rows = []
    eligible = []
    for species in candidates["species"]:
        info = ledger[species]
        coords = by_species.get(species, [])
        hits = sum(bool(polygon.covers(Point(lon, lat))) for lon, lat in coords)
        fraction = hits / len(coords) if coords else 0.0
        passes = (
            info.get("status") == "PASS_OCCURRENCE_GEOMETRY"
            and len(coords) >= 30
            and fraction >= threshold
        )
        row = {
            "species": species,
            "gbif_status": info.get("status"),
            "retained_occurrences": len(coords),
            "palearctic_occurrences": hits,
            "palearctic_fraction": fraction,
            "eligible": passes,
            "gbif_order": (info.get("taxonomy") or {}).get("order"),
            "gbif_family": (info.get("taxonomy") or {}).get("family"),
        }
        rows.append(row)
        if passes:
            eligible.append(species)

    ordered = role_order(eligible)
    cut = len(ordered) // 2
    sources = ordered[:cut]
    targets = ordered[cut:]
    minimum = int(contract["response_blind_panel_eligibility"]["minimum_species_to_construct_panel"])
    status = (
        "PASS_RESPONSE_BLIND_PALEARCTIC_INSECT_PANEL"
        if len(eligible) >= minimum and sources and targets
        else "NOT_EVALUABLE_PALEARCTIC_INSECT_PANEL_TOO_SMALL"
    )

    payload = {
        "schema": "ttf_palearctic_insect_lgm_panel_v0.1",
        "status": status,
        "candidate_species": len(candidates["species"]),
        "eligible_species": len(eligible),
        "source_species": sources,
        "target_species": targets,
        "source_count": len(sources),
        "target_count": len(targets),
        "directed_dyads": len(sources) * len(targets),
        "species_ledger": rows,
        "realm_feature_count": len(features),
        "realm_rule": {
            "field": realm_field,
            "value": realm_value,
            "minimum_fraction": threshold,
        },
        "response_firewall": candidates["response_firewall"],
        "next_step": (
            "fit fixed current/LGM SDMs and build response-blind R_LGM relation"
            if status == "PASS_RESPONSE_BLIND_PALEARCTIC_INSECT_PANEL"
            else "STOP_WITHOUT_SUBPANEL_GENETIC_RESPONSE"
        ),
    }
    args.output_panel.parent.mkdir(parents=True, exist_ok=True)
    args.output_panel.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": status,
        "eligible_species": len(eligible),
        "sources": len(sources),
        "targets": len(targets),
        "dyads": len(sources) * len(targets),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
