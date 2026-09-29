#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def _sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _point_on_segment(px, py, ax, ay, bx, by, eps=1e-10):
    cross = (px-ax)*(by-ay) - (py-ay)*(bx-ax)
    if abs(cross) > eps:
        return False
    return (
        min(ax,bx)-eps <= px <= max(ax,bx)+eps
        and min(ay,by)-eps <= py <= max(ay,by)+eps
    )


def _ring_contains(lon: float, lat: float, ring) -> bool:
    inside = False
    n = len(ring)
    if n < 4:
        return False
    for i in range(n - 1):
        ax, ay = map(float, ring[i])
        bx, by = map(float, ring[i+1])
        if _point_on_segment(lon, lat, ax, ay, bx, by):
            return True
        crosses = ((ay > lat) != (by > lat))
        if crosses:
            x = ax + (lat-ay) * (bx-ax) / (by-ay)
            if abs(x-lon) <= 1e-12:
                return True
            if x > lon:
                inside = not inside
    return inside


def _polygon_contains(lon: float, lat: float, polygon) -> bool:
    if not polygon:
        return False
    if not _ring_contains(lon, lat, polygon[0]):
        return False
    return not any(_ring_contains(lon, lat, hole) for hole in polygon[1:])


def geometry_contains(lon: float, lat: float, geometry: dict) -> bool:
    kind = geometry.get("type")
    coords = geometry.get("coordinates", [])
    if kind == "Polygon":
        return _polygon_contains(lon, lat, coords)
    if kind == "MultiPolygon":
        return any(_polygon_contains(lon, lat, polygon) for polygon in coords)
    raise RuntimeError(f"unsupported realm geometry type: {kind!r}")


def _load_palearctic_geometries(payload: dict) -> list[dict]:
    if payload.get("type") != "FeatureCollection":
        raise RuntimeError("realm input must be GeoJSON FeatureCollection")
    out = []
    for feature in payload.get("features", []):
        props = feature.get("properties", {})
        code = props.get("realmcode", props.get("REALMCODE"))
        if str(code).upper() == "PA":
            out.append(feature["geometry"])
    if not out:
        raise RuntimeError("GeoJSON contains no Palearctic (PA) feature")
    return out


def _split(species: list[str], contract_sha: str) -> tuple[list[str], list[str]]:
    keyed = sorted(
        (
            hashlib.sha256(
                f"{contract_sha}|split|{name}".encode("utf-8")
            ).hexdigest(),
            name,
        )
        for name in species
    )
    n_train = len(keyed) // 2
    return (
        sorted(name for _, name in keyed[:n_train]),
        sorted(name for _, name in keyed[n_train:]),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--localities", type=Path, required=True)
    ap.add_argument("--realm-geojson", type=Path, required=True)
    ap.add_argument("--contract", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    localities = json.loads(args.localities.read_text())
    contract = json.loads(args.contract.read_text())
    realm = json.loads(args.realm_geojson.read_text())
    if contract.get("schema") != "ttf_genetic_palearctic_lepidoptera_eligibility_v0.1":
        raise RuntimeError("unexpected eligibility contract")
    if contract.get("status") != "FROZEN_BEFORE_FORMAL_REALM_CLASSIFICATION":
        raise RuntimeError("eligibility contract was not frozen before formal classification")
    firewall = localities.get("outcome_firewall", {})
    score_keys = (
        "species_level_phase4_scores_used",
        "species_level_genetic_scores_used",
    )
    present_score_keys = [key for key in score_keys if key in firewall]
    if len(present_score_keys) != 1:
        raise RuntimeError(
            "response-blind locality input must expose exactly one recognized "
            "species-level genetic-score firewall key"
        )
    required_false = (
        present_score_keys
        + ["pairwise_genetic_distances_used", "nucleotide_identity_used"]
    )
    bad = {key: firewall.get(key) for key in required_false if firewall.get(key) is not False}
    if bad:
        raise RuntimeError(f"response-blind locality input firewall is open: {bad}")

    geoms = _load_palearctic_geometries(realm)
    minimum_fraction = float(contract["eligibility"]["minimum_palearctic_fraction"])
    rows = {}
    eligible = []
    for name, row in sorted(localities["species"].items()):
        tax = row["taxonomy"]
        if tax.get("class") != contract["eligibility"]["class_exact"]:
            continue
        if tax.get("order") != contract["eligibility"]["order_exact"]:
            continue
        pts = row["localities"]
        inside = sum(
            any(geometry_contains(float(lon), float(lat), geom) for geom in geoms)
            for lat, lon in pts
        )
        fraction = inside / len(pts)
        passed = fraction >= minimum_fraction
        rows[name] = {
            "family": tax.get("family"),
            "localities": len(pts),
            "palearctic_localities": inside,
            "palearctic_fraction": fraction,
            "eligible": passed,
        }
        if passed:
            eligible.append(name)

    contract_sha = _sha256_path(args.contract)
    train, evaluation = _split(eligible, contract_sha)
    minimum_species = int(contract["eligibility"]["minimum_eligible_species"])
    status = (
        "PASS_FORMAL_PALEARCTIC_LEPIDOPTERA_ELIGIBILITY"
        if len(eligible) >= minimum_species
        else "STOP_SUBGROUP_NOT_EVALUABLE_BEFORE_LGM_MODELING"
    )
    payload = {
        "schema": "ttf_genetic_palearctic_lepidoptera_eligibility_result_v0.1",
        "status": status,
        "eligibility_contract_sha256": contract_sha,
        "localities_sha256": _sha256_path(args.localities),
        "realm_geojson_sha256": _sha256_path(args.realm_geojson),
        "realm_source": contract["eligibility"]["realm_dataset"],
        "minimum_palearctic_fraction": minimum_fraction,
        "screened_lepidoptera_species": len(rows),
        "eligible_species_count": len(eligible),
        "eligible_species": sorted(eligible),
        "new_response_blind_split": {
            "train_species": train,
            "eval_species": evaluation,
            "train_count": len(train),
            "eval_count": len(evaluation),
        },
        "per_species": rows,
        "outcome_firewall": {
            "species_level_phase4_scores_used": False,
            "pairwise_genetic_distances_used": False,
            "nucleotide_identity_used": False,
        },
        "next_step": (
            "construct_frozen_LGM_refugial_predictor_then_TTF_Q_qualification"
            if status.startswith("PASS_")
            else "STOP_DO_NOT_OPEN_SUBGROUP_GENETIC_RESPONSE"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": status,
        "screened_lepidoptera": len(rows),
        "eligible": len(eligible),
        "train": len(train),
        "eval": len(evaluation),
        "genetic_response_used": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
