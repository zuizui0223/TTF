#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

from ttf.butterfly_resource_envelope import (
    descriptor_from_sources,
    select_resource_space_pilot,
    species_digest,
)
from ttf.lepidoptera_host_resource import build_insect_host_footprints


EXPECTED_S1_SCHEMA = "ttf_relational_prior_S1_species_exclusion_v0.1"
EXPECTED_S1_SPECIES_COUNT = 339
EXPECTED_S1_SPECIES_SHA256 = (
    "56d3e785134656b0393d448df12f9d3745bc5944812395e2c7db7c0e50688f4b"
)
EXPECTED_S1_DESIGN_SHA256 = (
    "f2f2088c294ecfac3c49659439fd7858330452bd8fbea090513a6b3d590c98f8"
)
EXPECTED_LEPTRAITS_SHA256 = (
    "6ec35b8a31e96c971aeaa228a48aae9f107c40c33695f0d470aa4382ca6d635b"
)
EXPECTED_INSECT_HOST_SHA256 = (
    "0a084fb5273e4780b03e015205c87e4545c69f0f8e517feb1a2a3879889508e9"
)
EXPECTED_HOST_DISTRIBUTION_SHA256 = (
    "c731906315f7452f83ce302c1d36c75ad94afc24ce7e24d461360312244ac558"
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def require_sha(path: Path, expected: str, label: str) -> None:
    observed = sha256_path(path)
    if observed != expected:
        raise RuntimeError(
            f"{label} SHA-256 drift: expected {expected}, observed {observed}"
        )


def load_s1_species(path: Path) -> tuple[str, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != EXPECTED_S1_SCHEMA:
        raise RuntimeError("unexpected S1 compact-manifest schema")
    if payload.get("source_design_sha256") != EXPECTED_S1_DESIGN_SHA256:
        raise RuntimeError("S1 source design SHA drift")
    names = tuple(map(str, payload.get("species", [])))
    if len(names) != EXPECTED_S1_SPECIES_COUNT or len(set(names)) != len(names):
        raise RuntimeError("S1 compact-manifest species count/uniqueness drift")
    if species_digest(names) != EXPECTED_S1_SPECIES_SHA256:
        raise RuntimeError("S1 species-list SHA drift")
    return names


def load_leptraits(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if "Species" not in set(reader.fieldnames or ()):
            raise RuntimeError("LepTraits schema drift: Species missing")
        for row in reader:
            species = str(row.get("Species", "")).strip()
            if species and species not in out:
                out[species] = {str(k): str(v or "") for k, v in row.items()}
    return out


def load_host_sidecars(
    interaction_path: Path,
    distribution_path: Path,
):
    pairs: list[tuple[str, str]] = []
    with interaction_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"insect_species", "accepted_plant_name_id"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("insect-host sidecar schema drift")
        for row in reader:
            insect = str(row.get("insect_species", "")).strip()
            host_id = str(row.get("accepted_plant_name_id", "")).strip()
            if insect and host_id:
                pairs.append((insect, host_id))

    native_units: dict[str, set[str]] = defaultdict(set)
    with distribution_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"accepted_plant_name_id", "area_code_l3"}
        if not required <= set(reader.fieldnames or ()):
            raise RuntimeError("host-distribution sidecar schema drift")
        for row in reader:
            host_id = str(row.get("accepted_plant_name_id", "")).strip()
            unit = str(row.get("area_code_l3", "")).strip()
            if host_id and unit:
                native_units[host_id].add(unit)
    return build_insect_host_footprints(pairs, native_units)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--s1-manifest", type=Path, required=True)
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--insect-host-csv", type=Path, required=True)
    ap.add_argument("--host-distribution-csv", type=Path, required=True)
    ap.add_argument("--output-table", type=Path, required=True)
    ap.add_argument("--output-footprints", type=Path, required=True)
    ap.add_argument("--output-pilot", type=Path, required=True)
    ap.add_argument("--pilot-size", type=int, default=10)
    args = ap.parse_args()

    require_sha(args.leptraits_csv, EXPECTED_LEPTRAITS_SHA256, "LepTraits")
    require_sha(
        args.insect_host_csv,
        EXPECTED_INSECT_HOST_SHA256,
        "HOSTS-WCVP insect-host sidecar",
    )
    require_sha(
        args.host_distribution_csv,
        EXPECTED_HOST_DISTRIBUTION_SHA256,
        "WCVP native WGSRPD3 sidecar",
    )

    s1_species = load_s1_species(args.s1_manifest)
    traits = load_leptraits(args.leptraits_csv)
    footprints, diagnostics = load_host_sidecars(
        args.insect_host_csv,
        args.host_distribution_csv,
    )

    missing_traits = [name for name in s1_species if name not in traits]
    if missing_traits:
        raise RuntimeError(
            "S1 species missing from exact LepTraits snapshot: "
            + ", ".join(missing_traits[:10])
        )

    descriptors = []
    for species in s1_species:
        descriptors.append(
            descriptor_from_sources(
                species,
                traits[species],
                diagnostics.get(
                    species,
                    {
                        "resolved_host_species": 0,
                        "hosts_with_primary_native_units": 0,
                        "primary_native_wgsrpd3_units": 0,
                    },
                ),
            )
        )

    pilot_species = select_resource_space_pilot(
        descriptors,
        pilot_size=args.pilot_size,
    )
    by_species = {d.species: d for d in descriptors}

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(asdict(descriptors[0]).keys())
    with args.output_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for descriptor in descriptors:
            writer.writerow(asdict(descriptor))

    footprint_payload = {
        "schema": "ttf_butterfly_resource_envelope_s1_footprints_v0.1",
        "status": "RECONSTRUCTED_FROM_PINNED_EXTERNAL_INPUTS",
        "s1_species_count": len(s1_species),
        "species_with_native_host_footprint": sum(
            1 for name in s1_species if name in footprints
        ),
        "species": {
            name: sorted(footprints[name])
            for name in s1_species
            if name in footprints
        },
        "input_sha256": {
            "s1_species_list": EXPECTED_S1_SPECIES_SHA256,
            "leptraits": EXPECTED_LEPTRAITS_SHA256,
            "insect_host_accepted": EXPECTED_INSECT_HOST_SHA256,
            "native_extant_nondoubtful_wgsrpd3": EXPECTED_HOST_DISTRIBUTION_SHA256,
        },
        "genetic_response_used": False,
        "gbif_butterfly_occurrences_used": False,
    }
    args.output_footprints.parent.mkdir(parents=True, exist_ok=True)
    args.output_footprints.write_text(
        json.dumps(footprint_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    pilot_payload = {
        "schema": "ttf_butterfly_resource_envelope_pilot_v0.1",
        "status": "EXPLORATORY_RESPONSE_BLIND_PILOT_SELECTED",
        "scientific_question": (
            "Within the geographic opportunity supplied by larval host plants, "
            "why do butterfly species fill some resource regions but not others?"
        ),
        "source_panel": {
            "name": "S1 major-complete butterfly panel",
            "species": len(s1_species),
            "species_list_sha256": EXPECTED_S1_SPECIES_SHA256,
            "source_design_sha256": EXPECTED_S1_DESIGN_SHA256,
        },
        "resource_eligible_species": sum(
            1 for d in descriptors if d.host_wgsrpd3_unit_count > 0
        ),
        "pilot_size": len(pilot_species),
        "pilot_species": list(pilot_species),
        "pilot_species_sha256": species_digest(pilot_species),
        "selection_rule": {
            "response_blind": True,
            "axes": [
                "z(log1p LepTraits NumberOfHostplantFamilies)",
                "z(log1p native host-resource WGSRPD3 unit count)",
            ],
            "algorithm": (
                "deterministic farthest-point sampling in the two-dimensional "
                "resource-specialization space; hash breaks exact ties"
            ),
            "uses_butterfly_gbif": False,
            "uses_genetic_response": False,
        },
        "pilot_descriptors": [asdict(by_species[name]) for name in pilot_species],
        "next_step": (
            "Acquire butterfly occurrence records only for these frozen 10 pilot "
            "species, map occurrences to WGSRPD3 level-3 units, and quantify "
            "whether host-available units contain enough occupied/unoccupied "
            "variation to justify a larger resource-envelope analysis."
        ),
        "input_sha256": footprint_payload["input_sha256"],
        "claim_boundary": (
            "Exploratory geometry/resolution pilot only; no confirmatory ecological "
            "effect or genetic-transfer conclusion is licensed."
        ),
    }
    args.output_pilot.parent.mkdir(parents=True, exist_ok=True)
    args.output_pilot.write_text(
        json.dumps(pilot_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "s1_species": len(s1_species),
                "resource_eligible_species": pilot_payload[
                    "resource_eligible_species"
                ],
                "pilot_species": list(pilot_species),
                "pilot_species_sha256": pilot_payload["pilot_species_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
