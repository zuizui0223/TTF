from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from typing import Mapping, Sequence

import numpy as np


PILOT_TAG = "butterfly-resource-envelope-pilot-v0.1"


@dataclass(frozen=True)
class ResourceEnvelopeDescriptor:
    species: str
    host_family_count: float
    host_wgsrpd3_unit_count: int
    wing_size_proxy: float
    voltinism: str
    canopy_affinity: str
    edge_affinity: str
    moisture_affinity: str
    disturbance_affinity: str
    resolved_host_species: int
    hosts_with_primary_native_units: int


def wing_size_proxy(row: Mapping[str, object]) -> float:
    keys = (
        "WS_L",
        "WS_U",
        "FW_L",
        "FW_U",
        "WS_L_Fem",
        "WS_U_Fem",
        "WS_L_Mal",
        "WS_U_Mal",
        "FW_L_Fem",
        "FW_U_Fem",
        "FW_L_Mal",
        "FW_U_Mal",
    )
    values: list[float] = []
    for key in keys:
        raw = str(row.get(key, "")).strip()
        if raw in {"", "NA"}:
            continue
        values.append(float(raw))
    if not values:
        raise ValueError("missing wing size")
    return float(np.mean(values))


def descriptor_from_sources(
    species: str,
    trait_row: Mapping[str, object],
    host_diagnostics: Mapping[str, int],
) -> ResourceEnvelopeDescriptor:
    host_family_raw = str(trait_row.get("NumberOfHostplantFamilies", "")).strip()
    if host_family_raw in {"", "NA"}:
        raise ValueError(f"{species}: missing NumberOfHostplantFamilies")
    host_family_count = float(host_family_raw)
    if not math.isfinite(host_family_count) or host_family_count < 0:
        raise ValueError(f"{species}: invalid NumberOfHostplantFamilies")

    voltinism = str(trait_row.get("Voltinism", "")).strip()
    habitats = {
        "canopy_affinity": str(trait_row.get("CanopyAffinity", "")).strip(),
        "edge_affinity": str(trait_row.get("EdgeAffinity", "")).strip(),
        "moisture_affinity": str(trait_row.get("MoistureAffinity", "")).strip(),
        "disturbance_affinity": str(trait_row.get("DisturbanceAffinity", "")).strip(),
    }
    if voltinism in {"", "NA"} or any(v in {"", "NA"} for v in habitats.values()):
        raise ValueError(f"{species}: incomplete major LepTraits axes")

    unit_count = int(host_diagnostics.get("primary_native_wgsrpd3_units", 0))
    return ResourceEnvelopeDescriptor(
        species=str(species),
        host_family_count=host_family_count,
        host_wgsrpd3_unit_count=unit_count,
        wing_size_proxy=wing_size_proxy(trait_row),
        voltinism=voltinism,
        canopy_affinity=habitats["canopy_affinity"],
        edge_affinity=habitats["edge_affinity"],
        moisture_affinity=habitats["moisture_affinity"],
        disturbance_affinity=habitats["disturbance_affinity"],
        resolved_host_species=int(host_diagnostics.get("resolved_host_species", 0)),
        hosts_with_primary_native_units=int(
            host_diagnostics.get("hosts_with_primary_native_units", 0)
        ),
    )


def _zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = float(np.std(values))
    if sd <= np.sqrt(np.finfo(float).eps):
        return np.zeros_like(values)
    return (values - float(np.mean(values))) / sd


def _tie_key(species: str, tag: str = PILOT_TAG) -> str:
    return hashlib.sha256(f"{tag}|{species}".encode("utf-8")).hexdigest()


def select_resource_space_pilot(
    descriptors: Sequence[ResourceEnvelopeDescriptor],
    *,
    pilot_size: int = 10,
    tag: str = PILOT_TAG,
) -> tuple[str, ...]:
    """Select a response-blind pilot spanning taxonomic and geographic host breadth.

    The selection uses only two pre-response axes:
      - log1p LepTraits host-family count
      - log1p native WGSRPD3 host-resource footprint size

    Both axes are population-standardized. The first species is the point
    farthest from the bivariate centroid. Remaining species are chosen by
    deterministic farthest-point sampling, maximizing minimum squared distance
    to the already selected set. Hash order resolves exact ties.
    """
    eligible = [
        d
        for d in descriptors
        if d.host_family_count > 0
        and d.host_wgsrpd3_unit_count > 0
        and math.isfinite(d.host_family_count)
    ]
    if pilot_size < 1:
        raise ValueError("pilot_size must be positive")
    if len(eligible) < pilot_size:
        raise ValueError(
            f"need at least {pilot_size} resource-eligible species; got {len(eligible)}"
        )

    eligible = sorted(eligible, key=lambda d: (_tie_key(d.species, tag), d.species))
    x = _zscore(np.log1p(np.asarray([d.host_family_count for d in eligible], float)))
    y = _zscore(
        np.log1p(np.asarray([d.host_wgsrpd3_unit_count for d in eligible], float))
    )
    points = np.column_stack([x, y])

    radius2 = np.sum(points * points, axis=1)
    first = min(
        range(len(eligible)),
        key=lambda i: (-float(radius2[i]), _tie_key(eligible[i].species, tag)),
    )
    selected = [first]
    remaining = set(range(len(eligible))) - {first}

    while len(selected) < pilot_size:
        best = min(
            remaining,
            key=lambda i: (
                -min(
                    float(np.sum((points[i] - points[j]) ** 2))
                    for j in selected
                ),
                _tie_key(eligible[i].species, tag),
            ),
        )
        selected.append(best)
        remaining.remove(best)

    return tuple(eligible[i].species for i in selected)


def species_digest(species: Sequence[str]) -> str:
    payload = "\n".join(map(str, species)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


__all__ = [
    "PILOT_TAG",
    "ResourceEnvelopeDescriptor",
    "descriptor_from_sources",
    "select_resource_space_pilot",
    "species_digest",
    "wing_size_proxy",
]
