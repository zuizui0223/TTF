from __future__ import annotations

from collections import defaultdict
import hashlib
from typing import Iterable, Mapping, Sequence


HOST_RESOURCE_PANEL_TAG = "lepidoptera-host-resource-v0.2"
SOURCE_ARCHIVE_SHA256 = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"


def panel_rank_key(
    species: str,
    *,
    source_sha256: str = SOURCE_ARCHIVE_SHA256,
    tag: str = HOST_RESOURCE_PANEL_TAG,
) -> tuple[str, str]:
    name = str(species).strip()
    if not name:
        raise ValueError("species must be non-empty")
    digest = hashlib.sha256(f"{tag}|{source_sha256}|{name}".encode("utf-8")).hexdigest()
    return digest, name


def select_and_split_species(
    species: Iterable[str],
    *,
    maximum_species: int = 500,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    labels = sorted({str(name).strip() for name in species if str(name).strip()})
    if maximum_species < 1:
        raise ValueError("maximum_species must be positive")
    ranked = sorted(labels, key=panel_rank_key)
    selected = tuple(ranked[: min(int(maximum_species), len(ranked))])
    cut = len(selected) // 2
    return selected, selected[:cut], selected[cut:]


def build_insect_host_footprints(
    insect_host_ids: Iterable[tuple[str, str]],
    host_native_units: Mapping[str, Iterable[str]],
) -> tuple[dict[str, frozenset[str]], dict[str, dict[str, int]]]:
    hosts_by_insect: dict[str, set[str]] = defaultdict(set)
    for insect, host_id in insect_host_ids:
        insect_name = str(insect).strip()
        host = str(host_id).strip()
        if insect_name and host:
            hosts_by_insect[insect_name].add(host)

    units_by_host = {
        str(host).strip(): frozenset(
            str(unit).strip() for unit in units if str(unit).strip()
        )
        for host, units in host_native_units.items()
        if str(host).strip()
    }

    footprints: dict[str, frozenset[str]] = {}
    diagnostics: dict[str, dict[str, int]] = {}
    for insect in sorted(hosts_by_insect):
        host_ids = hosts_by_insect[insect]
        with_native = {host for host in host_ids if units_by_host.get(host)}
        units: set[str] = set()
        for host in with_native:
            units.update(units_by_host[host])
        diagnostics[insect] = {
            "resolved_host_species": len(host_ids),
            "hosts_with_primary_native_units": len(with_native),
            "primary_native_wgsrpd3_units": len(units),
        }
        if units:
            footprints[insect] = frozenset(units)
    return footprints, diagnostics


def jaccard_units(a: Iterable[str], b: Iterable[str]) -> float:
    left = frozenset(map(str, a))
    right = frozenset(map(str, b))
    union = left | right
    if not union:
        raise ValueError("Jaccard is undefined for two empty footprints")
    return float(len(left & right) / len(union))


def species_list_sha256(species: Sequence[str]) -> str:
    payload = "\n".join(map(str, species)) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


__all__ = [
    "HOST_RESOURCE_PANEL_TAG",
    "SOURCE_ARCHIVE_SHA256",
    "build_insect_host_footprints",
    "jaccard_units",
    "panel_rank_key",
    "select_and_split_species",
    "species_list_sha256",
]
