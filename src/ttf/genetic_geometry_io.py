from __future__ import annotations

from dataclasses import dataclass
import csv
import hashlib
from pathlib import Path

import numpy as np

from .genetic_geometry import GeneticSamplingGeometry, prepare_genetic_sampling_geometry
from .heterogeneous_inference import density_scaled_k


_REQUIRED_COLUMNS = (
    "species",
    "locality_index",
    "x_km",
    "y_km",
    "z_km",
    "graph_k",
)


@dataclass(frozen=True)
class FrozenGeneticGeometryTable:
    """Audit-ready genetic geometry reconstructed from a frozen locality CSV."""

    geometries: dict[str, GeneticSamplingGeometry]
    csv_sha256: str
    species: tuple[str, ...]


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_frozen_genetic_geometry_csv(
    path: Path,
    *,
    expected_sha256: str | None = None,
    expected_neighbor_fraction: float | None = 0.15,
) -> FrozenGeneticGeometryTable:
    """Reconstruct species-local genetic geometry from a frozen locality CSV.

    The CSV is treated as authoritative for locality identity/order and graph ``k``.
    Reconstructed geometry must preserve that exact canonical coordinate order.
    When ``expected_neighbor_fraction`` is not ``None``, the stored ``graph_k`` must
    also equal the prospectively declared density-scaled rule for the locality count.
    No response or sequence information is read here.
    """
    source = Path(path)
    observed_sha = sha256_path(source)
    if expected_sha256 is not None and observed_sha != str(expected_sha256):
        raise RuntimeError("frozen genetic geometry CSV SHA256 drift")

    rows_by_species: dict[str, list[tuple[int, np.ndarray, int]]] = {}
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = tuple(reader.fieldnames or ())
        missing = [name for name in _REQUIRED_COLUMNS if name not in fieldnames]
        if missing:
            raise RuntimeError(f"frozen genetic geometry CSV missing columns: {missing}")
        for row in reader:
            species = str(row["species"]).strip()
            if not species:
                raise RuntimeError("empty species name in frozen genetic geometry CSV")
            try:
                locality_index = int(row["locality_index"])
                coordinates = np.asarray(
                    [float(row["x_km"]), float(row["y_km"]), float(row["z_km"])],
                    dtype=float,
                )
                graph_k = int(row["graph_k"])
            except (TypeError, ValueError) as exc:
                raise RuntimeError(f"invalid frozen geometry row for {species}") from exc
            if locality_index < 0 or graph_k < 1 or not np.isfinite(coordinates).all():
                raise RuntimeError(f"invalid frozen geometry values for {species}")
            rows_by_species.setdefault(species, []).append(
                (locality_index, coordinates, graph_k)
            )

    if not rows_by_species:
        raise RuntimeError("frozen genetic geometry CSV contains no species")

    geometries: dict[str, GeneticSamplingGeometry] = {}
    for species in sorted(rows_by_species):
        ordered = sorted(rows_by_species[species], key=lambda item: item[0])
        indices = [item[0] for item in ordered]
        if indices != list(range(len(ordered))):
            raise RuntimeError(f"noncanonical locality indices for {species}")
        coordinates = np.vstack([item[1] for item in ordered])
        if len(np.unique(coordinates, axis=0)) != len(coordinates):
            raise RuntimeError(f"duplicate frozen localities for {species}")
        k_values = {item[2] for item in ordered}
        if len(k_values) != 1:
            raise RuntimeError(f"graph_k varies within {species}")
        graph_k = int(next(iter(k_values)))
        if graph_k > len(coordinates) - 1:
            raise RuntimeError(f"graph_k exceeds available localities for {species}")
        if expected_neighbor_fraction is not None:
            expected_k = density_scaled_k(
                len(coordinates), fraction=float(expected_neighbor_fraction)
            )
            if graph_k != expected_k:
                raise RuntimeError(
                    f"density-scaled graph_k drift for {species}: {graph_k} != {expected_k}"
                )

        geometry = prepare_genetic_sampling_geometry(coordinates, k=graph_k)
        if geometry.graph_k != graph_k:
            raise RuntimeError(f"reconstructed graph_k drift for {species}")
        if not np.array_equal(np.asarray(geometry.coordinates), coordinates):
            raise RuntimeError(f"canonical coordinate order drift for {species}")
        geometries[species] = geometry

    species = tuple(sorted(geometries))
    return FrozenGeneticGeometryTable(
        geometries=geometries,
        csv_sha256=observed_sha,
        species=species,
    )


__all__ = [
    "FrozenGeneticGeometryTable",
    "load_frozen_genetic_geometry_csv",
    "sha256_path",
]
