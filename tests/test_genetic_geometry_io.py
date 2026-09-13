from __future__ import annotations

import csv
from pathlib import Path

import pytest

from ttf.genetic_geometry_io import load_frozen_genetic_geometry_csv, sha256_path


FIELDS = ["species", "locality_index", "x_km", "y_km", "z_km", "graph_k"]


def _write_geometry(path: Path, *, broken_index: bool = False, graph_k: int = 2) -> None:
    rows = []
    for species_index, species in enumerate(("sp_a", "sp_b", "sp_c", "sp_d")):
        for locality in range(4):
            index = locality
            if broken_index and species == "sp_a" and locality == 3:
                index = 4
            rows.append(
                {
                    "species": species,
                    "locality_index": index,
                    "x_km": species_index * 100.0 + locality,
                    "y_km": 0.0,
                    "z_km": 0.0,
                    "graph_k": graph_k,
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def test_frozen_geometry_loader_round_trips_canonical_csv(tmp_path: Path) -> None:
    path = tmp_path / "geometry.csv"
    _write_geometry(path)
    digest = sha256_path(path)
    table = load_frozen_genetic_geometry_csv(path, expected_sha256=digest)
    assert table.csv_sha256 == digest
    assert table.species == ("sp_a", "sp_b", "sp_c", "sp_d")
    assert all(geometry.n_localities == 4 for geometry in table.geometries.values())
    assert all(geometry.graph_k == 2 for geometry in table.geometries.values())


def test_frozen_geometry_loader_rejects_sha_drift(tmp_path: Path) -> None:
    path = tmp_path / "geometry.csv"
    _write_geometry(path)
    with pytest.raises(RuntimeError, match="SHA256 drift"):
        load_frozen_genetic_geometry_csv(path, expected_sha256="0" * 64)


def test_frozen_geometry_loader_rejects_noncanonical_indices(tmp_path: Path) -> None:
    path = tmp_path / "geometry.csv"
    _write_geometry(path, broken_index=True)
    with pytest.raises(RuntimeError, match="noncanonical locality indices"):
        load_frozen_genetic_geometry_csv(path)


def test_frozen_geometry_loader_rejects_density_scaled_k_drift(tmp_path: Path) -> None:
    path = tmp_path / "geometry.csv"
    _write_geometry(path, graph_k=3)
    with pytest.raises(RuntimeError, match="density-scaled graph_k drift"):
        load_frozen_genetic_geometry_csv(path)


def test_frozen_geometry_loader_can_validate_explicit_k_without_density_rule(tmp_path: Path) -> None:
    path = tmp_path / "geometry.csv"
    _write_geometry(path, graph_k=3)
    table = load_frozen_genetic_geometry_csv(path, expected_neighbor_fraction=None)
    assert all(geometry.graph_k == 3 for geometry in table.geometries.values())
