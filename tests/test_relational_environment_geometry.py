import importlib.util
from pathlib import Path

import numpy as np

from ttf.genetic_geometry import prepare_density_scaled_genetic_geometry

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "reconstruct_relational_environment_geometry.py"
_SPEC = importlib.util.spec_from_file_location("reconstruct_relational_environment_geometry", _SCRIPT)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)
edge_rows = _MODULE.edge_rows
locality_rows = _MODULE.locality_rows
verify_candidate = _MODULE.verify_candidate


def toy_geometry():
    coordinates = np.asarray([
        [6371.0, 0.0, 0.0],
        [6369.0, 100.0, 0.0],
        [6368.0, 0.0, 120.0],
        [6367.0, -110.0, 0.0],
        [6366.0, 0.0, -130.0],
        [6365.0, 80.0, 80.0],
        [6364.0, -80.0, 80.0],
        [6363.0, 80.0, -80.0],
        [6362.0, -80.0, -80.0],
        [6361.0, 160.0, 0.0],
        [6360.0, -160.0, 0.0],
        [6359.0, 0.0, 160.0],
    ])
    return prepare_density_scaled_genetic_geometry(coordinates, neighbor_fraction=0.15)


def test_frozen_candidate_geometry_contract_and_midpoints():
    geometry = toy_geometry()
    row = {
        "species": "Example species",
        "n_localities": str(geometry.n_localities),
        "graph_k": str(geometry.graph_k),
        "edges": str(geometry.n_edges),
        "min_endpoint_disjoint_training_edges": str(
            geometry.min_endpoint_disjoint_training_edges
        ),
    }
    verify_candidate(row, geometry)

    rows = list(edge_rows("Example species", geometry))
    assert len(rows) == geometry.n_edges
    first = rows[0]
    left, right = geometry.edge_nodes[0]
    expected = 0.5 * (geometry.coordinates[left] + geometry.coordinates[right])
    observed = np.asarray([
        float(first["mid_x_km"]),
        float(first["mid_y_km"]),
        float(first["mid_z_km"]),
    ])
    assert np.allclose(observed, expected)


def test_locality_rows_preserve_geometry_order():
    geometry = toy_geometry()
    latlon = np.column_stack((
        np.arange(geometry.n_localities, dtype=float),
        -np.arange(geometry.n_localities, dtype=float),
    ))
    rows = list(locality_rows("Example species", geometry, latlon))
    assert len(rows) == geometry.n_localities
    assert [row["locality_index"] for row in rows] == list(range(geometry.n_localities))
    assert float(rows[3]["latitude"]) == 3.0
    assert float(rows[3]["longitude"]) == -3.0
    assert np.isclose(float(rows[3]["x_km"]), geometry.coordinates[3, 0])
