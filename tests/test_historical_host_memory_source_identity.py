"""Response-blind source-identity firewall for the frozen 642-species host-memory panel."""
import hashlib
import json
from pathlib import Path

import pytest

from scripts.build_historical_host_memory_predictor import (
    FROZEN_INPUT_SHA256,
    verify_frozen_input_hashes,
    verify_host_pair_identity,
)


def test_host_pair_mapping_accepts_only_the_frozen_candidate_associations():
    candidates = [
        {"species": "Epinotia tedella", "accepted_host_ids": "380406.0;380741.0"},
        {"species": "Triodia sylvina", "accepted_host_ids": "2855039.0"},
    ]
    sidecar = [
        {"insect_species": "Epinotia tedella", "accepted_plant_name_id": "380741"},
        {"insect_species": "Epinotia tedella", "accepted_plant_name_id": "380406.0"},
        {"insect_species": "Triodia sylvina", "accepted_plant_name_id": "2855039"},
        {"insect_species": "Other insect", "accepted_plant_name_id": "123"},
    ]
    verify_host_pair_identity(candidates, sidecar)
    with pytest.raises(RuntimeError, match="mapping drift for Epinotia tedella"):
        verify_host_pair_identity(candidates, sidecar[1:])
    with pytest.raises(RuntimeError, match="mapping drift for Triodia sylvina"):
        verify_host_pair_identity(candidates, sidecar + [
            {"insect_species": "Triodia sylvina", "accepted_plant_name_id": "444"},
        ])
    with pytest.raises(RuntimeError, match="mapping drift for Triodia sylvina"):
        verify_host_pair_identity(candidates, sidecar[:2] + sidecar[3:])


def test_host_pair_mapping_rejects_invalid_candidate_ids_and_blank_sidecar():
    with pytest.raises(RuntimeError, match="invalid frozen host IDs"):
        verify_host_pair_identity([{"species": "A b", "accepted_host_ids": "3;3.0"}], [])
    with pytest.raises(RuntimeError, match="empty insect or accepted-host ID"):
        verify_host_pair_identity([], [
            {"insect_species": "A b", "accepted_plant_name_id": ""},
        ])


def test_exact_input_sha_gate_detects_same_shape_file_mutation(tmp_path, monkeypatch):
    from scripts import build_historical_host_memory_predictor as builder

    original = {}
    paths = {}
    for name in ("candidates", "localities", "edges"):
        p = tmp_path / (name + ".csv")
        p.write_text("species,edge_index\\nA b,1\\n", encoding="utf-8")
        paths[name] = p
        original[name] = hashlib.sha256(p.read_bytes()).hexdigest()
    monkeypatch.setattr(builder, "FROZEN_INPUT_SHA256", original)
    verify_frozen_input_hashes(paths["candidates"], paths["localities"], paths["edges"])
    paths["edges"].write_text("species,edge_index\\nA b,2\\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="edges SHA256 mismatch"):
        verify_frozen_input_hashes(paths["candidates"], paths["localities"], paths["edges"])


def test_pinned_hashes_match_preexisting_frozen_census_and_geometry_receipts():
    root = Path(__file__).resolve().parents[1]
    census = json.loads((root / "benchmarks/frozen/historical_host_memory_candidate_census_v0.1.json").read_text())
    geometry = json.loads((root / "benchmarks/frozen/historical_host_memory_geometry_receipt_v0.1.json").read_text())
    assert FROZEN_INPUT_SHA256 == {
        "candidates": census["candidate_table"]["sha256"],
        "localities": geometry["locality_csv"]["sha256"],
        "edges": geometry["edge_csv"]["sha256"],
    }
