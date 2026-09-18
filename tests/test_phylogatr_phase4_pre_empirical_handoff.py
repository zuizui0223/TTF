from __future__ import annotations

import json
from pathlib import Path


HANDOFF = Path(
    "benchmarks/frozen/genetic_phylogatr_phase4_pre_empirical_handoff_v0.1.json"
)
REGENERATION = Path(
    "docs/supporting/genetic_phylogatr_phase1_source_regeneration_v0.1.json"
)


def test_phase4_pre_empirical_handoff_is_closed_and_source_exact() -> None:
    handoff = json.loads(HANDOFF.read_text())
    regeneration = json.loads(REGENERATION.read_text())

    assert handoff["schema"] == "ttf_genetic_phylogatr_phase4_pre_empirical_handoff_v0.1"
    assert handoff["status"] == "PHASE4_AUTHORIZED_WAITING_FOR_EXACT_SOURCE_ARCHIVE_RECOVERY"
    assert handoff["source_main_sha"] == "c846bcf431a9aad52237a50009a91fd549b689ef"

    assert handoff["phase2"]["survivor_species"] == 211
    assert handoff["formal_phase3_gate_d"]["status"] == "PASS"
    assert handoff["formal_phase3_gate_d"]["max_private_null_wilson95_upper"] <= handoff["formal_phase3_gate_d"]["type1_ceiling"]
    assert handoff["formal_phase3_gate_d"]["shared_A2_wilson95_lower"] >= handoff["formal_phase3_gate_d"]["power_floor"]
    assert handoff["self_detectability"]["status"] == "PASS"
    assert handoff["fragility_diagnostic"]["status"] == "DIAGNOSTIC_FRAGILITY_CURVE_COMPLETE"
    assert handoff["fragility_diagnostic"]["diagnostic_metrics_gate_phase4"] is False
    assert handoff["phase4_authorization"]["status"] == "AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING"

    required = handoff["exact_source_archive_required"]
    frozen = regeneration["raw_source"]
    assert required["size_bytes"] == frozen["archive_size_bytes"] == 274988692
    assert required["sha256"] == frozen["raw_archive_sha256"]
    assert required["accept_only_byte_identical_phase1_manifest"] is True
    assert required["accept_only_byte_identical_phase1_geometry"] is True
    assert required["substitute_redownload_or_repack_allowed"] is False

    expected = regeneration["expected_phase1"]
    assert handoff["phase1_original"]["phase1_manifest_sha256"] == expected["phase1_manifest_sha256"]
    assert handoff["phase1_original"]["phase1_geometry_csv_sha256"] == expected["phase1_geometry_csv_sha256"]

    assert all(value is False for value in handoff["outcome_firewall"].values())
    assert handoff["repository_lineage"]["current_blocker_issue"] == 49
    assert handoff["repository_lineage"]["retired_archive_acquisition_issue"] == 21
    assert handoff["repository_lineage"]["retired_parallel_ncbi_source_pull_request"] == 37
