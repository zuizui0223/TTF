"""Check that old response-blind design overlap is disclosed, not data-driven filtered."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FREEZE=ROOT/"benchmarks/frozen"


def test_historical_host_cross_design_overlap_is_exact_and_non_authoritative():
    receipt=json.loads((FREEZE/"genetic_historical_host_cross_program_overlap_audit_v0.1.json").read_text())
    assert receipt["status"]=="FROZEN_RESPONSE_BLIND_CROSS_DESIGN_OVERLAP_DISCLOSURE"

    current=FREEZE/"genetic_historical_host_connectivity_panel_v0.2.csv"
    historical=FREEZE/"relational_historical_candidates_v0.1.csv"
    with current.open(newline="",encoding="utf8") as f:
        panel=[x["species"] for x in csv.DictReader(f)]
    with historical.open(newline="",encoding="utf8") as f:
        prior={x["species"] for x in csv.DictReader(f)}

    assert len(panel)==135
    assert len(panel)==len(set(panel))
    assert hashlib.sha256(current.read_bytes()).hexdigest()==receipt["fixed_primary_panel"]["sha256"]
    assert receipt["intersection_count"]==4
    assert sorted(set(panel)&prior)==receipt["intersecting_species"]
    assert receipt["fixed_primary_panel"]["backfill"] is False
    assert receipt["fixed_primary_panel"]["species_changed"] is False
    assert receipt["historical_C_geometry_design"]["genetic_response_opened"] is False
    assert all(x is False for x in receipt["response_firewall"].values())

    c=json.loads((FREEZE/"relational_program_final_state_v0.3.json").read_text())
    assert c["C_state"]=="NOT_EVALUABLE"
    assert c["final_state"]=="CLOSED_NO_EVALUABLE_TEST"
