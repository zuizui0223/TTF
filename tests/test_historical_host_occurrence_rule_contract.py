"""Ensure the frozen occurrence contract is parseable before any GBIF request."""
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def test_historical_host_occurrence_json_frozen_science():
    path=ROOT/"docs/supporting/genetic_historical_host_occurrence_rule_v0.1.json"
    rule=json.loads(path.read_text(encoding="utf-8"))
    assert rule["schema"]=="ttf_genetic_historical_host_occurrence_rule_v0.1"
    assert rule["status"]=="FROZEN_BEFORE_WCVP_PANEL_RESULT_OR_ANY_HOST_GBIF_RESULT"
    assert rule["native_range_geometry"]["primary_native_filter"]=={
        "introduced":0,"extinct":0,"location_doubtful":0
    }
    assert rule["gbif"]["page_size"]==300
    assert rule["gbif"]["maximum_pages"]==20
    assert rule["gbif"]["minimum_distance_km"]==10
    assert rule["gbif"]["minimum_retained_per_host"]==30
    assert rule["insect_host_completeness"]["required_fraction"]==0.8
    assert rule["insect_host_completeness"]["minimum_insects_after_occurrence_gate"]==90
    assert rule["selection_boundary"]==[
        "all_occurrence filtering and insect removal occur before palaeo host-resistance calculation and before any genetic nucleotide identity",
        "host or insect retention cannot depend on a genetic response",
        "request errors are technical incompleteness and must be resolved or the affected host remains unsupported; they cannot be relabelled ecological absence",
    ]
    assert all(v is False for v in rule["response_firewall"].values())
