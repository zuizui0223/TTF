"""Response-blind integrity checks for the prequery historical-host correction."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FROZEN=ROOT/"benchmarks"/"frozen"

V1_SHA="67e8d6ad838a6464e8557818af5324f107605050348ead8fb26d20e0f8afd0f4"
V2_SHA="a402483497bc7bdf9b7427bfc859ce9dcc9183eebf5964d4591f352847130484"
REMOVED={
    "Apoda y-inversum",
    "Chionodes electella",
    "Herpetogramma aeglealis",
    "Nisoniades rubescens",
    "Stigmella hybnerella",
}

def load_rows(name: str):
    path=FROZEN/name
    with path.open(newline="", encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    return hashlib.sha256(path.read_bytes()).hexdigest(),rows

def test_old_frozen_panel_restores_exact_receipt_and_host_id():
    h,rows=load_rows("genetic_historical_host_connectivity_panel_v0.1.csv")
    assert h==V1_SHA
    assert len(rows)==140
    row=next(x for x in rows if x["species"]=="Craniophora ligustri")
    assert row["accepted_host_ids"]=="353907;369664;6364"
    assert row["accepted_host_names"]=="Ligustrum vulgare;Fraxinus excelsior;Alnus glutinosa"

def test_v02_is_only_five_prequery_design_removals_and_no_backfill():
    h1,old=load_rows("genetic_historical_host_connectivity_panel_v0.1.csv")
    h2,new=load_rows("genetic_historical_host_connectivity_panel_v0.2.csv")
    assert (h1,h2)==(V1_SHA,V2_SHA)
    assert len(new)==135
    assert {x["species"] for x in old}-{x["species"] for x in new}==REMOVED
    assert [x["species"] for x in old if x["species"] not in REMOVED]==[x["species"] for x in new]
    o={x["species"]:x for x in old}
    for i,row in enumerate(new,1):
        before=o[row["species"]]
        assert int(row["host_panel_rank"])==i
        for key in ("n_hosts","accepted_host_ids","accepted_host_names"):
            assert row[key]==before[key]

def test_new_receipt_is_response_blind_and_invalidates_old_run_authority():
    audit=json.loads((FROZEN/"genetic_historical_host_prequery_integrity_audit_v0.2.json").read_text())
    assert audit["status"]=="RESPONSE_BLIND_PREQUERY_PANEL_CORRECTED_HOLD_EXTERNAL_EXECUTION"
    assert audit["full_design_freshness"]["complete_exclusion_union_species"]==2492
    assert set(audit["full_design_freshness"]["previous_panel_140_overlap_with_conditional_design"])==REMOVED
    assert audit["replacement_panel"]["species"]==135
    assert audit["replacement_panel"]["csv_sha256"]==V2_SHA
    assert audit["replacement_panel"]["no_species_added_or_backfilled"] is True
    assert audit["replacement_panel"]["requires_fresh_external_execution_authority"] is True
    assert audit["external_runs_invalidated_for_new_panel"]["both_stopped_before_any_external_data_query"] is True
    assert not any(audit["outcome_firewall"].values())
    assert audit["historical_biology_caution"]["present_host_association_does_not_prove_21ka_interaction"] is True
