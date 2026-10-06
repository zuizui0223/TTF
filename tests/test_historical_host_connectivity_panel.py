import csv
from pathlib import Path

from scripts.freeze_historical_host_connectivity_panel import build_panel


def _candidates(n=600):
    return [{"species":f"Insect {i:03d}"} for i in range(n)]


def test_historical_host_panel_keeps_only_narrow_resolved_native_hosts():
    candidates=_candidates()
    interactions=[]
    native=[]
    for i in range(120):
        # one to five hosts deterministically
        for j in range(1+(i%5)):
            hid=f"H{i:03d}_{j}"
            interactions.append({
                "insect_species":f"Insect {i:03d}",
                "accepted_plant_name_id":hid,
                "accepted_name":f"Plant {hid}",
            })
            native.append({"accepted_plant_name_id":hid})
    # one broad-host species should be excluded
    for j in range(6):
        hid=f"B_{j}"
        interactions.append({
            "insect_species":"Insect 120",
            "accepted_plant_name_id":hid,
            "accepted_name":f"Plant {hid}",
        })
        native.append({"accepted_plant_name_id":hid})
    panel,summary=build_panel(
        candidates,interactions,native,
        source_sha256="a"*64,
        minimum_species=100,
    )
    assert summary["status"]=="PASS_TO_HOST_OCCURRENCE_FEASIBILITY"
    assert len(panel)==120
    assert [row["host_panel_rank"] for row in panel]==list(range(1,121))
    assert summary["excluded_counts"]["host_breadth_gt_5"]==1
    assert summary["excluded_counts"]["no_wcvp_resolved_species_host"]==479


def test_historical_host_panel_fails_if_native_support_or_size_is_too_small():
    candidates=_candidates()
    interactions=[]
    native=[]
    for i in range(99):
        hid=f"H{i:03d}"
        interactions.append({
            "insect_species":f"Insect {i:03d}",
            "accepted_plant_name_id":hid,
            "accepted_name":f"Plant {hid}",
        })
        native.append({"accepted_plant_name_id":hid})
    # resolved but no native distribution
    interactions.append({
        "insect_species":"Insect 099",
        "accepted_plant_name_id":"NO_NATIVE",
        "accepted_name":"Plant missing",
    })
    panel,summary=build_panel(
        candidates,interactions,native,
        source_sha256="b"*64,
        minimum_species=100,
    )
    assert panel==[]
    assert summary["status"]=="NOT_EVALUABLE_HISTORICAL_HOST_CONNECTIVITY_HOST_DOMAIN_TOO_SMALL"
    assert summary["retained_species"]==99
    assert summary["excluded_counts"]["host_without_primary_native_wgsrpd3"]==1
