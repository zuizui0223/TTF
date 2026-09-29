import json
import runpy
from dataclasses import make_dataclass
from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]


def test_v03_nonlepidoptera_panel_is_disjoint_and_frozen():
    rule=json.loads(
        (ROOT/"docs/supporting/genetic_palearctic_lgm_subpanel_rule_v0.3.json").read_text()
    )
    census=json.loads(
        (ROOT/"benchmarks/frozen/genetic_palearctic_lgm_nonlepidoptera_census_v0.3.json").read_text()
    )
    panel=json.loads(
        (ROOT/"benchmarks/frozen/genetic_palearctic_lgm_panel_metadata_v0.3.json").read_text()
    )

    assert rule["status"].startswith("FROZEN_BEFORE_")
    assert census["counts"]["retained_species"] == 23
    assert census["counts"]["source_clusters"] == 11
    assert census["counts"]["target_clusters"] == 12
    assert census["counts"]["directed_source_target_dyads"] == 132

    rows=panel["rows"]
    assert len(rows)==23
    assert all(row["class"]=="Insecta" for row in rows)
    assert all(row["order"]!="Lepidoptera" for row in rows)
    assert {row["order"] for row in rows} == {"Coleoptera","Hymenoptera","Diptera"}

    source={row["species"] for row in rows if row["role"]=="source"}
    target={row["species"] for row in rows if row["role"]=="target"}
    assert len(source)==11
    assert len(target)==12
    assert source.isdisjoint(target)
    assert source | target == {row["species"] for row in rows}

    firewall=panel["response_firewall"]
    assert firewall["species_level_genetic_scores_used"] is False
    assert firewall["pairwise_v03_T_st_computed"] is False
    assert firewall["beta_LGM_computed"] is False


def test_v03_ttfq_signal_breadth_scales_with_realized_endpoint_count():
    ns=runpy.run_path(str(ROOT/"scripts/run_palearctic_lgm_ttf_q_v03.py"))
    Info=make_dataclass(
        "Info",
        [
            ("total_unique_variance_fraction",float),
            ("source_clusters",int),
            ("target_clusters",int),
        ],
    )
    Support=make_dataclass(
        "Support",
        [
            ("signal_effective_sources",float),
            ("signal_effective_targets",float),
            ("max_source_signal_share",float),
            ("max_target_signal_share",float),
        ],
    )
    rule={
        "minimum_total_unique_variance_fraction":0.15,
        "minimum_effective_source_fraction":0.50,
        "minimum_effective_target_fraction":0.50,
        "maximum_single_source_signal_share":0.15,
        "maximum_single_target_signal_share":0.15,
        "required_null_qualified_amplitudes":[0,1,2],
        "required_A2_evaluable_mde_max":0.10,
    }
    env={
        "0":{"calibration_pass":True,"evaluable_grid_mde":0.05},
        "1":{"calibration_pass":True,"evaluable_grid_mde":0.08},
        "2":{"calibration_pass":True,"evaluable_grid_mde":0.10},
        "3":{"calibration_pass":False,"evaluable_grid_mde":None},
    }
    gates,fs,ft=ns["evaluate_opening"](
        Info(0.20,11,12),
        Support(6.0,6.5,0.14,0.14),
        env,
        rule,
    )
    assert fs == 6.0/11.0
    assert ft == 6.5/12.0
    assert gates["overall_pass"] is True

    bad,_,_=ns["evaluate_opening"](
        Info(0.20,11,12),
        Support(5.0,6.5,0.14,0.14),
        env,
        rule,
    )
    assert bad["source_signal_breadth"] is False
    assert bad["overall_pass"] is False
