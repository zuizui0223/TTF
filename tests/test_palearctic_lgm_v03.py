import hashlib
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


def test_v03_role_hash_and_archive_provenance_are_canonical():
    rule=json.loads(
        (ROOT/"docs/supporting/genetic_palearctic_lgm_subpanel_rule_v0.3.json").read_text()
    )
    panel=json.loads(
        (ROOT/"benchmarks/frozen/genetic_palearctic_lgm_panel_metadata_v0.3.json").read_text()
    )
    correction=json.loads(
        (ROOT/"benchmarks/frozen/genetic_palearctic_lgm_archive_hash_provenance_correction_v0.1.json").read_text()
    )

    canonical=correction["canonical_exact_phylogatr_archive_sha256"]
    assert len(canonical)==64
    assert rule["role_assignment"]["exact_archive_sha256"]==canonical
    assert correction["status"]=="PROVENANCE_CORRECTION_NO_SCIENTIFIC_CHANGE"

    namespace=rule["role_assignment"]["namespace"]
    species=sorted(row["species"] for row in panel["rows"])
    ranked=sorted(
        species,
        key=lambda sp: (
            hashlib.sha256(f"{namespace}|{canonical}|{sp}".encode()).hexdigest(),
            sp,
        ),
    )
    expected_source=set(ranked[:11])
    expected_target=set(ranked[11:])
    observed_source={row["species"] for row in panel["rows"] if row["role"]=="source"}
    observed_target={row["species"] for row in panel["rows"] if row["role"]=="target"}
    assert observed_source==expected_source
    assert observed_target==expected_target


def test_v03_qualification_workflow_is_response_blind_and_complete():
    text=(ROOT/".github/workflows/palearctic-lgm-climate-ttfq.yml").read_text()
    required=[
        "36529160168",
        "genetic_palearctic_lgm_v03_climate_asset_binding_v0.1.json",
        "10791389409",
        "10791867875",
        "bind_palearctic_lgm_occurrences_v03.py",
        "prepare_palearctic_lgm_climate_inputs_v03.py",
        "build_palearctic_lgm_relation_v03.py",
        "run_palearctic_lgm_ttf_q_v03.py",
        "genetic_palearctic_lgm_ttf_q_execution_v0.3.json",
        "palearctic-lgm-v03-response-blind-qualification",
    ]
    for token in required:
        assert token in text
    forbidden=[
        "run_palearctic_lgm_genetic",
        "compute_palearctic_lgm_T_st",
        "pairwise_v03_T_st_computed: true",
        "beta_LGM_computed: true",
    ]
    for token in forbidden:
        assert token not in text


def test_v03_qualification_requires_frozen_aggregate_trigger():
    text=(ROOT/".github/workflows/palearctic-lgm-climate-ttfq.yml").read_text()
    token="benchmarks/frozen/genetic_palearctic_lgm_v03_qualification_trigger_v0.1.json"
    assert token in text
    assert "AUTHORIZE_RESPONSE_BLIND_V03_QUALIFICATION" in text
    assert "aggregate_artifact_id" in text
    assert "aggregate_artifact_digest_sha256" in text
