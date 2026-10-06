import copy

import pytest

from ttf.codistributed_execution import build_execution_receipt


def _fixtures():
    contract={
        "schema":"ttf_genetic_codistributed_recurrence_empirical_execution_contract_v0.1",
        "status":"FROZEN_BEFORE_AUTHORITATIVE_SURVIVOR_REQUALIFICATION_RESULT",
        "authoritative_survivor_qualification":{
            "workflow_run_id":37406111501,
            "required_result_status":"PASS_TO_ONE_SHOT_CONFIRMATORY_EMPIRICAL_OPENING",
        },
        "source":{"archive_sha256":"abc"},
        "response_firewall":{
            "confirmatory_nucleotide_identity_opened":False,
            "confirmatory_pairwise_genetic_distances_opened":False,
            "confirmatory_T_st_opened":False,
            "confirmatory_beta_G_opened":False,
        },
    }
    qualification={
        "schema":"ttf_genetic_codistributed_recurrence_survivor_fixed_tst_qualification_v0.1",
        "status":"PASS_TO_ONE_SHOT_CONFIRMATORY_EMPIRICAL_OPENING",
        "gate":{"overall_pass":True},
        "response_firewall":{
            "confirmatory_nucleotide_identity_opened":False,
            "confirmatory_pairwise_genetic_distances_opened":False,
            "confirmatory_T_st_opened":False,
            "confirmatory_beta_G_opened":False,
        },
    }
    authorization={
        "schema":"ttf_genetic_codistributed_recurrence_identity_opening_authorization_v0.1",
        "status":"AUTHORIZE_ONE_SHOT_CODISTRIBUTED_NUCLEOTIDE_IDENTITY_OPENING",
        "source":{"archive_sha256":"abc"},
        "qualification":{
            "sha256":"qualsha",
            "status":"PASS_TO_ONE_SHOT_CONFIRMATORY_EMPIRICAL_OPENING",
            "overall_pass":True,
        },
        "identity_opening":{
            "authorized_once":True,
            "nucleotide_identity_opened_at_authorization_time":False,
            "pairwise_genetic_distances_opened_at_authorization_time":False,
            "empirical_T_st_opened_at_authorization_time":False,
            "empirical_beta_G_opened_at_authorization_time":False,
        },
    }
    return contract,qualification,authorization


def test_execution_receipt_binds_authoritative_pass_and_closed_authorization():
    contract,qualification,authorization=_fixtures()
    p=build_execution_receipt(
        contract,qualification,authorization,
        contract_sha256="contractsha",
        qualification_sha256="qualsha",
        authorization_sha256="authsha",
        empirical_result_exists=False,
    )
    assert p["status"]=="AUTHORIZE_ONE_SHOT_EMPIRICAL_EXECUTION"
    assert p["survivor_qualification_workflow_run_id"]==37406111501
    assert p["survivor_qualification_sha256"]=="qualsha"
    assert p["identity_opening_authorization_sha256"]=="authsha"
    assert p["empirical_result_observed_before_execution"] is False
    assert all(v is False for v in p["response_firewall"].values())


def test_execution_receipt_refuses_scientific_fail_or_opened_state():
    contract,qualification,authorization=_fixtures()
    q=copy.deepcopy(qualification)
    q["gate"]["overall_pass"]=False
    with pytest.raises(RuntimeError,match="overall gate"):
        build_execution_receipt(
            contract,q,authorization,
            contract_sha256="contractsha",
            qualification_sha256="qualsha",
            authorization_sha256="authsha",
            empirical_result_exists=False,
        )

    a=copy.deepcopy(authorization)
    a["identity_opening"]["nucleotide_identity_opened_at_authorization_time"]=True
    with pytest.raises(RuntimeError,match="not closed"):
        build_execution_receipt(
            contract,qualification,a,
            contract_sha256="contractsha",
            qualification_sha256="qualsha",
            authorization_sha256="authsha",
            empirical_result_exists=False,
        )


def test_execution_receipt_refuses_wrong_artifact_or_existing_result():
    contract,qualification,authorization=_fixtures()
    a=copy.deepcopy(authorization)
    a["qualification"]["sha256"]="other"
    with pytest.raises(RuntimeError,match="different qualification"):
        build_execution_receipt(
            contract,qualification,a,
            contract_sha256="contractsha",
            qualification_sha256="qualsha",
            authorization_sha256="authsha",
            empirical_result_exists=False,
        )
    with pytest.raises(RuntimeError,match="already exists"):
        build_execution_receipt(
            contract,qualification,authorization,
            contract_sha256="contractsha",
            qualification_sha256="qualsha",
            authorization_sha256="authsha",
            empirical_result_exists=True,
        )
