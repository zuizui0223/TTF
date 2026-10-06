import copy

import pytest

from ttf.codistributed_post_result import (
    NULL_DECISION,
    POSITIVE_DECISION,
    audit_empirical_result,
)


AUTH="abc123"


def _result(beta=0.02,p=0.03):
    positive=beta>0 and p<=0.05
    decision=POSITIVE_DECISION if positive else NULL_DECISION
    return {
        "schema":"ttf_genetic_codistributed_recurrence_empirical_result_v0.1",
        "status":"ONE_SHOT_CONFIRMATORY_EMPIRICAL_RESULT_OPENED",
        "authorization_sha256":AUTH,
        "empirical_domain":{
            "species":326,"sources":166,"targets":160,"directed_dyads":26560,
        },
        "primary":{
            "coefficient":beta,
            "primary_envelope_p_value":p,
            "alpha":0.05,
            "positive":positive,
            "least_favourable_cell":"A3",
            "component_monte_carlo_p":{
                "A0p5":0.01,"A1":0.02,"A2":0.025,"A3":p,
            },
        },
        "decision":decision,
        "outcome_state":{
            "confirmatory_nucleotide_identity_opened":True,
            "confirmatory_pairwise_genetic_distances_computed_in_memory":True,
            "confirmatory_T_st_computed":True,
            "confirmatory_beta_G_computed":True,
            "serialized_nucleotide_identity":False,
            "serialized_edge_genetic_distance_vectors":False,
            "serialized_dyad_T_st":False,
        },
        "one_shot":{
            "post_result_retuning_allowed":False,
            "alternate_subgroup_result_allowed":False,
            "result_selection_rerun_allowed":False,
        },
    }


def test_post_result_audit_accepts_frozen_positive_branch():
    a=audit_empirical_result(_result(),expected_authorization_sha256=AUTH)
    assert a.positive is True
    assert a.decision==POSITIVE_DECISION


def test_post_result_audit_accepts_frozen_null_branch():
    a=audit_empirical_result(_result(beta=0.01,p=0.2),expected_authorization_sha256=AUTH)
    assert a.positive is False
    assert a.decision==NULL_DECISION


def test_post_result_audit_rejects_result_relabeling_or_serialization():
    r=_result()
    r["decision"]=NULL_DECISION
    with pytest.raises(RuntimeError,match="decision label"):
        audit_empirical_result(r,expected_authorization_sha256=AUTH)

    r=_result()
    r["outcome_state"]["serialized_dyad_T_st"]=True
    with pytest.raises(RuntimeError,match="serialization boundary"):
        audit_empirical_result(r,expected_authorization_sha256=AUTH)
