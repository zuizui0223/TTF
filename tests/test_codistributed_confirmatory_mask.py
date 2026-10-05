import json
from pathlib import Path
import sys

import pytest


ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

from freeze_codistributed_recurrence_confirmatory_mask import (
    authorize_formal,
    digest_labels,
    frozen_confirmatory_roles,
    read_candidates,
    read_excluded,
)


def test_confirmatory_roles_reproduce_frozen_design_digests():
    design=json.loads(
        (ROOT/"benchmarks/frozen/genetic_codistributed_recurrence_response_blind_design_v0.1.json").read_text()
    )
    rows=read_candidates(ROOT/"benchmarks/frozen/relational_fresh_candidate_1000_v0.1.csv")
    s1=read_excluded(ROOT/"benchmarks/frozen/relational_prior_S1_species_exclusion_v0.1.json")
    s2=read_excluded(ROOT/"benchmarks/frozen/relational_prior_S2_species_exclusion_v0.1.json")
    source,target=frozen_confirmatory_roles(rows,s1,s2,design)
    assert len(source)==190
    assert len(target)==191
    assert set(source).isdisjoint(target)
    assert digest_labels(source)==design["selection"]["confirmatory"]["source_digest_sha256"]
    assert digest_labels(target)==design["selection"]["confirmatory"]["target_digest_sha256"]


def test_confirmatory_mask_requires_formal_pass(tmp_path):
    fail={
        "schema":"ttf_genetic_codistributed_recurrence_fixed_tst_qualification_v0.3",
        "status":"NOT_EVALUABLE_CODISTRIBUTED_RECURRENCE_FIXED_TST",
        "gate":{"overall_pass":False},
        "response_firewall":{
            "development_nucleotide_identity_opened":False,
            "confirmatory_nucleotide_identity_opened":False,
            "empirical_T_st_opened":False,
            "empirical_beta_G_opened":False,
        },
    }
    path=tmp_path/"formal.json"
    path.write_text(json.dumps(fail))
    with pytest.raises(RuntimeError,match="does not authorize"):
        authorize_formal(path)

    fail["status"]="PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION"
    fail["gate"]["overall_pass"]=True
    fail["response_firewall"]["empirical_T_st_opened"]=True
    path.write_text(json.dumps(fail))
    with pytest.raises(RuntimeError,match="firewall"):
        authorize_formal(path)


def test_confirmatory_mask_accepts_only_closed_formal_pass(tmp_path):
    passed={
        "schema":"ttf_genetic_codistributed_recurrence_fixed_tst_qualification_v0.3",
        "status":"PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION",
        "gate":{"overall_pass":True},
        "response_firewall":{
            "development_nucleotide_identity_opened":False,
            "confirmatory_nucleotide_identity_opened":False,
            "empirical_T_st_opened":False,
            "empirical_beta_G_opened":False,
        },
    }
    path=tmp_path/"formal.json"
    path.write_text(json.dumps(passed))
    assert authorize_formal(path)["status"]=="PASS_TO_CONFIRMATORY_CHARACTER_MASK_PREPARATION"
