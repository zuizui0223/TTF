from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from pathlib import Path
from typing import Any


EXPECTED_CONTRACT="ttf_genetic_codistributed_recurrence_empirical_execution_contract_v0.1"
EXPECTED_QUALIFICATION="ttf_genetic_codistributed_recurrence_survivor_fixed_tst_qualification_v0.1"
EXPECTED_AUTHORIZATION="ttf_genetic_codistributed_recurrence_identity_opening_authorization_v0.1"


def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda:handle.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def require_closed_response_firewall(payload: dict[str,Any], *, label: str) -> None:
    state=payload.get("response_firewall")
    if not isinstance(state,dict) or not state:
        raise RuntimeError(f"{label} response firewall missing")
    if any(bool(value) for value in state.values()):
        raise RuntimeError(f"{label} response firewall is open")


def build_execution_receipt(
    contract: dict[str,Any],
    qualification: dict[str,Any],
    authorization: dict[str,Any],
    *,
    contract_sha256: str,
    qualification_sha256: str,
    authorization_sha256: str,
    empirical_result_exists: bool,
) -> dict[str,Any]:
    if contract.get("schema")!=EXPECTED_CONTRACT:
        raise RuntimeError("unexpected empirical execution contract")
    if contract.get("status")!="FROZEN_BEFORE_AUTHORITATIVE_SURVIVOR_REQUALIFICATION_RESULT":
        raise RuntimeError("empirical execution contract is not prospectively frozen")
    require_closed_response_firewall(contract,label="execution contract")

    qspec=contract.get("authoritative_survivor_qualification")
    if not isinstance(qspec,dict):
        raise RuntimeError("authoritative survivor qualification contract missing")
    run_id=int(qspec.get("workflow_run_id",-1))
    if run_id<1:
        raise RuntimeError("invalid authoritative survivor workflow run id")

    if qualification.get("schema")!=EXPECTED_QUALIFICATION:
        raise RuntimeError("unexpected survivor qualification schema")
    if qualification.get("status")!=qspec.get("required_result_status"):
        raise RuntimeError("survivor qualification did not reach the frozen required status")
    if qualification.get("gate",{}).get("overall_pass") is not True:
        raise RuntimeError("survivor qualification overall gate is not PASS")
    require_closed_response_firewall(qualification,label="survivor qualification")

    if authorization.get("schema")!=EXPECTED_AUTHORIZATION:
        raise RuntimeError("unexpected identity-opening authorization schema")
    if authorization.get("status")!="AUTHORIZE_ONE_SHOT_CODISTRIBUTED_NUCLEOTIDE_IDENTITY_OPENING":
        raise RuntimeError("identity opening is not authorized")
    aq=authorization.get("qualification")
    if not isinstance(aq,dict):
        raise RuntimeError("authorization qualification binding missing")
    if aq.get("sha256")!=qualification_sha256:
        raise RuntimeError("authorization is bound to a different qualification")
    if aq.get("status")!=qualification.get("status") or aq.get("overall_pass") is not True:
        raise RuntimeError("authorization qualification state drift")

    opening=authorization.get("identity_opening")
    if not isinstance(opening,dict) or opening.get("authorized_once") is not True:
        raise RuntimeError("authorization is not one-shot")
    for key in (
        "nucleotide_identity_opened_at_authorization_time",
        "pairwise_genetic_distances_opened_at_authorization_time",
        "empirical_T_st_opened_at_authorization_time",
        "empirical_beta_G_opened_at_authorization_time",
    ):
        if opening.get(key) is not False:
            raise RuntimeError(f"pre-opening authorization state is not closed: {key}")

    source=contract.get("source",{})
    if authorization.get("source",{}).get("archive_sha256")!=source.get("archive_sha256"):
        raise RuntimeError("authorization source archive differs from execution contract")

    if empirical_result_exists:
        raise RuntimeError("empirical result already exists; one-shot execution receipt refused")

    return {
        "schema":"ttf_genetic_codistributed_recurrence_empirical_execution_v0.1",
        "status":"AUTHORIZE_ONE_SHOT_EMPIRICAL_EXECUTION",
        "survivor_qualification_workflow_run_id":run_id,
        "survivor_qualification_sha256":qualification_sha256,
        "identity_opening_authorization_sha256":authorization_sha256,
        "empirical_execution_contract_sha256":contract_sha256,
        "source_archive_sha256":source.get("archive_sha256"),
        "empirical_result_observed_before_execution":False,
        "one_shot":{
            "github_event_required":"push",
            "github_run_attempt_required":1,
            "rerun_for_result_selection_allowed":False,
            "post_result_retuning_allowed":False,
        },
        "response_firewall":{
            "confirmatory_nucleotide_identity_opened":False,
            "confirmatory_pairwise_genetic_distances_opened":False,
            "confirmatory_T_st_opened":False,
            "confirmatory_beta_G_opened":False,
        },
    }


__all__=[
    "build_execution_receipt",
    "require_closed_response_firewall",
    "sha256_path",
]
