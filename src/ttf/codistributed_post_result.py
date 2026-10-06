from __future__ import annotations

from dataclasses import dataclass
from typing import Any


EXPECTED_SCHEMA="ttf_genetic_codistributed_recurrence_empirical_result_v0.1"
EXPECTED_STATUS="ONE_SHOT_CONFIRMATORY_EMPIRICAL_RESULT_OPENED"
POSITIVE_DECISION="CONFIRMATORY_POSITIVE_GEOGRAPHIC_RECURRENCE_CONDITIONAL_ON_SHARED_OPPORTUNITY"
NULL_DECISION="CONFIRMATORY_QUALIFIED_NULL_SHARED_GEOGRAPHIC_OPPORTUNITY_INSUFFICIENT"


@dataclass(frozen=True)
class CodistributedResultAudit:
    decision: str
    beta_G: float
    primary_p_value: float
    alpha: float
    positive: bool
    least_favourable_cell: str
    serialization_boundary_pass: bool
    one_shot_boundary_pass: bool
    domain_pass: bool


def audit_empirical_result(
    result: dict[str,Any],
    *,
    expected_authorization_sha256: str,
) -> CodistributedResultAudit:
    if result.get("schema")!=EXPECTED_SCHEMA:
        raise RuntimeError("unexpected codistributed empirical result schema")
    if result.get("status")!=EXPECTED_STATUS:
        raise RuntimeError("empirical result did not reach the frozen one-shot result state")
    if result.get("authorization_sha256")!=expected_authorization_sha256:
        raise RuntimeError("empirical result authorization hash drift")

    domain=result.get("empirical_domain",{})
    domain_pass=(
        int(domain.get("species",-1))==326
        and int(domain.get("sources",-1))==166
        and int(domain.get("targets",-1))==160
        and int(domain.get("directed_dyads",-1))==26560
    )
    if not domain_pass:
        raise RuntimeError("empirical result domain drift")

    primary=result.get("primary",{})
    beta=float(primary["coefficient"])
    p=float(primary["primary_envelope_p_value"])
    alpha=float(primary["alpha"])
    positive=bool(primary["positive"])
    if not (0.0<=p<=1.0 and alpha==0.05):
        raise RuntimeError("primary p-value/alpha drift")
    component={str(k):float(v) for k,v in primary["component_monte_carlo_p"].items()}
    if set(component)!={"A0p5","A1","A2","A3"}:
        raise RuntimeError("private-envelope component set drift")
    max_p=max(component.values())
    if abs(max_p-p)>1e-15:
        raise RuntimeError("primary p-value is not the least-favourable component maximum")
    least=str(primary["least_favourable_cell"])
    if least not in component or abs(component[least]-p)>1e-15:
        raise RuntimeError("least-favourable cell label drift")

    expected_positive=bool(beta>0.0 and p<=alpha)
    if positive!=expected_positive:
        raise RuntimeError("positive flag disagrees with frozen decision rule")
    expected_decision=POSITIVE_DECISION if expected_positive else NULL_DECISION
    if result.get("decision")!=expected_decision:
        raise RuntimeError("decision label disagrees with frozen interpretation branch")

    state=result.get("outcome_state",{})
    serialization_pass=(
        state.get("confirmatory_nucleotide_identity_opened") is True
        and state.get("confirmatory_pairwise_genetic_distances_computed_in_memory") is True
        and state.get("confirmatory_T_st_computed") is True
        and state.get("confirmatory_beta_G_computed") is True
        and state.get("serialized_nucleotide_identity") is False
        and state.get("serialized_edge_genetic_distance_vectors") is False
        and state.get("serialized_dyad_T_st") is False
    )
    if not serialization_pass:
        raise RuntimeError("empirical serialization boundary violated")

    one=result.get("one_shot",{})
    one_shot_pass=(
        one.get("post_result_retuning_allowed") is False
        and one.get("alternate_subgroup_result_allowed") is False
        and one.get("result_selection_rerun_allowed") is False
    )
    if not one_shot_pass:
        raise RuntimeError("one-shot boundary violated")

    return CodistributedResultAudit(
        decision=expected_decision,
        beta_G=beta,
        primary_p_value=p,
        alpha=alpha,
        positive=expected_positive,
        least_favourable_cell=least,
        serialization_boundary_pass=True,
        one_shot_boundary_pass=True,
        domain_pass=True,
    )


__all__=[
    "CodistributedResultAudit",
    "EXPECTED_SCHEMA",
    "EXPECTED_STATUS",
    "NULL_DECISION",
    "POSITIVE_DECISION",
    "audit_empirical_result",
]
