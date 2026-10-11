"""Unit tests use artificial taxon and mask-status tables; no raw sequence data."""
import copy

import pytest

from scripts.audit_historical_host_memory_mask_selection_bias import aggregate


def synthetic_sources():
    roles=[]
    candidates=[]
    diagnostics=[]
    for i in range(642):
        name=f"Taxon_{i:03d}"
        is_confirm=i<321
        if i<641:
            roles.append({"species":name,
                          "panel":"confirmatory" if is_confirm else "development"})
        if is_confirm:
            pass_mask=i<259
            edge_n=(43 if i<258 else 9268) if pass_mask else (100 if i<320 else 15130)
        else:
            edge_n=25
        candidates.append({"species":name,"family":"FamilyX" if i%2 else "FamilyY",
                           "edges":str(edge_n),"n_localities":"12"})
        diagnostics.append({"species":name,
                            "status":"complete" if i<641 else "nonfinite_edge_climate",
                            "unique_fraction_M_host":"0.24",
                            "predictor_condition_number":"4"})
    rows=[]
    for i in range(321):
        passed=i<259
        invalid=0 if passed or i==320 else (8 if i<319 else 33)
        rows.append({"species":f"Taxon_{i:03d}",
                     "status":"PASS_MASK" if passed else "NOT_EVALUABLE_CHARACTER_SUPPORT",
                     "survives":passed,
                     "reason_code":None if passed else
                         ("malformed_or_duplicate_frozen_alignment" if i==320 else
                          "one_or_more_frozen_edges_have_no_valid_cross_locality_pair"),
                     "invalid_edges":invalid})
    # 60 edge-support failures with 8 invalid edges each; one with 41;
    # one malformed record with 0: 480 + 41 = 521.
    rows[319]["invalid_edges"]=41
    # Original malformed alignment has no measurable edge-validity tally.
    rows[320]["invalid_edges"]=None
    mask={"schema":"ttf_historical_host_memory_independent_mask_replay_v0.1",
          "source_archive_sha256":"5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5",
          "genetic_response_not_opened":True,
          "genetic_distance_not_computed":True,
          "nucleotide_base_identity_not_persisted":True,
          "rows":rows}
    return mask,candidates,diagnostics,roles


def test_only_aggregate_output_and_original_graph_counts():
    p=aggregate(*synthetic_sources())
    assert p["mask_pass_species"]==259
    assert p["mask_fail_species"]==62
    assert p["original_confirmatory_edges"]==41592
    assert p["survivor_confirmatory_edges"]==20362
    assert p["invalid_original_edges_with_finite_masks"]==521
    assert p["group_medians"]["pass"]["host_history_unique_fraction"]==pytest.approx(0.24)
    assert len(p["original_graph_edge_count_bins"])==5
    assert p["interpretation"]["observed_genetic_distance_or_effect_unopened"] is True
    assert "rows" not in p and "species" not in p


def test_refuses_forged_mask_membership_and_unblinding():
    mask,candidates,diag,roles=synthetic_sources()
    mask["rows"][0]["species"]="Taxon_500"
    with pytest.raises(RuntimeError,match="does not equal original confirmatory"):
        aggregate(mask,candidates,diag,roles)
    mask,candidates,diag,roles=synthetic_sources()
    mask["genetic_distance_not_computed"]=False
    with pytest.raises(RuntimeError,match="not proven response blind"):
        aggregate(mask,candidates,diag,roles)


def test_refuses_inconsistent_original_edge_or_mask_failure_count():
    mask,candidates,diag,roles=synthetic_sources()
    candidates[0]["edges"]="99"
    with pytest.raises(RuntimeError,match="edge accounting drift"):
        aggregate(mask,candidates,diag,roles)
    mask,candidates,diag,roles=synthetic_sources()
    mask["rows"][258]["survives"]=False
    with pytest.raises(RuntimeError,match="inconsistent original mask"):
        aggregate(mask,candidates,diag,roles)


def test_result_is_invariant_to_row_order():
    a=synthetic_sources()
    b=copy.deepcopy(a)
    for part in (b[1],b[2],b[3]):
        part.reverse()
    b[0]["rows"].reverse()
    assert aggregate(*a)==aggregate(*b)
