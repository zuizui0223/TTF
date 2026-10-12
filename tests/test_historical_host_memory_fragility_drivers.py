"""Synthetic-data guard tests; never read real genetic sequences or private coefficients."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_historical_host_memory_fragility_drivers import summarize_inputs


def toy(n=20, with_invariant=False):
    species, rows, candidates = [], [], []
    for i in range(n):
        name = f"Toy species {i:02d}"
        beta = (-1) ** i * (0.03 + 0.04 * i)
        unique = 0.07 + (i % 9) * 0.085
        nloc = 12 + (i % 6)
        invariant = with_invariant and i in (0, 1, 2, 3)
        coefficient = 0 if invariant else beta * 0.006
        nflip = 0
        for loc in range(nloc):
            flipping = not invariant and (i % 3 == 0 or i % 7 == 0) and loc == 0
            b = -coefficient if flipping else coefficient
            nflip += bool(flipping)
            rows.append({"species": name, "locality_index": loc,
                         "remaining_edges": 18 + i,
                         "status": "defined", "slope": b, "reason": "", "post_ibd_sd": 0.003})
        species.append({
            "species": name, "family": f"Family{i%5}", "edges": 22+i, "localities": nloc,
            "standardized_beta": np.nan if invariant else beta,
            "unscaled_slope": coefficient,
            "unique_fraction": unique, "max_prediction_leverage": 0.15 + i / 100, "semipartial_correlation":
                np.nan if invariant else beta * np.sqrt(unique),
            "locality_deletions_attempted": nloc, "locality_deletions_defined": nloc,
            "locality_deletions_undefined": 0, "locality_deletions_reversing_sign": nflip
        })
        candidates.append({"species": name, "host_species": 1+i%3,
                           "native_units": 22 + i*2})
    return pd.DataFrame(species), pd.DataFrame(rows), pd.DataFrame(candidates)


def test_counts_and_quartiles_are_not_a_new_confirmatory_test():
    result = summarize_inputs(*toy(), expected_full_counts=False)
    assert result["panel"]["n_nonzero_defined"] == 20
    assert result["panel"]["n_baseline_invariant"] == 0
    assert result["panel"]["n_with_any_sign_reversal"] == 9
    assert sum(r["n"] for r in result["absolute_semipartial_quartiles"]) == 20
    assert all(row["n"] == 5 for row in result["absolute_semipartial_quartiles"])
    assert result["original_confirmatory_statistic_status"] == "NOT_EVALUABLE_EMPIRICAL_RESPONSE"
    assert result["new_p_values"] == 0
    assert "Toy species" not in str(result)


def test_four_invariant_species_are_never_imputed_or_falsely_called_stable():
    result = summarize_inputs(*toy(with_invariant=True), expected_full_counts=False)
    assert result["panel"]["n_all"] == 20
    assert result["panel"]["n_nonzero_defined"] == 16
    assert result["panel"]["n_baseline_invariant"] == 4
    assert sum(r["n"] for r in result["absolute_standardized_beta_quartiles"]) == 16


def test_inconsistent_locality_reversal_counts_are_rejected():
    a,b,c = toy()
    a.loc[0, "locality_deletions_reversing_sign"] += 1
    with pytest.raises(RuntimeError, match="sign reversals"):
        summarize_inputs(a,b,c,expected_full_counts=False)


def test_invalid_partial_correlation_identity_is_rejected():
    a,b,c = toy()
    a.loc[0, "semipartial_correlation"] = 0.99
    with pytest.raises(RuntimeError,match="semipartial"):
        summarize_inputs(a,b,c,expected_full_counts=False)


def test_duplicate_species_and_locality_are_rejected():
    a,b,c = toy()
    with pytest.raises(RuntimeError,match="duplicate"):
        summarize_inputs(pd.concat([a,a.iloc[[0]]]),b,c,expected_full_counts=False)
    with pytest.raises(RuntimeError,match="duplicate"):
        summarize_inputs(a,pd.concat([b,b.iloc[[0]]]),c,expected_full_counts=False)


def test_shuffled_private_rows_do_not_change_public_report():
    a,b,c=toy()
    result=summarize_inputs(a,b,c,expected_full_counts=False)
    alt=summarize_inputs(a.sample(frac=1,random_state=1),
                         b.sample(frac=1,random_state=2),
                         c.sample(frac=1,random_state=3),expected_full_counts=False)
    assert result == alt


def test_one_invalid_locality_is_counted_but_not_zero_filled():
    a,b,c=toy()
    mask=(b.species=="Toy species 00")&(b.locality_index==2)
    b.loc[mask,"status"]="undefined"
    b.loc[mask,"slope"]=np.nan
    b.loc[mask,"reason"]="insufficient_endpoint_disjoint_edges"
    a.loc[0,"locality_deletions_defined"] -= 1
    a.loc[0,"locality_deletions_undefined"] += 1
    result=summarize_inputs(a,b,c,expected_full_counts=False)
    assert result["panel"]["n_deletions_undefined"] == 1
    assert result["panel"]["n_deletions_defined"] == int(b.status.eq("defined").sum())


def test_no_unblinding_of_development_response_via_extra_candidate_species():
    a,b,c=toy()
    extra=pd.DataFrame([{"species":"Development-only species", "host_species":2,"native_units":99}])
    result=summarize_inputs(a,b,pd.concat([c,extra]),expected_full_counts=False)
    assert result["panel"]["n_all"] == len(a)
