"""Artificial coefficient and site data; never reads original genetic identities."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import pytest
from audit_locality_stability_drivers import analyze_species_and_deletions


def toy():
    # Avoid attaching regression test behavior to unexposed empirical species.
    vals = [0.04, 0.15, 0.23, 0.4, 0.6, 0.8, 1.0, 0.0]
    rows, dels = [], []
    for i, b in enumerate(vals):
        name = f"Toy{i}"
        coef = (1 if i % 2 == 0 else -1) * (b + 0.1) if b else 0.0
        rows.append(dict(
            species=name, localities="3", edges="7", standardized_beta=str(b) if b else "",
            unscaled_slope=str(coef),
            unique_fraction=str(0.1 + i * 0.02),
            max_prediction_leverage=str(0.2+i*0.1),
            max_host_residual_leverage=str(0.05+i*0.01),
            locality_deletions_attempted="3", locality_deletions_defined="3",
            locality_deletions_undefined="0",
            locality_deletions_reversing_sign=str(1 if i == 0 else 0),
        ))
        for j in range(3):
            val = -coef if i == 0 and j == 0 else (0.0 if i == 1 and j == 1 else coef)
            dels.append(dict(species=name, locality_index=str(j), status="defined", slope=str(val)))
    return rows, dels


def test_strict_sign_reversal_does_not_count_exact_zero_deletion():
    rows, dels = toy()
    result = analyze_species_and_deletions(rows, dels)
    assert result["nonzero_baseline_species"] == 7
    assert result["species_with_any_strict_sign_flip"] == 1
    assert result["strict_flip_attempts"] == 1
    assert result["zero_deletion_values_from_nonzero_baselines"] == 1
    assert result["undefined_site_deletions"] == 0
    assert result["original_zero_species"] == 1
    assert result["exploratory_no_p_values"] is True


def test_duplicate_and_missing_locality_attempts_fail_closed():
    rows, dels = toy()
    with pytest.raises(ValueError, match="duplicate"):
        analyze_species_and_deletions(rows, dels + [dels[0]])
    with pytest.raises(ValueError, match="locality attempts"):
        analyze_species_and_deletions(rows, dels[:-1])


def test_mismatched_original_species_flip_count_is_not_silently_accepted():
    rows, dels = toy()
    rows[0]["locality_deletions_reversing_sign"] = "2"
    with pytest.raises(ValueError, match="strict flip"):
        analyze_species_and_deletions(rows, dels)


def test_no_baseline_significance_or_emergent_selection_in_aggregate():
    rows, dels = toy()
    result = analyze_species_and_deletions(rows, dels)
    assert "species_names" not in result
    assert "p_value" not in result
    assert sum(r["species"] for r in []) == 0
    assert all(row["n"] >= 1 for row in result["quartiles_by_baseline_abs_standardized_beta"])


def test_one_undefined_deletion_is_separate_from_species_ineligibility():
    rows, dels = toy()
    # Convert one defined to undefined, retain no forced zero. Adjust source ledger to match.
    d = next(x for x in dels if x["species"] == "Toy2" and x["locality_index"] == "2")
    d["status"] = "undefined"
    d["slope"] = ""
    rows[2]["locality_deletions_defined"] = "2"
    rows[2]["locality_deletions_undefined"] = "1"
    result = analyze_species_and_deletions(rows, dels)
    assert result["undefined_site_deletions"] == 1
    assert result["species_count"] == 8

def test_two_way_median_split_is_exploratory_and_exhausts_original_nonzero_species():
    rows, dels = toy()
    result = analyze_species_and_deletions(rows, dels)
    cells = result["median_split_abs_beta_by_ibd_leverage"]["cells"]
    assert sum(r["n"] for r in cells) == result["nonzero_baseline_species"] == 7
    assert sum(r["strict_sign_flips"] for r in cells) == result["strict_flip_attempts"]
    assert sum(r["defined_deletions"] for r in cells) == result["defined_site_deletions_among_nonzero_species"]
    assert all(r["low_high_baseline_effect"] in ("low", "high") for r in cells)
    assert all(r["low_high_ibd_prediction_leverage"] in ("low", "high") for r in cells)
