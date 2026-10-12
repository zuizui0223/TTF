#!/usr/bin/env python3
"""Post-outcome, read-only locality-deletion influence drivers for frozen TTF.

Only consumes already-opened derived coefficient and deletion sensitivity CSVs.
NO sequence identities, pairwise genetic distance, novel response tests, or
empirical p-values are read or computed. The 259-species primary result remains
NOT_EVALUABLE_EMPIRICAL_RESPONSE and is never modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

STABILITY_SHA256 = "b447b360edc3d73ae62fe4c668f2cadf9c13f94a0292d8aee5e96828c6724f23"
DELETIONS_SHA256 = "ba6f0f5f1ec9b9da4e697a0a271f0156ac2ad80d385da0a793bba8e200f97686"

DRIVERS = {
    "baseline_abs_standardized_beta": "abs_standardized",
    "frozen_localities": "localities",
    "host_history_unique_fraction": "unique_fraction",
    "max_host_residual_leverage": "max_host_residual_leverage",
    "max_ibd_prediction_leverage": "max_prediction_leverage",
}


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda: f.read(1 << 20), b""):
            h.update(data)
    return h.hexdigest()


def finite_float(obj, field):
    try:
        v = float(obj[field])
    except (ValueError, TypeError, KeyError) as exc:
        raise ValueError(f"cannot parse finite {field}") from exc
    if not math.isfinite(v):
        raise ValueError(f"non-finite {field}")
    return v


def ties_average_ranks(values):
    order = sorted(range(len(values)), key=lambda i: (values[i], i))
    result = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j) / 2
        for k in range(i, j):
            result[order[k]] = avg
        i = j
    return result


def descriptive_spearman(x, y):
    if len(x) != len(y) or len(x) < 3:
        return None
    a = ties_average_ranks(x)
    b = ties_average_ranks(y)
    am = statistics.mean(a)
    bm = statistics.mean(b)
    c = sum((u-am) * (v-bm) for u, v in zip(a, b))
    aa = sum((u-am)**2 for u in a)
    bb = sum((v-bm)**2 for v in b)
    return c/math.sqrt(aa*bb) if aa*bb > 0 else None


def quartiles_of_species(species, metric):
    ordered = sorted(species, key=lambda row: (row[metric], row["name"]))
    n = len(ordered)
    bounds = [round(k*n/4) for k in range(5)]
    groups = []
    for q in range(4):
        subset = ordered[bounds[q]:bounds[q+1]]
        if not subset:
            continue
        flips = sum(x["flip_count"] for x in subset)
        attempts = sum(x["defined_deletions"] for x in subset)
        groups.append({
            "quartile": q+1,
            "n": len(subset),
            "median_original_metric": statistics.median(row[metric] for row in subset),
            "species_with_any_flip": sum(row["flip_count"] > 0 for row in subset),
            "fraction_species_with_any_flip": sum(row["flip_count"] > 0 for row in subset)/len(subset),
            "strict_sign_flips": flips,
            "defined_site_deletions": attempts,
            "fraction_deletions_flipping": flips/attempts if attempts else None,
        })
    return groups


def descriptive_median_cross_tab(nonzero):
    """Unadjusted 2x2 description; neither a causal adjustment nor a test."""
    m_beta = statistics.median(r["abs_standardized"] for r in nonzero)
    m_ibd = statistics.median(r["max_prediction_leverage"] for r in nonzero)
    bins = defaultdict(list)
    for r in nonzero:
        effect = "low" if r["abs_standardized"] < m_beta else "high"
        ibd = "low" if r["max_prediction_leverage"] < m_ibd else "high"
        bins[(effect, ibd)].append(r)
    out = []
    for effect in ("low", "high"):
        for ibd in ("low", "high"):
            vals=bins.get((effect,ibd),[])
            if not vals:
                continue
            flips=sum(r["flip_count"] for r in vals)
            n_defined=sum(r["defined_deletions"] for r in vals)
            out.append({
                "low_high_baseline_effect": effect,
                "low_high_ibd_prediction_leverage": ibd,
                "n": len(vals),
                "species_with_any_flip": sum(r["flip_count"]>0 for r in vals),
                "fraction_species_with_any_flip": sum(r["flip_count"]>0 for r in vals)/len(vals),
                "strict_sign_flips": flips,
                "defined_deletions": n_defined,
                "fraction_deletions_flipping": flips/n_defined if n_defined else None,
                "median_effect_absolute": statistics.median(r["abs_standardized"] for r in vals),
                "median_ibd_prediction_leverage": statistics.median(r["max_prediction_leverage"] for r in vals),
            })
    return {
        "baseline_abs_standardized_beta_median":m_beta,
        "max_ibd_prediction_leverage_median":m_ibd,
        "split_rule":"strictly below median = low; equal/above = high; ties not broken by outcome",
        "cells":out,
        "claim_scope":"Descriptive joint stratification only; groups are post-outcome, not adjusted effects, and no p-values computed."
    }


def analyze_species_and_deletions(species_rows, deletion_rows):
    if not species_rows:
        raise ValueError("no original species")
    by_species = {}
    for row in species_rows:
        name = str(row.get("species", "")).strip()
        if not name or name in by_species:
            raise ValueError("duplicate or missing species in stability table")
        by_species[name] = row
    deletions = defaultdict(dict)
    for row in deletion_rows:
        sp = row.get("species", "")
        if sp not in by_species:
            raise ValueError("unknown species in locality deletion")
        j = int(row.get("locality_index", -1))
        if j in deletions[sp]:
            raise ValueError("duplicate locality deletion")
        deletions[sp][j] = row

    processed = []
    defined_global = []
    baseline_zero_count = 0
    zero_deletions_from_nonzero = 0
    undefined_total = 0
    for name, row in by_species.items():
        n_sites = int(row["localities"])
        attempts = deletions.get(name, {})
        if n_sites < 2 or set(attempts) != set(range(n_sites)):
            raise ValueError(f"frozen locality attempts mismatch: {name}")
        if int(row["locality_deletions_attempted"]) != n_sites:
            raise ValueError("attempted-site count mismatch")
        orig = finite_float(row, "unscaled_slope")
        beta_txt = str(row.get("standardized_beta", "")).strip()
        if orig == 0:
            baseline_zero_count += 1
            if beta_txt:
                raise ValueError("original zero-response species has a standardized beta")
        elif not beta_txt:
            raise ValueError("nonzero response has no standardized beta")
        b_abs = abs(float(beta_txt)) if beta_txt else 0.0
        vals = []
        undefined = 0
        flips = 0
        for j in range(n_sites):
            r = attempts[j]
            status = r["status"]
            if status == "defined":
                val = finite_float(r, "slope")
                if orig == 0 and abs(val) > 1e-13:
                    raise ValueError("original zero species became nonzero without definition")
                if orig != 0 and val == 0:
                    zero_deletions_from_nonzero += 1
                if orig * val < 0:  # A zero deletion is NOT an opposite-sign coefficient.
                    flips += 1
                vals.append(val)
                defined_global.append((orig, val))
            elif status == "undefined":
                undefined += 1
            else:
                raise ValueError("unknown locality deletion status")
        if (len(vals) != int(row["locality_deletions_defined"])
                or undefined != int(row["locality_deletions_undefined"])):
            raise ValueError(f"defined/undefined locality count mismatch: {name}")
        if flips != int(row["locality_deletions_reversing_sign"]):
            raise ValueError(f"strict flip count does not reproduce frozen ledger: {name}")
        undefined_total += undefined
        processed.append({
            "name": name,
            "baseline": orig,
            "abs_standardized": b_abs,
            "localities": n_sites,
            "edges": int(row["edges"]),
            "unique_fraction": finite_float(row, "unique_fraction"),
            "max_host_residual_leverage": finite_float(row, "max_host_residual_leverage"),
            "max_prediction_leverage": finite_float(row, "max_prediction_leverage"),
            "flip_count": flips,
            "defined_deletions": len(vals),
            "undefined_deletions": undefined,
            "flip_fraction": flips/len(vals) if vals else None,
        })

    nonzero = [r for r in processed if r["baseline"] != 0]
    if not nonzero:
        raise ValueError("no defined nonzero baseline coefficients")
    flip_total = sum(x["flip_count"] for x in nonzero)
    n_defined = sum(x["defined_deletions"] for x in nonzero)
    result = {
        "schema": "ttf_hhm259_post_outcome_locality_sign_stability_drivers_v0.1",
        "status": "EXPLORATORY_POST_OUTCOME_NO_CONFIRMATORY_REPLACEMENT",
        "original_confirmatory_status": "NOT_EVALUABLE_EMPIRICAL_RESPONSE",
        "species_count": len(processed),
        "original_zero_species": baseline_zero_count,
        "nonzero_baseline_species": len(nonzero),
        "attempted_site_deletions": sum(x["localities"] for x in processed),
        "defined_site_deletions": sum(x["defined_deletions"] for x in processed),
        "undefined_site_deletions": undefined_total,
        "species_with_any_strict_sign_flip": sum(r["flip_count"] > 0 for r in nonzero),
        "strict_flip_attempts": flip_total,
        "defined_site_deletions_among_nonzero_species": n_defined,
        "fraction_defined_nonzero_site_deletions_flipping": flip_total/n_defined,
        "median_per_species_flip_fraction_nonzero": statistics.median(r["flip_fraction"] for r in nonzero),
        "zero_deletion_values_from_nonzero_baselines": zero_deletions_from_nonzero,
        "original_unscaled_259_equal_species_mean": statistics.mean(r["baseline"] for r in processed),
        "single_site_replacement_mean_min": min(
            sum(r["baseline"] for r in processed)/len(processed) +
            (modified-baseline)/len(processed)
            for baseline, modified in defined_global
        ),
        "single_site_replacement_mean_max": max(
            sum(r["baseline"] for r in processed)/len(processed) +
            (modified-baseline)/len(processed)
            for baseline, modified in defined_global
        ),
        "median_split_abs_beta_by_ibd_leverage": descriptive_median_cross_tab(nonzero),
        "quartiles_by_baseline_abs_standardized_beta": quartiles_of_species(nonzero, "abs_standardized"),
        "quartiles_by_frozen_localities": quartiles_of_species(nonzero, "localities"),
        "quartiles_by_host_history_unique_fraction": quartiles_of_species(nonzero, "unique_fraction"),
        "quartiles_by_max_host_residual_leverage": quartiles_of_species(nonzero, "max_host_residual_leverage"),
        "quartiles_by_max_ibd_prediction_leverage": quartiles_of_species(nonzero, "max_prediction_leverage"),
        "descriptive_rank_correlations_of_flip_fraction": {
            key: descriptive_spearman([r[field] for r in nonzero],
                                      [r["flip_fraction"] for r in nonzero])
            for key, field in DRIVERS.items()
        },
        "interpretation_limits": {
            "exploratory_no_p_values": True,
            "no_new_genetic_hypothesis_tests": True,
            "no_species_or_edge_drop": True,
            "no_alternative_confirmatory_statistic": True,
            "no_development_sequence_opening": True,
            "no_nucleotide_identity_or_per_edge_distance_used_in_this_derived_audit": True,
            "quartiles_not_prospective_decision_rules": True,
            "any_flip_frequency_depends_on_number_of_localities": True,
            "baseline_near_zero_is_mathematically_easier_to_reverse": True,
            "single_site_replacement_is_not_an_independent_holdout_test": True,
        },
    }
    result["exploratory_no_p_values"] = True
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--species-stability", required=True, type=Path)
    p.add_argument("--site-deletions", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    for path, digest in ((args.species_stability, STABILITY_SHA256),
                         (args.site_deletions, DELETIONS_SHA256)):
        if file_sha256(path) != digest:
            raise RuntimeError(f"immutable post-outcome derived input SHA-256 mismatch: {path.name}")
    with args.species_stability.open(newline="", encoding="utf-8") as f:
        s = list(csv.DictReader(f))
    with args.site_deletions.open(newline="", encoding="utf-8") as f:
        d = list(csv.DictReader(f))
    result = analyze_species_and_deletions(s, d)
    if (result["species_count"],result["nonzero_baseline_species"],
        result["attempted_site_deletions"],result["defined_site_deletions"],
        result["undefined_site_deletions"],result["strict_flip_attempts"]) != (259,255,6374,6373,1,420):
        raise RuntimeError("historical frozen lineage or strict sign-reversal count drift")
    result["input_sha256"] = {"species_stability": STABILITY_SHA256,
                              "site_deletions": DELETIONS_SHA256}
    result["origin"] = "The two SHA-pinned PRIVATE per-species/leave-one-locality-out CSVs saved in TTF_hhm259_locality_stability_20261012.zip; only aggregate statistics are persisted."
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2)+"\n")
    print(json.dumps({k:result[k] for k in (
        "species_count","species_with_any_strict_sign_flip","strict_flip_attempts",
        "single_site_replacement_mean_min","single_site_replacement_mean_max")},
        sort_keys=True))


if __name__ == "__main__":
    main()
