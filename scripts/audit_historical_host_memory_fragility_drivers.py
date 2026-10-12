#!/usr/bin/env python3
"""Owner-private inputs -> public aggregates describing *post-outcome* COI fragility.

No sequence or edge-genetic-distance input is accepted. Original confirmatory
one-shot p-value remains undefined. This is descriptive, not a new hypothesis
test, and does not redefine a species' genetic eligibility.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

PINS = {
    "species_stability": "b447b360edc3d73ae62fe4c668f2cadf9c13f94a0292d8aee5e96828c6724f23",
    "locality_deletions": "ba6f0f5f1ec9b9da4e697a0a271f0156ac2ad80d385da0a793bba8e200f97686",
    "frozen_candidates": "c36cbb2ba0cbf7d0222645a04538c78236cfda392dd0a3d11dd443f12347d35b",
    "primary_result": "5d994a76f7f20a778f0394f545c08a729182ce32298503c128d95c00da29f965",
}


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _rate(v, n):
    return float(v / n) if n else None


def _quartiles(s: pd.DataFrame, column: str, label: str) -> list[dict]:
    bins = pd.qcut(s[column], 4, labels=False, duplicates="raise")
    rows=[]
    for j in range(4):
        g=s.loc[bins == j]
        if len(g)==0:
            raise RuntimeError("empty original predictor quartile")
        flips=int(g.flips.sum());ops=int(g.defined.sum())
        rows.append({
            "quartile": j + 1,
            "n": int(len(g)),
            "median_abs_signal":float(g[column].median()),
            "n_with_any_reversal":int((g.flips > 0).sum()),
            "species_any_reversal_fraction":_rate(int((g.flips > 0).sum()),len(g)),
            "defined_locality_deletions":ops,
            "reversal_events":flips,
            "locality_deletion_reversal_fraction":_rate(flips,ops),
            "signal_type":label,
        })
    return rows


def _cv(s: pd.DataFrame) -> dict:
    # Each species is one outcome row, even though it has many deletion trials.
    # Held-out families prevent directly training on the test family's species.
    groups=s.family.to_numpy()
    labels=(s.flips.to_numpy() > 0).astype(int)
    if len(np.unique(groups)) < 5 or len(np.unique(labels)) < 2:
        raise RuntimeError("insufficient independent families or outcome variation")
    feats={
        "constant_prevalence": [],
        "absolute_standardized_beta_and_locality_count": ["log_abs_beta", "log_localities"],
        "absolute_semipartial_and_locality_count": ["log_abs_semi", "log_localities"],
        "semipartial_plus_recorded_host_count": ["log_abs_semi", "log_localities", "host_species"],
    }
    results={}
    for label, cols in feats.items():
        probabilities=np.empty(len(s),dtype=float)
        for train,test in GroupKFold(n_splits=5).split(s,labels,groups):
            if len(cols)==0:
                p=np.repeat(labels[train].mean(),len(test))
            else:
                if len(np.unique(labels[train]))<2:
                    raise RuntimeError("family held-out training fold has only one outcome class")
                mdl=make_pipeline(StandardScaler(),LogisticRegression(C=1.0,max_iter=1000,solver="lbfgs"))
                mdl.fit(s.iloc[train][cols],labels[train])
                p=mdl.predict_proba(s.iloc[test][cols])[:,1]
            probabilities[test]=p
        results[label]={
            "brier_score":float(brier_score_loss(labels,probabilities)),
            "roc_auc":float(roc_auc_score(labels,probabilities)),
            "n_species":int(len(s)),
        }
    return {"cv":"five_fold_family_held_out_GroupKFold", "no_hyperparameter_search":True,
            "species_weighting":"one species = one binary event (at least one flip)",
            "output_not_a_genetic_effect_p_value":True,"models":results}


def summarize_inputs(
    private_species:pd.DataFrame,
    private_deletions:pd.DataFrame,
    frozen_candidates:pd.DataFrame,
    *,expected_full_counts:bool=True,
)->dict:
    s=private_species.copy().sort_values("species").reset_index(drop=True)
    d=private_deletions.copy().sort_values(["species","locality_index"]).reset_index(drop=True)
    c=frozen_candidates.copy()
    for table,cols in [(s,["species"]),(d,["species","locality_index"]),(c,["species"])]:
        if table.duplicated(cols).any():
            raise RuntimeError("duplicate frozen species/locality identity")
    if set(d.species)!=set(s.species) or not set(s.species).issubset(c.species):
        raise RuntimeError("private locality deletions do not match qualified species")
    c=c.set_index("species")
    s["host_species"]=s.species.map(c.host_species)
    s["native_units"]=s.species.map(c.native_units)
    if s[["host_species","native_units"]].isna().any().any():
        raise RuntimeError("missing frozen recorded-host metadata")
    if ((s.host_species < 1) | (s.host_species > 3) | (s.native_units < 1)).any():
        raise RuntimeError("original frozen host-resource metadata outside eligibility")
    if not np.isfinite(s["unique_fraction"]).all() or not s.unique_fraction.between(0,1,inclusive="right").all():
        raise RuntimeError("invalid frozen host-memory unique information fraction")
    if not np.isfinite(s.unscaled_slope).all():
        raise RuntimeError("nonfinite original unscaled coefficient")
    counts=d.groupby("species",sort=False).agg(
        attempted=("status","size"),defined=("status",lambda x:(x=="defined").sum()),
        undefined=("status",lambda x:(x=="undefined").sum()),
    ).reindex(s.species)
    if set(d.status)!={"defined"} and not set(d.status).issubset({"defined","undefined"}):
        raise RuntimeError("illegal original locality-withholding status")
    if (counts.attempted.to_numpy()!=s.locality_deletions_attempted.to_numpy()).any() or (counts.attempted.to_numpy()!=s.localities.to_numpy()).any():
        raise RuntimeError("original locality-deletion opportunity counts drift")
    if (counts.defined.to_numpy()!=s.locality_deletions_defined.to_numpy()).any() or (counts.undefined.to_numpy()!=s.locality_deletions_undefined.to_numpy()).any():
        raise RuntimeError("undefined original locality deletions not retained")
    expected_reversals={}
    slopes=s.set_index("species").unscaled_slope
    for sp,g in d.groupby("species"):
        valid=g.loc[g.status=="defined"]
        if len(valid) and (not np.isfinite(valid.slope).all()):
            raise RuntimeError("nonfinite defined locality fit")
        if len(g.loc[g.status=="undefined"]) and g.loc[g.status=="undefined","slope"].notna().any():
            raise RuntimeError("undefined locality fit carries an invented coefficient")
        origin=float(slopes.loc[sp]);v=valid.slope.to_numpy()
        expected_reversals[sp] = int(np.sum(np.sign(v)==-np.sign(origin))) if abs(origin)>1e-12 else 0
    if any(expected_reversals[sp]!=int(r.locality_deletions_reversing_sign)
           for sp,r in s.set_index("species").iterrows()):
        raise RuntimeError("private locality sign reversals do not reproduce source counters")
    s["flips"]=s.locality_deletions_reversing_sign.astype(int)
    s["defined"]=s.locality_deletions_defined.astype(int)
    invariant=s.unscaled_slope.abs()<=1e-12
    if (s.loc[invariant,"standardized_beta"].notna().any() or
        s.loc[~invariant,"standardized_beta"].isna().any()):
        raise RuntimeError("four invariant COI responses must remain standardized-undefined")
    if (s.loc[invariant,"semipartial_correlation"].notna().any() or
        s.loc[~invariant,"semipartial_correlation"].isna().any()):
        raise RuntimeError("invalid semipartial identity for invariant COI species")
    nonzero=s.loc[~invariant].copy()
    expected_r=nonzero.standardized_beta.to_numpy()*np.sqrt(nonzero.unique_fraction.to_numpy())
    max_error=float(np.max(np.abs(expected_r-nonzero.semipartial_correlation.to_numpy())))
    if max_error>1e-10:
        raise RuntimeError("standardized beta / unique host fraction / semipartial identity mismatch")
    if np.any(np.abs(nonzero.semipartial_correlation)>1+1e-10):
        raise RuntimeError("semipartial magnitude exceeds correlation bound")
    if expected_full_counts:
        if (len(s),len(nonzero),int(invariant.sum()),len(d),int(counts.defined.sum()),int(counts.undefined.sum()),int(s.flips.sum())) != (259,255,4,6374,6373,1,420):
            raise RuntimeError("259-species post-outcome frozen locality-deletion census drift")
    nonzero["abs_beta"]=nonzero.standardized_beta.abs()
    nonzero["abs_semi"]=nonzero.semipartial_correlation.abs()
    nonzero["log_abs_beta"]=np.log10(nonzero.abs_beta)
    nonzero["log_abs_semi"]=np.log10(nonzero.abs_semi)
    nonzero["log_localities"]=np.log(nonzero.localities)
    if not np.isfinite(nonzero[["log_abs_beta","log_abs_semi","log_localities"]]).all().all():
        raise RuntimeError("undefined absolute effect log transforms")
    bins_beta=_quartiles(nonzero,"abs_beta","absolute_standardized_beta")
    bins_semi=_quartiles(nonzero,"abs_semi","absolute_host_history_semipartial_correlation")
    host_summary=[]
    for k,g in s.groupby("host_species"):
        valid=g[~(g.unscaled_slope.abs()<=1e-12)]
        host_summary.append({"number_of_frozen_recorded_host_species":int(k),
                             "total_species":int(len(g)),
                             "finite_standardized_species":int(len(valid)),
                             "positive_standardized_coefficients":int((valid.standardized_beta>0).sum()),
                             "negative_standardized_coefficients":int((valid.standardized_beta<0).sum()),
                             "species_any_locality_flip":int((valid.flips>0).sum()),
                             "median_host_history_unique_fraction":float(g.unique_fraction.median())})
    corr={}
    for name,col in [("abs_standardized_beta", "abs_beta"),
                     ("abs_semipartial", "abs_semi"),
                     ("host_history_unique_fraction", "unique_fraction"),
                     ("original_edges", "edges"),
                     ("max_ibd_prediction_leverage", "max_prediction_leverage")]:
        # Rank correlation is computed without ever requesting an inferential p-value.
        ranks_x=rankdata(nonzero[col].to_numpy(),method="average")
        ranks_y=rankdata((nonzero.flips/nonzero.defined).to_numpy(),method="average")
        corr[name]=float(np.corrcoef(ranks_x,ranks_y)[0,1])
    return {
        "schema":"ttf_hhm259_post_outcome_fragility_driver_aggregates_v0.1",
        "status":"EXPLORATORY_AFTER_EMPIRICAL_OUTCOME_EXPOSURE",
        "original_confirmatory_statistic_status":"NOT_EVALUABLE_EMPIRICAL_RESPONSE",
        "panel":{
            "n_all":int(len(s)),"n_nonzero_defined":int(len(nonzero)),
            "n_baseline_invariant":int(invariant.sum()),
            "n_with_any_sign_reversal":int((nonzero.flips>0).sum()),
            "n_deletions_attempted":int(len(d)),
            "n_deletions_defined":int(counts.defined.sum()),
            "n_deletions_undefined":int(counts.undefined.sum()),
            "n_defined_deletions_among_nonzero_species":int(nonzero.defined.sum()),
            "n_reversal_events":int(nonzero.flips.sum()),
            "fraction_species_with_any_reversal":_rate(int((nonzero.flips>0).sum()),len(nonzero)),
            "fraction_nonzero_species_deletions_reversing":_rate(int(nonzero.flips.sum()),int(nonzero.defined.sum())),
        },
        "absolute_standardized_beta_quartiles":bins_beta,
        "absolute_semipartial_quartiles":bins_semi,
        "spearman_descriptive_against_species_locality_flip_fraction":corr,
        "recorded_host_metadata_by_count":host_summary,
        "family_held_out_internal_prediction":_cv(nonzero),
        "identity_verification":{"r_semi_equals_beta_times_sqrt_unique_fraction_max_absolute_error":max_error},
        "limitations":{
            "the_absolute_coefficient_and_semi_partial_are_calculated_from_the_same_opened_COI_response":True,
            "CV_score_is_not_out_of_sample_validation_of_an_ecological_host_memory_mechanism":True,
            "host_species_is_recorded_subset_not_total_diet_breadth":True,
            "no_p_values_or_confidence_intervals":True,
            "original_primary_estimate_still_undefined":True,
            "no_new_confirmatory_panel":True,
            "no_raw_nucleotide_sequences_or_edge_p_distances_read_by_this_script":True,
            "no_public_species_names_or_coefficients":True,
        },
        "new_p_values":0,
    }


def main():
    ap=argparse.ArgumentParser()
    for key in ("species_stability","locality_deletions","frozen_candidates","primary_result"):
        ap.add_argument("--"+key.replace("_","-"),type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    for key,pin in PINS.items():
        if sha256(getattr(args,key))!=pin:
            raise RuntimeError("private/source hash mismatch for "+key)
    terminal=json.loads(args.primary_result.read_text())
    if (terminal.get("status")!="NOT_EVALUABLE_EMPIRICAL_RESPONSE" or
        terminal.get("n_species_frozen")!=259 or terminal.get("undefined_species")!=4 or
        terminal.get("monte_carlo_p") is not None):
        raise RuntimeError("original one-shot empirical terminal result was changed")
    s=pd.read_csv(args.species_stability)
    d=pd.read_csv(args.locality_deletions)
    c=pd.read_csv(args.frozen_candidates)
    result=summarize_inputs(s,d,c)
    result["source_sha256"]=dict(PINS)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"n_nonzero":result["panel"]["n_nonzero_defined"],
                      "n_flip_species":result["panel"]["n_with_any_sign_reversal"],
                      "n_flip_events":result["panel"]["n_reversal_events"]},sort_keys=True))


if __name__=="__main__":
    main()
