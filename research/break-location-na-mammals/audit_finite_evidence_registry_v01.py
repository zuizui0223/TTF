"""Read-only validation of frozen TTF Quaternary information-vs-outcome synthesis.

Offline mode verifies arithmetic and claim boundaries from the immutable
cross-program registry. Online mode additionally checks *exact published Git
blobs* and their scalar contents from the original source branches. It never
fetches underlying genetic sequences, licensed occurrence coordinates, or
runs any empirical model. No pooled effect or meta-p value is estimated.
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import os
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
REGISTRY = ROOT / "cross_program_finite_evidence_v01.json"


def close(a, b, tag, tolerance=1e-11):
    if not math.isfinite(a) or not math.isfinite(b) or abs(a-b) > tolerance:
        raise AssertionError(f"{tag}: {a!r} != {b!r}")


def interval_has_zero(ci, label):
    if not len(ci) == 2 or not ci[0] < 0 < ci[1]:
        raise AssertionError(f"{label}: interval not bracketing zero: {ci}")


def get_nested(d, *keys):
    for key in keys:
        d = d[key]
    return d


def source_entries(registry):
    p = registry["response_blind_predictor"]
    return [
        ("TTFQ_C1", p["source"]),
        ("TTFQ_C2", p["geography_source"]),
        *[(entry["id"], entry["source"]) for entry in registry["genetic_endpoints"]],
        ("Stage23_jittered_climate", registry["predictor_positional_sensitivity"]["source"]),
    ]


def check_registry(r):
    if r["status"] != "POST_RESULT_DESCRIPTIVE_SYNTHESIS_ONLY_ALL_ORIGINAL_CONFIRMATORY_ENDPOINTS_CLOSED":
        raise AssertionError("Not a read-only synthesis")
    b = r["logic_boundary"]
    for key in ("predictor_nonredundancy_is_not_genetic_variance_explained",
                "source_datasets_overlap_not_fully_audited",
                "no_pooled_p_values_effect_estimates_or_species_denominators",
                "no_new_hypothesis_test_or_posthoc_subgroup",
                "no_model_reruns",
                "separate_TTFQ_scientific_freeze_unchanged"):
        if b.get(key) is not True:
            raise AssertionError(f"Mandatory interpretation firewall missing: {key}")
    p = r["response_blind_predictor"]
    if (p["n_species"],p["development_directed_dyads"],p["confirmatory_directed_dyads"]) != (653,15625,40602):
        raise AssertionError("TTF-Q population or dyad geometry changed")
    for role in ("development","confirmatory"):
        pre = p["historical_unique_fraction_before_geography"][role]
        post = p["historical_unique_fraction_after_geography"][role]
        fraction = p["geographic_adjustment_retention_fraction"][role]
        if not 0 < post <= pre <= 1:
            raise AssertionError("History unique predictor variance decomposition invalid")
        close(post/pre, fraction, f"{role} geographic retention fraction")
    if p["outcome_genetic_response_evaluated"] is not False:
        raise AssertionError("TTF-Q predictor qualification is not a genetic response test")
    ep={e["id"]:e for e in r["genetic_endpoints"]}
    if set(ep) != {"COI_transfer_211","GenDivRange_He_holdout",
                   "CaliPopGen_Fst_holdout","Jensen_mammal_break_location_holdout"}:
        raise AssertionError("Focal closed biological response family changed")
    for e in ep.values():
        if e["frozen"] is not True:
            raise AssertionError("Previously opened outcome was silently unfrozen")
        if e["n_species"] <= 0:
            raise AssertionError("Invalid biological species count")
        if e["training_species"] + e["heldout_species"] != e["n_species"]:
            raise AssertionError(f"Species split invalid for {e['id']}")
        source=e["source"]
        if not source.get("git_blob_sha") or len(source["git_blob_sha"])!=40:
            raise AssertionError("Source blob SHA missing")
    coi=ep["COI_transfer_211"]
    if coi["qualified_positive_transfer"] is not False:
        raise AssertionError("Closed COI transfer result changed")
    if not (0 < coi["profiled_private_p"] < 1 and 0 < coi["secondary_within_species_null_calibrated_p"] < 1):
        raise AssertionError("COI p-value outside [0,1]")
    if coi["secondary_raw_self_statistic"] >= 0:
        raise AssertionError("Self-structure raw statistic must not be claimed positive")
    for name in ("GenDivRange_He_holdout","CaliPopGen_Fst_holdout"):
        e=ep[name]
        close((e["baseline_mse"]-e["augmented_mse"])/e["baseline_mse"],
              e["gain"],name+" relative heldout MSE gain")
        interval_has_zero(e["uncertainty_interval_95"],name)
        if e["qualified_positive_additional_prediction"] is not False:
            raise AssertionError(name+" was not a positive endpoint")
    e=ep["Jensen_mammal_break_location_holdout"]
    if (e["case_site_records"],e["surveyed_background_records"],e["training_species"],e["heldout_species"]) != (10013,66687,44,36):
        raise AssertionError("Frozen Jensen response geography changed")
    for role in ("original","jitter50km"):
        a=e[role]
        close(a["augmented_auc"]-a["baseline_auc"],a["delta_auc"],f"Jensen {role} delta AUC")
        interval_has_zero(a["bootstrap_interval_9833"],f"Jensen {role}")
        if not (0 <= a["augmented_auc"] <= 1 and 0 <= a["baseline_auc"] <= 1):
            raise AssertionError("AUC outside valid range")
    if (e["qualified_positive_additional_prediction"] is not False or
        e["one_shot_action_disabled_after_success"] is not True):
        raise AssertionError("Jensen terminal first-shot protection failed")
    s=r["predictor_positional_sensitivity"]
    if s["is_not_genetic_effect"] is not True:
        raise AssertionError("Climatic predictor positional stability must not be a genetic effect")
    if s["n_species"]!=80 or s["case_sites"]!=10013:
        raise AssertionError("Stage23 source scope changed")
    for label,src in source_entries(r):
        if len(src.get("git_blob_sha",""))!=40:
            raise AssertionError("Unverified source blob: "+label)
    return ep


def fetch_json_source(source):
    repo=source["repository"]
    path=source["file"]
    ref=source["ref"]
    url=("https://api.github.com/repos/"+quote(repo,safe="/")+"/contents/"+
         quote(path,safe="/")+"?ref="+quote(ref,safe=""))
    headers={"Accept":"application/vnd.github+json",
             "X-GitHub-Api-Version":"2022-11-28",
             "User-Agent":"TTF-read-only-finite-evidence-validation"}
    if os.getenv("GH_TOKEN"):
        headers["Authorization"]="Bearer "+os.environ["GH_TOKEN"]
    req=Request(url,headers=headers)
    with urlopen(req,timeout=30) as response:
        record=json.load(response)
    if record.get("sha") != source["git_blob_sha"]:
        raise AssertionError(f"Referenced source SHA drift: {source['file']}")
    body=base64.b64decode(record["content"])
    return json.loads(body)


def cross_check(r, docs):
    p=r["response_blind_predictor"]
    c1,c2=docs["TTFQ_C1"],docs["TTFQ_C2"]
    if c1["design"]["admissible_species"]!=p["n_species"]:
        raise AssertionError("TTF-Q source species mismatch")
    for role in ("development","confirmatory"):
        close(c1[role]["total_unique_variance_fraction"],
              p["historical_unique_fraction_before_geography"][role], f"C1 {role}")
        close(c2[role]["C2_total_unique_variance_fraction"],
              p["historical_unique_fraction_after_geography"][role], f"C2 {role}")
        close(c2[role]["geography_retained_fraction"],
              p["geographic_adjustment_retention_fraction"][role], f"C2 retained {role}")
    ep={e["id"]:e for e in r["genetic_endpoints"]}
    coi=docs["COI_transfer_211"]["frozen_science"]
    close(coi["primary_place_beyond_ibd"]["statistic"],ep["COI_transfer_211"]["statistic"],"COI original")
    close(coi["primary_place_beyond_ibd"]["profiled_private_p_value"],
          ep["COI_transfer_211"]["profiled_private_p"],"COI p")
    close(coi["within_species_self"]["statistic"],ep["COI_transfer_211"]["secondary_raw_self_statistic"],"COI self")
    close(coi["within_species_self"]["upper_tail_p_value"],
          ep["COI_transfer_211"]["secondary_within_species_null_calibrated_p"],"COI self p")
    he=docs["GenDivRange_He_holdout"]["empirical"]
    e=ep["GenDivRange_He_holdout"]
    for a,b in (("baseline_pair_MSE","baseline_mse"),
                ("historical_added_pair_MSE","augmented_mse"),
                ("relative_MSE_gain","gain")):
        close(he[a],e[b], "He "+a)
    if he["species_cluster_bootstrap_95"] != e["uncertainty_interval_95"]:
        raise AssertionError("He frozen CI mismatch")
    fs=docs["CaliPopGen_Fst_holdout"]["empirical"]
    e=ep["CaliPopGen_Fst_holdout"]
    for a,b in (("baseline_group_equal_MSE","baseline_mse"),
                ("plus_history_group_equal_MSE","augmented_mse"),
                ("relative_heldout_MSE_gain","gain")):
        close(fs[a],e[b], "Fst "+a)
    if fs["species_cluster_bootstrap95"] != e["uncertainty_interval_95"]:
        raise AssertionError("Fst frozen CI mismatch")
    j=docs["Jensen_mammal_break_location_holdout"]
    e=ep["Jensen_mammal_break_location_holdout"]
    for src_name,out_name in (("primary_original_sites","original"),
                              ("prespecified_50km_jitter","jitter50km")):
        p1=j[src_name];s1=e[out_name]
        for a,b in (("baseline_mean_heldout_auc","baseline_auc"),
                    ("baseline_plus_historical_mean_heldout_auc","augmented_auc"),
                    ("historical_added_mean_delta_auc","delta_auc")):
            close(p1[a],s1[b], "Jensen "+src_name+" "+a)
        if p1["bootstrap9833_two_sided_ci"] != s1["bootstrap_interval_9833"]:
            raise AssertionError("Jensen CI mismatch: "+src_name)
    s=docs["Stage23_jittered_climate"]["history_predictor"]
    r23=r["predictor_positional_sensitivity"]
    close(s["pooled_site_spearman"],r23["pooled_original_vs_50km_spearman"],"Stage23 pooled rho")
    close(s["species_median_within_species_spearman"],r23["median_of_species_spearman"],"Stage23 cluster rho")
    return True


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--online",action="store_true",
                        help="Verify all source contents and immutable blobs through GitHub REST GET")
    args=parser.parse_args()
    r=json.loads(REGISTRY.read_text())
    check_registry(r)
    if args.online:
        docs={k:fetch_json_source(src) for k,src in source_entries(r)}
        cross_check(r,docs)
    print(json.dumps({"status":"PASS_FROZEN_CROSS_PROGRAM_SCALAR_AND_CLAIM_BOUNDARIES",
                      "sources_checked_live":len(source_entries(r)) if args.online else 0,
                      "frozen_biological_endpoints":len(r["genetic_endpoints"]),
                      "new_genetic_outcome_opened":False,
                      "pooled_inference_performed":False},sort_keys=True))


if __name__=="__main__":
    main()
