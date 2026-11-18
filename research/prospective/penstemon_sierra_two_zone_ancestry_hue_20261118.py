"""Source-pinned exploratory Penstemon Sierra cross-zone ancestry-to-hue transfer.

Design committed BEFORE full original data opening:
PENSTEMON_SIERRA_ANCESTRY_HUE_ZONE_HOLDOUT_PLAN_20261118.md
Original author, published paper, no new sample collection. TTF-Q original
frozen scientific claims unchanged; 2 zones != 2 species replication.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import math
import statistics
import urllib.request
from pathlib import Path

AUTHOR_REPO="benstemon/sierra-hybrid-zone"
SOURCE_COMMIT="0ba4842824bcd9ddc824eedea303ad1c176fc82b"
SOURCE_PATH="1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv"
SOURCE_BLOB="7cd5abd0f2a5ad6de1563d3ddd2c0d857b037de8"
ZONES=("hybrids_virginialakes","hybrids_gemlakes")
MIN_HYBRIDS_PER_ZONE=20
EXPECTED_SOURCE_ROWS=421


def git_blob_sha(data:bytes)->str:
    return hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()


def fnum(value:str,field:str)->float|None:
    value=(value or "").strip()
    if not value or value.lower() in {"na","nan","none","null"}:return None
    val=float(value)
    if not math.isfinite(val):raise ValueError(f"nonfinite original source {field}")
    return val


def parse_source(raw:bytes)->dict:
    if git_blob_sha(raw)!=SOURCE_BLOB:
        raise ValueError("original Git blob SHA mismatch; data drift must stop")
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    fields=reader.fieldnames or []
    expected={"newID","pop","class","subpopulation","modified_hue","K1","K2","region"}
    if len(fields)!=len(set(fields)) or not expected.issubset(fields):
        raise ValueError("missing or duplicated original input headers")
    rows=[]
    ids=set()
    counts=collections.Counter()
    missing=collections.Counter()
    for row in reader:
        sid=(row.get("newID") or "").strip()
        if not sid or sid in ids:raise ValueError("missing/duplicate original individual newID")
        ids.add(sid)
        label=(row.get("class") or "").strip()
        zone=(row.get("subpopulation") or "").strip()
        counts[(label,zone)]+=1
        if label!="hybrid" or zone not in ZONES:continue
        x=fnum(row.get("K1"),"K1")
        y=fnum(row.get("modified_hue"),"modified_hue")
        if x is None:missing[(zone,"K1")]+=1
        if y is None:missing[(zone,"modified_hue")]+=1
        if x is None or y is None:continue
        if not 0<=x<=1:
            raise ValueError("original admixture ancestry K1 outside [0,1]")
        rows.append({"newID":sid,"zone":zone,"x":x,"y":y})
    if len(ids)!=EXPECTED_SOURCE_ROWS:
        raise ValueError(f"original source row count changed: {len(ids)}")
    eligible=collections.Counter(r["zone"] for r in rows)
    if any(eligible[z]<MIN_HYBRIDS_PER_ZONE for z in ZONES):
        raise ValueError("insufficient hybrids in one of the preregistered two zones")
    return {"rows":rows,"source_count":len(ids),
      "complete_case_hybrid_counts":dict(sorted(eligible.items())),
      "original_class_by_subpopulation_counts":{
          f"{k[0]}:{k[1]}":v for k,v in sorted(counts.items())},
      "eligible_hybrid_missingness":{
          f"{k[0]}:{k[1]}":v for k,v in sorted(missing.items())},
      "excluded_nonhybrids_and_unassigned":len(ids)-sum(v for (cl,zone),v in counts.items()
                                  if cl=="hybrid" and zone in ZONES),
      "phenotypes_and_genome_ancestry_same_newID":True,
    }


def ols_fit(xs:list[float],ys:list[float])->tuple[float,float]:
    if len(xs)!=len(ys) or len(xs)<MIN_HYBRIDS_PER_ZONE:
        raise ValueError("training sample count or grain failure")
    mx=sum(xs)/len(xs)
    my=sum(ys)/len(ys)
    denom=sum((x-mx)**2 for x in xs)
    if denom<=1e-12:raise ValueError("no ancestry K1 variation in training zone")
    slope=sum((x-mx)*(y-my) for x,y in zip(xs,ys))/denom
    return my-slope*mx,slope


def pearson(xs:list[float],ys:list[float])->float|None:
    if len(xs)!=len(ys) or len(xs)<2:return None
    x0=sum(xs)/len(xs)
    y0=sum(ys)/len(ys)
    denom=sum((x-x0)**2 for x in xs)*sum((y-y0)**2 for y in ys)
    if denom<=0:return None
    return sum((x-x0)*(y-y0) for x,y in zip(xs,ys))/math.sqrt(denom)


def one_fold(training:list[dict],test:list[dict],tr_name:str,te_name:str)->dict:
    fit_a,fit_b=ols_fit([r["x"] for r in training],[r["y"] for r in training])
    training_mean=sum(r["y"] for r in training)/len(training)
    actual=[r["y"] for r in test]
    pred=[fit_a+fit_b*r["x"] for r in test]
    sse=sum((y-p)**2 for y,p in zip(actual,pred))
    b_sse=sum((y-training_mean)**2 for y in actual)
    if b_sse<=0:
        raise ValueError("zero denominator in heldout training-mean benchmark")
    return {
      "training_zone":tr_name,"heldout_zone":te_name,
      "train_n_individuals":len(training),"heldout_n_individuals":len(test),
      "fitted_training_intercept":fit_a,"fitted_training_ancestry_slope":fit_b,
      "training_hue_mean_benchmark":training_mean,
      "heldout_rmse":math.sqrt(sse/len(test)),
      "heldout_mae":sum(abs(y-p) for y,p in zip(actual,pred))/len(test),
      "heldout_training_mean_baseline_rmse":math.sqrt(b_sse/len(test)),
      "heldout_training_mean_baseline_mae":sum(abs(y-training_mean) for y in actual)/len(test),
      "heldout_sse":sse,"heldout_baseline_sse":b_sse,
      "heldout_sse_skill_vs_training_mean":1-sse/b_sse,
      "heldout_pearson_r_pred_vs_observed":pearson(pred,actual),
    }


def analyze(raw:bytes)->dict:
    source=parse_source(raw)
    rows=source.pop("rows")
    byzone={z:[r for r in rows if r["zone"]==z] for z in ZONES}
    folds=[
       one_fold(byzone[tr],byzone[te],tr,te)
       for tr,te in [(ZONES[0],ZONES[1]),(ZONES[1],ZONES[0])]
    ]
    return {
      "status":"EXPLORATORY_SIERRA_TWO_ZONE_ANCESTRY_HUE_TRANSFER",
      "published_study":"Stone et al. 2026 Evolution doi:10.1093/evolut/qpag199",
      "original_source_repo":AUTHOR_REPO,
      "original_source_commit":SOURCE_COMMIT,
      "original_source_path":SOURCE_PATH,
      "original_source_blob":SOURCE_BLOB,
      "source_audit":source,
      "heldout_unit":"entire hybrid-zone group",
      "distinct_zone_count":2,
      "response":"original authored modified_hue",
      "predictor":"original pooled ADMIXTURE K1",
      "baseline":"training-zone mean hue only",
      "model":"one training-zone univariate OLS with intercept and no tuning",
      "folds":folds,
      "source_only_descriptive_mean_skill":sum(x["heldout_sse_skill_vs_training_mean"] for x in folds)/2,
      "limitations":[
         "one species pair across only two source zones",
         "source pooled ancestry K1 uses source test zone genotypes transductively",
         "no independent neutral-history term",
         "not causal locus or reproductive barrier inference",
         "source author paper already reported hue versus genomic barrier decoupling",
         "pre-analysis source peek included four rows, not fully blinded",
      ],
      "no_pooled_individual_as_independent_zone_pvalue":True,
      "not_TTF_Q_confirmatory":True,
      "original_TTF_main_frozen_science_unmodified":True,
    }


def get_original_bytes()->bytes:
    url=f"https://raw.githubusercontent.com/{AUTHOR_REPO}/{SOURCE_COMMIT}/{SOURCE_PATH}"
    with urllib.request.urlopen(url,timeout=40) as response:
        return response.read(700_001)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    data=analyze(a.input.read_bytes() if a.input else get_original_bytes())
    out=json.dumps(data,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if a.output:a.output.write_text(out)
    else:print(out)


if __name__=="__main__":main()
