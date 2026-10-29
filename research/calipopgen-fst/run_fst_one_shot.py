"""Authorized one-shot CaliPopGen Fst; never run before locked predictor/calibration gates."""
import argparse,collections,hashlib,io,json,math,re,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from inspect_xlsx_schema import exact_source_download
from climate_feature_gate import read_climate,features
from qualify_fst_blind import frozen_binomial,colnum,S
from model_fst_holdout import evaluate

KEEP={"CitationID","Pop1ScientificName","Pop2ScientificName","Pop1LatitudeDD",
  "Pop1LongitudeDD","Pop2LatitudeDD","Pop2LongitudeDD","GenDist",
  "GenDistMetric","MarkerType","GeneTarget","sepAnalyses","HabitatType","specialComparisonType","CoordError"}
def response_rows(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        strings=[]
        with z.open("xl/sharedStrings.xml") as f:
            for _,e in ET.iterparse(f,events=("end",)):
                if e.tag==S+"si":
                    strings.append("".join(v.text or "" for v in e.iter(S+"t")))
                    e.clear()
        with z.open("xl/worksheets/sheet1.xml") as f:
            columns=None
            for _,row in ET.iterparse(f,events=("end",)):
                if row.tag!=S+"row":continue
                cells={}
                for cell in row.findall(S+"c"):
                    i=colnum(cell.attrib.get("r",""))
                    if columns is not None and columns.get(i) not in KEEP:continue
                    v=cell.find(S+"v")
                    if v is None:cells[i]=None;continue
                    if columns is not None and columns.get(i)=="GenDist":
                        cells[i]=(cell.get("t"),v.text) # shared-string index only
                    elif cell.get("t")=="s":cells[i]=strings[int(v.text)]
                    elif cell.get("t")=="inlineStr":
                        cells[i]="".join(t.text or "" for t in cell.iter(S+"t"))
                    else:cells[i]=v.text
                if columns is None:
                    columns={i:str(t or "") for i,t in cells.items()}
                    if not KEEP.issubset(set(columns.values())):raise ValueError("schema missing fields")
                else:
                    out={name:cells.get(i) for i,name in columns.items() if name in KEEP}
                    if (str(out.get("GenDistMetric") or "").strip().lower()=="fst"
                        and str(out.get("HabitatType") or "").strip().lower()=="terrestrial"
                        and str(out.get("specialComparisonType") or "").strip().lower()=="regular"):
                        encoded=out.get("GenDist")
                        if isinstance(encoded,tuple):
                            t,value=encoded
                            out["GenDist"]=strings[int(value)] if t=="s" and value is not None else value
                    else:
                        out["GenDist"]=None
                    yield out
                row.clear()

def source_response_lookup(raw,qualified_groups):
    allowed={}
    for g in qualified_groups:
        key=(g["study_id"],g["species"],g["marker"],g["gene_target"],g["sep_analyses"])
        pairs={tuple(x) for x in g["climate_valid_pair_coordinates"]}
        allowed[key]=pairs
    result={k:{} for k in allowed}
    counts=collections.Counter()
    for r in response_rows(raw):
        counts["original_rows"]+=1
        metric=re.sub(r"[^a-z]","",str(r.get("GenDistMetric") or "").lower())
        if metric!="fst":continue
        if str(r.get("specialComparisonType") or "").strip().lower()!="regular":continue
        if str(r.get("HabitatType") or "").strip().lower()!="terrestrial":continue
        sp=frozen_binomial(r.get("Pop1ScientificName"))
        if sp!=frozen_binomial(r.get("Pop2ScientificName")):continue
        key=(str(r.get("CitationID") or "").strip(),sp,
          str(r.get("MarkerType") or "").strip(),
          str(r.get("GeneTarget") or "").strip(),
          str(r.get("sepAnalyses") or "").strip())
        if key not in allowed:continue
        try:
            a=(float(r["Pop1LatitudeDD"]),float(r["Pop1LongitudeDD"]))
            b=(float(r["Pop2LatitudeDD"]),float(r["Pop2LongitudeDD"]))
        except (TypeError,ValueError):continue
        aa,bb=sorted((a,b))
        pair=(*aa,*bb)
        if pair not in allowed[key]:continue
        if pair in result[key]:
            counts["duplicate_rejected_first_row_locked"]+=1
            continue
        # The first original source row determines the pair record.
        counts["eligible_source_pair_first_occurrence"]+=1
        try:
            val=float(r["GenDist"])
            if not math.isfinite(val) or not -1<=val<=1:raise ValueError()
            counts["numeric_Fst_values_accepted"]+=1
        except (TypeError,ValueError):
            val=None;counts["invalid_or_non_scalar_Fst"]+=1
        result[key][pair]=val
    if counts["original_rows"]!=14073:raise ValueError("source record count drift")
    return result,dict(counts)

def analyze(args):
    auth=json.loads(args.authorization.read_text())
    info=json.loads(args.information_gate.read_text())
    if auth.get("status")!="AUTHORIZE_FRESH_CALIPOPGEN_FST_ONCE":
        raise ValueError("Fst numerical response not authorized")
    if info["status"]!="PASS_PREDICTOR_INFORMATION_BEFORE_FST" or info["Fst_numeric_values_opened"]:
        raise ValueError("historical predictor-information prequalification absent")
    report,cohorts,by_site=read_climate(args.climate_dir)
    if (report["climate_value_matrix_sha256"]!=info["climate_input_sha256"]
        or report["climate_value_matrix_sha256"]!=auth["climate_matrix_sha256"]):
        raise ValueError("wrong immutable climate cohort")
    sd=np.asarray(info["training_current_bioclim_sd"],float)
    if sd.shape!=(4,) or (sd<=0).any():raise ValueError("bad climate normalization")
    # Fst numeric values opened ONLY after all of the above checks pass.
    lock=json.loads(args.source_lock.read_text())
    raw=exact_source_download(lock)
    responses,qc=source_response_lookup(raw,cohorts)
    model_groups=[]
    for g in cohorts:
        key=(g["study_id"],g["species"],g["marker"],g["gene_target"],g["sep_analyses"])
        v=responses[key]
        x=features(g,by_site,sd)
        y=np.asarray([v.get(tuple(p)) for p in g["climate_valid_pair_coordinates"]],object)
        valid=np.array([z is not None for z in y],bool)
        if valid.sum()<5:continue
        X=x[valid].astype(float);Y=np.asarray(y[valid],float)
        X-=X.mean(axis=0);Y-=Y.mean()
        model_groups.append({"species":g["species"],"role":g["role"],"x":X,"y":Y})
    tr={g["species"] for g in model_groups if g["role"]==0}
    te={g["species"] for g in model_groups if g["role"]==1}
    per_role=[sum(g["role"]==r for g in model_groups) for r in (0,1)]
    if len(tr|te)<30 or min(len(tr),len(te))<10 or min(per_role)<12 or len(model_groups)<35 or tr&te:
        out={"status":"NOT_EVALUABLE_POST_OPENING_SUPPORT","qc":qc,
             "species_by_role":[len(tr),len(te)],"groups_by_role":per_role,
             "genetic_Fst_numeric_opened":True,"positive_claim_authorized":False}
    else:
        out=evaluate(model_groups,bootstrap_draws=1999,seed=20261029)
        out["status"]="FROZEN_CALIPOPGEN_FST_EMPIRICAL_ENDPOINT"
        out["source_qc"]=qc
        out["opened_source"]="CaliPopGen Dataset 2 Figshare 15015927.v1, original MD5"
        out["genetic_Fst_numeric_opened"]=True
        out["closed_He_COI_TTF_responses_reopened"]=False
        out["interpretation_scope"]="Within-study Fst spatial contrasts in California, not exact genetic break positions or universal history causation."
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,sort_keys=True))
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--authorization",type=Path,required=True)
    p.add_argument("--information-gate",type=Path,required=True)
    p.add_argument("--source-lock",type=Path,required=True)
    p.add_argument("--climate-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    analyze(p.parse_args())
