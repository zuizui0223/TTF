"""Response-blind source-only historical climate innovation after current climate & geography.

INPUTS: genuine WorldClim v1.4 *matched* 1960–1990 BIO1/BIO12 cells and
~22ka MPI-ESM-P BIO1/BIO12 at 95 deduplicated Albano grid cells.
No plant fitness, genetics, traits, 211-species TTF or 44/50 Bontrager outcomes.

EXPLORATORY GCM: MPI-ESM-P is a *sensitivity* model in prior frozen contract;
the prespecified PRIMARY CCSM4 is not yet available, so this alone cannot
pass the full historical predictor gate.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

HERE=Path(__file__).resolve().parent
CURRENT=HERE/"source_registry/albano2025_WC14_modern_BIO1_BIO12_unlabeled_cells_v01.json"
LGM=HERE/"source_registry/albano2025_MPI_LGM_BIO1_BIO12_unlabeled_cells_v01.json"
RES=2.5/60

def cell_lon_lat(row,col):
    return -180+(col+0.5)*RES,90-(row+0.5)*RES

def prepare(current,lgm):
    valid_lgm_schemas={
       "ttf_Albano2025_MPI_ESM_P_source_only_2_5m_LGM_cells_v01":"MPI-ESM-P",
       "ttf_Albano2025_CCSM4_PRIMARY_source_only_2_5m_LGM_cells_v01":"CCSM4",
    }
    if (current.get("schema")!="ttf_Albano2025_WorldClim14_modern_1960_1990_actual_cells_v01"
        or lgm.get("schema") not in valid_lgm_schemas):
        raise ValueError("Original matched-modern versus frozen LGM receipt schemas missing")
    if current.get("baseline_period")!="1960-1990" or lgm.get("unlabeled_grid_cell_count")!=95:
        raise ValueError("Original climate period/cell-family mismatch")
    if current["phenotype_fitness_genetic_responses_used"] or lgm["fitness_genetics_adaptation_responses_opened"]:
        raise ValueError("Original biological response unexpectedly opened")
    old={(x["row"],x["col"]):x for x in current["cells"]}
    past={(x["row"],x["col"]):x for x in lgm["cell_profiles"]}
    if len(old)!=95 or set(old)!=set(past):raise ValueError("Modern vs LGM cell join incomplete or duplicate")
    dataset=[]
    for key in sorted(old):
        lon,lat=cell_lon_lat(*key)
        continent="NAM" if -125<lon<-50 else "EUR" if -20<lon<45 else None
        if continent is None:raise ValueError("Source pixel outside two original continents")
        p=old[key];q=past[key]
        if any(p[f"BIO{z}_raw"] is None or q[f"BIO{z}_raw"] is None
               for z in (1,12)):raise ValueError("Missing original matched climate values")
        dataset.append({"row":key[0],"col":key[1],"lat":lat,"lon":lon,"origin":continent,
                        "modern_bio1":float(p["BIO1_raw"])/10,
                        "modern_bio12":float(p["BIO12_raw"]),
                        "delta_bio1":(float(q["BIO1_raw"])-float(p["BIO1_raw"]))/10,
                        "delta_bio12":float(q["BIO12_raw"])-float(p["BIO12_raw"])})
    if (sum(x["origin"]=="NAM" for x in dataset),
        sum(x["origin"]=="EUR" for x in dataset))!=(47,48):
        raise ValueError("Expected 47 NAM and 48 distinct EUR 2.5m cells")
    for first,cont in [(0,"NAM"),(5,"EUR")]:
        items=sorted((x for x in dataset if x["origin"]==cont),
                     key=lambda x:(x["lat"],x["row"],x["col"]))
        base,extra=divmod(len(items),5)
        left=0
        for i in range(5):
            num=base+(i<extra)
            for d in items[left:left+num]:d["fold"]=first+i
            left+=num
    if {x["fold"] for x in dataset}!=set(range(10)):raise ValueError("Source spatial fold setup failed")
    return dataset

def rank_and_prediction(data,fold,target):
    import numpy as np
    rows=np.array([[float(x["origin"]=="EUR"),x["lat"],x["lon"],
                     x["modern_bio1"],x["modern_bio12"]] for x in data],dtype=float)
    y=np.array([x[target] for x in data],dtype=float)
    train=np.array([x["fold"]!=fold for x in data])
    test=~train
    means=rows[train,1:].mean(axis=0)
    sd=rows[train,1:].std(axis=0)
    if np.any(sd==0):raise ValueError("Zero modern/geography training variation")
    xtrain=np.column_stack((np.ones(train.sum()),rows[train,0],(rows[train,1:]-means)/sd))
    xtest=np.column_stack((np.ones(test.sum()),rows[test,0],(rows[test,1:]-means)/sd))
    if np.linalg.matrix_rank(xtrain)!=6:raise ValueError("Fold current+geometry baseline rank failed")
    coef=np.linalg.lstsq(xtrain,y[train],rcond=None)[0]
    return y[test],xtest@coef

def evaluate(current,lgm):
    import numpy as np
    rows=prepare(current,lgm)
    primary=(lgm["schema"]=="ttf_Albano2025_CCSM4_PRIMARY_source_only_2_5m_LGM_cells_v01")
    folds={str(i):sum(x["fold"]==i for x in rows) for i in range(10)}
    out={"schema":"ttf_Albano2025_MPI_source_only_95grid_historical_innovation_v01",
         "status":"PRIMARY_GCM_SOURCE_ONLY_HISTORY_PREDICTOR_DIAGNOSTIC_NOT_BIOLOGICAL_EFFECT"
                   if primary else "SINGLE_SENSITIVITY_GCM_HISTORY_PREDICTOR_DIAGNOSTIC_NOT_BIOLOGICAL_EFFECT",
         "source_grid_count":95,
         "source_individual_population_count":96,
         "source_continent_cell_counts":{"NAM":47,"EUR":48},
         "modern_period":"1960-1990","LGM_period":"approx_22ka",
         "GCM":"CCSM4" if primary else "MPI-ESM-P",
         "primary_CCSM4_raster_not_yet_tested":not primary,
         "current_environmental_and_geographical_controls":["intercept","continent","latitude","longitude","modern_BIO1","modern_BIO12"],
         "spatial_holdout_folds":folds,
         "source_data_only":True,
         "fitness_genotype_trait_response_tested":False,
         "independent_LGM_effect_on_fitness_tested":False,
         "historical_predictor_exploratory_after_GCM_source_availability":not primary,
         "climate_metrics":{}}
    for label,target,units in [("BIO1","delta_bio1","degrees_C"),("BIO12","delta_bio12","mm")]:
        y=np.array([x[target] for x in rows])
        prediction=np.full(len(rows),np.nan)
        for i in range(10):
            truth,pred=rank_and_prediction(rows,i,target)
            test=[k for k,x in enumerate(rows) if x["fold"]==i]
            if len(test)!=len(pred):raise ValueError("Fold prediction count mismatch")
            prediction[test]=pred
        if not np.isfinite(prediction).all():raise ValueError("Incomplete blocked predictions")
        residual=y-prediction
        denom=np.sum((y-y.mean())**2)
        if denom==0:raise ValueError("No historical change variability")
        metrics={
           "units":units,
           "mean_historical_change_LGM_minus_modern":float(y.mean()),
           "sd_historical_change":float(y.std(ddof=1)),
           "blocked_OOF_R2":float(1-np.sum(residual**2)/denom),
           "blocked_OOF_RMSE":float(np.sqrt(np.mean(residual**2))),
           "blocked_OOF_residual_sd":float(np.std(residual,ddof=1)),
           "OOF_residual_variance_over_original_change_variance":float(np.var(residual,ddof=1)/np.var(y,ddof=1)),
           "continent_summaries":{},
        }
        for cont in ("NAM","EUR"):
            idx=np.array([x["origin"]==cont for x in rows])
            metrics["continent_summaries"][cont]={
                 "n":int(idx.sum()),
                 "mean_change":float(y[idx].mean()),
                 "change_sd":float(y[idx].std(ddof=1)),
                 "blocked_OOF_residual_sd":float(residual[idx].std(ddof=1)),
            }
        out["climate_metrics"][label]=metrics
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--modern",type=Path,default=CURRENT)
    p.add_argument("--lgm",type=Path,default=LGM)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    o=evaluate(json.loads(a.modern.read_text("utf-8")),json.loads(a.lgm.read_text("utf-8")))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(o,sort_keys=True))
if __name__=="__main__":main()
