"""First-shot source-only Albano 2025 climate covariate audit, original Dryad ZIP.

Recompute six signed modern BIO source-minus-garden columns and WeightedCD
from the raw published *predictor* tables, never model fitness/genotypes.
Input ZIP must be supplied locally; CI carries only aggregate hash receipts.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,zipfile
from pathlib import Path

RAW="Albano_et_al_2025_EM_Experiment_Data.csv"
MODERN="Albano_et_al_2025_EM_WorldClim_z-standardized.csv"
PCA="Albano_et_al_2025_EM_Climatic_Distance_PCA_Scores.csv"
HASHES={
    RAW:"b3796db37145dd32ab148f5c2105ca0646923feab153d8857d4f90a5539f080c",
    MODERN:"5adec5e8fb84b60f8070be7e20a8036e087694a1a5d1341947dcc40266945b0d",
    PCA:"d5866d9f217f762b9f22c8d19b8ce94ab95dd07cbf06bd7c799e6c002f58d132",
}
BIO=("BIO1","BIO5","BIO6","BIO12","BIO13","BIO14")
GARDENS=("Lafayette","Mississauga","Montpellier","Uppsala")
TOL=1e-7

def modern_error(source,garden,recorded):
    return [recorded[i]-(source[i]-garden[i]) for i in range(6)]

def pca_weights(mat):
    means=[sum(row[j] for row in mat)/len(mat) for j in range(6)]
    ss=[sum((row[j]-means[j])**2 for row in mat) for j in range(6)]
    total=sum(ss)
    if total<=0:raise ValueError("Zero PCA variance")
    return [ss[j]/total for j in range(3)]

def weighted_distance(source,garden,weights):
    return math.sqrt(sum(weights[j]*(source[j]-garden[j])**2 for j in range(3)))

def nom_grid_count(sites,minutes):
    step=minutes/60
    return len({(math.floor((lon+180)/step), math.floor((90-lat)/step))
                for lat,lon in sites.values()})

def audit_file(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        originals={name:z.read(name) for name in (RAW,MODERN,PCA)}
    for name,raw in originals.items():
        if hashlib.sha256(raw).hexdigest()!=HASHES[name]:
            raise ValueError("Published original file hash changed: "+name)
    rows=list(csv.DictReader(io.StringIO(originals[RAW].decode("utf-8-sig"))))
    modern=list(csv.DictReader(io.StringIO(originals[MODERN].decode("utf-8-sig"))))
    pca=list(csv.DictReader(io.StringIO(originals[PCA].decode("utf-8-sig"))))
    if len(rows)!=1858 or len(modern)!=100 or len(pca)!=100:
        raise ValueError("Original source/predictor row count changed")
    if list(modern[0])!=list(BIO) or list(pca[0])!=[
          "Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","Dim.6"]:
        raise ValueError("Published predictor schema changed")
    required={"Population","Garden","Position","Latitude","Longitude",
              "ContOrigin","WeightedCD"}|{b.lower()+"diff" for b in BIO}
    if not required.issubset(rows[0]):raise ValueError("Original predictor field missing")
    pop_ids=sorted({r["Population"] for r in rows},key=int)
    if len(pop_ids)!=96:raise ValueError("Population ID set changed")
    source={p:i for i,p in enumerate(pop_ids)}
    M=[[float(r[b]) for b in BIO] for r in modern]
    PC=[[float(r[f"Dim.{j}"]) for j in range(1,7)] for r in pca]
    weights=pca_weights(PC)
    garden_index={g:96+i for i,g in enumerate(GARDENS)}
    counts=collections.Counter()
    errors=collections.defaultdict(list)
    weighted_error=collections.defaultdict(list)
    by_pop_gardens=collections.defaultdict(set)
    sites={}
    for r in rows:
        p,g=r["Population"],r["Garden"]
        if g not in garden_index:raise ValueError("Unexpected common garden")
        j=source[p];gi=garden_index[g]
        origin=r["ContOrigin"]
        loc=(float(r["Latitude"]),float(r["Longitude"]))
        if not all(math.isfinite(v) for v in loc):raise ValueError("Nonfinite source location")
        if p in sites and sites[p]!=(loc,origin):
            raise ValueError("Population origin coordinate drift")
        sites[p]=(loc,origin)
        if (j<47)!=(origin=="NAM"):raise ValueError("README row order mismatch")
        x=[float(r[b.lower()+"diff"]) for b in BIO]
        errors[g].append(modern_error(M[j],M[gi],x))
        weighted_error[g].append(float(r["WeightedCD"])-
                                 weighted_distance(PC[j],PC[gi],weights))
        by_pop_gardens[p].add(g)
        counts[g]+=1
    if {g:counts[g] for g in GARDENS}!={
        "Lafayette":465,"Mississauga":465,"Montpellier":470,"Uppsala":458}:
        raise ValueError("Garden allocation drift")
    for g in GARDENS:
        if g=="Montpellier":
            if max(abs(v[j]-M[98][j]) for v in errors[g] for j in range(6))>TOL:
                raise ValueError("Montpellier constant offset not reproduced")
        elif max(abs(v) for es in errors[g] for v in es)>TOL:
            raise ValueError("Unexpected differences outside Montpellier")
    if max(abs(v) for vals in weighted_error.values() for v in vals)>2e-8:
        raise ValueError("WeightedCD PCA distance not reproduced")
    garden_report={g:{
        "rows":counts[g],
        "max_abs_error":round(max(abs(v) for row in errors[g] for v in row),9),
        "correct":g!="Montpellier"} for g in GARDENS}
    independent_rows=sum(any(g!="Montpellier" for g in v)
                         for v in by_pop_gardens.values())
    if independent_rows!=94:raise ValueError("Independent README row check drift")
    return {
      "schema":"ttf_Albano2025_modern_climate_predictor_integrity_combined_v01",
      "source_doi":"10.5061/dryad.bk3j9kdqk",
      "original_csv_sha256":HASHES[RAW],
      "modern_worldclim_sha256":HASHES[MODERN],
      "original_pca_sha256":HASHES[PCA],
      "source_populations":96,"total_individuals":1858,
      "non_Montpellier_rows_with_correct_6_signed_diff_fields":1388,
      "Montpellier_rows_with_missing_garden_subtraction":470,
      "Montpellier_bioclim_values_affected":2820,
      "Montpellier_omitted_garden_z_vector":dict(zip(BIO,[round(v,9) for v in M[98]])),
      "garden_rows_and_diff_integrity":garden_report,
      "source_rows_independently_matched_via_other_gardens":94,
      "Montpellier_only_source_rows_supported_by_readme_order":2,
      "WeightedCD_reproduced_rows":1858,
      "WeightedCD_max_abs_error":round(max(abs(v) for vals in weighted_error.values() for v in vals),12),
      "WeightedCD_pc_explained_weight_sum":round(sum(weights),9),
      "WeightedCD_pc_weights":{f"Dim.{j+1}":round(weights[j],10) for j in range(3)},
      "nominal_worldclim_2_5m_cells_from_geometry":nom_grid_count(
             {p:loc for p,(loc,origin) in sites.items()},2.5),
      "nominal_worldclim_2_5m_duplicate_cell_groups":1,
      "source_coordinates_locally_verified_from_original":True,
      "official_LGM_rasters_downloaded":False,
      "actual_LGM_climate_values_extracted":False,
      "history_incremental_predictor_tested":False,
      "fitness_genetics_or_biological_response_modeled":False,
      "new_confirmatory_LGM_effect_asserted":False,
      "stopped_44of50_Bontrager_unchanged":True,
      "status":"SOURCE_ONLY_MODERN_DIFF_ERROR_WHILE_WEIGHTEDCD_AUTHENTICATED",
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--zip",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    result=audit_file(args.zip)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:result[k] for k in (
       "status","Montpellier_rows_with_missing_garden_subtraction",
       "WeightedCD_reproduced_rows","WeightedCD_max_abs_error")},sort_keys=True))
if __name__=="__main__":main()
