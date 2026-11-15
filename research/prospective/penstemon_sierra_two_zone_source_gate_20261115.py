"""Original 2026 Sierra Penstemon two-hybrid-zone matched phenotype/genome source gate.

Author: benstemon/sierra-hybrid-zone, pinned original public Git blob data.
Only original specimen ID and categorical site/pop/class keys and MISSINGNESS
of prespecified hue/admixture predictors are checked. NO actual floral
values, K1/K2 ancestry magnitudes, GWAS, SNP cline effects decoded.
"""
from __future__ import annotations
import collections
import csv
import hashlib
import io
import json
import urllib.request

REPO="benstemon/sierra-hybrid-zone"
COMMIT="0ba4842824bcd9ddc824eedea303ad1c176fc82b"
ORIGINAL={
 "joint_phenotype_admixture":{
   "path":"1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv",
   "blob":"7cd5abd0f2a5ad6de1563d3ddd2c0d857b037de8",
 },
 "PCA_phenotype_admixture":{
   "path":"3.admixture-pca/data/PCA_values_phenotypes_admixture-sierras.csv",
   "blob":"d02506158716404a57f8f371c0077dc9989b5cb0",
 }
}
REQUIRED=("newID","modified_hue","K1","K2","pop","subpopulation","class")
CATEGORICAL=("subpopulation","class","pop","region")
MISSING_TOKENS={"","NA","NaN","nan","NULL","N/A"}

def fetch_source(original):
    url=f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{original['path']}"
    with urllib.request.urlopen(url,timeout=40) as response:
        raw=response.read(800_001)
    if len(raw)>800_000:raise ValueError("Author source file too large")
    sha=hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()
    if sha!=original["blob"]:
        raise ValueError("Original author Git source SHA mismatch")
    return raw

def audit(raw):
    rd=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    fields=rd.fieldnames or []
    if len(fields)!=len(set(fields)) or not set(REQUIRED)<=set(fields):
        raise ValueError(f"Source missing expected columns: {set(REQUIRED)-set(fields)}")
    ids=set()
    repeated=[]
    counts=collections.Counter()
    cat={name:collections.Counter() for name in CATEGORICAL}
    missing=collections.Counter()
    for row in rd:
        counts["rows"]+=1
        id=(row["newID"] or "").strip()
        if not id:raise ValueError("Empty source specimen identifier")
        if id in ids:repeated.append(id)
        ids.add(id)
        for name in CATEGORICAL:
            cat[name][(row.get(name) or "").strip()]+=1
        for name in ("modified_hue","chroma","avg_floral_tube_length",
                     "avg_tube_width_at_mouth","K1","K2"):
            if name in fields and (row.get(name) or "").strip() in MISSING_TOKENS:
                missing[name]+=1
        if all((row.get(x) or "").strip() not in MISSING_TOKENS for x in ("modified_hue","K1")):
            counts["hue_ancestry_nonmissing_rows"]+=1
        # No numeric cell content read or returned.
    return {
       "rows":counts["rows"],"unique_specimen_IDs":len(ids),
       "duplicate_original_ID_count":len(repeated),
       "duplicates_only_if_present":repeated[:30],
       "header_columns":fields,
       "categorical_source_counts":{
           key:dict(sorted(value.items())) for key,value in cat.items()},
       "missingness_by_source_field":dict(sorted(missing.items())),
       "hue_K1_both_present_structural_rows":counts["hue_ancestry_nonmissing_rows"],
       "no_flower_or_genomic_values_read":True,
       "ids_for_internal_crossfile_set":ids,
    }

def main():
    sources={}
    identifier_sets={}
    for key,info in ORIGINAL.items():
        data=audit(fetch_source(info))
        identifier_sets[key]=data.pop("ids_for_internal_crossfile_set")
        data.update({"source_path":info["path"],"original_git_blob":info["blob"]})
        sources[key]=data
    a,b=identifier_sets.values()
    print(json.dumps({
       "status":"PENSTEMON_SIERRA_TWO_ZONE_SOURCE_KEYS_ONLY",
       "source_commit":COMMIT,
       "source_files":sources,
       "source_specimen_ID_intersection":len(a&b),
       "joint_only_ID_count":len(a-b),
       "PCA_only_ID_count":len(b-a),
       "same_individual_genome_and_phenotype_columns_present_but_not_validated":True,
       "not_a_pheno_genomic_transfer_score":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
