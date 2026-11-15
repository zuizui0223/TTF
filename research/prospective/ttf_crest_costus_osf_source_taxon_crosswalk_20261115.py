"""Source-only key geometry for OSF Costus scent sample/mean/phylogeny labels.

Never reads flower scent compounds, numerical traits, pollinator outcome labels,
sample IDs, or RNA expression. Only source taxon strings explicitly identified
by the 2025 authors. No automated taxonomic synonym merging.
"""
from __future__ import annotations
import csv
import io
import json
import re
import sys
from pathlib import Path

THIS=str(Path(__file__).resolve().parent)
if THIS not in sys.path:sys.path.insert(0,THIS)
from ttf_crest_costus_scent_phylogeny_source_join_20261115 import FILES, fetch_file


SOURCE_KEYS={
 "2025_floral_samples":["Species"],
 "2025_scent_samples":["Species"],
 "2025_species_mean_scent":[""],
 "2025_scent_phylo_info":["Name","Name_in_info"],
}

def norm(x:str)->str:
    return re.sub(r"\s+"," ",x.strip()).lower()

def inspect(raw:bytes,fields:list[str])->dict:
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    hd=reader.fieldnames or []
    if len(set(hd))!=len(hd) or not set(fields).issubset(hd):
        raise ValueError(f"original author key fields differ: required={fields}; provided={hd[:6]}")
    result={f:[] for f in fields}
    n=0
    for row in reader:
        n+=1
        for field in fields:
            k=norm(row.get(field) or "")
            if k:
                result[field].append(k)
    if n==0:raise ValueError("no published source rows")
    return {
      "source_row_count":n,
      "source_key_columns":fields,
      "distinct_taxon_keys_per_source_column":{k:len(set(v)) for k,v in result.items()},
      "keys":{k:sorted(set(v)) for k,v in result.items()},
      "no_numeric_value_read":True,
    }

def main():
    table={}
    for label,fields in SOURCE_KEYS.items():
        raw,problem=fetch_file(FILES[label])
        if problem:
            raise ValueError(f"original OSF source not accessible for {label}: {problem}")
        val=inspect(raw,fields)
        val["verified_sha256"]=FILES[label]["sha256"]
        table[label]=val
    flower=set(table["2025_floral_samples"]["keys"]["Species"])
    allsample=set(table["2025_scent_samples"]["keys"]["Species"])
    means=set(table["2025_species_mean_scent"]["keys"][""])
    phylo_tree=set(table["2025_scent_phylo_info"]["keys"]["Name"])
    phylo_external=set(table["2025_scent_phylo_info"]["keys"]["Name_in_info"])
    out={
      "status":"ORIGINAL_COSTUS_SAMPLE_MEAN_PHYLOGENY_ID_KEYS_NO_OUTCOMES",
      "datasets":table,
      "exact_key_overlap":{
         "flower_vs_all_scent_samples":sorted(flower&allsample),
         "flower_vs_species_mean_row_names":sorted(flower&means),
         "flower_vs_scent_phylo_Name_in_info":sorted(flower&phylo_external),
         "species_mean_vs_scent_phylo_Name_in_info":sorted(means&phylo_external),
         "species_mean_vs_scent_phylo_Name":sorted(means&phylo_tree),
      },
      "unmatched_exact_keys":{
         "flower_missing_from_means":sorted(flower-means),
         "means_missing_from_flowers":sorted(means-flower),
         "flower_missing_from_phylo_Name_in_info":sorted(flower-phylo_external),
         "phylo_Name_in_info_missing_from_flowers":sorted(phylo_external-flower),
      },
      "observed_pollinator_labels_not_authenticated":True,
      "no_scent_values_opened":True,
      "source_tree_tip_ancestry_not_yet_parsed":True,
      "no_species_level_transfer_test_run":True,
    }
    print(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
