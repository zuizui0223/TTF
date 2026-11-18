"""Source-only 52 Penstemon Habroanthus BioSample taxon assignment gate.

Prospectively documented in PENSTEMON_52_BIOSAMPLE_SOURCE_TAXON_GATE_20261115.md
BEFORE opening the remaining 47 NCBI BioSample source taxon summaries.
Uses ENA original 52 SAMN IDs -> NCBI ESearch numeric UIDs -> ESummary.
No genotypes, source sampledata, RNA expression, phenotypes, RI, or DNA reads.
"""
from __future__ import annotations
import collections
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

HERE=str(Path(__file__).resolve().parent)
if HERE not in sys.path:sys.path.insert(0,HERE)
from penstemon_habroanthus_ncbi_biosample_source_pilot_20261115 import (
    STUDY, NCBI_SEARCH, NCBI, MAX_BYTES, source_samples,
    numeric_uids_from_esearch, source_metadata_aggregate,
)

N_SOURCE=52
SOURCE_COMPARISONS={
    "barbatus_to_virgatus":("Penstemon barbatus","Penstemon virgatus"),
    "cardinalis_to_glaber":("Penstemon cardinalis","Penstemon glaber"),
    "eatonii_to_laevis":("Penstemon eatonii","Penstemon laevis"),
    "labrosus_to_putus":("Penstemon labrosus","Penstemon putus"),
    "labrosus_to_deaveri":("Penstemon labrosus","Penstemon deaveri"),
}


def source_coverage(summary:dict,accessions:list[str])->dict:
    if summary.get("status")!="NCBI_BIOSAMPLE_SOURCE_SUMMARY_ONLY":
        return {"status":"NCBI_SUMMARY_SOURCE_NOT_EVALUABLE",
                "original_source_accessions":len(accessions)}
    if len(accessions)!=N_SOURCE or len(set(accessions))!=N_SOURCE:
        raise ValueError("Original ENA BioProject IDs not 52 unique BioSamples")
    records=summary.get("records") or []
    ids=[r.get("source_biosample_accession") for r in records]
    if len(records)!=N_SOURCE or len(set(ids))!=N_SOURCE or set(ids)!=set(accessions):
        return {"status":"NCBI_SOURCE_SAMPLE_ACCESSIONS_NOT_FULLY_JOINED",
                "expected_source_sample_count":N_SOURCE,
                "observed_source_summary_rows":len(records),
                "unmatched_source_sample_count":len(set(accessions)-set(ids)),
                "not_a_biological_null":True}
    labels=collections.Counter((r.get("scientific_name") or "").strip() for r in records)
    blank=labels.pop("",0)
    contrast={}
    for key,(derived,bee) in SOURCE_COMPARISONS.items():
        n_derived=labels.get(derived,0)
        n_bee=labels.get(bee,0)
        contrast[key]={"hummingbird_source_taxon":derived,
                       "bee_source_taxon":bee,
                       "hummingbird_genome_samples":n_derived,
                       "bee_genome_samples":n_bee,
                       "source_pair_covered":bool(n_derived>0 and n_bee>0)}
    return {
       "status":"ALL_52_ORIGINAL_BIOSAMPLES_TAXON_SOURCE_MATCHED" if blank==0
          else "BIOPROJECT_SOURCE_COMPLETE_TAXON_UNRESOLVED",
       "source_genome_sample_count":N_SOURCE,
       "source_distinct_taxon_labels":len(labels),
       "source_taxon_sample_counts":dict(sorted(labels.items())),
       "blank_scientific_names":blank,
       "four_origin_sister_source_coverage":contrast,
       "number_of_four_primary_sister_source_pairs_covered":sum(
           contrast[k]["source_pair_covered"] for k in
           ("barbatus_to_virgatus","cardinalis_to_glaber","eatonii_to_laevis","labrosus_to_putus")),
       "source_biosample_taxon_crosswalk":sorted(
           ({"biosample_accession":r["source_biosample_accession"],
             "scientific_name":r.get("scientific_name")} for r in records),
           key=lambda r:r["biosample_accession"]),
       "sample_attribute_values_never_exported":True,
       "not_a_molecular_phenotype_or_species_transfer_result":True,
    }


def request_json(url:str):
    req=urllib.request.Request(url,headers={
        "User-Agent":"TTF-Habroanthus-52-original-BioSample-taxonomy/1.0",
        "Accept":"application/json",
    })
    with urllib.request.urlopen(req,timeout=55) as resp:
        raw=resp.read(MAX_BYTES+1)
    if len(raw)>MAX_BYTES:raise ValueError("NCBI public source response too large")
    return raw


def main():
    accs=source_samples()
    if len(accs)!=N_SOURCE or len(set(accs))!=N_SOURCE:
        raise ValueError("Original BioProject 52 sample source drift")
    search_term=" OR ".join(f"{acc}[Accession]" for acc in accs)
    q=urllib.parse.urlencode({
      "db":"biosample","term":search_term,"retmode":"json","retmax":"100",
      "tool":"TTF_CREST_source_taxonomy_gate",
    })
    try:sr=request_json(NCBI_SEARCH+"?"+q)
    except urllib.error.HTTPError as e:
        if e.code in (400,401,403,404,429,503):
            print(json.dumps({"status":"NCBI_FULL_52_ESEARCH_ACCESS_UNAVAILABLE",
              "http_status":e.code,"source_ena_samples":len(accs),
              "source_taxon_alignment":"NOT_EVALUABLE"},indent=2))
            return
        raise
    ids=numeric_uids_from_esearch(sr,accs)
    if ids["status"]!="NCBI_SOURCE_UIDS_RESOLVED":
        print(json.dumps({"status":"NCBI_FULL_52_ESEARCH_NOT_EVALUABLE",
            "original_source_samples":len(accs),"uid_result":ids,
            "not_a_biological_null":True},indent=2,sort_keys=True))
        return
    q2=urllib.parse.urlencode({
      "db":"biosample","id":",".join(ids["ids"]),"retmode":"json",
      "tool":"TTF_CREST_source_taxonomy_gate",
    })
    try:res=request_json(NCBI+"?"+q2)
    except urllib.error.HTTPError as e:
        if e.code in (400,401,403,404,429,503):
            print(json.dumps({"status":"NCBI_FULL_52_ESUMMARY_ACCESS_UNAVAILABLE",
               "http_status":e.code,"numeric_UIDs_resolved":len(ids["ids"]),
               "original_source_samples":len(accs)},indent=2))
            return
        raise
    metadata=source_metadata_aggregate(res,accs)
    report=source_coverage(metadata,accs)
    print(json.dumps({
      "status":"PENSTEMON_HABROANTHUS_52_ORIGINAL_BIOSAMPLE_SOURCE_TAXON_GATE",
      "source_ena_bioproject":STUDY,
      "source_input_accession_count":len(accs),
      "source_esearch_numeric_uid_count":len(ids["ids"]),
      "source_matched_taxonomy":report,
      "no_floral_or_genotype_values_read":True,
      "not_a_species_transfer_or_ancestral_trait_prediction":True,
      "older_2023_tree_45_specimens_not_equated_to_2026_sra":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
