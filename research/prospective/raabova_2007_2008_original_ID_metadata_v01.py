"""Response-blind original 2007-vs-2008 population identifier form audit.

Confirm the exact original article identity first; only tabulate source/garden
ID code tokens, never fit fitness or infer geographic coordinates. Record
layout/abstention rather than post-result alias repair or any response.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re
from collections import Counter
from pathlib import Path
from bontrager_hargreaves_finegrain_overlap_v01 import original_bytes, norm_taxon

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTRACT=HERE/"raabova_2007_2008_original_ID_metadata_contract_v01.json"
CODE=re.compile(r"(?<![A-Za-z0-9])([KS][123])(?![A-Za-z0-9])")

def classify_code(value:str,allowed:list[str])->str|None:
    matches=set(CODE.findall(value.strip()))
    if not matches:return None
    if len(matches)!=1:return "AMBIGUOUS_MULTIPLE_CODES"
    value=matches.pop()
    return value if value in allowed else None

def census(raw:bytes,contract:dict)->dict:
    src=contract["original_data"]
    if hashlib.sha256(raw).hexdigest()!=src["original_SHA256"]:
        raise ValueError("Original Hargreaves file SHA mismatch")
    if contract["status"]!="FROZEN_BEFORE_RAABOVA_PAIR_ORIGINAL_SOURCE_ID_RECORD_INSPECTION":
        raise ValueError("Original-source code census not authorized")
    allowed=src["columns"]
    if allowed!=["study","taxon","functional_group","source","site"] or "fitness" in allowed:
        raise ValueError("Metadata allowlist changed")
    paper=contract["validated_primary_papers"]
    if paper["sudova_2010"]["table1_reference_50_is"]!=paper["different_2008"]["doi"]:
        raise ValueError("Sudova original reference #50 is not 2008")
    if paper["different_2008"]["doi"]==paper["target_2007"]["doi"]:
        raise ValueError("Raabova 2007 and 2008 must be distinct")
    reader=csv.reader(io.StringIO(raw.decode("cp1252",errors="strict"),newline=""),strict=True)
    header=next(reader,[])
    if len(header)!=47 or len(set(header))!=47 or not set(allowed).issubset(header):
        raise ValueError("Original 47 field source metadata missing")
    at={col:header.index(col) for col in allowed}
    target={paper["target_2007"]["study_label"],paper["different_2008"]["study_label"]}
    allowed_codes=contract["exact_id_rule"]["allowed_known_2007_population_codes"]
    if allowed_codes!=["K1","K2","K3","S1","S2","S3"]:
        raise ValueError("Original 2007 region-code set changed")
    rows={s:0 for s in target}
    taxon_counts={s:Counter() for s in target}
    ids_source={s:set() for s in target}
    ids_garden={s:set() for s in target}
    physical=0
    for row in reader:
        if not row:continue
        physical+=1
        if physical>200000 or len(row)!=47:
            raise ValueError("Source row cap/width drift")
        meta={k:row[idx].strip() for k,idx in at.items()}
        st=meta["study"]
        if st not in target:continue
        if meta["functional_group"] not in {"annual","herbaceous_perennial","woody_perennial"}:continue
        rows[st]+=1
        tax=meta["taxon"].casefold().strip()
        if norm_taxon(meta["taxon"])=="aster amellus":
            taxon_counts[st]["strict_binomial_aster_amellus"]+=1
        elif "aster" in tax and "amellus" in tax:
            taxon_counts[st]["aster_amellus_nonbinomial_with_extra_text"]+=1
        else:
            taxon_counts[st]["other_taxon_text"]+=1
        for k,d in (("source",ids_source),("site",ids_garden)):
            if meta[k]:d[st].add(meta[k])

    if physical!=11088:
        raise ValueError("Original published source row count drift")
    y2007=paper["target_2007"]["study_label"]
    y2008=paper["different_2008"]["study_label"]
    def result(st):
        x=ids_source[st];y=ids_garden[st]
        csrc=Counter(classify_code(s,allowed_codes) or "NONE" for s in x)
        cdst=Counter(classify_code(s,allowed_codes) or "NONE" for s in y)
        return {
            "original_study_label":st,
            "original_plant_rows":rows[st],
            "unique_original_source_identifiers":len(x),
            "unique_original_destination_identifiers":len(y),
            "taxon_text_categories":dict(sorted(taxon_counts[st].items())),
            "literal_original_2007_population_code_source_identifier_counts":dict(sorted(csrc.items())),
            "literal_original_2007_population_code_destination_identifier_counts":dict(sorted(cdst.items())),
            "source_identifiers_with_2007_code_tokens":sum(v for k,v in csrc.items() if k in allowed_codes),
            "destination_identifiers_with_2007_code_tokens":sum(v for k,v in cdst.items() if k in allowed_codes),
            "no_raw_original_identifiers_exported":True,
        }
    return {
        "schema":"ttf_Raabova2007_2008_literal_original_ID_form_source_only_v01",
        "status":"SOURCE_ONLY_PRIMARY_STUDY_LABEL_FORM_NOT_GEOGRAPHIC_VERIFICATION",
        "source_original_SHA256":src["original_SHA256"],
        "original_total_source_rows":physical,
        "raabova2007":result(y2007),
        "raabova2008":result(y2008),
        "sudova_2010_Table1_references_2008_not_2007":True,
        "independent_pop_site_2007_georeference_certified":False,
        "exact_study_2007_locality_code_to_Hargreaves_ID_mapping_certified":False,
        "fitness_or_historical_climate_values_opened":False,
        "previous_invalid_2007_Sudova_georeference_benchmark_reopened":False,
        "Bontrager44of50_stopped_gate_unchanged":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    c=json.loads(CONTRACT.read_text(encoding="utf-8"))
    source=json.loads((ROOT/c["original_data"]["original_metadata_contract"]).read_text(encoding="utf-8"))["original"]["hargreaves"]
    raw=original_bytes(source)
    out=census(raw,c)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print("RAABOVA_PRIMARY_LITERAL_CODE_AUDIT=",json.dumps({
        "raabova2007":out["raabova2007"],
        "raabova2008":out["raabova2008"],
        "status":out["status"]
    },sort_keys=True))

if __name__=="__main__":main()
