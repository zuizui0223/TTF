"""Source-only transplant original ZIP CSV HEADERS and README metadata excerpts.

Input original ZIP passes both publisher Dryad SHA256 and UBC mirror MD5
before inspecting. Never decode any ZIP CSV data row (row 2+); no
fitness, genotype, local adaptation sign, magnitude or site values used.
"""
from __future__ import annotations
import argparse,csv,io,json,zipfile
from pathlib import Path
from urllib.request import Request,urlopen
from bontrager_mirror_zip_v01 import verify_original_zip,ORIGINAL_URL

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"bontrager_member_header_contract_v01.json"
ORIGINAL=HERE/"bontrager_mirror_zip_contract_v01.json"

def csv_headers_from_archive(raw,selected):
    items=[]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names=set(z.namelist())
        for key in selected:
            if key not in names or not key.endswith(".csv"):
                raise ValueError("Original predeclared ZIP CSV missing")
            with z.open(key) as f: first=f.readline(16384)
            try:
                headers=next(csv.reader([first.decode("utf-8-sig").rstrip("\r\n")],strict=True))
            except (UnicodeDecodeError,csv.Error) as err:
                raise ValueError("Original source CSV header cannot be decoded") from err
            if not headers or len(headers)>300 or len(set(headers))!=len(headers):
                raise ValueError("Original source CSV header absent or ambiguous")
            items.append({"archive_member":key,"headers":headers,"column_count":len(headers),
                          "source_rows_read":False,"fitness_values_read":False})
    return items

def codebook_excerpts(raw,name,max_bytes):
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if name not in z.namelist():raise ValueError("Original README missing")
        with z.open(name) as f:txt=f.read(max_bytes).decode("utf-8-sig","replace")
    keys=("source","location","latitude","longitude","fitness","study","site",
          "climate","transplant","population","file","data")
    lines=[]
    for i,s in enumerate(txt.splitlines(),1):
        if s.strip() and any(k in s.lower() for k in keys):
            lines.append({"line":i,"source_codebook_excerpt":s.strip()[:180]})
        if len(lines)>=30:break
    return lines

def run():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    c=json.loads(CONTRACT.read_text());lock=json.loads(ORIGINAL.read_text())
    if (c.get("status")!="FROZEN_AFTER_ORIGINAL_ZIP_SHA_MATCH_BEFORE_MEMBER_HEADER_ACCESS"
        or lock.get("status")!="FROZEN_BEFORE_UBC_BOREALIS_ORIGINAL_ZIP_CONTENT"
        or c["zip_original_sha256"]!=lock["dryad_sha256"]
        or c["zip_original_md5"]!=lock["mirror_md5"]):
        raise ValueError("Original pre-response ZIP identity contract drifted")
    output={"schema":"ttf_Bontrager_transplant_original_source_CSV_headers_v01",
            "original_Dryad_sha256":c["zip_original_sha256"],
            "raw_fitness_rows_read":False,"rowwise_location_values_read":False,
            "historical_climate_association_fit":False,
            "original_Haldane_terminal_unchanged":True}
    try:
        q=Request(ORIGINAL_URL,headers={"User-Agent":"TTF-Bontrager-original-ZIP-column-names-only"})
        with urlopen(q,timeout=90) as r:raw=r.read(lock["max_download_bytes"]+1)
        verify_original_zip(raw,lock)
        output["source_member_columns"]=csv_headers_from_archive(raw,c["selected_source_csv_members"])
        output["source_README_relevant_codebook_excerpts"]=codebook_excerpts(
            raw,c["readme_member"],c["allow_readme_first_bytes"])
        output["status"]="SOURCE_CSV_HEADERS_AND_CODEBOOK_VERIFIED_NO_FITNESS_ROWS"
    except Exception as err:
        output.update(status="SOURCE_ZIP_HEADER_STOP_NOT_BIOLOGICAL_NULL",
                      error_type=type(err).__name__,error_message=str(err)[:300])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print(output["status"],len(output.get("source_member_columns",[])))

if __name__=="__main__":run()
