"""Verify official XLSX bytes and inspect header/row count only; do NOT inspect GenDist."""
import argparse
import hashlib
import io
import json
import urllib.request
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

S="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def exact_source_download(lock):
    assert lock["status"]=="FROZEN_BEFORE_NUMERICAL_FST_OPENING"
    id=int(lock["xlsx_file_id"])
    url="https://ndownloader.figshare.com/files/"+str(id)
    if url!=lock["download_url"]:raise ValueError("wrong source URL")
    n=int(lock["size_bytes"])
    req=urllib.request.Request(url,headers={"User-Agent":"TTF-CaliPopGen-source-headers-v01"})
    with urllib.request.urlopen(req,timeout=120) as response:
        raw=response.read(n+1)
    if len(raw)!=n:raise ValueError("Source byte size mismatch")
    if hashlib.md5(raw).hexdigest()!=lock["md5"]:
        raise ValueError("Original Figshare XLSX MD5 mismatch")
    return raw

def column_index(ref):
    n=0
    for char in ref:
        if not char.isalpha():break
        n=n*26+ord(char.upper())-64
    return n-1

def schema_only(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        if "xl/worksheets/sheet1.xml" not in z.namelist():
            raise ValueError("Unexpected first worksheet")
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            with z.open("xl/sharedStrings.xml") as f:
                for _,el in ET.iterparse(f,events=("end",)):
                    if el.tag==S+"si":
                        strings.append("".join(t.text or "" for t in el.iter(S+"t")))
                        el.clear()
        columns={}
        rows=0
        with z.open("xl/worksheets/sheet1.xml") as f:
            for _,row in ET.iterparse(f,events=("end",)):
                if row.tag!=S+"row":continue
                rows+=1
                if rows==1:
                    for cell in row.findall(S+"c"):
                        i=column_index(cell.attrib.get("r",""))
                        value=cell.find(S+"v")
                        inline=cell.find(S+"is")
                        if cell.attrib.get("t")=="s" and value is not None:
                            text=strings[int(value.text)]
                        elif inline is not None:
                            text="".join(t.text or "" for t in inline.iter(S+"t"))
                        else:text=value.text if value is not None else ""
                        columns[i]=text
                row.clear()
        if rows<=1:raise ValueError("No data rows")
        ordered=[columns.get(i,"") for i in range(max(columns)+1)]
        if not any("gendist" in str(c).lower() for c in ordered):
            raise ValueError("Expected genetic differentiation metric column absent")
        return {"rows_in_first_worksheet":rows-1,"columns":ordered,
                "column_count":len(ordered),"genetic_numeric_values_opened":False}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--lock",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    lock=json.loads(args.lock.read_text())
    blob=exact_source_download(lock)
    schema=schema_only(blob)
    result={"schema":"ttf_calipopgen_fst_xlsx_header_census_v01",
            "status":"SOURCE_BYTES_VERIFIED_HEADERS_ONLY",
            "source_article":int(lock["source_article"]),
            "source_file_id":int(lock["xlsx_file_id"]),
            "source_md5":hashlib.md5(blob).hexdigest(),**schema}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
