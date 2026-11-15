"""Nature Communications 2024 Timema XLSX: schema/row counts only.

Does not decode numeric cells, reproductive-isolation or Fst observations.
Original XLSX bytes verified with SHA256 from prior publisher metadata.
"""
from __future__ import annotations
import io
import json
import hashlib
import re
import zipfile
import urllib.request
import xml.etree.ElementTree as ET
from posixpath import normpath

BASE="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-46294-9/MediaObjects/"
FILES={
 "source_data":("41467_2024_46294_MOESM6_ESM.xlsx",
  "c309402e1b1fa58ab0fb0763086a155aff028cdaa04b804c2c91796c5bed61e2"),
 "supplement_1_5":("41467_2024_46294_MOESM4_ESM.xlsx",
  "7e2b6210c182859d7b4d77ad8b0ca0858a1e437c0106a9c42e6c3787c24a33db")
}
NS="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
PKG="http://schemas.openxmlformats.org/package/2006/relationships"
OFFICE="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
HINTS=re.compile(r"pair|species|taxon|population|host|cross|genetic|distance|isolation|site|location|plant|nucdna|gb[s]?|fst|sample",re.I)

def safe_string(cell:ET.Element, strings:list[str])->str|None:
    # Only text-like metadata; no numeric cells are ever returned.
    t=cell.get("t")
    if t=="s":
        v=cell.find(f"{{{NS}}}v")
        if v is None or v.text is None: return None
        idx=int(v.text)
        return strings[idx] if 0<=idx<len(strings) else None
    if t=="inlineStr":
        return "".join(x.text or "" for x in cell.findall(f".//{{{NS}}}t"))
    return None

def parse_workbook(raw:bytes,expected_sha:str)->dict:
    if hashlib.sha256(raw).hexdigest()!=expected_sha:
        raise ValueError("publisher XLSX original SHA256 mismatch")
    if len(raw)>1_000_000 or not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("invalid or overlarge original XLSX")
    with zipfile.ZipFile(io.BytesIO(raw)) as book:
        root=ET.fromstring(book.read("xl/workbook.xml"))
        rels=ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
        relmap={r.get("Id"):r.get("Target") for r in rels.findall(f"{{{PKG}}}Relationship")}
        ss=[]
        if "xl/sharedStrings.xml" in book.namelist():
            shared=ET.fromstring(book.read("xl/sharedStrings.xml"))
            for si in shared.findall(f"{{{NS}}}si"):
                ss.append("".join(t.text or "" for t in si.findall(f".//{{{NS}}}t")))
        result=[]
        for sh in root.findall(f".//{{{NS}}}sheets/{{{NS}}}sheet"):
            target=relmap.get(sh.get(f"{{{OFFICE}}}id"))
            if not target:raise ValueError("sheet relation unresolved")
            path=normpath(target.lstrip("/") if target.startswith("/xl/") else "xl/"+target)
            if target.startswith("/xl/"):path=target.lstrip("/")
            if path not in book.namelist():raise ValueError("sheet XML absent")
            tree=ET.fromstring(book.read(path))
            dimension=tree.find(f"{{{NS}}}dimension")
            sheetdata=tree.find(f"{{{NS}}}sheetData")
            rows=list(sheetdata.findall(f"{{{NS}}}row")) if sheetdata is not None else []
            hints=[]
            for row in rows[:4]:
                for cell in row.findall(f"{{{NS}}}c"):
                    label=safe_string(cell,ss)
                    if label and len(label)<100 and HINTS.search(label):
                        # Short header hints only; may include numerical text in
                        # unstructured labels but *never* worksheet numeric values.
                        hints.append({"cell":cell.get("r"),"header_hint":label})
            result.append({"sheet":sh.get("name"),"declared_range":dimension.get("ref") if dimension is not None else None,
                           "row_element_count":len(rows),"first_four_rows_header_hints":hints[:30],
                           "original_numeric_cells_not_decoded":True})
    return {"worksheet_count":len(result),"sheets":result}


def main():
    reports=[]
    for label,(file,sha) in FILES.items():
        req=urllib.request.Request(BASE+file,headers={"User-Agent":"TTF-source-schema/1.0"})
        with urllib.request.urlopen(req,timeout=35) as resp:
            raw=resp.read(1_000_001)
        out=parse_workbook(raw,sha)
        out.update({"source":label,"file":file,"source_sha256":sha})
        reports.append(out)
    print(json.dumps({"status":"PUBLISHER_SHEET_GRAIN_ONLY_NO_BIOLOGICAL_RESPONSES",
                      "sources":reports,"dryad_original_csv_access":"AUTH_401",
                      "species_disjoint_TTF_claim":"NOT_QUALIFIED",
                      "read_numeric_outcomes":False},indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
