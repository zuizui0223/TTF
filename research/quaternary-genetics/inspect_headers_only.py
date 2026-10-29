"""Inspect first-row headers of SHA-verified GenDivRange XLSX. No genetic values."""
import argparse, hashlib, io, json, urllib.request, zipfile, xml.etree.ElementTree as ET
from pathlib import Path
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def examine(item):
    url="https://ndownloader.figshare.com/files/"+str(int(item["id"]))
    with urllib.request.urlopen(url,timeout=90) as r:raw=r.read(int(item["size_bytes"])+1)
    if len(raw)!=int(item["size_bytes"]):raise ValueError("size mismatch")
    if hashlib.md5(raw).hexdigest()!=item["md5_expected"]:raise ValueError("MD5 mismatch")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            with z.open("xl/sharedStrings.xml") as f:
                for _,el in ET.iterparse(f,events=("end",)):
                    if el.tag==NS+"si":
                        shared.append("".join(t.text or "" for t in el.iter(NS+"t")))
                        el.clear()
        with z.open("xl/worksheets/sheet1.xml") as f:
            for _,row in ET.iterparse(f,events=("end",)):
                if row.tag==NS+"row":break
            else:raise ValueError("empty workbook")
            names=[]
            for c in row.findall(NS+"c"):
                v=c.find(NS+"v")
                if c.get("t")=="s":names.append(shared[int(v.text)])
                elif c.get("t")=="inlineStr":names.append("".join(t.text or "" for t in c.iter(NS+"t")))
                else:names.append(v.text if v is not None else "")
    return {"filename":item["filename"],"md5":hashlib.md5(raw).hexdigest(),
            "headers":names,"response_values_opened":False}
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--metadata",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();meta=json.loads(args.metadata.read_text())
    if meta["status"]!="FILE_LIST_ONLY_NO_GENETIC_RESPONSE":raise ValueError("not blind")
    info=[examine(f) for f in meta["files"] if f["ext"]==".xlsx"]
    out={"schema":"ttf_fresh_genetics_header_only_v01","files":info,"he_numeric_values_opened":False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"HEADERS_ONLY","files":[{"name":x["filename"],"columns":x["headers"]} for x in info]}))
