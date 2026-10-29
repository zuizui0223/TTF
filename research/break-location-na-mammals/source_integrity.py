"""Verify published Jensen2024 source archives, inspect only ZIP inventory and XLSX headers."""
import argparse,collections,hashlib,io,json,urllib.request,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def fetch(entry,dest):
    u="https://figshare.com/ndownloader/files/"+str(entry["id"])+"?private_link="+entry["share"]
    req=urllib.request.Request(u,headers={"User-Agent":"TTF-break-source-integrity-v01"})
    n=0;h=hashlib.sha256();m=hashlib.md5()
    with urllib.request.urlopen(req,timeout=120) as r,dest.open("wb") as f:
        while True:
            chunk=r.read(1048576)
            if not chunk:break
            n+=len(chunk)
            if n>entry["bytes"]:raise ValueError("larger than original receipt")
            f.write(chunk);h.update(chunk);m.update(chunk)
    if n!=entry["bytes"] or m.hexdigest()!=entry["md5"]:
        raise ValueError("original file size/MD5 mismatch")
    return {"bytes":n,"md5":m.hexdigest(),"sha256":h.hexdigest()}
def zipinventory(dest):
    with zipfile.ZipFile(dest) as z:
        files=[f for f in z.infolist() if not f.is_dir()]
        names=[f.filename for f in files]
        if any(k.startswith("/") or ".." in Path(k).parts for k in names):
            raise ValueError("unsafe ZIP path")
        return {"member_count":len(files),"uncompressed_bytes":sum(f.file_size for f in files),
          "extensions":dict(collections.Counter(Path(x).suffix.lower() for x in names)),
          "names":names[:500],"coordinates_or_break_features_read":False}
def headers(dest):
    with zipfile.ZipFile(dest) as z:
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            with z.open("xl/sharedStrings.xml") as f:
                for _,el in ET.iterparse(f,events=("end",)):
                    if el.tag==NS+"si":
                        strings.append("".join(t.text or "" for t in el.iter(NS+"t")));el.clear()
        with z.open("xl/worksheets/sheet1.xml") as f:
            for _,row in ET.iterparse(f,events=("end",)):
                if row.tag==NS+"row":break
            else:raise ValueError("no header row")
            names=[]
            for c in row.findall(NS+"c"):
                v=c.find(NS+"v")
                if c.get("t")=="s" and v is not None:names.append(strings[int(v.text)])
                elif c.get("t")=="inlineStr":names.append("".join(t.text or "" for t in c.iter(NS+"t")))
                else:names.append(v.text if v is not None else "")
        return {"header_names":names,"data_rows_read":False}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();m=json.loads(a.manifest.read_text())
    if m["status"]!="FROZEN_BEFORE_DOWNLOAD_OR_BREAK_LOCATION_OPENING":
        raise ValueError("unfrozen metadata")
    a.work.mkdir(parents=True,exist_ok=True);a.out.parent.mkdir(parents=True,exist_ok=True)
    outputs=[]
    for e in m["files"]:
        dest=a.work/(str(e["id"])+Path(e["name"]).suffix)
        receipt=fetch(e,dest)
        ext=dest.suffix.lower()
        structure=zipinventory(dest) if ext==".zip" else headers(dest) if ext==".xlsx" else {"not_executed":True}
        outputs.append({"role":e["role"],"article_id":e["article_id"],"file_id":e["id"],
          "file_name":e["name"],**receipt,**structure})
        dest.unlink()
        print("VERIFY",e["role"],receipt["bytes"],"archive_members",structure.get("member_count"),flush=True)
    result={"schema":"ttf_jensen_original_source_qualified_inventory_v01",
       "status":"CHECKSUMS_AND_ARCHIVE_INVENTORY_COMPLETE_NO_RESPONSE_OPENING",
       "files":outputs,"break_coordinates_opened":False,"model_tested":False}
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"files":len(outputs)}))
if __name__=="__main__":main()
