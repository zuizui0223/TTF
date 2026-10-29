"""Conservative response-blind geographical availability census; no He values."""
import argparse,collections,io,json,math,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
from inspect_headers_only import fetch_verified
S="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def ci(ref):
    v=0
    for c in ref:
        if not c.isalpha():break
        v=v*26+ord(c.upper())-64
    return v-1
def census(blob):
    groups=collections.defaultdict(set);total=present=geo=0
    with zipfile.ZipFile(io.BytesIO(blob)) as z,z.open("xl/worksheets/sheet1.xml") as f:
        for _,row in ET.iterparse(f,events=("end",)):
            if row.tag!=S+"row":continue
            if row.get("r")=="1":row.clear();continue
            total+=1;vals={}
            for cell in row.findall(S+"c"):
                k=ci(cell.get("r",""))
                if k not in (0,1,10,11,21):continue
                v=cell.find(S+"v")
                if k==21:
                    # Genetic-response numeric value remains unread.
                    vals[k]=(cell.get("t") in ("n",None)) and v is not None and bool(v.text)
                else:vals[k]=v.text if v is not None else None
            if not vals.get(21):row.clear();continue
            present+=1
            try:
                lat=float(vals[10]);lon=float(vals[11])
                assert math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180
            except (KeyError,TypeError,ValueError,AssertionError):
                row.clear();continue
            key=(vals.get(0),vals.get(1))
            if None not in key:groups[key].add((lat,lon));geo+=1
            row.clear()
    eligible=[(sp,study) for (sp,study),pts in groups.items() if
        len(pts)>=5 and (max(p[0] for p in pts)-min(p[0] for p in pts))*111>=100]
    return dict(status="CONSERVATIVE_GEOGRAPHIC_ONLY_PRECHECK",population_rows=total,
        he_presence_mask_rows=present,valid_locations=geo,
        study_species_groups=len(groups),eligible_groups_lower_bound=len(eligible),
        eligible_species_id_lower_bound=len({sp for sp,study in eligible}),
        genetic_response_values_opened=False)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--metadata",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();m=json.loads(a.metadata.read_text())
    item=next(v for v in m["files"] if v["filename"]=="pop_tab_v2025-03-31.xlsx")
    out=census(fetch_verified(item))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out))
