"""Jensen2024 verified common-name to scientific-name join, without opening AUC or break coordinates."""
import argparse,collections,hashlib,json,re,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
from source_integrity import fetch
S="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
NAMESPACE="ttf-jensen-na-mammal-break-v01|"
def normalize(x):return re.sub(r"[^a-z0-9]+","",str(x or "").lower())
def excel_common_scientific(path):
    rows=[]
    with zipfile.ZipFile(path) as zip:
        idx=[]
        with zip.open("xl/worksheets/sheet1.xml") as f:
            for _,r in ET.iterparse(f,events=("end",)):
                if r.tag!=S+"row":continue
                values={}
                for cell in r.findall(S+"c"):
                    ref=cell.get("r","")
                    if not ref.startswith(("A","B")):continue
                    v=cell.find(S+"v")
                    if v is None:continue
                    values[ref[0]]=(cell.get("t"),v.text)
                if r.get("r")!="1":
                    rows.append((values.get("A"),values.get("B")))
                r.clear()
        used={int(v) for pair in rows for cell in pair if cell and cell[0]=="s" for v in (cell[1],)}
        strings={}
        with zip.open("xl/sharedStrings.xml") as f:
            j=0
            for _,node in ET.iterparse(f,events=("end",)):
                if node.tag!=S+"si":continue
                if j in used:
                    strings[j]="".join(t.text or "" for t in node.iter(S+"t"))
                j+=1;node.clear()
    def decode(entry):
        if entry is None:return ""
        t,v=entry
        return strings[int(v)] if t=="s" else str(v)
    return [(decode(a).strip(),decode(b).strip()) for a,b in rows]
def join(public_species,grid):
    mapping=collections.defaultdict(set)
    for common,latin in public_species:
        if common and latin:mapping[normalize(common)].add(latin)
    records=[];missing=[];ambiguous=[]
    for v in grid["full_sampling_counts"]:
        if not v.get("eligible"):continue
        names=mapping.get(normalize(v["species_common_name"]),set())
        if len(names)==0:missing.append(v["species_common_name"]);continue
        if len(names)>1:
            ambiguous.append(v["species_common_name"]);continue
        latin=next(iter(names))
        split=hashlib.sha256((NAMESPACE+latin.lower()).encode()).digest()[0]%2
        records.append({"common_name":v["species_common_name"],
         "scientific_name":latin,"split_role":split,
         "case_points":v["positive_break_grid_points"],
         "background_points":v["background_surveyed_grid_points"],
         "synthetic_estimator_min10_cases":v["positive_break_grid_points"]>=10
             and v["background_surveyed_grid_points"]>=10})
    stricter=[x for x in records if x["synthetic_estimator_min10_cases"]]
    roles=[{r["scientific_name"] for r in stricter if r["split_role"]==i} for i in (0,1)]
    if len({r["scientific_name"] for r in records})!=len(records):
        raise ValueError("different source common names map to same scientific taxon")
    status="PASS_LATIN_SPECIES_HELDOUT_SUPPORT" if len(stricter)>=50 and min(map(len,roles))>=20 else "NOT_EVALUABLE_LATIN_SPECIES_SPLIT"
    return {"schema":"ttf_jensen2024_species_identity_and_holdout_gate_v01",
       "status":status,"eligible_common_name_species":sum(x["eligible"] for x in grid["full_sampling_counts"]),
       "exact_latin_name_matched":len(records),
       "unmatched_names":missing,"ambiguous_names":ambiguous,
       "qualified_estimator_species_min10":len(stricter),
       "heldout_species_by_role":[len(v) for v in roles],
       "cross_role_species_overlap":len(roles[0]&roles[1]),
       "qualified_species_records":stricter,
       "species_metadata_auc_values_opened":False,
       "genetic_break_locations_exported":False,
       "empirical_model_run":False}
def main():
    p=argparse.ArgumentParser()
    for k in ("manifest","grid","out","work"):p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    m=json.loads(a.manifest.read_text())
    g=json.loads(a.grid.read_text())
    if g["status"]!="PASS_SAMPLED_BREAK_BACKGROUND_SUPPORT":raise ValueError("not spatial qualified")
    e=next(x for x in m["files"] if x["role"]=="species_metadata")
    a.work.mkdir(parents=True,exist_ok=True)
    file=a.work/"original_species.xlsx";fetch(e,file)
    result=join(excel_common_scientific(file),g)
    file.unlink(missing_ok=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="qualified_species_records"},sort_keys=True))
if __name__=="__main__":main()
