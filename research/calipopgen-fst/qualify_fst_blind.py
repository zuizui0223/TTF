"""CaliPopGen Fst eligibility: read geographical metadata + GenDist presence, NEVER values."""
import argparse,collections,hashlib,io,json,math,re,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
from inspect_xlsx_schema import exact_source_download

S="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
NEEDED={"CitationID","SpeciesID","Pop1ScientificName","Pop2ScientificName",
 "Pop1LatitudeDD","Pop1LongitudeDD","Pop2LatitudeDD","Pop2LongitudeDD",
 "CoordError","GenDist","GenDistMetric","MarkerType","GeneTarget",
 "sepAnalyses","specialComparisonType","HabitatType"}

def colnum(ref):
    n=0
    for ch in ref:
        if not ch.isalpha():break
        n=n*26+ord(ch.upper())-64
    return n-1
def sample_metadata_rows(xlsx):
    with zipfile.ZipFile(io.BytesIO(xlsx)) as z:
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            with z.open("xl/sharedStrings.xml") as f:
                for _,node in ET.iterparse(f,events=("end",)):
                    if node.tag==S+"si":
                        strings.append("".join(x.text or "" for x in node.iter(S+"t")))
                        node.clear()
        with z.open("xl/worksheets/sheet1.xml") as f:
            columns=None
            for _,row in ET.iterparse(f,events=("end",)):
                if row.tag!=S+"row":continue
                record={}
                for cell in row.findall(S+"c"):
                    j=colnum(cell.get("r",""))
                    if columns is not None and columns.get(j) not in NEEDED:continue
                    value=cell.find(S+"v")
                    typ=cell.get("t")
                    name=columns.get(j) if columns is not None else None
                    if name=="GenDist":
                        # Outcome values NEVER decoded. Data may be stored as shared strings;
                        # test cell presence without decoding its genetic value.
                        record[j]=(value is not None and typ!="e")
                        continue
                    if typ=="s" and value is not None:
                        txt=strings[int(value.text)]
                    elif typ=="inlineStr":
                        txt="".join(e.text or "" for e in cell.iter(S+"t"))
                    else:txt=value.text if value is not None else None
                    record[j]=txt
                if columns is None:
                    columns={j:str(value or "").strip() for j,value in record.items()}
                    if not NEEDED.issubset(set(columns.values())):
                        raise ValueError("CaliPopGen schema changed: "+str(NEEDED-set(columns.values())))
                else:yield {name:record.get(j) for j,name in columns.items() if name in NEEDED}
                row.clear()

def as_float(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError,OverflowError):return None
def gc_km(a,b):
    p1,l1,p2,l2=map(math.radians,(*a,*b))
    z=math.sin((p2-p1)/2)**2+math.cos(p1)*math.cos(p2)*math.sin((l2-l1)/2)**2
    return 12742.0176*math.asin(min(1.,math.sqrt(z)))
def frozen_binomial(v):
    words=str(v or "").strip().split()
    if len(words)<2:return None
    if words[1].lower() in ("sp.","sp","spp.","cf.","aff.") or not words[0][0].isalpha():return None
    return (words[0]+" "+words[1]).lower()

def evaluate(xlsx):
    c=collections.Counter()
    groups=collections.defaultdict(lambda:{"pairs":set(),"sites":set(),"missing_coord_error":0,"known_coord_error":0})
    markers=collections.Counter()
    metriclabels=collections.Counter()
    habitatlabels=collections.Counter()
    speciallabels=collections.Counter()
    for r in sample_metadata_rows(xlsx):
        c["total_rows"]+=1
        metric=str(r.get("GenDistMetric") or "").strip()
        metriclabels[metric]+=1
        if re.sub(r"[^a-z]","",metric.lower())!="fst":
            c["non_Fst_metric"]+=1;continue
        c["rows_fst_label"]+=1
        if not r.get("GenDist"):
            c["no_numeric_type_response_cell"]+=1;continue
        habitatlabels[str(r.get("HabitatType") or "").strip()]+=1
        if str(r.get("HabitatType") or "").strip().lower()!="terrestrial":
            c["not_terrestrial"]+=1;continue
        speciallabels[str(r.get("specialComparisonType") or "").strip()]+=1
        if str(r.get("specialComparisonType") or "").strip().lower() not in ("","none","na","nan"):
            c["special_temporal_hybrid_or_interspecific"]+=1;continue
        sp=frozen_binomial(r.get("Pop1ScientificName"))
        target=frozen_binomial(r.get("Pop2ScientificName"))
        if not sp or not target or sp!=target:
            c["missing_taxon_or_interspecific"]+=1;continue
        a=tuple(as_float(r.get(k)) for k in ("Pop1LatitudeDD","Pop1LongitudeDD"))
        b=tuple(as_float(r.get(k)) for k in ("Pop2LatitudeDD","Pop2LongitudeDD"))
        if any(x is None for x in (*a,*b)):
            c["missing_geolocation"]+=1;continue
        if any(not(-90<=lat<=90 and -180<=lon<=180) or (lat==0 and lon==0) for lat,lon in (a,b)):
            c["invalid_geolocation"]+=1;continue
        if a==b:
            c["zero_distance_pair"]+=1;continue
        err=as_float(r.get("CoordError"))
        if err is not None and (err<0 or err>25):
            c["coordinate_error_above25km"]+=1;continue
        study=str(r.get("CitationID") or "").strip()
        if not study:
            c["missing_study"]+=1;continue
        marker=str(r.get("MarkerType") or "").strip()
        gene=str(r.get("GeneTarget") or "").strip()
        analysis=str(r.get("sepAnalyses") or "").strip()
        if not marker:
            c["missing_marker"]+=1;continue
        if gene.lower()=="none":gene=""
        key=(study,sp,marker,gene,analysis)
        gp=groups[key]
        pair=tuple(sorted((a,b)))
        if pair in gp["pairs"]:
            c["repeat_unordered_pair_in_group"]+=1;continue
        gp["pairs"].add(pair)
        gp["sites"].update((a,b))
        gp["missing_coord_error" if err is None else "known_coord_error"]+=1
        markers[marker]+=1
        c["prequalified_pair_rows"]+=1
    accepted=[]
    for (study,sp,marker,gene,analysis),g in sorted(groups.items()):
        p=g["sites"]
        if len(g["pairs"])<5 or len(p)<4:
            c["groups_too_few_pairs_or_sites"]+=1;continue
        # Reject groups restricted to a tiny sampled region before response opening.
        if max(gc_km(a,b) for a,b in g["pairs"])<100:
            c["groups_extent_under100km"]+=1;continue
        parity=hashlib.sha256(("calipopgen-fst-ttf-fresh-v01|"+sp).encode()).digest()[0]%2
        accepted.append({"study_id":study,"species":sp,"marker":marker,
          "gene_target":gene,"sep_analyses":analysis,"role":parity,
          "pair_count":len(g["pairs"]),"location_count":len(p),
          "missing_coord_error_pairs":g["missing_coord_error"],
          "known_coord_error_pairs":g["known_coord_error"]})
    roles=[{g["species"] for g in accepted if g["role"]==i} for i in (0,1)]
    per_role_groups=[sum(g["role"]==i for g in accepted) for i in (0,1)]
    eligible_species=len(roles[0]|roles[1])
    passed=eligible_species>=30 and len(accepted)>=35 and min(map(len,roles))>=10 and min(per_role_groups)>=12
    return {"schema":"ttf_calipopgen_fst_response_blind_support_v01",
       "status":"PASS_RESPONSE_BLIND_GEOGRAPHIC_SUPPORT" if passed else "NOT_EVALUABLE_RESPONSE_BLIND_SUPPORT",
       "source_rows_scanned":c["total_rows"],"counts":dict(c),
       "raw_metric_labels":dict(metriclabels),"raw_habitat_labels_in_Fst":dict(habitatlabels),
       "raw_special_comparison_labels_in_terrestrial_Fst":dict(speciallabels),
       "prequalified_markers":dict(markers),
       "eligible_species":eligible_species,"eligible_study_species_marker_groups":len(accepted),
       "species_by_hash_role":[len(x) for x in roles],
       "study_groups_by_hash_role":per_role_groups,
       "eligible_groups":accepted,
       "Fst_numeric_values_opened":False,
       "old_TTF_or_GenDivRange_responses_reopened":False}

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--lock",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    xlsx=exact_source_download(json.loads(args.lock.read_text()))
    result=evaluate(xlsx)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="eligible_groups"},sort_keys=True))
