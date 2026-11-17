"""Exact-rational rank audit of frozen Ellis source/garden comparison topology.

Source-only: no original fitness, AFLP, coordinates or historical climate.
"""
from fractions import Fraction
from itertools import product
import hashlib, json, argparse
from pathlib import Path

HERE=Path(__file__).resolve().parent
RECEIPT=HERE/"source_registry/ellis2007_24_to_18_origin_garden_arm_audit_v02.json"
SHA="28cea27d59affd846c9ee490f20daedd36220e6e75041f429663a208dcd239ce"
CLASSES=("local_same_source_and_garden","same_basin_other_quartz",
         "other_basin_same_quartz","other_basin_other_quartz_not_in_paper_design")

def rank(matrix):
    a=[[Fraction(x) for x in row] for row in matrix]
    if not a:return 0
    if any(len(row)!=len(a[0]) for row in a):raise ValueError("ragged")
    r=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        v=a[r][col]
        a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][col]:
                v=a[i][col]
                a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

def category(s,g):
    if s==g:return CLASSES[0]
    if s[0]==g[0]:return CLASSES[1]
    if s[1]==g[1]:return CLASSES[2]
    return CLASSES[3]

def analyze(data):
    if data.get("schema")!="ttf_Ellis2007_primary_24_source_garden_vs_Hargreaves_18_source_only_edges_diagnosis_v02":
        raise ValueError("wrong frozen receipt")
    for key in ("no_fitness_growth_survival_reproduction_or_AFLP_response_rows_read",
                "no_LGM_21ka_0BP_climate_effect_estimated",
                "old_Bontrager_44of50_terminal_stop_unchanged"):
        if data.get(key) is not True:raise ValueError("frozen outcome guard changed")
    rows=data["published_origin_site_pairs_by_class"]
    if len(rows)!=4 or len({r["class"] for r in rows})!=4:raise ValueError("duplicate classes")
    by={r["class"]:r for r in rows}
    if set(by)!=set(CLASSES):raise ValueError("unknown class")
    if ([by[k]["paper_expected"] for k in CLASSES]!=[6,6,12,0]
        or [by[k]["original_Hargreaves_extract_observed"] for k in CLASSES]!=[6,0,12,0]):
        raise ValueError("class topology drift")
    nodes=list(product(range(3),range(2)))
    primary=[(s,g) for s in nodes for g in nodes if category(s,g) in CLASSES[:3]]
    secondary=[(s,g) for s in nodes for g in nodes if category(s,g) in (CLASSES[0],CLASSES[2])]
    if (len(primary),len(secondary))!=(data["original_primary_expected_directed_source_garden_pairs"],
                                        data["observed_Hargreaves_distinct_directed_source_garden_pairs"]):
        raise ValueError("total drift")
    def row(s,g):
        return [1,s[1],g[1],int(s[0]==1),int(s[0]==2),int(g[0]==1),int(g[0]==2)]
    p=rank([row(s,g) for s,g in primary])
    q=rank([row(s,g) for s,g in secondary])
    if (p,q)!=(7,6):raise AssertionError("rank changed")
    if not all(s[1]==g[1] for s,g in secondary):raise AssertionError("quartz join changed")
    return {
      "schema":"ttf_Ellis2007_source_only_design_rank_v01",
      "input_sha256":SHA,
      "status":"SOURCE_ONLY_ALGEBRAIC_IDENTIFIABILITY_NOT_BIOLOGICAL_EFFECT",
      "design_columns":["intercept","source_quartz","garden_quartz","source_basin_2","source_basin_3","garden_basin_2","garden_basin_3"],
      "original_directed_arms":24,"secondary_directed_arms":18,
      "original_design_rank":p,"secondary_design_rank":q,
      "original_full_column_rank":True,"secondary_full_column_rank":False,
      "source_vs_garden_quartz_separately_identifiable_original":True,
      "source_vs_garden_quartz_separately_identifiable_secondary":False,
      "secondary_source_quartz_equals_garden_quartz_on_all_observed_arms":True,
      "scope":"linear additive fixed-effect design with separate source/garden basin and quartz, not fitted fitness effects",
      "fitness_genetics_climate_or_LGM_response_opened":False,
      "main_or_frozen_endpoints_modified":False
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    raw=RECEIPT.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SHA:raise ValueError("receipt SHA drift")
    result=analyze(json.loads(raw))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))

if __name__=="__main__":main()
