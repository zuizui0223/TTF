"""Read exact published Jensen analysis script to inventory original barrier covariates.
No mapped break coordinates, feature outcomes, or phenotypic effects are read.
"""
import argparse, collections, hashlib, json, re
from pathlib import Path
from source_integrity import fetch

MATCH = re.compile(r"mountain|water|river|lake|ecoreg|ecotone|rugged|terrain|bio[0-9]|worldclim|chelsa|raste?r|geod|dist|st_read|read_?sf|readRDS|^library\(|^source\(|^load\(|^read\.",re.I)
SOURCE_OPS = re.compile(r"st_read\s*\(|readRDS\s*\(|raster\s*\(|rast\s*\(|read_?sf\s*\(|read\.csv\s*\(|load\s*\(|read\.table\s*\(",re.I)
ASSET_PATTERN = re.compile(r"""(?i)([^\s"'()]+?\.(?:tif|tiff|shp|gpkg|rds|rda|csv|txt|gdb|zip))""")

def scan(text):
    lines=text.splitlines()
    excerpts=[{"line":i,"text":s[:300]} for i,s in enumerate(lines,1) if MATCH.search(s)]
    reads=[{"line":i,"text":s[:300]} for i,s in enumerate(lines,1) if SOURCE_OPS.search(s)]
    assets=sorted(set(v for s in lines for v in ASSET_PATTERN.findall(s)))
    return {"total_lines":len(lines),"input_reader_context":[{"line":i,"text":s[:400]} for i,s in enumerate(lines,1) if 35<=i<=112], "candidate_environmental_lines":excerpts,
            "file_read_calls":reads,"mentioned_file_tokens":assets,
            "response_geometry_or_genetic_outcomes_opened":False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--work",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    meta=json.loads(args.manifest.read_text())
    src=next(f for f in meta["files"] if f["role"]=="original_analysis_code")
    args.work.mkdir(parents=True,exist_ok=True)
    file=args.work/"original_verified_analysis.R"
    receipt=fetch(src,file)
    body=file.read_text(encoding="utf-8-sig",errors="replace")
    analysis=scan(body)
    file.unlink()
    result={"schema":"ttf_jensen_original_barrier_sources_and_script_audit_v01",
            "status":"ORIGINAL_CODE_SHA256_VERIFIED; NO_EMPIRICAL_BREAK_EFFECT",
            "verified_md5":receipt["md5"],"verified_sha256":receipt["sha256"],
            "source_file_id":src["id"],**analysis}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"total_lines":analysis["total_lines"],
        "geographic_source_matches":len(analysis["candidate_environmental_lines"]),
        "read_calls":len(analysis["file_read_calls"]),
        "file_tokens":analysis["mentioned_file_tokens"][:40],
        "genetic_response_opened":False},sort_keys=True))
if __name__=="__main__":main()
