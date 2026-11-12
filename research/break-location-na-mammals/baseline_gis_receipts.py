"""Verify present availability, ZIP structure, and SHA for external barrier GIS.
No genetic break geometries, species labels or outcome values are loaded.
The publisher's source ZIPs are NOT uploaded or redistributed.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,urllib.request,urllib.error,zipfile
from pathlib import Path

HOSTS={"dmap-prod-oms-edc.s3.us-east-1.amazonaws.com",
       "naturalearth.s3.amazonaws.com",
       "datadownload-production.s3.amazonaws.com"}
def obtain(item,dir):
    from urllib.parse import urlparse
    url=item["url"]; host=urlparse(url).hostname
    if urlparse(url).scheme!="https" or host not in HOSTS:
        raise ValueError("Unexpected provider host")
    destination=dir/(item["term"]+".zip")
    max_bytes=350_000_000
    try:
        request=urllib.request.Request(url,headers={"User-Agent":"TTF-independent-barrier-GIS-preflight/1.0"})
        total=0;hash=hashlib.sha256()
        with urllib.request.urlopen(request,timeout=100) as response,destination.open("wb") as f:
            http=int(response.status)
            for _ in range(max_bytes//(1024*1024)+2):
                block=response.read(1024*1024)
                if not block:break
                total+=len(block)
                if total>max_bytes:raise ValueError("asset too large for preflight")
                f.write(block);hash.update(block)
            else:raise ValueError("stream exceeded bounded size")
        if not zipfile.is_zipfile(destination):raise ValueError("source not a ZIP")
        with zipfile.ZipFile(destination) as z:
            members=[x for x in z.infolist() if not x.is_dir()]
            names=[f.filename for f in members]
            if any(n.startswith("/") or ".." in Path(n).parts for n in names):
                raise ValueError("unsafe ZIP source")
            extensions=dict(collections.Counter(Path(n).suffix.lower() for n in names))
            result={"term":item["term"],"provider":item["provider"],
                "source_url":url,"status":"RETRIEVED_ORIGINAL_PROVIDER_CANDIDATE",
                "http_status":http,"download_bytes":total,
                "sha256":hash.hexdigest(),
                "member_count":len(members),
                "extensions":extensions,"member_names":names[:120],
                "source_rights":item["rights"],
                "genetic_break_or_background_coordinates_opened":False}
    except Exception as exc:
        result={"term":item["term"],"source_url":url,"status":"SOURCE_RETRIEVAL_FAILED",
                "error_type":type(exc).__name__,
                "error":str(exc)[:250],"genetic_break_or_background_coordinates_opened":False}
    finally:
        destination.unlink(missing_ok=True)
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    data=json.loads(args.manifest.read_text())
    if data["status"]!="PRE_RESPONSE_CANDIDATE_URLS_REQUIRING_SOURCE_SHA_AND_SCHEMA":
        raise ValueError("Unfrozen source candidates")
    args.work.mkdir(parents=True,exist_ok=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    receipts=[]
    for e in data["source_candidates"]:
        q=obtain(e,args.work)
        receipts.append(q)
        print("BARRIER GIS SOURCE",q["term"],q["status"],q.get("download_bytes"),flush=True)
    complete=all(x["status"]=="RETRIEVED_ORIGINAL_PROVIDER_CANDIDATE" for x in receipts)
    result={"schema":"ttf_jensen_barrier_gis_source_receipts_v01",
            "status":"ALL_FIVE_PROVIDER_CANDIDATES_VERIFIED" if complete else "PARTIAL_SOURCE_AVAILABILITY",
            "verified_candidate_count":sum(x["status"]=="RETRIEVED_ORIGINAL_PROVIDER_CANDIDATE" for x in receipts),
            "candidate_count":len(receipts),
            "all_required_original_baseline_variables_ready":False,
            "earthenv_tri_1km_source_pending":True,
            "extra_waterway_source_pending":True,
            "receipts":receipts,
            "genetic_break_location_response_tested":False}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],
        "verified":result["verified_candidate_count"],
        "needs_TRI":True,"outcome_tested":False}))
if __name__=="__main__":main()
