"""New Penstemon 239-species nuclear phylogeny Dryad source-only audit.

Only official version/file metadata and, when public file bytes are accessible,
small RAxML/ASTRAL archive *member names*. No sequence alignments, tree
topologies, flower traits, model results or genetic outcomes are decoded.
"""
from __future__ import annotations
import hashlib
import io
import json
import urllib.error
import urllib.parse
import urllib.request
import zipfile

DOI="10.5061/dryad.g79cnp66z"
EXPECTED_SOURCE_DATE="2026-11-02"
TREE_FILES={"2.RAxML.zip","3.Astral.zip"}
MAX_FILE_BYTES=1_500_000
BASE="https://datadryad.org/api/v2"


def json_api(url):
    req=urllib.request.Request(url,headers={
       "Accept":"application/json","User-Agent":"TTF-Penstemon-public-source-audit/1.0",
    })
    with urllib.request.urlopen(req,timeout=50) as response:
        raw=response.read(2_000_001)
    if len(raw)>2_000_000:raise ValueError("Oversized source metadata JSON")
    return json.loads(raw)


def embedded(payload,needle):
    for key,value in (payload.get("_embedded") or {}).items():
        if needle in key.lower() and isinstance(value,list):return value
    raise ValueError(f"source HAL collection {needle} missing")


def href(entry):
    x=(entry.get("_links") or {}).get("self")
    return x if isinstance(x,str) else x.get("href") if isinstance(x,dict) else ""


def obj_id(entry):
    if isinstance(entry.get("id"),int):return entry["id"]
    tail=href(entry).rsplit("/",1)[-1]
    return int(tail) if tail.isdigit() else None


def version_rows(d):
    result=[]
    for x in embedded(d,"versions"):
        ident=obj_id(x)
        if ident is None:raise ValueError("source version id missing")
        result.append({"id":ident,"version":x.get("versionNumber"),
                       "published":x.get("publicationDate"),"updated":x.get("lastModificationDate")})
    if not result:raise ValueError("no source versions")
    return result


def file_rows(d):
    result=[]
    for x in embedded(d,"files"):
        ident=obj_id(x)
        if ident is None:raise ValueError("source file id missing")
        result.append({"id":ident,"name":x.get("path") or x.get("fileName"),
                       "size":x.get("size"),"digest_type":x.get("digestType"),
                       "sha256":x.get("digest"),"self":href(x)})
    if not result:raise ValueError("no source files")
    return result


def archive_member_access(f):
    size=int(f.get("size") or -1)
    if size<100 or size>MAX_FILE_BYTES:
        return {"status":"SKIPPED_UNSAFE_OR_TOO_LARGE"}
    url=f"{BASE}/files/{f['id']}/download"
    req=urllib.request.Request(url,headers={"User-Agent":"TTF-Penstemon-source-only/1.0"})
    try:
        with urllib.request.urlopen(req,timeout=45) as response:
            raw=response.read(size+1)
    except urllib.error.HTTPError as e:
        if e.code in (401,403,404):
            return {"status":"SOURCE_FILE_ACCESS_BLOCKED","http_status":e.code,"url":url}
        raise
    if len(raw)!=size or hashlib.sha256(raw).hexdigest()!=f["sha256"]:
        raise ValueError(f"source file {f['name']} size/SHA256 mismatch")
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("original tree archive not ZIP")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        infos=[{"path":x.filename,"size":x.file_size} for x in archive.infolist() if not x.is_dir()]
        if len(infos)>500:raise ValueError("too many archived source members")
    return {"status":"ORIGINAL_TREE_ARCHIVE_BYTES_VERIFIED_MEMBERS_ONLY",
            "member_count":len(infos),"member_paths":infos,
            "SHA256_verified":True,"data_sequences_or_tree_nodes_read":False}


def main():
    encoded=urllib.parse.quote("doi:"+DOI,safe="")
    v=version_rows(json_api(f"{BASE}/datasets/{encoded}/versions"))
    selected=max(v,key=lambda x:(int(x["version"] or 0),x["id"]))
    files=file_rows(json_api(f"{BASE}/versions/{selected['id']}/files"))
    present={f["name"] for f in files}
    if not TREE_FILES.issubset(present):
        raise ValueError(f"original RAxML/ASTRAL archives absent: {TREE_FILES-present}")
    result={
      "status":"PENSTEMON_239SPP_NUCLEAR_SOURCE_METADATA_AND_TREE_ACCESS_ONLY",
      "doi":DOI,"public_release_date_paper":EXPECTED_SOURCE_DATE,
      "version":selected,"versions":v,
      "files":files,
      "small_tree_archive_access":{
        f["name"]:archive_member_access(f) for f in files if f["name"] in TREE_FILES
      },
      "no_genetic_or_floral_response_values_read":True,
      "no_sequences_or_tree_topologies_decoded":True,
      "not_a_phenotype_phylogeny_join":True,
    }
    print(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
