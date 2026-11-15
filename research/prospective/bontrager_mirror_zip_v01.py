"""Cross-check publisher Dryad and UBC mirror original ZIP: ZIP NAMES ONLY.

Metadata/byte provenance only; no extracted CSV rows, fitness values,
origin coordinates, source study IDs or genetic responses inspected.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from pathlib import Path
from urllib.request import urlopen,Request

META_URL="https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP2%2FDU3ZTQ"
ORIGINAL_URL="https://borealisdata.ca/api/access/datafile/150984?format=original"

def verify_original_zip(raw,lock):
    if (len(raw)!=lock["expected_file_bytes"] or
        hashlib.sha256(raw).hexdigest()!=lock["dryad_sha256"] or
        hashlib.md5(raw).hexdigest()!=lock["mirror_md5"]):
        raise ValueError("Original archive not identical to both published source digests")
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("Original source is not a ZIP")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:names=z.namelist()
    if len(names)>2000 or any(n.startswith("/") or ".." in Path(n).parts for n in names):
        raise ValueError("Unsafe original archive member names")
    return {"member_names":names,"member_count":len(names),
            "sha256_verified_against_Dryad":lock["dryad_sha256"],
            "md5_verified_against_Borealis":lock["mirror_md5"],
            "zip_member_contents_read":False,
            "fitness_or_site_records_read":False}

def run():
    p=argparse.ArgumentParser()
    p.add_argument("--contract",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args();lock=json.loads(args.contract.read_text())
    if lock.get("status")!="FROZEN_BEFORE_UBC_BOREALIS_ORIGINAL_ZIP_CONTENT":
        raise ValueError("Mirror source identity not frozen")
    receipt={"schema":"ttf_Bontrager_original_Dryad_Borealis_zip_identity_members_v01",
             "fitness_or_site_records_read":False,"biological_model_fitted":False,
             "independent_new_studies_qualified":False,"old_Haldane_opened":False}
    try:
        q=Request(META_URL,headers={"Accept":"application/json",
                                    "User-Agent":"TTF-original-transplant-Dataverse-mirror-source-only"})
        with urlopen(q,timeout=45) as r:obj=json.load(r)
        meta=obj["data"];version=meta["latestVersion"]
        if (obj.get("status")!="OK" or meta.get("id")!=lock["mirror_dataset_id"]
            or version.get("datasetPersistentId")!=lock["mirror_dataset_doi"]
            or (version.get("versionNumber"),version.get("versionMinorNumber"))!=(1,0)):
            raise ValueError("Exact UBC source version/DOI mismatch")
        entries=[z for z in version["files"] if z.get("dataFile",{}).get("id")==lock["mirror_source_file_id"]]
        if len(entries)!=1 or entries[0].get("restricted"):
            raise ValueError("Original mirrored archive unavailable or restricted")
        f=entries[0]["dataFile"]
        if (f.get("filename")!=lock["file_name"] or
            f.get("filesize")!=lock["expected_file_bytes"] or
            f.get("md5")!=lock["mirror_md5"]):
            raise ValueError("Original mirrored archive manifest changed")
        q=Request(ORIGINAL_URL,headers={"User-Agent":"TTF-original-ZIP-central-directory-only"})
        with urlopen(q,timeout=90) as r:raw=r.read(lock["max_download_bytes"]+1)
        if len(raw)>lock["max_download_bytes"]:
            raise ValueError("Source download exceeds capped size")
        receipt.update(verify_original_zip(raw,lock))
        receipt["status"]="IDENTICAL_DRYAD_BOREALIS_ORIGINAL_ZIP_MEMBER_INDEX_VERIFIED"
    except Exception as e:
        receipt.update(status="MIRROR_ORIGINAL_ZIP_UNAVAILABLE_NOT_BIOLOGICAL_NULL",
                       error_type=type(e).__name__,error_message=str(e)[:280])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(receipt["status"],receipt.get("member_count"))
if __name__=="__main__":run()
