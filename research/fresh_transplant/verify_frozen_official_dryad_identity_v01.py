"""Fail-closed comparison of live Dryad API file metadata to frozen first-shot receipt.

Authenticates DOI, published version ID, file path, file ID, byte size and
repository-DECLARED SHA256. Does not hash downloaded original file bytes,
read source GPS rows, or use fitness/genetic outcomes.
"""
from __future__ import annotations
import argparse,json,re
from pathlib import Path

HERE=Path(__file__).resolve().parent
RECEIPT=HERE/"source_registry/fresh_dryad_original_file_identity_v01.json"

def verify(receipt,live):
    if receipt.get("schema")!="ttf_fresh_Dryad_official_original_file_identity_frozen_v01":
        raise ValueError("Not a frozen original source receipt")
    if live.get("schema")!="ttf_fresh_official_Dryad_DOI_file_manifest_v01":
        raise ValueError("Live official manifest schema changed")
    if not receipt["stopped_bontrager_44of50_unchanged"]:
        raise ValueError("Frozen prior stop changed")
    if (receipt["source_GPS_rows_downloaded_or_verified"]!=0 or
        receipt["original_fitness_or_genetic_responses_read"] or
        receipt["history_21ka_contrast_tested"]):
        raise ValueError("Original source-only receipt boundary changed")
    records=receipt["datasets"]
    observed=live["studies"]
    if len(records)!=3 or len(observed)!=3:
        raise ValueError("Study family changed")
    checked=0
    for frozen,now in zip(records,observed):
        if (frozen["study_id"],frozen["doi"],frozen["version_id"])!=(
            now["study_id"],now["dataset_doi"],now["latest_public_version_id"]):
            raise ValueError("Original published study/version identity drift")
        if now["status"]!="ALL_PINNED_FILES_REGISTERED_IN_OFFICIAL_METADATA":
            raise ValueError("Original repository file manifest not fully accessible")
        got=now["files_metadata"]
        want=frozen["files"]
        if len(got)!=len(want):raise ValueError("Original source file family length drift")
        for expected,actual in zip(want,got):
            if not re.fullmatch(r"[0-9a-f]{64}",expected["sha256_from_repository_metadata"]):
                raise ValueError("Malformed repository-reported SHA256")
            if actual["status"]!="ORIGINAL_FILE_METADATA_MATCHES_FROZEN_SOURCE":
                raise ValueError("Original file identity not qualified")
            if (expected["path"],expected["id"],expected["size_bytes"],
                expected["sha256_from_repository_metadata"])!=(
                    actual["expected_original_path"],actual["official_file_record_id"],
                    actual["official_file_size_bytes"],actual["official_file_digest"]):
                raise ValueError("Official file identity/digest drift")
            if actual["official_file_digest_type"].lower() not in ("sha-256","sha256"):
                raise ValueError("Unexpected digest format")
            if not actual["file_id_matches_frozen_page"]:
                raise ValueError("Original dataset public file link changed")
            checked+=1
    if (live["fitness_or_genetics_response_opened"] or
        live["source_GPS_independently_certified"]!=0 or
        live["LGM21ka_contrast_qualified"]):
        raise ValueError("Data response unexpectedly promoted")
    return {
      "status":"OFFICIAL_ORIGINAL_METADATA_MATCHED_FROZEN_FIRSTSHOT",
      "original_dataset_count":3,
      "original_file_count":checked,
      "frozen_receipt_provenance":"OFFICIAL_METADATA_DIGEST_ONLY_NOT_FILE_BYTE_HASH",
      "new_independent_GPS_populations_certified":0,
      "new_21ka_LGM_fitness_effect_tested":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--observed",required=True,type=Path)
    p.add_argument("--receipt",default=RECEIPT,type=Path)
    a=p.parse_args()
    result=verify(json.loads(a.receipt.read_text("utf-8")),
                  json.loads(a.observed.read_text("utf-8")))
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":
    main()
