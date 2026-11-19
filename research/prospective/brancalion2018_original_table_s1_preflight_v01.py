"""Brancalion 2018 Table S1 public-original supplement source-only preflight.

Official-source routes only. If acquired, read only OOXML text from the
environmental Table S1 to discover whether the archive carries exact source
GPS/maternal IDs. No fitness, genotype, LGM or other response is analyzed.
A protected-area coordinate does NOT qualify as original source-tree GPS.
"""
from __future__ import annotations
import argparse, hashlib, io, json, re, tarfile, xml.etree.ElementTree as ET, zipfile
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlparse
from urllib.request import Request,urlopen

DOI="10.1002/ece3.4248"
PMCID="PMC6106193"
SUPPLEMENT="ece34248-sup-0003-TableS1.docx"
MAX_DOCX=2_000_000
MAX_ARCHIVE=40_000_000
PUBLISHER="https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fece3.4248&file="+SUPPLEMENT
PMC_DIRECT=[
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC6106193/bin/"+SUPPLEMENT,
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/6106193/bin/"+SUPPLEMENT,
]
PMC_OA="https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id="+PMCID
ALLOWED={"onlinelibrary.wiley.com","pmc.ncbi.nlm.nih.gov",
         "www.ncbi.nlm.nih.gov","ftp.ncbi.nlm.nih.gov"}
SITE_NAMES=("Caetetus","Carlos Botelho","Ilha do Cardoso")
TOKEN_RE=re.compile(r"\b(?:latitude|longitude|coordinates?|GPS|mother|maternal|source population|seed lot)\b",re.I)
COORDINATE_RE=re.compile(r"(?<!\d)[+-]?\d{1,2}\.\d{3,}\s*[,;/ ]\s*[+-]?\d{1,3}\.\d{3,}(?!\d)")


def fetch(url,max_bytes):
    if urlparse(url).scheme!="https" or urlparse(url).hostname not in ALLOWED:
        raise ValueError("Source URL is not an approved official HTTPS route")
    req=Request(url,headers={"Accept":"application/octet-stream, application/xml, */*",
                            "User-Agent":"TTF-original-source-study-audit/2026-11-19"})
    with urlopen(req,timeout=22) as res:
        if urlparse(res.geturl()).hostname not in ALLOWED:
            raise ValueError("Unexpected redirect away from an official source")
        data=res.read(max_bytes+1)
        if len(data)>max_bytes:raise ValueError("Official file exceeds source-only size cap")
        return data


def oai_archive_url(xml_bytes):
    root=ET.fromstring(xml_bytes)
    record=root.find(".//record")
    if record is None or record.get("id")!=PMCID:
        return None
    for node in record.findall(".//link"):
        if node.get("format") in ("tgz","tar.gz"):
            url=node.get("href","")
            if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
                return "https://"+url[len("ftp://"):]
            if url.startswith("https://ftp.ncbi.nlm.nih.gov/"):
                return url
    return None


def docx_metadata(raw):
    if len(raw)>MAX_DOCX:raise ValueError("DOCX too large")
    if not raw.startswith(b"PK\x03\x04"):raise ValueError("Not a DOCX ZIP")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        members=set(z.namelist())
        if "[Content_Types].xml" not in members or "word/document.xml" not in members:
            raise ValueError("Not an Office Open XML document")
        entry=z.getinfo("word/document.xml")
        if entry.file_size>400_000:raise ValueError("OOXML inflation size too large")
        xml=z.read(entry)
    root=ET.fromstring(xml)
    ns="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    text=" ".join((node.text or "") for node in root.iter(ns+"t"))
    # Only disclose high-level source-field inventory, not the table contents.
    # The existence of numbers/coordinates is not source-GPS certification.
    return {
        "sha256":hashlib.sha256(raw).hexdigest(),
        "docx_bytes":len(raw),
        "original_site_names_mentioned":[s for s in SITE_NAMES if s.casefold() in text.casefold()],
        "geographic_or_parentage_field_keyword_count":len(TOKEN_RE.findall(text)),
        "coordinate_shaped_text_count":len(COORDINATE_RE.findall(text)),
        "document_text_extracted_to_receipt":False,
        "fitness_or_genetics_or_LGM_response_values_exported":False,
    }


def select_from_oa_archive(raw):
    with tarfile.open(fileobj=io.BytesIO(raw),mode="r:gz") as z:
        candidates=[m for m in z.getmembers() if m.isfile()
                    and (m.name.endswith(SUPPLEMENT)
                         or (m.name.lower().endswith(".docx")
                             and ("s001" in m.name.lower() or "tables1" in m.name.lower())))]
        if len(candidates)!=1:
            raise ValueError("OA archive does not uniquely identify original Table S1")
        file=candidates[0]
        if file.size>MAX_DOCX:raise ValueError("Supplement larger than source-only cap")
        return z.extractfile(file).read(MAX_DOCX+1)


def probe(fetcher=fetch):
    attempts=[]
    verified=None
    for route,url in [("PUBLISHER_OFFICIAL",PUBLISHER)]+[
        ("PMC_DIRECT_"+str(i),u) for i,u in enumerate(PMC_DIRECT,1)]:
        attempt={"route":route}
        try:
            raw=fetcher(url,MAX_DOCX)
            verified=docx_metadata(raw)
            attempt["status"]="ORIGINAL_TABLE_S1_DOCX_AUTHENTIC_STORAGE_CONFIRMED"
            attempts.append(attempt)
            break
        except (ValueError,HTTPError,URLError,TimeoutError,ET.ParseError,
                zipfile.BadZipFile) as e:
            attempt.update({"status":"ACCESS_OR_STORAGE_NOT_CONFIRMED",
                            "error_type":type(e).__name__,
                            "http_code":getattr(e,"code",None)})
        attempts.append(attempt)
    if verified is None:
        attempt={"route":"PMC_OFFICIAL_OA_PACKAGE"}
        try:
            archive_url=oai_archive_url(fetcher(PMC_OA,100_000))
            if archive_url is None:raise ValueError("No matching PMC OA package record")
            raw=select_from_oa_archive(fetcher(archive_url,MAX_ARCHIVE))
            verified=docx_metadata(raw)
            attempt["status"]="ORIGINAL_TABLE_S1_DOCX_AUTHENTIC_STORAGE_CONFIRMED"
        except (ValueError,HTTPError,URLError,TimeoutError,ET.ParseError,
                tarfile.TarError,zipfile.BadZipFile) as e:
            attempt.update({"status":"ACCESS_OR_STORAGE_NOT_CONFIRMED",
                            "error_type":type(e).__name__,
                            "http_code":getattr(e,"code",None)})
        attempts.append(attempt)
    return {
      "schema":"ttf_Brancalion2018_original_TableS1_source_only_preflight_v01",
      "original_doi":DOI,"pmcid":PMCID,"published_supplement_locator":SUPPLEMENT,
      "source_only_table_scope":"author environmental site Table S1, not Table S2 GO or fitness",
      "status":"ORIGINAL_TABLE_S1_STORAGE_CONFIRMED_SOURCE_GPS_NOT_CERTIFIED" if verified
               else "ORIGINAL_TABLE_S1_BYTES_NOT_RECOVERED",
      "attempts":attempts,"verified_docx_metadata":verified,
      "independent_original_maternal_tree_GPS_certified":0,
      "independent_source_population_GPS_certified":0,
      "cross_assay_maternal_ID_join_completed":False,
      "historical_climate_contrast_qualified":False,
      "fitness_SNP_LGM_adaptation_plasticity_speciation_response_opened":False,
      "prior_Bontrager_44_of_50_stop_unchanged":True,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    result=probe()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("BRANCALION2018_TABLE_S1_PREFLIGHT="+json.dumps(result,sort_keys=True,ensure_ascii=False))

if __name__=="__main__":
    main()
