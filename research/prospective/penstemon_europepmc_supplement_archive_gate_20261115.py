"""Official Europe PMC/NCBI supplemental source alternative index.

The publication Stone et al. 2026 (PMC13419788) has a supplemental Table S1,
but direct PMC /bin paths returned HTTP404 and bioRxiv HTTP429.
Query official alternative archival APIs only; list source attachment metadata,
SHA256 and original XML xlink hrefs; do NOT decode voucher, floral or SNP values.
"""
from __future__ import annotations
import hashlib
import io
import json
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ARTICLE_ID="PMC13419788"
SOURCE_ENDPOINTS={
 "EUROPE_PMC_SUPPLEMENT_ZIP":
  f"https://www.ebi.ac.uk/europepmc/webservices/rest/{ARTICLE_ID}/supplementaryFiles",
 "EUROPE_PMC_JATS_XML":
  f"https://www.ebi.ac.uk/europepmc/webservices/rest/{ARTICLE_ID}/fullTextXML",
 "NCBI_PMC_OA_MANIFEST":
  f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={ARTICLE_ID}",
}
MAX_BYTES=30_000_000
MAX_INDEX_ROWS=100


def xml_support_metadata(raw:bytes)->dict:
    if len(raw)>MAX_BYTES:raise ValueError("official XML exceeds declared size cap")
    root=ET.fromstring(raw)
    refs=[]
    for elem in root.iter():
        leaf=elem.tag.rsplit("}",1)[-1].lower()
        if "supplement" in leaf or leaf in ("link","file"):
            attrs={k.rsplit("}",1)[-1]:v for k,v in elem.attrib.items()}
            # Only source links and ids; do not export text/tables from paper.
            limited={k:attrs[k] for k in ("id","href","xlink:href","locator","format","type") if k in attrs}
            if limited and len(refs)<MAX_INDEX_ROWS:
                refs.append({"element":leaf,"refs":limited})
    return {
      "status":"ORIGINAL_XML_SOURCE_LINKS_ONLY",
      "root_tag":root.tag.rsplit("}",1)[-1],
      "source_bytes":len(raw),
      "sha256":hashlib.sha256(raw).hexdigest(),
      "reference_count":len(refs),
      "source_link_records":refs,
      "no_original_sample_or_phenotypic_values_read":True,
    }


def zip_support_metadata(raw:bytes)->dict:
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("official supplemental source package is not ZIP")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names=[
            {"path":i.filename,"uncompressed_size":i.file_size,
             "is_directory":i.is_dir()}
            for i in z.infolist()
        ]
    if len(names)>MAX_INDEX_ROWS:
        return {"status":"TOO_MANY_ZIP_MEMBERS_SOURCE_HOLD","count":len(names)}
    return {
       "status":"OFFICIAL_SUPPLEMENT_ZIP_MEMBER_INDEX_ONLY",
       "bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),
       "member_count":len(names),
       "member_metadata":names,
       "no_original_pdf_or_tabular_contents_read":True,
    }


def fetch_index(url,expected):
    req=urllib.request.Request(url,headers={
      "User-Agent":"TTF-Habroanthus-PMC-original-supplement-index/1.0",
      "Accept":"application/zip,application/xml,text/xml,*/*",
    })
    try:
        with urllib.request.urlopen(req,timeout=50) as response:
            raw=response.read(MAX_BYTES+1)
            kind=response.headers.get("Content-Type","")
    except urllib.error.HTTPError as e:
        if e.code in (400,401,403,404,406,429,503):
            return {"status":"OFFICIAL_SOURCE_HTTP_BLOCKED","http_status":e.code}
        raise
    except urllib.error.URLError as e:
        return {"status":"OFFICIAL_SOURCE_CONNECTION_UNAVAILABLE",
                "error_type":type(e).__name__}
    if len(raw)>MAX_BYTES:return {"status":"OFFICIAL_SOURCE_TOO_LARGE"}
    if expected=="zip" and zipfile.is_zipfile(io.BytesIO(raw)):
        result=zip_support_metadata(raw)
    elif raw.lstrip().startswith(b"<"):
        result=xml_support_metadata(raw)
    else:
        result={"status":"UNRECOGNIZED_OFFICIAL_SOURCE_CONTENT",
                "bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}
    result["reported_content_type"]=kind
    return result


def main():
    result={
       name:fetch_index(u,"zip" if name.endswith("ZIP") else "xml")
       for name,u in SOURCE_ENDPOINTS.items()
    }
    print(json.dumps({
      "status":"PENSTEMON_ALTERNATIVE_OFFICIAL_SUPPLEMENT_SOURCE_INDEX_ONLY",
      "paper_pmc":ARTICLE_ID,
      "candidate_statuses":result,
      "no_34_phenotyped_IDs_assigned":True,
      "no_flower_or_genetic_outcomes_read":True,
    },indent=2,sort_keys=True,ensure_ascii=False))


if __name__=="__main__":main()
