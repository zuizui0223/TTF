"""ZIP central-directory metadata via bounded official WorldClim HTTP ranges.

No climate raster values or user source coordinates are read in this step.
A server that ignores Range is rejected before an oversized download.
"""
from __future__ import annotations
import argparse, json, re, struct
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"lgm_paired_climate_source_only_contract_v01.json"
HOST="geodata.ucdavis.edu"
MAX_TAIL=65536
MAX_CENTRAL=500000
CENTRAL=struct.Struct("<4s6H3I5H2I")
EOCD=struct.Struct("<4s4H2IH")

def read_range(url,start,end):
    if urlparse(url).scheme!="https" or urlparse(url).hostname!=HOST:
        raise ValueError("Source must be official HTTPS WorldClim archive")
    spec=(f"bytes={start}-{end}" if start is not None else f"bytes=-{end}")
    req=Request(url,headers={"Range":spec,"Accept":"application/octet-stream",
                             "User-Agent":"TTF-WorldClim-ZIP-index-only/2026-11-20"})
    with urlopen(req,timeout=25) as res:
        if res.status!=206 or urlparse(res.geturl()).hostname!=HOST:
            raise ValueError("Official server did not honor ranged HTTPS access")
        expected=(end-start+1) if start is not None else end
        raw=res.read(expected+1)
        if len(raw)!=expected:
            raise ValueError("Incorrect HTTP Range size / excess payload")
        content_range=res.headers.get("Content-Range","")
        match=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",content_range)
        if not match or int(match.group(2))-int(match.group(1))+1!=len(raw):
            raise ValueError("HTTP Content-Range not internally consistent")
        if start is not None and (int(match.group(1))!=start or int(match.group(2))!=end):
            raise ValueError("HTTP Content-Range does not match requested bytes")
        return raw,int(match.group(3))

def parse_zip_central(tail,total,fetch_cd):
    idx=tail.rfind(b"PK\x05\x06")
    if idx<0 or len(tail)-idx<EOCD.size:raise ValueError("ZIP EOCD not found in last 64KiB")
    sig,disk,disk_cd,n_disk,n_total,sz,offset,comment=EOCD.unpack_from(tail,idx)
    if (disk or disk_cd or n_disk!=n_total or n_total==65535 or
        sz==4294967295 or offset==4294967295 or sz>MAX_CENTRAL or
        offset+sz>total or idx+EOCD.size+comment!=len(tail)):
        raise ValueError("ZIP unsupported multi-volume/ZIP64/inconsistent central directory")
    cd=fetch_cd(offset,offset+sz-1)
    out=[];pos=0
    for i in range(n_total):
        if pos+CENTRAL.size>len(cd):raise ValueError("Truncated ZIP central directory")
        v=CENTRAL.unpack_from(cd,pos)
        if v[0]!=b"PK\x01\x02":raise ValueError("Central record signature invalid")
        name_len,extra_len,comment_len=v[10:13]
        name_bytes=cd[pos+CENTRAL.size:pos+CENTRAL.size+name_len]
        name=name_bytes.decode("utf-8" if v[3]&0x800 else "cp437")
        if len(out)>=100:raise ValueError("Unexpectedly large WorldClim archive index")
        if name.startswith("/") or ".." in Path(name).parts:
            raise ValueError("Untrusted ZIP path")
        out.append({
            "filename":name,
            "method":v[4],
            "CRC32_hex":f"{v[7]:08x}",
            "compressed_bytes":v[8],
            "uncompressed_bytes":v[9],
            "local_header_offset":v[16],
        })
        pos+=CENTRAL.size+name_len+extra_len+comment_len
    if pos!=len(cd):raise ValueError("Trailing ZIP central-directory bytes")
    return out

def source_urls(contract):
    if contract.get("schema")!="ttf_Albano2025_pre_response_LGM_climate_exposure_contract_v01":
        raise ValueError("Frozen 21ka source-only contract changed")
    expected={
      "WC1.4_MODERN":"https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/bio_2-5m_bil.zip",
      "CCSM4":"https://geodata.ucdavis.edu/climate/cmip5/lgm/cclgmbi_2-5m.zip",
      "MIROC-ESM":"https://geodata.ucdavis.edu/climate/cmip5/lgm/mrlgmbi_2-5m.zip",
      "MPI-ESM-P":"https://geodata.ucdavis.edu/climate/cmip5/lgm/melgmbi_2-5m.zip",
    }
    found={"WC1.4_MODERN":contract["baseline_url"]}
    for m in contract["models"]:found[m["code"]]=m["official_url"]
    if found!=expected:raise ValueError("GCM/source list drift")
    return found

def inspect_one(url,ranger=read_range):
    last,total=ranger(url,None,MAX_TAIL)
    entries=parse_zip_central(last,total,
           lambda start,end:ranger(url,start,end)[0])
    return {"status":"OFFICIAL_ZIP_CENTRAL_DIRECTORY_VERIFIED",
            "official_archive_total_bytes":total,
            "archive_file_count":len(entries),
            "files":entries,
            "raster_data_values_read":False}

def evaluate(contract,inspector=inspect_one):
    archives={}
    for code,url in source_urls(contract).items():
        try:
            archives[code]=inspector(url)
        except (HTTPError,URLError,TimeoutError,ValueError,OSError,UnicodeError,
                struct.error) as exc:
            archives[code]={"status":"OFFICIAL_ZIP_INDEX_ACCESS_NOT_CONFIRMED",
                            "error_type":type(exc).__name__,
                            "http_code":getattr(exc,"code",None),
                            "reason":str(exc)[:160] if isinstance(exc,ValueError) else None}
    return {
       "schema":"ttf_Albano2025_official_WorldClim_lgm_ZIP_directory_metadata_v01",
       "status":"CLIMATE_ZIP_INDEX_ONLY_NO_RASTER_OR_LGM_PREDICTOR_RESULT",
       "archives":archives,
       "actual_96_source_coordinates_submitted_to_ci":False,
       "all_LGM_climate_values_extracted":False,
       "fitness_or_genetic_responses_read":False,
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    a=p.parse_args();x=evaluate(json.loads(CONTRACT.read_text("utf-8")))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(x,sort_keys=True))

if __name__=="__main__":main()
