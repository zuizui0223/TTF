"""Source-only WorldClim 1.4 1960-1990 BIO1/BIO12 actual BIL source-cell extraction.

Only four OFFICIAL ZIP members (two data BIL + HDR pairs) are fetched
by prevalidated byte ranges with official central-directory CRC checks.
Unlabeled 95 original 2.5' grid cells; no original site GPS, species
response, or genetic values in the GitHub repository.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re,struct,tempfile,time,zlib
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen

HERE=Path(__file__).resolve().parent
CELLS=HERE/"source_registry/albano2025_unlabeled_worldclim_2_5m_cells_v01.json"
MANIFEST=HERE/"source_registry/worldclim_14_current_bio1_bio12_BIL_member_index_v01.json"
URL="https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/bio_2-5m_bil.zip"
ARCHIVE_BYTES=129319755
HEADER=struct.Struct("<4s5H3I2H")
RES=2.5/60

def cell_center(row,col):
    if not(0<=row<3600 and 0<=col<8640):raise ValueError("Grid location outside WorldClim geography")
    return (-180+(col+0.5)*RES,90-(row+0.5)*RES)

def request(url,start,end):
    if url!=URL:raise ValueError("Only official paired WorldClim baseline supported")
    n=end-start+1
    if n<1 or n>12_000_000:raise ValueError("Member request outside bounded source")
    errors=[]
    for i in range(3):
        try:
            req=Request(url,headers={"Range":f"bytes={start}-{end}","Accept":"application/octet-stream",
                                     "User-Agent":"TTF-WC14-BIO-LGM-paired-source-only/2026-11-20"})
            with urlopen(req,timeout=95) as r:
                if r.status!=206 or r.geturl().split("/")[2]!="geodata.ucdavis.edu":
                    raise ValueError("Unapproved official archive response")
                if r.headers.get("Content-Range")!=f"bytes {start}-{end}/{ARCHIVE_BYTES}":
                    raise ValueError("Original archive Content-Range identity drift")
                result=r.read(n+1)
                if len(result)!=n:raise ValueError("Member HTTP size mismatch")
                return result
        except (HTTPError,URLError,OSError,TimeoutError) as exc:
            errors.append(type(exc).__name__)
            if i<2:time.sleep(3)
    raise ValueError("Official baseline member network unavailable: "+",".join(errors))

def read_member(entry,fetcher=request):
    n=entry["compressed_bytes"]
    name=entry["filename"]
    offset=entry["local_header_offset"]
    data=fetcher(URL,offset,offset+30+len(name.encode())+n-1)
    if len(data)!=30+len(name)+n:raise ValueError("Incomplete ZIP source member")
    head=HEADER.unpack_from(data)
    sig,ver,flags,method,mtime,mdate,crc,comp,uncomp,fnlen,extra=head
    if(sig!=b"PK\x03\x04" or method!=8 or extra!=0 or fnlen!=len(name) or
       data[HEADER.size:HEADER.size+fnlen].decode()!=name):
        raise ValueError("Original ZIP local member identity mismatch")
    if comp and comp!=n:raise ValueError("Original compressed member size drift")
    if uncomp and uncomp!=entry["uncompressed_bytes"]:
        raise ValueError("Original uncompressed member size drift")
    encoded=data[HEADER.size+fnlen:]
    decoder=zlib.decompressobj(-15)
    raw=decoder.decompress(encoded,65_000_000)+decoder.flush()
    if not decoder.eof or decoder.unused_data or len(raw)!=entry["uncompressed_bytes"]:
        raise ValueError("Original deflate member not fully verified")
    if f"{zlib.crc32(raw)&0xffffffff:08x}"!=entry["CRC32_hex"]:
        raise ValueError("Official ZIP member CRC32 check failed")
    return raw

def extract_values(bil,hdr,points):
    import rasterio
    with tempfile.TemporaryDirectory(prefix="ttf_wc14_sourceonly_") as d:
        p=Path(d)/"data.bil";p.write_bytes(bil)
        p.with_suffix(".hdr").write_bytes(hdr)
        with rasterio.open(p) as src:
            if src.count!=1 or (src.crs is not None and not src.crs.is_geographic):
                raise ValueError("Baseline BIL geographic raster metadata invalid")
            if not(math.isclose(abs(src.transform.a),RES,abs_tol=2e-5)
                   and math.isclose(abs(src.transform.e),RES,abs_tol=2e-5)):
                raise ValueError("Baseline BIL raster resolution incompatible")
            if (abs(src.bounds.left+180)>0.1 or abs(src.bounds.top-90)>0.1):
                raise ValueError("Baseline BIL raster origin not global WGS84")
            xy=[cell_center(row,col) for row,col in points]
            samples=list(src.sample(xy,masked=True))
            out=[]
            for z in samples:
                v=z[0]
                if getattr(v,"mask",False) is not False and bool(v.mask):
                    out.append(None)
                elif not math.isfinite(float(v)) or (src.nodata is not None and float(v)==src.nodata):
                    out.append(None)
                else:
                    out.append(float(v))
            metadata={
              "raster_crs_in_embedded_header":str(src.crs) if src.crs else None,
              "raster_geographic_origin":["-180 longitude","+90 latitude"],
              "raster_width":src.width,"raster_height":src.height,
              "raster_dtype":src.dtypes[0],"raster_nodata":src.nodata,
              "valid_values":sum(x is not None for x in out),
              "sampled_cells":len(out),
            }
    return out,metadata

def evaluate(grid,manifest,fetcher=request,sampler=extract_values):
    if grid.get("schema")!="ttf_Albano2025_nominal_WorldClim_2_5m_unlabeled_grid_query_v01":
        raise ValueError("No original source-only cell query")
    cells=[(p["row"],p["col"]) for p in grid["cells"]]
    if len(cells)!=95 or len(set(cells))!=95 or grid["source_population_count"]!=96:
        raise ValueError("Original 95-cell/96-source design drift")
    if (manifest.get("schema")!="ttf_Albano2025_WC14_modern_1960_1990_BIL_frozen_member_manifest_v01"
        or manifest["original_archive_total_bytes"]!=ARCHIVE_BYTES
        or manifest["original_archive_url"]!=URL
        or manifest["original_archive_member_count"]!=38):
        raise ValueError("Original paired baseline archive identity drift")
    files={v["filename"]:v for v in manifest["members"]}
    if set(files)!={"bio1.bil","bio1.hdr","bio12.bil","bio12.hdr"}:
        raise ValueError("Four required BIL/HDR members not frozen")
    sampled={};original={}
    for bio in ("bio1","bio12"):
        rawbil=read_member(files[bio+".bil"],fetcher)
        rawhdr=read_member(files[bio+".hdr"],fetcher)
        values,metadata=sampler(rawbil,rawhdr,cells)
        if len(values)!=95:raise ValueError("Source-cell sample count drift")
        sampled[bio.upper()]=values
        original[bio.upper()]={
          "bil_sha256":hashlib.sha256(rawbil).hexdigest(),
          "hdr_sha256":hashlib.sha256(rawhdr).hexdigest(),
          "metadata":metadata,
        }
    return {
      "schema":"ttf_Albano2025_WorldClim14_modern_1960_1990_actual_cells_v01",
      "status":"MATCHED_BASELINE_RASTER_CELLS_EXTRACTED_HISTORICAL_PREDICTOR_UNTESTED",
      "official_archive_url":URL,
      "baseline_period":"1960-1990",
      "cell_count":95,
      "variable_units":{"BIO1":"degC times 10","BIO12":"precipitation millimeters"},
      "cells":[{"row":row,"col":col,
                "BIO1_raw":sampled["BIO1"][i],
                "BIO12_raw":sampled["BIO12"][i]}
               for i,(row,col) in enumerate(cells)],
      "member_provenance":original,
      "original_exact_GPS_or_population_IDs_shared":False,
      "LGM_22ka_effect_tested":False,
      "phenotype_fitness_genetic_responses_used":False,
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    g=json.loads(CELLS.read_text("utf-8"))
    m=json.loads(MANIFEST.read_text("utf-8"))
    result=evaluate(g,m)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
