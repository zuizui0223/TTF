"""First-shot actual MPI-ESM-P LGM BIO1/BIO12 raster extraction on 95 UNLABELED cells.

Use HTTP Range to retrieve ONLY two CRC-pinned original WorldClim
GeoTIFF ZIP members (~28.5MB compressed), without full archive or
individual plant data. Raster coordinates are coarse 2.5 arcmin cells
derived from the authentic source-only CSV; NO original lat/lon, source
population ID, garden, phenotype, fitness or genotype are exported.
This is an LGM climate extraction, not a climate-effect test.
"""
from __future__ import annotations
import argparse, hashlib, json, math, re, struct, time, zlib
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

HERE=Path(__file__).resolve().parent
CELLS=HERE/"source_registry/albano2025_unlabeled_worldclim_2_5m_cells_v01.json"
URL="https://geodata.ucdavis.edu/climate/cmip5/lgm/melgmbi_2-5m.zip"
HOST="geodata.ucdavis.edu"
ARCHIVE_BYTES=223064411
ENTRIES={
  "BIO1":{"name":"melgmbi1.tif","local_header_offset":0,
          "compressed_bytes":11503947,"uncompressed_bytes":12404728,
          "crc32":"7687a926"},
  "BIO12":{"name":"melgmbi12.tif","local_header_offset":35777679,
           "compressed_bytes":16949518,"uncompressed_bytes":19155477,
           "crc32":"1f15f5db"},
}
LOCAL=struct.Struct("<4s5H3I2H")
RES=2.5/60

def cell_center(row,col):
    if not (0<=row<4320 and 0<=col<8640):raise ValueError("Grid cell outside WorldClim")
    return -180+(col+0.5)*RES,90-(row+0.5)*RES

def member_span(entry):
    return entry["local_header_offset"],entry["local_header_offset"]+30+len(entry["name"].encode())+entry["compressed_bytes"]-1

def remote_member_bytes(entry,fetcher=None):
    if fetcher is None:fetcher=range_fetch
    start,end=member_span(entry)
    data=fetcher(URL,start,end)
    need=end-start+1
    if len(data)!=need:raise ValueError("Partial/oversized original member")
    v=LOCAL.unpack_from(data)
    sig,vers,flags,method,mtime,mdate,crc,comp,uncomp,n,extra=v
    if (sig!=b"PK\x03\x04" or method!=8 or
        data[LOCAL.size:LOCAL.size+n].decode("utf-8")!=entry["name"] or extra!=0 or
        LOCAL.size+n+entry["compressed_bytes"]!=len(data)):
        raise ValueError("Official ZIP local header/member geometry drift")
    if comp and comp!=entry["compressed_bytes"]:raise ValueError("Compressed member size drift")
    if uncomp and uncomp!=entry["uncompressed_bytes"]:raise ValueError("Uncompressed member size drift")
    encoded=data[LOCAL.size+n:]
    decoder=zlib.decompressobj(-15)
    raster=decoder.decompress(encoded,30_000_000)+decoder.flush()
    if not decoder.eof or decoder.unused_data:
        raise ValueError("Original ZIP member deflate truncated or malformed")
    if len(raster)!=entry["uncompressed_bytes"]:
        raise ValueError("Official uncompressed TIFF byte count mismatch")
    crc_actual=zlib.crc32(raster)&0xffffffff
    if f"{crc_actual:08x}"!=entry["crc32"]:
        raise ValueError("Official central directory CRC32 drift")
    return raster

def range_fetch(url,start,end):
    if url!=URL:raise ValueError("Only previously validated official MPI LGM archive")
    length=end-start+1
    if length<1 or length>18_000_000:raise ValueError("Unbounded original member request")
    errors=[]
    for attempt in range(3):
        try:
            req=Request(url,headers={"Range":f"bytes={start}-{end}",
                                     "Accept":"application/octet-stream",
                                     "User-Agent":"TTF-MPI-LGM-unlabeled-source-cell/2026-11-20"})
            with urlopen(req,timeout=95) as res:
                if res.status!=206 or res.geturl().split("/")[2]!=HOST:
                    raise ValueError("Official ranged request not verified")
                info=res.headers.get("Content-Range","")
                if info!=f"bytes {start}-{end}/{ARCHIVE_BYTES}":
                    raise ValueError("Official total archive byte range mismatch")
                blob=res.read(length+1)
                if len(blob)!=length:
                    raise ValueError("Official member range returned wrong byte count")
                return blob
        except (HTTPError,URLError,OSError,TimeoutError) as e:
            errors.append(type(e).__name__)
            if attempt<2:time.sleep(3)
    raise ValueError("Official member range transport failed after retries: "+",".join(errors))

def extract_cell_values(tiff,coords):
    import rasterio
    from rasterio.io import MemoryFile
    with MemoryFile(tiff) as mem:
        with mem.open() as src:
            if src.count!=1 or src.crs is None or not src.crs.is_geographic:
                raise ValueError("Original GeoTIFF missing 1-band geographic CRS")
            if not(math.isclose(abs(src.transform.a),RES,rel_tol=0,abs_tol=2e-5)
                   and math.isclose(abs(src.transform.e),RES,rel_tol=0,abs_tol=2e-5)):
                raise ValueError("Official GeoTIFF grid resolution inconsistent")
            xy=[cell_center(r,c) for r,c in coords]
            sampled=list(src.sample(xy,masked=True))
            values=[]
            for sample in sampled:
                item=sample[0]
                if getattr(item,"mask",False) is not False:
                    if bool(item.mask):values.append(None);continue
                val=float(item)
                if not math.isfinite(val) or (src.nodata is not None and val==src.nodata):
                    values.append(None)
                else:
                    values.append(val)
            info={
              "raster_width":src.width,"raster_height":src.height,
              "raster_crs":src.crs.to_string(),"raster_dtype":src.dtypes[0],
              "raster_resolution":list(src.res),"raster_nodata":src.nodata,
              "sampled_nominal_cells":len(coords),
              "valid_cells":sum(x is not None for x in values),
            }
    return values,info

def evaluate(grid,fetcher=None,sampler=extract_cell_values):
    if grid.get("schema")!="ttf_Albano2025_nominal_WorldClim_2_5m_unlabeled_grid_query_v01":
        raise ValueError("Unrecognized source-only cell family")
    if grid.get("source_csv_sha256")!="b3796db37145dd32ab148f5c2105ca0646923feab153d8857d4f90a5539f080c":
        raise ValueError("Original source provenance hash changed")
    cells=[(x["row"],x["col"]) for x in grid["cells"]]
    if grid["source_population_count"]!=96 or len(cells)!=95 or len(set(cells))!=95:
        raise ValueError("Expected 95 unique original grid cells")
    for row,col in cells:cell_center(row,col)
    profiles={}
    for code,entry in ENTRIES.items():
        raster=remote_member_bytes(entry,fetcher=fetcher)
        value,header=sampler(raster,cells)
        if len(value)!=95:raise ValueError("Raster sample count drift")
        profiles[code]={"original_tiff_sha256":hashlib.sha256(raster).hexdigest(),
                       "original_tiff_size":len(raster),"raster":header,
                       "values":value}
    return {
       "schema":"ttf_Albano2025_MPI_ESM_P_source_only_2_5m_LGM_cells_v01",
       "status":"ACTUAL_LGM_GCM_RASTER_VALUES_SAMPLED_NO_INDEPENDENT_PREDICTOR_EFFECT_YET",
       "official_source_URL":URL,"source_Dryad_original_hash":grid["source_csv_sha256"],
       "source_population_count":96,"unlabeled_grid_cell_count":95,
       "grid_coordinate_system":"WorldClim 1.4 global 2.5 arcmin cell row/col, center-point lookup",
       "climate_age":"Last Glacial Maximum about 22 ka, not exactly 21 ka",
       "variable_units":{"BIO1":"degrees C times 10","BIO12":"precipitation mm"},
       "cell_profiles":[{"row":r,"col":c,
           "BIO1_raw":profiles["BIO1"]["values"][i],
           "BIO12_raw":profiles["BIO12"]["values"][i]} for i,(r,c) in enumerate(cells)],
       "geotiff_provenance":{k:{x:v[x] for x in ("original_tiff_sha256","original_tiff_size","raster")} for k,v in profiles.items()},
       "original_exact_GPS_coordinates_or_population_IDs_exposed":False,
       "modern_matched_baseline_rasters_extracted":False,
       "fitness_genetics_adaptation_responses_opened":False,
       "LGM_independent_predictor_gain_tested":False,
    }

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    a=p.parse_args();grid=json.loads(CELLS.read_text("utf-8"))
    result=evaluate(grid)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))
if __name__=="__main__":main()
