#!/usr/bin/env python3
"""Cloud-runner exact CHELSA-TraCE21k 0BP sampling, without taxa or DNA data.

Input is *unlabelled* frozen C1 coordinates in exact occurrence-row order.
Outputs only 4 climate columns plus byte/provenance hashes; source coordinates
are not placed in workflow artifacts. Verify all original COG SHA-256 receipts.
"""
from __future__ import annotations
import argparse,base64,hashlib,json,lzma,os,time
from pathlib import Path
import numpy as np

EXPECTED_COORDINATES_SHA256 = 'ae5d5b1e192500e7b12b42789b4e8ee7931f439e8c46905b267b0554ef25cef6'
EXPECTED_COORDINATES_N = 96607
VARIABLES = ('bio01','bio07','bio12','bio15')


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<22),b''):h.update(block)
    return h.hexdigest()


def decode_coordinates(path: Path) -> np.ndarray:
    encoded=path.read_text(encoding='ascii')
    raw=lzma.decompress(base64.b64decode(encoded,validate=False))
    observed=hashlib.sha256(raw).hexdigest()
    if observed!=EXPECTED_COORDINATES_SHA256:
        raise ValueError(f'coordinate SHA-256 mismatch: {observed}')
    a=np.frombuffer(raw,dtype='<f8').reshape(-1,2).copy()
    if a.shape != (EXPECTED_COORDINATES_N, 2):
        raise ValueError(f'coordinate count mismatch: {a.shape}')
    if not np.isfinite(a).all() or np.any(np.abs(a[:,0])>180) or np.any(np.abs(a[:,1])>90):
        raise ValueError('nonfinite or impossible coordinates')
    return a


def download_exact(item: dict, target: Path, attempts: int=3) -> Path:
    import requests
    out=target/item['logical_basename']
    if out.exists() and out.stat().st_size==item['size_bytes'] and sha256(out)==item['sha256']:
        print(f'REUSE VERIFIED {out.name}',flush=True)
        return out
    partial=target/(out.name+'.part')
    for attempt in range(attempts):
        pos=partial.stat().st_size if partial.exists() else 0
        if pos>item['size_bytes']:
            partial.unlink();pos=0
        req_headers={'Range':f'bytes={pos}-'} if pos else {}
        try:
            with requests.get(item['resolved_url'],stream=True,timeout=(25,120),headers=req_headers) as r:
                r.raise_for_status()
                if pos and r.status_code==206:
                    if not r.headers.get('Content-Range','').startswith(f'bytes {pos}-'):
                        raise RuntimeError('invalid HTTP byte-range response')
                    mode='ab'
                elif r.status_code==200:
                    mode='wb'
                else:
                    raise RuntimeError(f'bad status {r.status_code}')
                with partial.open(mode) as f:
                    for chunk in r.iter_content(chunk_size=4<<20):
                        if chunk:f.write(chunk)
            if partial.stat().st_size != item['size_bytes']:
                raise IOError(f'partial {partial.stat().st_size} bytes (required {item["size_bytes"]})')
            observed=sha256(partial)
            if observed != item['sha256']:
                raise ValueError(f'CORRUPT {out.name}: SHA {observed}')
            partial.replace(out)
            print(f'FETCH VERIFIED {out.name}, {out.stat().st_size:,} bytes',flush=True)
            return out
        except (requests.RequestException,OSError) as e:
            print(f'TRANSIENT download problem, attempt {attempt+1}/{attempts}: {e}',flush=True)
            if attempt+1 == attempts:raise
            time.sleep(5*(attempt+1))
    raise RuntimeError('download failed')


def collect(manifest_path: Path, encoded_coords: Path, work: Path, output: Path):
    import rasterio
    m=json.loads(manifest_path.read_text())
    items=m['assets']
    if len(items)!=4 or [a['logical_basename'].split('_')[2] for a in items]!=list(VARIABLES):
        raise ValueError('must have exactly four ordered 0BP assets')
    coords=decode_coordinates(encoded_coords)
    print(f'COORDINATES VERIFIED: {len(coords):,}, {EXPECTED_COORDINATES_SHA256}',flush=True)
    work.mkdir(parents=True,exist_ok=True)
    output.mkdir(parents=True,exist_ok=True)
    values=np.empty((len(coords),len(items)),dtype='<f8')
    raster_details=[]
    for j,item in enumerate(items):
        path=download_exact(item,work)
        with rasterio.open(path) as ds:
            if ds.count!=1 or ds.crs is None or not ds.crs.is_geographic:
                raise ValueError(f'bad raster georeferencing or band count: {path}')
            sample=ds.sample([tuple(x) for x in coords],indexes=1,masked=True)
            for i,v in enumerate(sample):
                if np.ma.is_masked(v[0]) or not np.isfinite(v[0]) or (ds.nodata is not None and v[0]==ds.nodata):
                    raise ValueError(f'nodata/nonfinite raster value {item["logical_basename"]} at row {i}')
                values[i,j]=float(v[0])
            raster_details.append({
                'variable':VARIABLES[j], 'sha256':item['sha256'],
                'size_bytes':int(path.stat().st_size),
                'raster_dtype':str(ds.dtypes[0]), 'crs':str(ds.crs),
                'shape':[int(ds.height),int(ds.width)],
                'transform':str(ds.transform)
            })
        print(f'REAL SAMPLES {VARIABLES[j]}: {len(coords)}; range {float(values[:,j].min()):.5f}..{float(values[:,j].max()):.5f}',flush=True)
    np.save(output/'trace0_climate_values.npy',values,allow_pickle=False)
    manifest={
      'schema':'ttf_q_trace0_exact_samples_v0.1',
      'status':'SAMPLED_ALL_FOUR_VERIFIED_COG_ASSETS',
      'units':'Original raw provider raster values, no unit conversions',
      'coordinate_sha256':EXPECTED_COORDINATES_SHA256,
      'coordinate_rows':EXPECTED_COORDINATES_N,
      'variables':list(VARIABLES),
      'matrix_shape':list(values.shape),
      'matrix_sha256':sha256(output/'trace0_climate_values.npy'),
      'assets':raster_details,
      'includes_occurrence_coordinates':False,
      'includes_species_or_genetic_responses':False,
    }
    (output/'trace0_exact_samples_receipt.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print('SUCCESS: exact 0BP point samples, no genetic outcome opened',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--coordinates',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    collect(args.manifest,args.coordinates,args.work,args.output)
