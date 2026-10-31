"""Stage20: fixed 80-species predictors; no genetic/climate fit or effect testing.

All restricted source geometries and site coordinates stay in the ephemeral work
folder. Only counts, hashes, software versions, and masked-feature QA are exported.
"""
from __future__ import annotations
import argparse,gc,hashlib,io,json,math,os,platform,sys,time,urllib.request,zipfile
from pathlib import Path
import numpy as np
import rasterio,shapefile,shapely,scipy,pyproj
from rasterio.windows import Window
from scipy.ndimage import binary_erosion
from scipy.spatial import cKDTree
from shapely.geometry import shape
from shapely.ops import transform
from shapely.validation import make_valid
from pyproj import CRS,Transformer
from site_feature_core import dense_vertices,sphere_xyz,vertex_distances,raster_values,complete_rows

HERE=Path(__file__).resolve().parent

def load(name):return json.loads((HERE/name).read_text())
def canonical_sha(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def obtain(url,path,size,sha256):
    if path.exists() and path.stat().st_size==size:
        h=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
        if h==sha256:return
    req=urllib.request.Request(url,headers={'User-Agent':'TTF-frozen80-input-extraction-v1'})
    h=hashlib.sha256();n=0
    with urllib.request.urlopen(req,timeout=180) as r,path.open('wb') as f:
        while True:
            chunk=r.read(2**20)
            if not chunk:break
            n+=len(chunk)
            if n>size:raise ValueError('asset larger than source receipt')
            f.write(chunk);h.update(chunk)
    if n!=size or h.hexdigest()!=sha256:
        path.unlink(missing_ok=True);raise ValueError('asset original size/SHA256 mismatch')

def sample_species(line,poly,latin):
    from frozen_background_sample import select_background
    if not poly.is_valid:poly=make_valid(poly)
    center=poly.representative_point()
    crs=CRS.from_proj4(f'+proj=aeqd +lat_0={center.y} +lon_0={center.x} +datum=WGS84 +units=m +no_defs')
    forward=Transformer.from_crs(4326,crs,always_xy=True)
    inverse=Transformer.from_crs(crs,4326,always_xy=True)
    line=transform(forward.transform,line);poly=transform(forward.transform,poly)
    d=np.arange(5000.,line.length,10000.)
    pts=shapely.line_interpolate_point(line,d)
    pts=pts[shapely.covers(poly,pts)]
    case_xy=shapely.get_coordinates(pts)
    lo,la=inverse.transform(case_xy[:,0],case_xy[:,1]);cases=np.column_stack([lo,la])
    xmin,ymin,xmax,ymax=poly.bounds
    xx=np.arange(math.ceil(xmin/20000)*20000,xmax,20000)
    yy=np.arange(math.ceil(ymin/20000)*20000,ymax,20000)
    if len(xx)*len(yy)>400000:raise ValueError('grid safety limit changed')
    gx,gy=np.meshgrid(xx,yy);gx=gx.ravel();gy=gy.ravel()
    inside=shapely.contains_xy(poly,gx,gy)
    near=shapely.intersects_xy(line.buffer(50000.),gx,gy)
    valid=inside&~near
    lo,la=inverse.transform(gx[valid],gy[valid])
    background=select_background(latin,zip(lo,la))
    return cases,np.asarray(background,float),int(valid.sum())

def rebuild_sites(panel,work):
    from source_integrity import fetch
    from source_overlap_gate import members,load_group,norm
    assets={a['role']:a for a in load('original_files_v02.json')['files']}
    b=work/'breaks.zip';r=work/'ranges.zip'
    fetch(assets['pleistocene_break_lines'],b);fetch(assets['surveyed_or_clipped_ranges'],r)
    kept=[s for s in panel['qualified_species_records'] if s['synthetic_estimator_min10_cases']]
    if len(kept)!=80:raise ValueError('not frozen 80-species input')
    allxy=[];labels=[];species=[];records=[]
    with zipfile.ZipFile(b) as bz,zipfile.ZipFile(r) as rz:
        bb=members(bz,'break');rr=members(rz,'range')
        for i,s in enumerate(kept):
            key=norm(s['common_name'])
            line,_=load_group(bz,bb[key][1]);poly,_=load_group(rz,rr[key][1])
            cases,bg,original_bg=sample_species(line,poly,s['scientific_name'])
            if len(cases)!=s['case_points'] or original_bg!=s['background_points']:
                raise ValueError('frozen source point counts changed for '+s['scientific_name'])
            if len(bg)!=min(1000,original_bg):raise ValueError('background cap mismatch')
            xy=np.vstack([cases,bg]);sphere_xyz(xy)
            allxy.append(xy);labels.extend([1]*len(cases)+[0]*len(bg));species.extend([i]*len(xy))
            records.append({**s,'capped_background_points':len(bg),'coordinate_sha256':hashlib.sha256(xy.astype('<f8').tobytes()).hexdigest()})
            print('SITE_DESIGN_VERIFIED',i+1,s['scientific_name'],len(cases),len(bg),flush=True)
    xy=np.vstack(allxy);y=np.asarray(labels,'u1');sp=np.asarray(species,'i4')
    if len(xy)!=76700 or y.sum()!=10013:raise ValueError('frozen global sample counts changed')
    return xy,y,sp,records

def vector_geometries(path,role,contract):
    with zipfile.ZipFile(path) as z:
        stem=next(n[:-4] for n in z.namelist() if n.lower().endswith('.shp'))
        rdr=shapefile.Reader(**{e:io.BytesIO(z.read(stem+'.'+e)) for e in ('shp','shx','dbf')})
        crs=CRS.from_wkt(z.read(stem+'.prj').decode('utf8','replace'))
        trans=Transformer.from_crs(crs,4326,always_xy=True)
        selected=[];names=set();n_all=0
        lakes={n.lower() for n in contract['source_gis'][1]['lake_names']}
        for sr in rdr.iterShapeRecords():
            n_all+=1;a=sr.record.as_dict()
            if role.startswith('ecoregion'):
                field='NA_L1CODE' if role.endswith('1') else 'NA_L2CODE'
                if float(a[field])==0:continue
            elif role=='naturalearth_rivers':
                if float(a['scalerank'])>3:continue
            elif role=='naturalearth_lakes':
                name=str(a.get('name') or a.get('name_en') or '').strip().lower()
                if name not in lakes:continue
                names.add(name)
            geom=transform(trans.transform,shape(sr.shape.__geo_interface__))
            if not geom.is_empty:selected.append(geom)
        if role=='naturalearth_lakes' and names!=lakes:raise ValueError('frozen eight lake source names differ')
    return selected,{'source_features':n_all,'retained_features':len(selected)}

def mountain_distances(path,xy):
    """Native 30 arcsec mountain cells; global boundary-center spherical distance."""
    xyz=[];boundary_count=0
    with rasterio.open(path) as ds:
        if ds.crs!=CRS.from_epsg(4326) or ds.count!=1:raise ValueError('wrong native mountain grid')
        for start in range(0,ds.height,512):
            end=min(ds.height,start+512);a=max(0,start-1);b=min(ds.height,end+1)
            image=ds.read(1,window=Window(0,a,ds.width,b-a))
            mountain=(image>=1)&(image<=7)
            edge=mountain&~binary_erosion(mountain,structure=np.ones((3,3),bool),border_value=0)
            edge=edge[start-a:end-a]
            rr,cc=np.nonzero(edge)
            if len(rr):
                lo,la=rasterio.transform.xy(ds.transform,rr+start,cc,offset='center')
                pts=np.column_stack([lo,la]);pts[:,0]=(pts[:,0]+180)%360-180
                xyz.append(sphere_xyz(pts));boundary_count+=len(rr)
    tree=cKDTree(np.vstack(xyz));del xyz;gc.collect()
    d=vertex_distances(xy,tree)
    classes=raster_values(path,xy,mask_nodata=False)
    inside=(classes>=1)&(classes<=7);d[inside]=0
    d[~np.isfinite(classes)]=np.nan
    return d,{'boundary_cell_centers':boundary_count,'inside_mountain_sites':int(inside.sum()),'distance_method':'global spherical nearest native boundary-cell center; <=0.66km cell-center uncertainty plus sphere/ellipsoid error','outside_grid_sites':int((~np.isfinite(classes)).sum())}

def extract(args):
    work=args.work;work.mkdir(parents=True,exist_ok=True)
    panel=json.loads(args.panel.read_text());contract=load('frozen_sampling_barriers_history_v01.json')
    xy,y,sp,records=rebuild_sites(panel,work)
    fields=[];names=[];receipts=[];details={}
    climate=json.loads(args.climate_assets.read_text())
    if len(climate['assets'])!=8:raise ValueError('expected eight original C1 rasters')
    for i,a in enumerate(climate['assets']):
        path=work/f'climate_{i}.tif';obtain(a['url'],path,a['size_bytes'],a['sha256'])
        val=raster_values(path,xy);fields.append(val);names.append(f"{a['variable']}_{a['time_index']}")
        receipts.append({'role':names[-1],'bytes':a['size_bytes'],'sha256':a['sha256']})
        path.unlink();print('CLIMATE_INPUT_VERIFIED',names[-1],int(np.isfinite(val).sum()),flush=True)
    gis=load('original_barrier_baseline_source_candidates_v01.json')
    sha={a['term']:a for a in load('barrier_gis_publisher_SHA_lock_v01.json')['asset_hashes']}
    vectors={}
    for a in gis['source_candidates']:
        term=a['term'];h=sha[term];path=work/(term+'.zip')
        obtain(a['url'],path,h['bytes'],h['sha256']);receipts.append({'role':term,'bytes':h['bytes'],'sha256':h['sha256']})
        if term=='unep_mountains':
            with zipfile.ZipFile(path) as z:
                if any(n.startswith('/') or '..' in Path(n).parts for n in z.namelist()):raise ValueError('unsafe source ZIP')
                z.extractall(work/'mountains')
            grid=work/'mountains/Mountain_UNEP-WCMC_2002/Mountain_UNEP-WCMC_2002/geog/mtn_02g'
            val,detail=mountain_distances(grid,xy);fields.append(val);names.append('mountain_km');details[term]=detail
        else:
            vectors[term],details[term]=vector_geometries(path,term,contract)
        path.unlink();print('BARRIER_SOURCE_VERIFIED',term,flush=True)
    for role,geometries in [('water_km',vectors['naturalearth_rivers']+vectors['naturalearth_lakes']),('ecoregion1_km',vectors['ecoregion_level1']),('ecoregion2_km',vectors['ecoregion_level2'])]:
        vertices=dense_vertices(geometries);val=vertex_distances(xy,vertices)
        if role=='water_km':
            for lake in vectors['naturalearth_lakes']:
                val[shapely.intersects_xy(lake,xy[:,0],xy[:,1])]=0
        fields.append(val);names.append(role);details[role]={'dense_vertices':len(vertices),'max_spherical_discretization_km':.79}
        del vertices;gc.collect();print('BARRIER_VALUES_EXTRACTED',role,len(val),flush=True)
    tri=contract['source_gis'][-1];path=work/'tri.tif'
    obtain(tri['url'],path,tri['bytes'],tri['sha256']);fields.append(raster_values(path,xy));names.append('TRI')
    receipts.append({'role':'TRI','bytes':tri['bytes'],'sha256':tri['sha256']});path.unlink()
    x=np.column_stack(fields);valid=complete_rows(x)
    climate_valid=complete_rows(x[:,:8]);rows=[]
    for i,s in enumerate(records):
        m=sp==i;case=m&(y==1);bg=m&(y==0)
        cp=int((case&valid).sum());bp=int((bg&valid).sum())
        rows.append({'species':s['scientific_name'],'role':s['split_role'],'before_cases':int(case.sum()),'before_background':int(bg.sum()),'climate_complete_cases':int((case&climate_valid).sum()),'climate_complete_background':int((bg&climate_valid).sum()),'all_features_complete_cases':cp,'all_features_complete_background':bp,'eligible_after_fields':bool(cp>=10 and bp>=20 and bp>=cp),'coordinate_sha256':s['coordinate_sha256']})
    eligible=[s for s in rows if s['eligible_after_fields']];roles=[sum(s['role']==i for s in eligible) for i in (0,1)]
    passed=len(eligible)>=50 and min(roles)>=20
    # Ephemeral array only. Workflow artifact allow-list excludes this directory.
    np.savez_compressed(work/'ephemeral_site_features.npz',xy=xy,case=y,species_index=sp,features=x,feature_names=np.asarray(names))
    report={'schema':'ttf_jensen80_stage20_actual_predictor_join_v01','status':'PASS_INPUT_SUPPORT_NOT_EMPIRICAL_EFFECT' if passed else 'NOT_EVALUABLE_INPUT_SUPPORT','prequalified_species':80,'prequalified_roles':[44,36],'sites':len(xy),'case_sites':int(y.sum()),'background_sites':int((y==0).sum()),'actual_site_counts_match_frozen_design':True,'features':names,'complete_counts_by_feature':{n:int(np.isfinite(x[:,i]).sum()) for i,n in enumerate(names)},'all_eight_climates_complete':int(climate_valid.sum()),'all_thirteen_features_complete':int(valid.sum()),'eligible_species_after_all_fields':len(eligible),'eligible_species_per_role':roles,'source_receipts':receipts,'barrier_details':details,'species_support':rows,'coordinates_raw_sha256':hashlib.sha256(xy.astype('<f8').tobytes()).hexdigest(),'feature_matrix_raw_sha256':hashlib.sha256(x.astype('<f8').tobytes()).hexdigest(),'class_labels_sha256':hashlib.sha256(y.tobytes()).hexdigest(),'restrictions':{'raw_coordinates_exported':False,'IUCN_UNEP_sources_exported':False,'author_baseline_exactly_reproduced':False,'St_Lawrence_Alaska_manual_waterways_present':False},'distance_approximation':{'line_densification_degrees':.01,'line_sampling_error_upper_km':.79,'mountain_cell_center_error_upper_km':.66,'sphere_not_WGS84_ellipsoid':True},'analysis_boundary':{'only_missingness_and_source_support_inspected':True,'genetic_history_AUC_calculated':False,'genetic_model_fitted':False,'prior_closed_endpoints_reopened':False},'software':{'python':platform.python_version(),'numpy':np.__version__,'rasterio':rasterio.__version__,'shapely':shapely.__version__,'pyproj':pyproj.__version__,'scipy':scipy.__version__}}
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('panel','climate-assets','work','output'):p.add_argument('--'+k,type=Path,required=True)
    args=p.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    result=extract(args);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','sites','all_eight_climates_complete','all_thirteen_features_complete','eligible_species_after_all_fields','eligible_species_per_role')},sort_keys=True),flush=True)
