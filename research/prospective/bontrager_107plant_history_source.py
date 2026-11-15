"""New Bontrager plant reciprocal transplant HISTORY CONTRASTS — SOURCE ONLY.

No fitness value, genotype, climate–fitness association, historical range or
other TTF response is read or modeled. Original Dryad/UBC archive plus all
eight original CHELSA 21ka/0BP COG sizes and SHA256 must match.

Rows in a common experimental garden share paper/species/site/fitness-type.
All pairwise origin-source comparisons and weighted historical-predictor
partial residuals are purely environmental design-qualification measures.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,itertools,json,math,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from urllib.request import Request,urlopen

import numpy as np
from bontrager_mirror_zip_v01 import verify_original_zip,ORIGINAL_URL
from bontrager_transplant_geo_census_v01 import (
    metadata_census,original_lonlat,binomial,connected_components
)
from proceed_plant_chelsa_source_support import verified_raster
from site_feature_core import raster_values

HERE=Path(__file__).resolve().parent
CONTRACT=HERE/"bontrager_107plant_history_contrast_contract_v01.json"
GEO=HERE/"bontrager_transplant_geo_contract_v01.json"
ZIP_LOCK=HERE/"bontrager_mirror_zip_contract_v01.json"
ASSETS=HERE/"chelsa_trace21k_fixed8_assets_v01.json"

def sha_float64(a):
    return hashlib.sha256(np.asarray(a,dtype="<f8").tobytes(order="C")).hexdigest()

def haversine_km(a,b):
    lo1,la1=a;lo2,la2=b
    x=math.sin(math.radians(la2-la1)/2)**2
    x+=math.cos(math.radians(la1))*math.cos(math.radians(la2))*math.sin(
       math.radians((lo2-lo1+180)%360-180)/2)**2
    return 6371.0088*2*math.asin(math.sqrt(max(0,min(x,1))))

def roles_for_graph(graph):
    blocks=connected_components(graph)
    blocks.sort(key=lambda b:(-len(b[0]),
         hashlib.sha256("|".join(sorted(b[0])).encode()).hexdigest()))
    role={};chosen={"train":set(),"heldout":set()}
    for taxa,papers in blocks:
        dst="train" if len(chosen["train"])<=len(chosen["heldout"]) else "heldout"
        for name in taxa:
            if name in role:raise ValueError("Species belongs to two disconnected components")
            role[name]=dst
        chosen[dst].update(taxa)
    digest=hashlib.sha256(json.dumps(
        {k:sorted(v) for k,v in chosen.items()},
        sort_keys=True,separators=(",",":")).encode()).hexdigest()
    return role,digest

def extract_source_cells(csv_bytes,base,allowed_groups):
    """Return only original site metadata. Index biological response NEVER."""
    rdr=csv.reader(io.TextIOWrapper(io.BytesIO(csv_bytes),encoding="utf-8-sig",newline=""),
                   delimiter=",",strict=True)
    fields=next(rdr)
    keys=base["allowed_metadata_only_columns"]
    if len(fields)!=47 or len(set(fields))!=47 or not set(keys).issubset(fields):
        raise ValueError("47-field original fitness CSV source metadata changed")
    ix={k:fields.index(k) for k in keys}
    candidates=[];sourced=defaultdict(set);gardens=defaultdict(set)
    for values in rdr:
        if not values:continue
        if len(values)!=47:raise ValueError("Original source CSV row width changed")
        m={k:values[i].strip() for k,i in ix.items()}
        if m["functional_group"] not in allowed_groups:continue
        name=binomial(m["taxon"]);study=m["study"];src=m["source"];site=m["site"];ft=m["fitness_type"]
        if not(name and study and src and site and ft):continue
        sc=original_lonlat(m["long_source"],m["lat_source"])
        gc=original_lonlat(m["long_site"],m["lat_site"])
        if sc is None or gc is None:continue
        sourced[(study,name,src)].add(sc)
        gardens[(study,name,site)].add(gc)
        candidates.append((study,name,src,site,ft,sc,gc))
    conflicted_sources={k for k,v in sourced.items() if len(v)!=1}
    conflicted_gardens={k for k,v in gardens.items() if len(v)!=1}
    cells=defaultdict(dict);destinations={}
    for study,sp,src,site,ft,sc,gc in candidates:
        if ((study,sp,src) in conflicted_sources or
            (study,sp,site) in conflicted_gardens):continue
        k=(study,sp,site,ft)
        if src in cells[k] and cells[k][src]!=sc:
            raise ValueError("Inconsistent source population coordinate within cell")
        cells[k][src]=sc
        destinations[k]=gc
    good={k:v for k,v in cells.items() if len(v)>=2 and len(set(v.values()))>=2}
    graph=defaultdict(set)
    for key in good:graph[key[0]].add(key[1])
    return good,{k:destinations[k] for k in good},graph

def species_equal_training_scales(cells,roles,site_features,index):
    unique_units=defaultdict(dict)
    for (study,sp,site,ft),sources in cells.items():
        if roles[sp]!="train":continue
        for name,xy in sources.items():
            key=(study,sp,name)
            if key in unique_units[sp] and unique_units[sp][key]!=xy:
                raise ValueError("Original source coordinate contradiction")
            unique_units[sp][key]=xy
    species=sorted(unique_units)
    if len(species)!=54:raise ValueError("Original 54 training plant species drifted")
    positions=[];weights=[]
    for sp in species:
        locs=list(unique_units[sp].values())
        for loc in locs:
            positions.append(site_features[index[loc],4:8])
            weights.append(1/(len(species)*len(locs)))
    arr=np.asarray(positions,float);w=np.asarray(weights,float)
    if not np.isclose(w.sum(),1.):raise ValueError("Training source-site species weights invalid")
    mu=np.sum(arr*w[:,None],axis=0)
    sd=np.sqrt(np.maximum(np.sum(((arr-mu)**2)*w[:,None],axis=0),0))
    if not np.isfinite(sd).all() or (sd<1e-10).any():
        raise ValueError("No training variation in required original present climate")
    return sd,unique_units

def h_source_summary(features,sd):
    x=np.asarray(features,float)
    if x.ndim!=2 or x.shape[1]!=8 or not np.isfinite(x).all():
        raise ValueError("Eight finite original CHELSA variables required")
    return np.linalg.norm((x[:,:4]-x[:,4:8])/sd,axis=1)

def matched_source_only(cells,destinations,roles,sites,features,contract):
    index={xy:i for i,xy in enumerate(sites)}
    if len(index)!=len(sites) or len(features)!=len(sites):
        raise ValueError("Exact source site unique mapping failed")
    source_sd,units=species_equal_training_scales(cells,roles,features,index)
    h=h_source_summary(features,source_sd)
    hw=[];weights=[]
    for sp,locs in units.items():
        for xy in locs.values():
            hw.append(h[index[xy]]);weights.append(1/(len(units)*len(locs)))
    hv=np.asarray(hw);w=np.asarray(weights)
    hm=float(np.sum(hv*w))
    hsd=float(np.sqrt(np.sum(w*(hv-hm)**2)))
    if hsd<=1e-12:raise ValueError("Training historical predictor has zero variation")
    paired=set();matched=set();per_type=Counter();usable=set()
    conditional={};cellweights=defaultdict(list)
    count_cells=Counter()
    all_pair=0;passmodern=0;passhistory=0;passboth=0
    by_species_cells=Counter(k[1] for k in cells)
    grouprows=[]
    for key in sorted(cells):
        study,sp,site,ft=key
        sources=cells[key];dest=destinations[key]
        # Historical and modern environmental values are source metadata only.
        if not np.isfinite(features[index[dest],4:8]).all():
            continue
        available=[(name,xy) for name,xy in sorted(sources.items())
                   if np.isfinite(features[index[xy]]).all()]
        if len(available)<2 or len(set(xy for _,xy in available))<2:continue
        count_cells[roles[sp]]+=1
        usable.add(sp)
        if roles[sp]=="train":
            bw=1/(54*by_species_cells[sp]*len(available))
            for name,xy in available:
                i=index[xy]
                g=features[i,4:8]
                row=np.concatenate([g,np.array([xy[1],xy[0],haversine_km(xy,dest)])])
                grouprows.append((key,sp,float(bw),float(h[i]),row))
        for (srca,xya),(srcb,xyb) in itertools.combinations(available,2):
            if xya==xyb:continue
            all_pair+=1
            ia=index[xya];ib=index[xyb]
            delta0=float(np.linalg.norm(
              (features[ia,4:8]-features[ib,4:8])/source_sd))
            deltah=abs(float(h[ia]-h[ib]))/hsd
            a=delta0<=contract["matching"]["source_only_modern_matched_threshold"]
            b=deltah>=contract["matching"]["source_only_H_contrasted_threshold"]
            if a:passmodern+=1
            if b:passhistory+=1
            if a and b:
                passboth+=1
                matched.add((study,sp,ft,site))
                per_type[ft]+=1
                paired.add((study,sp,ft,min(srca,srcb),max(srca,srcb)))
            if all_pair>1_000_000:
                raise ValueError("Source-only repeated origin combinations unbounded")
    # Nonredundancy *within* each same-garden cell; fix all cell means.
    groups=defaultdict(list)
    for key,sp,bw,hi,row in grouprows:
        groups[key].append((sp,bw,hi,row))
    ys=[];xs=[];ws=[];spec=[]
    for key,rs in groups.items():
        mat=np.stack([v[3] for v in rs])
        hh=np.asarray([v[2] for v in rs])
        mm=mat.mean(axis=0);mh=hh.mean()
        for sp,bw,hi,row in rs:
            xs.append(row-mm);ys.append(hi-mh);ws.append(bw);spec.append(sp)
    if not xs:raise ValueError("No qualified training source garden cell")
    xx=np.asarray(xs); yy=np.asarray(ys);ww=np.asarray(ws)
    variance=float(np.sum(ww*yy*yy))
    if variance<=1e-12:raise ValueError("No within-garden historical magnitude contrast")
    coef=np.linalg.lstsq(xx*np.sqrt(ww[:,None]),yy*np.sqrt(ww),rcond=1e-12)[0]
    residual=yy-xx@coef
    num=float(np.sum(ww*residual**2))
    frac=num/variance
    if not (0<=frac<=1.000001):raise ValueError("Source history residual fraction invalid")
    cnt_species={sp for study,sp,ft,site in matched}
    cnt_studies={study for study,sp,ft,site in matched}
    matched_role={k:len({sp for st,sp,ft,si in matched if roles[sp]==k})
                  for k in ("train","heldout")}
    pre=contract["matching"]["minimal_candidate_panel"]
    passed=(len(cnt_species)>=pre["min_species_with_matched_history_contrasts"]
            and len(cnt_studies)>=pre["min_studies_with_matched_history_contrasts"]
            and len({sp.split()[0] for sp in cnt_species})>=pre["min_genera_with_matched_history_contrasts"]
            and min(matched_role.values())>=pre["min_species_per_original_disjoint_role"])
    return {
        "status":("PASS_SOURCE_ONLY_MODERN_MATCHED_HISTORY_CONTRAST_GEO"
                  if passed else pre["if_insufficient"]),
        "origin_pair_comparisons_within_same_garden":all_pair,
        "origin_pairs_modern_climate_distance_le_1":passmodern,
        "origin_pairs_history_magnitude_delta_ge_half_SD":passhistory,
        "origin_pairs_both_criteria":passboth,
        "unique_matched_study_species_origins_fitness_type_pairs_deduped_across_gardens":len(paired),
        "matched_study_species_garden_fitness_type_cells":len(matched),
        "matched_species":len(cnt_species),"matched_studies":len(cnt_studies),
        "matched_genera":len({sp.split()[0] for sp in cnt_species}),
        "matched_species_per_original_role":matched_role,
        "matched_pairs_by_fitness_type":dict(sorted(per_type.items())),
        "climate_finite_comparable_cells_by_role":dict(sorted(count_cells.items())),
        "source_training_0BP_four_bioclim_SD":source_sd.tolist(),
        "source_training_history_magnitude_SD":hsd,
        "within_garden_history_residual_variance_ratio_after_7_modern_geographic_controls":frac,
        "within_garden_source_H_unadjusted_weighted_ss":variance,
        "within_garden_source_H_residual_weighted_ss":num,
        "within_garden_covariate_matrix_rank":int(np.linalg.matrix_rank(xx)),
        "within_garden_covariate_dim":7,
        "within_garden_source_comparison_rows_in_train":len(ys),
        "historical_information_fraction_is_NOT_fitness_R2_or_adaptation_effect":True,
        "pairwise_pseudoreplication_not_resolved_by_pair_count":True
    }

def run():
    ap=argparse.ArgumentParser()
    ap.add_argument("--work",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    c=json.loads(CONTRACT.read_text());orig=json.loads(GEO.read_text())
    lock=json.loads(ZIP_LOCK.read_text());assets=json.loads(ASSETS.read_text())
    if c.get("status")!="FROZEN_BEFORE_107PLANT_ORIGINAL_CHELSA_8BAND_SOURCE_SAMPLING":
        raise ValueError("No original 107-plant pre-climate source authority")
    out={"schema":"ttf_Bontrager_107plant_source_only_modern_matched_LGM_history_contrast_v01",
        "fitness_sign_magnitude_sample_rows_NEVER_accessed":True,
        "no_phenotype_fitness_or_response_fit":True,
        "no_reuse_of_PROCEED_56plant_stopped_Haldane_response":True,
        "no_new_phylogeographic_or_speciation_inference":True,
        "source_only_status_not_confirmatory_evolutionary_prediction":True,
        "raw_sites_and_climate_not_exported":True}
    args.work.mkdir(parents=True,exist_ok=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    try:
        q=Request(ORIGINAL_URL,headers={"User-Agent":"TTF-Bontrager-climate-source-only-v01"})
        with urlopen(q,timeout=90) as resp:raw=resp.read(lock["max_download_bytes"]+1)
        verify_original_zip(raw,lock)
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            with z.open(c["original_source"]["original_member"]) as f:
                studycsv=f.read(25_000_000)
        if len(studycsv)>=25_000_000:raise ValueError("Original source CSV over cap")
        groups=set(c["original_source"]["original_functional_groups"])
        check=metadata_census(studycsv,orig,allowed_functional_groups=groups)
        authority=c["original_source"]
        if (check["qualified_binomial_taxa"]!=authority["source_nominal_plant_taxa"] or
           check["qualified_original_source_studies"]!=authority["source_original_studies"] or
           check["shared_garden_fitness_type_cells_with_multiple_source_IDs_and_coordinates"]!=authority["source_original_garden_fitness_type_cells_with_two_original_origins"] or
           check["study_species_component_split"]["source_only_membership_sha256"]!=authority["source_study_species_disjoint_membership_sha256"]):
            raise ValueError("Original plant source-only geometry not reproduced")
        cells,gardens,graph=extract_source_cells(studycsv,orig,groups)
        roles,digest=roles_for_graph(graph)
        if digest!=authority["source_study_species_disjoint_membership_sha256"] or len(cells)!=711:
            raise ValueError("Exact original plant cells or roles differ from source-only freeze")
        sites=sorted({v for values in cells.values() for v in values.values()} |
                     set(gardens.values()))
        xy=np.asarray(sites,float)
        if len(sites)!=len(set(sites)) or xy.shape[1]!=2 or not np.isfinite(xy).all():
            raise ValueError("Source locations malformed")
        fields=assets["assets"]
        if len(fields)!=8:
            raise ValueError("Original eight SHA-locked CHELSA raster assets changed")
        names=[f"{f['variable']}_{f['time_index']}" for f in fields]
        if names!=["bio01_-190","bio07_-190","bio12_-190","bio15_-190",
                  "bio01_20","bio07_20","bio12_20","bio15_20"]:
            raise ValueError("Original 21ka/0BP CHELSA variable order changed")
        extracted=[];receipts=[]
        for i,field in enumerate(fields):
            dest=args.work/f"Bontrager107_CHELSA_{i}.tif"
            verified_raster(field,dest)
            values=raster_values(dest,xy)
            dest.unlink()
            extracted.append(values)
            receipts.append({"band":names[i],"sha256":field["sha256"],
                             "bytes":field["size_bytes"],
                             "finite_unique_source_or_garden_sites":int(np.isfinite(values).sum())})
            print("BONTRAGER_CHELSEA_SOURCE_ORIGINAL_SHA_OK",names[i],flush=True)
        matrix=np.column_stack(extracted)
        if matrix.shape!=(len(sites),8):raise ValueError("CHELSA original source field geometry mismatch")
        out.update(matched_source_only(cells,gardens,roles,sites,matrix,c))
        out.update({
           "original_plant_taxa_source_qualified":len(roles),
           "original_plant_garden_cells_source_qualified":len(cells),
           "unique_source_and_destination_original_geographies":len(sites),
           "unique_fully_finite_eight_band_sites":int(np.isfinite(matrix).all(axis=1).sum()),
           "original_site_lonlat_sha256":sha_float64(xy),
           "original_8_bioclim_numeric_matrix_sha256":sha_float64(matrix),
           "original_study_species_disjoint_membership_sha256":digest,
           "original_source_zip_sha256":lock["dryad_sha256"],
           "original_source_bioclim_sha_receipts":receipts
        })
    except Exception as e:
        out.update(status="NOT_EVALUABLE_SOURCE_CHELSA_OR_MATCHED_HISTORY_METADATA",
                   error_type=type(e).__name__,error_message=str(e)[:440])
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("BONTRAGER_107PLANT_HISTORY_GEOMETRY",out["status"],
          out.get("matched_species"),out.get("matched_studies"),
          out.get("matched_study_species_garden_fitness_type_cells"))
if __name__=="__main__":run()
