"""Post-terminal numeric identity diagnosis; never reopen a biological outcome.

The existing source-only geometry is reconstructed once. Separate processes
compare runtime CPU dispatch, exact hashes and aggregate byte differences.
No alternate source hash, rounding tolerance, fitted response or random world
is accepted. Input arrays remain private temporary files, never artifacts.
"""
from __future__ import annotations
import argparse
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import numpy as np

HERE=Path(__file__).resolve().parent
ORIGINAL_COMMIT='a4664cc1fccf7a16fed1b85f239dc0f123389a87'
CLOSED='TERMINAL_NOT_EVALUABLE_H_SYNTHETIC_NULL_CEILING_FAIL_NO_BIOLOGICAL_OUTCOME_OPENING'
MODULES=(
 'proceed70_haldane_source_null_calibration.py',
 'proceed70_haldane_rate_model_core.py',
 'proceed70_haldane_presence_mask.py',
 'proceed70_plant_history_information.py',
 'proceed_plant_chelsa_source_support.py',
 'proceed_v62_plant_genetic_source_census.py',
 'proceed_v62_metadata_row_census.py',
 'proceed_v62_genetic_geometry_intersection.py',
 'proceed_v62_borealis_header_only.py',
 'site_feature_core.py',
)


def array_hash(a):
    return hashlib.sha256(np.asarray(a,dtype='<f8').tobytes(order='C')).hexdigest()


def numeric_difference(a,b):
    a=np.asarray(a,dtype='<f8');b=np.asarray(b,dtype='<f8')
    if a.shape!=b.shape or a.ndim not in (1,2):
        raise ValueError('aligned one- or two-dimensional arrays required')
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise ValueError('nonfinite numeric source values')
    aa=np.ascontiguousarray(a);bb=np.ascontiguousarray(b)
    different=aa.view('<u8')!=bb.view('<u8')
    delta=np.abs(a-b)
    columns=np.flatnonzero(different.any(axis=0)) if a.ndim==2 else ([] if not different.any() else [0])
    return {
      'shape':list(a.shape),'exact_byte_identity':bool(not different.any()),
      'bitwise_different_elements':int(different.sum()),
      'different_columns':[int(i) for i in columns],
      'max_absolute_difference':float(delta.max(initial=0.)),
      'left_sha256':array_hash(a),'right_sha256':array_hash(b),
      'raw_values_exported':False,
    }


def exact_pair_matches(x,h,expected_x,expected_h):
    return x==expected_x and h==expected_h


def runtime_metadata():
    from numpy._core import _multiarray_umath as um
    buffer=io.StringIO()
    with contextlib.redirect_stdout(buffer):np.show_runtime()
    cpu=[]
    p=Path('/proc/cpuinfo')
    if p.exists():
        cpu=sorted({s.split(':',1)[1].strip() for s in p.read_text().splitlines()
                    if s.startswith('model name')})
    try:
        from numpy.lib.introspect import opt_func_info
        dispatch=opt_func_info(func_name='sin|cos|arctan2|arcsin|sqrt|log1p',signature='float64')
    except (ImportError,TypeError,AttributeError) as exc:dispatch={'unavailable':type(exc).__name__}
    return {'python':platform.python_version(),'numpy':np.__version__,
            'machine':platform.machine(),'cpu_models':cpu,
            'compiled_dispatch_targets':list(um.__cpu_dispatch__),
            'disabled_cpu_features':os.environ.get('NPY_DISABLE_CPU_FEATURES',''),
            'numpy_runtime':buffer.getvalue(), 'ufunc_dispatch':dispatch,
            'openblas_threads':os.environ.get('OPENBLAS_NUM_THREADS')}


def code_identity():
    root=HERE.parents[1]
    results={}
    for name in MODULES:
        rel='research/prospective/'+name
        old=subprocess.run(['git','show',ORIGINAL_COMMIT+':'+rel],cwd=root,
                           capture_output=True,check=True).stdout
        current=(HERE/name).read_bytes()
        results[name]={'original_sha256':hashlib.sha256(old).hexdigest(),
                       'current_sha256':hashlib.sha256(current).hexdigest(),
                       'identical':old==current}
    return results


def prepare_sources(work,receipt):
    from proceed_v62_borealis_header_only import validated_get_catalog,verify_frozen_catalog
    from proceed_v62_metadata_row_census import fetch_original_datafile_bytes
    from proceed_v62_genetic_geometry_intersection import SOURCE_CONTRACT,NATIVE_LOCK
    from proceed_plant_chelsa_source_support import source_geometry,verified_raster
    from site_feature_core import raster_values
    lock=receipt['source_pre_registered']
    catalog=validated_get_catalog('https://borealisdata.ca/api/datasets/:persistentId/?persistentId=doi%3A10.5683%2FSP3%2FNXSL3Q')
    verify_frozen_catalog(catalog,json.loads(SOURCE_CONTRACT.read_text())['current_source'])
    raw,md5,sha=fetch_original_datafile_bytes(1024206)
    native=json.loads(NATIVE_LOCK.read_text())['native_original_transport']
    if sha!=lock['source_proceed_sha256'] or md5!=native['observed_md5'] or len(raw)!=native['observed_size_bytes']:
        raise ValueError('original native PROCEED bytes changed')
    design=json.loads((HERE/'proceed_plant_chelsa_source_support_contract_v01.json').read_text())
    cohort,_=source_geometry(raw,design)
    coords=sorted({v for _,_,a,b in cohort for v in (a,b)})
    xy=np.asarray(coords,float)
    if array_hash(xy)!='aed6486ee9aa49b0d454426c4d5964ebd79a8a96db2aac452c80fd5b7f52af32':
        raise ValueError('original coordinate identity failed')
    assets=json.loads((HERE/'chelsa_trace21k_fixed8_assets_v01.json').read_text())['assets']
    columns=[]
    for i,asset in enumerate(assets):
        path=work/f'climate_{i}.tif';verified_raster(asset,path)
        columns.append(raster_values(path,xy));path.unlink()
        print('SOURCE_RASTER_SHA_VERIFIED',i+1,flush=True)
    matrix=np.column_stack(columns)
    if array_hash(matrix)!=lock['source_CHELSA8_raw_numeric_SHA256']:
        raise ValueError('original climate matrix identity failed')
    # Ephemeral input is never uploaded; original response fields are not decoded.
    source=work/'source.csv';source.write_bytes(raw);source.chmod(0o600)
    np.savez(work/'inputs.npz',coords=xy,climate=matrix)
    return {'native_source_sha256':sha,'coordinate_sha256':array_hash(xy),
            'climate_matrix_sha256':array_hash(matrix),'original_sites':len(xy)}


def worker(work,label):
    from proceed70_haldane_source_null_calibration import source_only_geometry,source_history_concentration
    design=json.loads((HERE/'proceed_plant_chelsa_source_support_contract_v01.json').read_text())
    with np.load(work/'inputs.npz',allow_pickle=False) as source:
        xy=source['coords'];climate=source['climate']
    coords=[tuple(row) for row in xy]
    X,H,sp,st,role,sd,_=source_only_geometry((work/'source.csv').read_bytes(),climate,coords,design)
    before=(array_hash(X),array_hash(H))
    # This original helper preceded hashes in the frozen run; test mutation.
    info=source_history_concentration(H,X,sp,role)
    after=(array_hash(X),array_hash(H))
    np.savez(work/(label+'.npz'),X=X,H=H,sp=sp,st=st,role=role)
    report={'runtime':runtime_metadata(),'X_sha256':before[0],'H_sha256':before[1],
            'X_column_sha256':[array_hash(X[:,i]) for i in range(X.shape[1])],
            'array_shapes':{'X':list(X.shape),'H':list(H.shape)},
            'post_information_call_hashes':{'X':after[0],'H':after[1]},
            'information_call_mutates_inputs':before!=after,
            'source_current_climate_SD_sha256':array_hash(sd),
            'source_roles_sha256':hashlib.sha256(json.dumps(
                 [sp.tolist(),st.tolist(),role.tolist()],separators=(',',':')).encode()).hexdigest()}
    (work/(label+'.json')).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')


def main(output):
    receipt=json.loads((HERE/'proceed70_haldane_synthetic_terminal_receipt_v01.json').read_text())
    if receipt['status']!=CLOSED:raise ValueError('terminal closure required')
    lock=receipt['source_pre_registered']
    expected_x=lock['source_21_baseline_numeric_SHA256'];expected_h=lock['source_one_history_numeric_SHA256']
    report={'schema':'ttf_postterminal_numeric_provenance_v01',
            'terminal_status_unchanged':CLOSED,'real_phenotype_rate_values_read':False,
            'additional_null_or_positive_worlds_drawn':0,'scientific_thresholds_changed':False,
            'source_hashes_replaced_or_tolerance_relaxed':False,
            'expected_X_sha256':expected_x,'expected_H_sha256':expected_h,
            'original_source_commit':ORIGINAL_COMMIT}
    output.parent.mkdir(parents=True,exist_ok=True)
    try:
        report['source_code_identity']=code_identity()
        if not all(v['identical'] for v in report['source_code_identity'].values()):
            raise ValueError('source pipeline changed since original run; inspect diff first')
        from numpy._core import _multiarray_umath as um
        dispatch=list(um.__cpu_dispatch__)
        wide=[s for s in dispatch if 'AVX512' in s or s=='X86_V4']
        modes=[('native',''),('native_repeat',''),('without_avx512',','.join(wide)),
               ('baseline_dispatch',','.join(dispatch))]
        with tempfile.TemporaryDirectory(prefix='ttf-private-numeric-') as directory:
            work=Path(directory)
            report['source_identity']=prepare_sources(work,receipt)
            report['modes']={};report['comparisons_to_native']={}
            for label,disabled in modes:
                env=dict(os.environ);env['NPY_DISABLE_CPU_FEATURES']=disabled
                proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',str(work),
                                     '--label',label],env=env,capture_output=True,text=True,timeout=90)
                if proc.returncode:
                    report['modes'][label]={'worker_failed':True,'stderr_tail':proc.stderr[-1500:]}
                    continue
                result=json.loads((work/(label+'.json')).read_text())
                result['matches_frozen_X']=result['X_sha256']==expected_x
                result['matches_frozen_H']=result['H_sha256']==expected_h
                result['exact_reference_pair']=exact_pair_matches(result['X_sha256'],result['H_sha256'],expected_x,expected_h)
                report['modes'][label]=result
                print('NUMERIC_REPLAY',label,result['matches_frozen_X'],result['matches_frozen_H'],flush=True)
            if not (work/'native.npz').exists():raise ValueError('native worker failed')
            with np.load(work/'native.npz',allow_pickle=False) as base:
                for label,_ in modes[1:]:
                    path=work/(label+'.npz')
                    if path.exists():
                        with np.load(path,allow_pickle=False) as candidate:
                            report['comparisons_to_native'][label]={name:numeric_difference(base[name],candidate[name]) for name in ('X','H')}
            matched=[label for label,data in report['modes'].items() if data.get('exact_reference_pair')]
            report['exact_reference_modes']=matched
            report['status']='EXACT_ORIGINAL_NUMERIC_REPLAY_FOUND' if matched else 'DERIVED_NUMERIC_IDENTITY_UNRESOLVED_TERMINAL_UNCHANGED'
            # A separate post-terminal analytical diagnostic is allowed ONLY
            # on exact original input bytes, never as a replacement gate.
            if matched:
                from terminal_analytic_audit import analyse
                with np.load(work/(matched[0]+'.npz'),allow_pickle=False) as a:
                    report['closed_form_diagnosis_on_exact_reference']=analyse(a['X'],a['H'],a['sp'],a['st'],a['role'])
                report['analytical_diagnostic_not_biological_result']=True
    except Exception as exc:
        report['status']='NUMERIC_PROVENANCE_ERROR_TERMINAL_UNCHANGED'
        report['error']=type(exc).__name__+': '+str(exc)[:500]
        output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
        raise
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('NUMERIC_PROVENANCE_STATUS',report['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path)
    p.add_argument('--worker',type=Path);p.add_argument('--label')
    args=p.parse_args()
    if args.worker:worker(args.worker,args.label)
    elif args.output:main(args.output)
    else:p.error('--output or --worker required')
