import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).parents[1]/'scripts'))
import acquire_lgm_external_references as a


def feature(name):
    assert hasattr(a,name), f'missing transport recovery feature: {name}'
    return getattr(a,name)


def result(name, status='INSUFFICIENT_IN_RETRIEVED_REFERENCE'):
    return ({'species':name,'status':status,'external_reference_points_50km':0},[],[{'url':'https://example.org','status':200,'response_sha256':'a'*64}])


def test_species_checkpoint_written_before_later_worker_fails(tmp_path):
    run=feature('checkpointed_acquisition')
    calls=[]
    def acquire(name,support):
        calls.append(name)
        if name=='B': raise RuntimeError('simulated interruption')
        return result(name)
    with pytest.raises(RuntimeError,match='simulated interruption'):
        run(['A','B','C'],{x:None for x in 'ABC'},tmp_path,{'source':'fixed'},acquire,workers=1)
    assert calls==['A','B']
    assert len(list(tmp_path.glob('*.json')))==1
    calls.clear()
    def resume(name,support):calls.append(name);return result(name)
    rows,pending=run(['A','B','C'],{x:None for x in 'ABC'},tmp_path,{'source':'fixed'},resume,workers=1)
    assert calls==['B','C']
    assert [x[0]['species'] for x in rows]==['A','B','C'] and pending==[]


def test_no_automatic_retry_of_recorded_transport_error(tmp_path):
    run=feature('checkpointed_acquisition')
    calls=[]
    def acquire(name,support):calls.append(name);return result(name,'REQUEST_ERROR')
    first,_=run(['A'],{'A':None},tmp_path,{'source':'fixed'},acquire,workers=1)
    second,_=run(['A'],{'A':None},tmp_path,{'source':'fixed'},acquire,workers=1)
    assert first==second and calls==['A']


def test_bounded_batch_retains_pending_as_pending(tmp_path):
    run=feature('checkpointed_acquisition')
    rows,pending=run(['A','B','C'],dict.fromkeys('ABC'),tmp_path,{},lambda n,s:result(n),workers=1,max_new_species=2)
    assert [x[0]['species'] for x in rows]==['A','B']
    assert pending==['C']
    rows,pending=run(['A','B','C'],dict.fromkeys('ABC'),tmp_path,{},lambda n,s:result(n),workers=1,max_new_species=2)
    assert len(rows)==3 and not pending


def test_reject_checkpoint_under_changed_input_binding(tmp_path):
    run=feature('checkpointed_acquisition')
    run(['A'],{'A':None},tmp_path,{'source':'old'},lambda n,s:result(n),workers=1)
    with pytest.raises(ValueError,match='binding'):
        run(['A'],{'A':None},tmp_path,{'source':'new'},lambda n,s:result(n),workers=1)


def test_tampered_result_rejected(tmp_path):
    run=feature('checkpointed_acquisition')
    run(['A'],{'A':None},tmp_path,{},lambda n,s:result(n),workers=1)
    p=next(tmp_path.glob('*.json')); j=json.loads(p.read_text());j['result'][0]['status']='faked'
    p.write_text(json.dumps(j))
    with pytest.raises(ValueError,match='digest'):
        run(['A'],{'A':None},tmp_path,{},lambda n,s:result(n),workers=1)


def test_refuse_species_identity_mismatch(tmp_path):
    run=feature('checkpointed_acquisition')
    with pytest.raises(ValueError,match='species'):
        run(['A'],{'A':None},tmp_path,{},lambda n,s:result('WRONG'),workers=1)
    assert not list(tmp_path.glob('*.json'))


def test_clean_resume_bytes_stable(tmp_path):
    run=feature('checkpointed_acquisition')
    names=['A','B','C'];support=dict.fromkeys(names)
    run(names,support,tmp_path,{},lambda n,s:result(n),workers=2)
    before={p.name:p.read_bytes() for p in tmp_path.glob('*.json')}
    def forbidden(*args):raise AssertionError('network called on complete checkpoint')
    rows,pending=run(names,support,tmp_path,{},forbidden,workers=2)
    assert len(rows)==3 and pending==[]
    assert before=={p.name:p.read_bytes() for p in tmp_path.glob('*.json')}


def test_reject_duplicate_species(tmp_path):
    run=feature('checkpointed_acquisition')
    with pytest.raises(ValueError,match='duplicate'):
        run(['A','A'],{'A':None},tmp_path,{},lambda n,s:result(n),workers=1)


def test_opaque_names_do_not_escape_checkpoint_directory(tmp_path):
    run=feature('checkpointed_acquisition')
    run(['../species'],{'../species':None},tmp_path,{},lambda n,s:result(n),workers=1)
    assert len(list(tmp_path.glob('*.json')))==1


def test_all_panel_export_has_unattempted_status_not_zero_failure():
    export=feature('panel_summary_rows')
    rows=export(['A','B'],[result('A')])
    assert rows[1]=={'species':'B','status':'NOT_ATTEMPTED_IN_CURRENT_SNAPSHOT'}
    assert 'external_reference_points_50km' not in rows[1]
