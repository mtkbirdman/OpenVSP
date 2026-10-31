"""NCPU runtime settings: API request, Run provenance and resumption safety."""
import json

import pytest

from openvsppy.vv.VVRunner import resolve_ncpu, run_case, run_cases, check_case
from openvsppy.vv.__main__ import positive_ncpu, main
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_ncpu_resolution_default_case_cli_override():
    assert resolve_ncpu({}) == 1
    assert resolve_ncpu({'vspaero': {'NCPU': 2}}) == 2
    assert resolve_ncpu({'vspaero': {'NCPU': 2}}, 4) == 4
    assert positive_ncpu('16') == 16


@pytest.mark.parametrize('invalid', [0, -1, 1.5, '4', None, True, False])
def test_ncpu_case_rejects_non_positive_or_non_integer(invalid):
    with pytest.raises(ValueError, match='NCPU'):
        resolve_ncpu({'vspaero': {'NCPU': invalid}})


@pytest.mark.parametrize('invalid', ['0', '-2', 'NaN', '2.5'])
def test_ncpu_cli_validation(invalid):
    with pytest.raises(SystemExit) as result:
        main(['run', '--case', 'GRID-01-MEDIUM', '--ncpu', invalid])
    assert result.value.code == 2


def test_case_json_ncpu_preflight_rejected(study, configured_case):
    root,db,_ = configured_case
    db.add_setting('VSP3505-VV-001', {
        'setting_id':'NCPU-INVALID','analysis_mode':'STEADY_SWEEP','mesh_level':'MEDIUM',
        'requested':{'thin_set_index':3,'thick_set_index':4,'vspaero':{'NCPU':0}}
    })
    db.resolve_case('GRID-01-MEDIUM','G3','FC-001','NCPU-INVALID')
    assert any('NCPU' in msg for msg in check_case(db,'GRID-01-MEDIUM',root/'study.json',ROOT/'README.md'))


def test_run_records_requested_ncpu_without_modifying_case(study, configured_case, monkeypatch):
    root,db,_=configured_case
    def unavailable(*args, **kwargs):
        raise FileNotFoundError('offline worker deliberately unavailable')
    monkeypatch.setattr('openvsppy.vv.VVRunner.subprocess.Popen',unavailable)
    input_json=db.row('numerical_setting','setting_id','SET1')['requested_json']
    records=[]
    for ncpu in (1,4):
        result=run_case(db.path,'GRID-01-MEDIUM',root/'study.json',
                        binary_path=ROOT/'README.md', ncpu=ncpu)
        assert result['status']=='EXECUTION_FAILED'
        row=db.row('run','run_id',result['run_id'])
        env=json.loads(row['environment_json'])
        snapshot=json.loads(row['requested_snapshot_json'])
        job=json.loads((root/row['run_dir']/'input/job.json').read_text())
        assert env['ncpu']==ncpu
        assert snapshot['runtime']['ncpu']==ncpu
        assert job['runtime']['ncpu']==ncpu
        assert job['setting']['requested']==json.loads(input_json)
        records.append(row['fingerprint_sha256'])
    assert records[0] != records[1]
    assert db.row('numerical_setting','setting_id','SET1')['requested_json']==input_json
    assert all(v=='VERIFIED' for _,v in db.verify_artifacts(root))


def test_resumption_never_reuses_other_ncpu(study,configured_case,monkeypatch):
    root,db,_=configured_case
    binary=ROOT/'README.md'
    runid, rundir=db.register_run('GRID-01-MEDIUM',root,binary,ROOT/'src/openvsppy/vv/VVRunner.py',
        {'test':'prior-succeeded'}, ncpu=1)
    # Insert a prior successful linked Point and artifact so it WOULD be reusable
    # if NCPU were not considered.  The test must not accept another NCPU.
    sample=rundir/'dummy.txt';sample.write_text('evidence')
    db.register_artifact(root,'VSP3505-VV-001',sample,run_id=runid)
    # Work with a stubbed candidate: the matcher should only consider NCPU=4.
    def fake_run(*args,**kwargs):
        return {'status':'BLOCKED','run_id':'NOT_REUSED','ncpu':kwargs['ncpu']}
    monkeypatch.setattr('openvsppy.vv.VVRunner.run_case',fake_run)
    with db.transaction() as c:
        c.execute("UPDATE run SET status='SUCCEEDED' WHERE run_id=?",(runid,))
        c.execute("INSERT INTO case_point (point_id,run_id,ordinal,target_condition_id,observed_condition_json,point_status) VALUES(?,?,?,?,?,?)", ('PT-FAKE',runid,0,'FC-001','{}','MATCHED'))
        c.execute("INSERT INTO case_run_link VALUES(?,?,?,?)",('GRID-01-MEDIUM','PT-FAKE','2026-10-31','MOCK'))
    same=run_cases(db.path,root/'study.json',['GRID-01-MEDIUM'],
                   binary_path=binary,ncpu=1,resume=True)
    assert same[0]['status']=='REUSED'
    results=run_cases(db.path,root/'study.json',['GRID-01-MEDIUM'],
                      binary_path=binary,ncpu=4,resume=True)
    assert results[0]['status']=='BLOCKED'
    assert results[0]['ncpu']==4


def test_cli_threads_reach_run_batch_and_execute_all(study,configured_case,monkeypatch,capsys):
    root,db,_=configured_case
    from openvsppy.vv import __main__ as cli
    forwarded=[]
    def fake_run(*args, **kwargs):
        forwarded.append(('run',kwargs.get('ncpu')))
        return {'status':'SUCCEEDED','run_id':'TEST'}
    def fake_batch(*args, **kwargs):
        forwarded.append(('batch',kwargs.get('ncpu')))
        return [{'status':'SUCCEEDED','case_id':'GRID-01-MEDIUM'}]
    def fake_all(*args, **kwargs):
        forwarded.append(('execute-all',kwargs.get('ncpu')))
        return {'status':'REVIEW_REPORT_GENERATED'}
    monkeypatch.setattr(cli,'run_case',fake_run)
    monkeypatch.setattr(cli,'run_cases',fake_batch)
    monkeypatch.setattr(cli,'execute_study',fake_all)
    common=['--study',str(root/'study.json'),'--db',str(db.path)]
    assert main([*common,'run','--case','GRID-01-MEDIUM','--ncpu','4'])==0
    assert main([*common,'run-test','--test','GRID-01','--binary',str(ROOT/'README.md'),'--ncpu','4'])==0
    assert main([*common,'execute-all','--config',str(root/'manifest.json'),'--ncpu','4'])==0
    assert forwarded==[('run',4),('batch',4),('execute-all',4)]
