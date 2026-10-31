"""2026.10.31.11: end-to-end safety and native failure visibility."""
from pathlib import Path
import json
import subprocess

import pytest

from openvsppy.vv.VVNative import (
    classify_control_surfaces, inspect_native_mesh, diagnose_solver_log)
from openvsppy.vv.VVRunner import select_case_timeout, wait_for_worker


def test_control_group_missing_entire_geometry_is_detected():
    surfaces=['HTailGeom_Surf0_ELEVATOR','HTailGeom_Surf1_ELEVATOR']
    present, absent=classify_control_surfaces(surfaces, {'WingGeom','FuselageGeom'},
                                              {'WingGeom','FuselageGeom','HTailGeom'})
    assert not present and absent==surfaces
    present, absent=classify_control_surfaces(surfaces, {'HTailGeom','WingGeom'},
                                              {'WingGeom','FuselageGeom','HTailGeom'})
    assert present==surfaces and not absent


def test_partial_control_group_does_not_silently_pass():
    surfaces=['WingGeom_Surf0_AILERON','WingTipGeom_Surf0_AILERON']
    present,absent=classify_control_surfaces(surfaces,{'WingGeom'}, {'WingGeom','WingTipGeom'})
    assert len(present)==len(absent)==1
    with pytest.raises(ValueError,match='parent geometry'):
        classify_control_surfaces(['Unknown_Surf0_ELEVATOR'], {'WingGeom'}, {'WingGeom'})


def test_mesh_preflight_detects_observed_nref_pathology(tmp_path):
    geom=tmp_path/'m.vspgeom'
    key=tmp_path/'m.vkey'
    key.write_text('# VSPGEOM v3 Tag Key File\nm.vspgeom\n6\n')
    geom.write_text('# vspgeom v3\n1\n5498 7 4\n')
    result=inspect_native_mesh(geom,key)
    assert result['status']=='INVALID'
    assert '5498 nodes' in result['issues'][0]
    geom.write_text('# vspgeom v3\n1\n5501 5243 4\n')
    result=inspect_native_mesh(geom,key)
    assert result['status']=='PASS' and result['vkey_parts']==6


def test_observed_failure_stages(tmp_path):
    log=tmp_path/'vspaero.log'
    log.write_text('Updating wake grids...\nDone\n')
    d=diagnose_solver_log(log)
    assert d['classification']=='SOLVER_STOPPED_AFTER_WAKE_GRID_UNEXPLAINED'
    log.write_text('wtf ! Loop 1 is not associated with any component!')
    assert diagnose_solver_log(log)['classification']=='INVALID_MESH_COMPONENT'
    log.write_text('Could not find control surface: HTailGeom_Surf0_ELEVATOR')
    assert diagnose_solver_log(log)['classification']=='INVALID_CONTROL_SURFACE'


def test_case_timeout_overrides_and_rejects_bad_values():
    assert select_case_timeout('STA-01-ADJOINT','ADJOINT',3600,
        {'ADJOINT':10800})==10800
    assert select_case_timeout('STA-01-ADJOINT','ADJOINT',3600,
        {'ADJOINT':10800},{'STA-01-ADJOINT':7200})==7200
    assert select_case_timeout('GRID-01-COARSE','STEADY_SWEEP',3600,
        {'ADJOINT':10800})==3600
    for invalid in (True,0,-1,'7200'):
        with pytest.raises(ValueError,match='Invalid timeout'):
            select_case_timeout('X','ADJOINT',3600,{'ADJOINT':invalid})


def test_adjoint_progress_with_cr_delimited_native_log(tmp_path):
    log=tmp_path/'vspaero_api_stdout.log'
    events=[]
    class FakeProc:
        def __init__(self):self.calls=0
        def wait(self,timeout):
            self.calls+=1
            if self.calls==1:
                log.write_bytes(b'Wake Iter: 20 / 20\rSolving adjoint linear systems ...\r'
                                b'Inviscid Force in x...\rGMRES Iter: 0 ...\r'
                                b'GMRES Iter: 12 ...\r')
                raise subprocess.TimeoutExpired(['solver'],timeout)
            return 0
    assert wait_for_worker(FakeProc(),['solver'],1000,log,verbose=True,
        progress_callback=events.append,poll_seconds=.01,heartbeat_seconds=2)==0
    assert any(e['event']=='WAKE' and e['iteration']==20 for e in events)
    assert any(e['event']=='ADJOINT_STAGE' and e['stage']=='Inviscid Force in x' for e in events)
    assert any(e['event']=='ADJOINT_GMRES' and e['iteration']==12 for e in events)


def test_manifest_default_only_changes_runtime_timeout():
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((root/'vv_study/execution.reviewed.json').read_text())
    assert manifest['timeouts_seconds_by_mode']['ADJOINT']==10800


def test_native_worker_exception_keeps_preflight_reason():
    from openvsppy.vv.VVNative import classify_worker_error
    fallback={'classification':'NATIVE_OUTPUT_MISSING_UNEXPLAINED'}
    assert classify_worker_error('RuntimeError: VSPAERO mesh preflight failed: too few loops',
        fallback)=='MESH_PREFLIGHT_FAILED'
    assert classify_worker_error('RuntimeError: Nonzero deflection for omitted Control Group: X',
        fallback)=='OMITTED_CONTROL_GROUP_WITH_DEFLECTION'


def test_case_timeout_reaches_worker_without_changing_case(study,configured_case,monkeypatch):
    from openvsppy.vv.VVRunner import run_cases
    root,db,_=configured_case
    calls=[]
    def fake_run(*args,**kwargs):
        calls.append(kwargs['timeout_seconds'])
        return {'status':'EXECUTION_FAILED','run_id':'RUNTIME_TEST'}
    monkeypatch.setattr('openvsppy.vv.VVRunner.run_case',fake_run)
    result=run_cases(db.path,root/'study.json',['GRID-01-MEDIUM'],
        binary_path=Path(__file__).resolve().parents[1]/'README.md',
        resume=False,timeout_seconds=3600,
        timeouts_by_mode={'ADJOINT':10800},
        timeouts_by_case={'GRID-01-MEDIUM':14400})
    assert calls==[14400]
    assert result[0]['status']=='EXECUTION_FAILED'
