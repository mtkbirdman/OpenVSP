"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path
import shutil

import pytest

from openvsppy.vv.VVDatabase import VVDatabase
from openvsppy.vv.VVDatabase import execution_fingerprint
from openvsppy.vv.VVDatabase import postprocess_fingerprint
from openvsppy.vv.VVPostprocess import import_polar
from openvsppy.vv.VVPostprocess import import_stab
from openvsppy.vv.VVPostprocess import insert_dataset
from openvsppy.vv.VVRunner import check_case

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'tests/fixtures/SampleGlider_native_3.50.5.stab'
MAPPING = {'ConGrp_1':'AILERON_GROUP','ConGrp_2':'ELEVATOR_GROUP','ConGrp_3':'RUDDER_GROUP'}


def test_atomic_success_and_case_link_rollback(tmp_path):
    """A failed Case-link INSERT must not leave a completed Run."""
    root=tmp_path/'study';root.mkdir()
    db=VVDatabase(root/'vv.sqlite')
    db.initialize(ROOT/'vv_study/study.json', ROOT/'vv_study/01_database/schema.sql')
    sid='VSP3505-VV-001'
    db.add_geometry(sid,{'geometry_id':'A','stage':'G3','vsp3_source':str(ROOT/'examples/models/G103A/G103A.vsp3')})
    db.add_condition(sid,{'condition_id':'C','mach':0.,'re_cref':2e6,'alpha_deg':0.,'beta_deg':0.,
        'controls_deg':{},'trim_mode':'FIXED_CONTROL','sref':18.,'bref':27.,'cref':.625,
        'moment_x':1.4,'moment_y':0.,'moment_z':0.,'reference_length_unit':'m',
        'axes_convention_id':'NATIVE_UNVERIFIED'})
    db.add_setting(sid,{'setting_id':'SET','analysis_mode':'STEADY','mesh_level':'BASE',
        'requested':{'thin_set_index':4,'thick_set_index':3}})
    db.resolve_case('STA-01-STEADY','A','C','SET')
    rid,folder=db.register_run('STA-01-STEADY',root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'case':'atomic'})
    for state in ('PRECHECKED','RUNNING','CAPTURED'):db.transition(rid,state)
    file=folder/'raw'/'model_effective.stab';file.parent.mkdir();shutil.copyfile(FIXTURE,file)
    with db.transaction():
        point=import_stab(db,study_root=root,run_id=rid,path=file,control_mapping=MAPPING)['point_id']
        db.transition(rid,'PARSED')
    with pytest.raises(Exception):
        with db.transaction():
            db.transition(rid,'SUCCEEDED')
            db.conn.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',('INVALID-CASE',point,'2026-10-31','TEST'))
    assert db.row('run','run_id',rid)['status']=='PARSED'
    assert db.conn.execute('SELECT COUNT(*) FROM case_run_link').fetchone()[0]==0
    with db.transaction():
        db.transition(rid,'SUCCEEDED')
        db.conn.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',('STA-01-STEADY',point,'2026-10-31','TEST'))
    assert db.row('run','run_id',rid)['status']=='SUCCEEDED'
    assert db.conn.execute('SELECT COUNT(*) FROM case_run_link').fetchone()[0]==1
    db.close()


def test_new_schema_requires_new_database_and_preserves_study(study):
    root, db = study
    assert db.conn.execute('SELECT count(*) FROM planned_case').fetchone()[0] == 68
    assert db.conn.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'").fetchone()[0] == 19
    before = db.conn.execute('SELECT * FROM study').fetchall()
    with pytest.raises(ValueError, match='Already initialized'):
        db.initialize(root/'study.json', root/'01_database/schema.sql')
    assert db.conn.execute('SELECT * FROM study').fetchall() == before


def test_source_fingerprints_change_on_worker_edit(tmp_path):
    # Use explicit manifest, test that mutating an included file changes digest.
    from openvsppy.vv.VVDatabase import source_fingerprint
    (tmp_path/'worker.py').write_text('version 1')
    first=source_fingerprint(tmp_path,('worker.py',))
    (tmp_path/'worker.py').write_text('version 2')
    assert source_fingerprint(tmp_path,('worker.py',))!=first
    assert len(execution_fingerprint(ROOT))==64
    assert len(postprocess_fingerprint(ROOT))==64


def test_catalog_and_constraints(study):
    _,db=study
    assert db.conn.execute('SELECT count(*) FROM test_definition').fetchone()[0]==25
    assert db.conn.execute('SELECT count(*) FROM planned_case').fetchone()[0]==68
    assert db.conn.execute('PRAGMA foreign_keys').fetchone()[0]==1
    assert db.row('planned_case','case_id','GRID-01-MEDIUM')['execute_policy']=='UNRESOLVED'
    with pytest.raises(ValueError):db.resolve_case('GRID-01-MEDIUM','none','none','none')


def test_geometry_reference_hash_and_case(study, configured_case):
    root,db,cond=configured_case
    errors=check_case(db,'GRID-01-MEDIUM',root/'study.json')
    assert any('executable' in x for x in errors)
    assert not check_case(db,'GRID-01-MEDIUM',root/'study.json',binary_path=ROOT/'README.md')
    with pytest.raises(ValueError):db.add_condition('VSP3505-VV-001',{**cond,'mach':float('nan'),'condition_id':'NAN'})


def test_run_states_and_native_polar_import(study, configured_case):
    root,db,_=configured_case
    runid,run_dir=db.register_run('GRID-01-MEDIUM',root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'test':1})
    with pytest.raises(ValueError):db.transition(runid,'SUCCEEDED')
    db.transition(runid,'PRECHECKED')
    db.transition(runid,'RUNNING')
    db.transition(runid,'CAPTURED',effective={'test':1})
    raw=run_dir/'raw';raw.mkdir()
    polar=raw/'sample.polar';shutil.copy2(ROOT/'examples/models/CRM-HS/CRM-HS.polar',polar)
    p=import_polar(db,study_root=root,run_id=runid,path=polar)
    assert p['matched'] and p['qoi_count']>=40
    assert db.conn.execute('SELECT count(*) FROM quantity WHERE dataset_id=?',(p['dataset_id'],)).fetchone()[0]>=40
    db.transition(runid,'PARSED');db.transition(runid,'SUCCEEDED')
    assert db.row('run','run_id',runid)['status']=='SUCCEEDED'
    assert db.verify_artifacts(root)==[(str(polar.relative_to(root)).replace('\\','/'),'VERIFIED')]
    # If source file is later modified its audit status becomes CHANGED_OR_MISSING.
    polar.write_text('modified',encoding='utf-8')
    assert db.verify_artifacts(root)[0][1]=='CHANGED_OR_MISSING'


def test_external_cd0_provenance(study, configured_case):
    root,db,_=configured_case
    p=insert_dataset(db,study_id='VSP3505-VV-001',kind='VSPAERO',source='test',point_id='missing',method='MOCK') if False else None
    # The DB requires a real VSPAERO point. Ensure that reference-only CD0 cannot
    # be mislabelled as a VSPAERO run or treated as aerodynamic verification.
    with pytest.raises(Exception):
        insert_dataset(db,study_id='VSP3505-VV-001',kind='VSPAERO',source='test',point_id='not_a_point',method='TEST')


def test_case_redefinition_and_missing_reference(study, configured_case):
    root,db,cond=configured_case
    with pytest.raises(ValueError,match='cannot be reused'):
        db.add_condition('VSP3505-VV-001',{**cond,'alpha_deg':5.})
    with pytest.raises(ValueError,match='cannot be reused'):
        db.add_setting('VSP3505-VV-001',{'setting_id':'SET1','analysis_mode':'STEADY_SWEEP',
                            'mesh_level':'FINE','requested':{'thin_set_index':3,'thick_set_index':4}})
    with pytest.raises(ValueError,match='immutable'):
        db.add_geometry('VSP3505-VV-001',{'geometry_id':'G3','stage':'G3',
             'vsp3_source':str(ROOT/'examples/models/G103A/G103A.vsp3'),'representation':{'set_indices':[3,4]}})
