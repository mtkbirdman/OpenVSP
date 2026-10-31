"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path
import json
import math
import shutil

import pytest

from openvsppy.vv.VVDatabase import VVDatabase
from openvsppy.vv.VVDatabase import digest
from openvsppy.vv.VVPostprocess import central_difference
from openvsppy.vv.VVPostprocess import compare_quantities
from openvsppy.vv.VVPostprocess import create_drag_dataset
from openvsppy.vv.VVPostprocess import drag_a
from openvsppy.vv.VVPostprocess import drag_b
from openvsppy.vv.VVPostprocess import evaluate_control_sweep
from openvsppy.vv.VVPostprocess import import_polar
from openvsppy.vv.VVPostprocess import import_reference_csv
from openvsppy.vv.VVPostprocess import import_stab
from openvsppy.vv.VVPostprocess import insert_dataset
from openvsppy.vv.VVPostprocess import insert_quantity
from openvsppy.vv.VVPostprocess import reprocess_run
from openvsppy.vv.VVPostprocess import richardson_gci
from openvsppy.vv.VVPostprocess import screen_convergence
from openvsppy.vv.VVPostprocess import store_assessment
from openvsppy.vv.VVPostprocess import validation_uncertainty

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'tests/fixtures/SampleGlider_native_3.50.5.stab'
MAPPING = {'ConGrp_1':'AILERON_GROUP','ConGrp_2':'ELEVATOR_GROUP','ConGrp_3':'RUDDER_GROUP'}


def test_nonfinite_optional_polar_is_recorded_in_sqlite(study):
    root, db = study
    sid = 'VSP3505-VV-001'
    db.add_geometry(sid, {'geometry_id': 'g', 'stage': 'G3',
                          'vsp3_source': str(ROOT / 'examples/models/CRM-HS/CRM-HS.vsp3')})
    db.add_condition(sid, {'condition_id': 'f', 'mach': .6, 're_cref': 1e7,
                           'alpha_deg': 3., 'beta_deg': 0., 'controls_deg': {},
                           'trim_mode': 'FIXED_CONTROL', 'sref': 10., 'bref': 12., 'cref': 2.,
                           'moment_x': 0., 'moment_y': 0., 'moment_z': 0.,
                           'reference_length_unit': 'm', 'axes_convention_id': 'VSPAERO_NATIVE'})
    db.add_setting(sid, {'setting_id': 's', 'analysis_mode': 'STEADY_SWEEP',
                         'mesh_level': 'MEDIUM', 'requested': {'thin_set_index': 4,
                                                              'thick_set_index': 3}})
    db.resolve_case('GRID-01-MEDIUM', 'g', 'f', 's')
    run_id, run_dir = db.register_run('GRID-01-MEDIUM', root, ROOT / 'README.md',
                                      ROOT / 'src/openvsppy/vv/VVRunner.py', {'test': 'nonfinite'})
    (run_dir / 'raw').mkdir()
    path = run_dir / 'raw' / 'native.polar'
    lines = (ROOT / 'examples/models/CRM-HS/CRM-HS.polar').read_text().splitlines()
    header = next(x.split() for x in lines if x.split()[:3] == ['Beta', 'Mach', 'AoA'])
    i = next(i for i, x in enumerate(lines) if len(x.split()) == len(header)
             and x.split()[0] == '0.000000000000')
    row = lines[i].split()
    row[header.index('L/D')] = 'nan'
    lines[i] = ' '.join(row)
    path.write_text('\n'.join(lines))
    result = import_polar(db, study_root=root, run_id=run_id, path=path)
    assert result['matched']
    assert result['nonfinite_columns'] == ['L/D']
    status = {x['qoi']: (x['value_status'], x['native_value']) for x in db.conn.execute(
        'SELECT qoi,value_status,native_value FROM quantity WHERE dataset_id=?',
        (result['dataset_id'],))}
    assert status['L/D'] == ('NOT_AVAILABLE', None)
    assert status['CLi'][0] == 'AVAILABLE' and status['CLi'][1] is not None


def test_stability_db_import_and_idempotent_reprocess(tmp_path):
    study=tmp_path/'vv_study';study.mkdir()
    db=VVDatabase(study/'vv.sqlite')
    db.initialize(ROOT/'vv_study/study.json',ROOT/'vv_study/01_database/schema.sql')
    sid='VSP3505-VV-001'
    db.add_geometry(sid,{'geometry_id':'SAMPLE','stage':'G3',
                         'vsp3_source':str(ROOT/'examples/models/G103A/G103A.vsp3')})
    cond={'condition_id':'SAMPLE-FC','mach':0.,'re_cref':2e6,'alpha_deg':0.,'beta_deg':0.,
        'controls_deg':{},'trim_mode':'FIXED_CONTROL','sref':18.,'bref':27.,'cref':.625,
        'moment_x':1.4,'moment_y':0.,'moment_z':0.,'reference_length_unit':'m',
        'axes_convention_id':'NATIVE_UNVERIFIED'}
    db.add_condition(sid,cond)
    db.add_setting(sid,{'setting_id':'SAMPLE-SET','analysis_mode':'STEADY',
       'mesh_level':'BASE','requested':{'thin_set_index':4,'thick_set_index':3,
       'control_group_names':{'aileron':'AILERON_GROUP','elevator':'ELEVATOR_GROUP','rudder':'RUDDER_GROUP'}}})
    db.resolve_case('STA-01-STEADY','SAMPLE','SAMPLE-FC','SAMPLE-SET')
    run,folder=db.register_run('STA-01-STEADY',study,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'test':'synthetic'})
    for phase in ('PRECHECKED','RUNNING','CAPTURED'):db.transition(run,phase)
    (folder/'raw').mkdir()
    file=folder/'raw'/'model_effective.stab';shutil.copyfile(FIXTURE,file)
    with pytest.raises(ValueError,match='not verified'):
        import_stab(db,study_root=study,run_id=run,path=file)
    result=import_stab(db,study_root=study,run_id=run,path=file,control_mapping=MAPPING)
    assert result['quantity_count']==132
    assert db.conn.execute("SELECT COUNT(*) FROM quantity WHERE dataset_id=? AND value_status='NATIVE_ONLY'",(result['dataset_id'],)).fetchone()[0]==132
    assert db.conn.execute("SELECT COUNT(*) FROM quantity WHERE dataset_id=? AND canonical_value IS NOT NULL",(result['dataset_id'],)).fetchone()[0]==0
    db.transition(run,'PARSED')
    with db.transaction():
        db.transition(run,'SUCCEEDED')
        db.conn.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',('STA-01-STEADY',result['point_id'],'2026-10-31','ONE_POINT_MATCH'))
    before=digest(file)
    new=reprocess_run(db,study_root=study,run_id=run)
    assert new['status']=='REPROCESSED'
    assert reprocess_run(db,study_root=study,run_id=run)['status']=='REUSED_POSTPROCESS'
    assert before==digest(file)
    assert db.conn.execute("SELECT COUNT(*) FROM quantity WHERE dataset_id=? AND value_status='NATIVE_ONLY'",(new['dataset_id'],)).fetchone()[0]==132
    assert db.conn.execute('PRAGMA foreign_key_check').fetchall()==[]
    db.close()


def test_gci_only_for_monotonic_uniform_refinement():
    valid=richardson_gci([1.25,1.0625,1.015625],[0.4,0.2,0.1])
    assert valid['status']=='CANDIDATE_ONLY'
    assert valid['order']==pytest.approx(2.)
    assert valid['extrapolated']==pytest.approx(1.)
    assert valid['asymptotic_verified'] is False
    assert richardson_gci([1.,1.2,1.1],[.4,.2,.1])['status']=='NOT_ASSESSED'
    assert richardson_gci([1.25,1.0625,1.015625],[.4,.3,.1])['status']=='NOT_ASSESSED'


def test_validation_uncertainty_requires_review():
    assert validation_uncertainty(numerical=1,input_=2,experimental=2)['status']=='NOT_ASSESSED'
    result=validation_uncertainty(numerical=1,input_=2,experimental=2,independent=True,standard=True)
    assert result['standard_uncertainty']==pytest.approx(3)
    assert validation_uncertainty(numerical=None,input_=2,experimental=2,independent=True,standard=True)['status']=='NOT_ASSESSED'


def test_control_effectiveness_per_rad_without_acceptance():
    q=evaluate_control_sweep({-5:-.5,-2:-.2,0:0.,2:.2,5:.5,10:1.},'CMiy')
    assert q['symmetric_per_rad']['2']==pytest.approx(.4/(4*math.pi/180))
    assert q['slope_difference_per_rad']==pytest.approx(0)
    assert q['status']=='INSUFFICIENT_EVIDENCE'


def test_assessment_remains_unapproved(study):
    _,db=study
    row=store_assessment(db,study_id='VSP3505-VV-001',test_id='ACC-01',qoi='ANY',evidence=[],reason='test')
    record=db.conn.execute('SELECT * FROM assessment WHERE assessment_id=?',(row,)).fetchone()
    assert record['engineering_status']=='NOT_EVALUATED'
    assert record['reviewer'] is None and record['approved_utc'] is None
    with pytest.raises(ValueError,match='cannot approve'):
        store_assessment(db,study_id='VSP3505-VV-001',test_id='ACC-01',qoi='ANY',evidence=[],reason='test',engineering='PASS')


def test_reprocessing_existing_raw_never_changes_native(study):
    root,db=study
    geom=ROOT/'examples/models/CRM-HS/CRM-HS.vsp3'
    db.add_geometry('VSP3505-VV-001',{'geometry_id':'G','stage':'G3','vsp3_source':str(geom)})
    fc={'condition_id':'F','mach':.6,'re_cref':1e7,'alpha_deg':3.,'beta_deg':0.,
        'controls_deg':{},'trim_mode':'FIXED_CONTROL','sref':10.,'bref':12.,'cref':2.,
        'moment_x':0.,'moment_y':0.,'moment_z':0.,'reference_length_unit':'m',
        'axes_convention_id':'VSPAERO_NATIVE'}
    db.add_condition('VSP3505-VV-001',fc)
    db.add_setting('VSP3505-VV-001',{'setting_id':'S','analysis_mode':'STEADY_SWEEP',
        'mesh_level':'MEDIUM','requested':{'thin_set_index':4,'thick_set_index':3}})
    db.resolve_case('GRID-01-MEDIUM','G','F','S')
    rid,rd=db.register_run('GRID-01-MEDIUM',root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'test':'reprocess'})
    db.transition(rid,'PRECHECKED');db.transition(rid,'RUNNING');db.transition(rid,'CAPTURED')
    (rd/'raw').mkdir()
    shutil.copy2(ROOT/'examples/models/CRM-HS/CRM-HS.polar',rd/'raw'/'orig.polar')
    imp=import_polar(db,study_root=root,run_id=rid,path=rd/'raw'/'orig.polar')
    db.transition(rid,'PARSED');db.transition(rid,'SUCCEEDED')
    original=digest(rd/'raw'/'orig.polar')
    result=reprocess_run(db,study_root=root,run_id=rid)
    assert result['status']=='REPROCESSED'
    assert result['dataset_id']!=imp['dataset_id']
    assert reprocess_run(db,study_root=root,run_id=rid)['status']=='REUSED_POSTPROCESS'
    assert digest(rd/'raw'/'orig.polar')==original
    assert db.conn.execute("SELECT count(*) FROM dataset WHERE point_id=?",(imp['point_id'],)).fetchone()[0]==2


def test_drag_and_central_difference():
    assert drag_a(.015,.02)==pytest.approx(.035)
    assert drag_b(.015,.007,.005,.020)==pytest.approx(.037)
    assert central_difference(.2,-.2,2.)==pytest.approx(.4/(2*math.radians(2)))
    assert screen_convergence([.5,.51,.512],.005)['provisional_screen']=='WITHIN_PROVISIONAL_LIMIT'
    assert screen_convergence([.5,.51,.512],.005)['engineering_verdict']=='INCONCLUSIVE'


def test_invalid_control_group_mapping(study, configured_case):
    root,db,_=configured_case
    runid,rd=db.register_run('GRID-01-MEDIUM',root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'test':1})
    raw=rd/'raw';raw.mkdir()
    shutil.copy2(ROOT/'examples/models/CRM-HS/CRM-HS.polar',raw/'example.polar')
    import_polar(db,study_root=root,run_id=runid,path=raw/'example.polar')
    shutil.copy2(ROOT/'tests/fixtures/SampleGlider_native_3.50.5.stab',raw/'example.stab')
    with pytest.raises(ValueError,match='not verified'):
        import_stab(db,study_root=root,run_id=runid,path=raw/'example.stab')
    # A reviewed control mapping does not excuse an unrelated native Mach,
    # AoA or reference geometry: this fixture comes from a different aircraft.
    with pytest.raises(ValueError,match='Native stability condition mismatch'):
        import_stab(db,study_root=root,run_id=runid,path=raw/'example.stab',
                    control_mapping={'ConGrp_1':'AILERON_GROUP','ConGrp_2':'ELEVATOR_GROUP','ConGrp_3':'RUDDER_GROUP'})


def test_reference_import_and_comparison(study,tmp_path, configured_case):
    root,db,_=configured_case
    src=tmp_path/'ref.csv'
    src.write_text('qoi,value,unit\nCLi,0.4,1\n',encoding='utf-8')
    meta={'kind':'WIND_TUNNEL','source_id':'wind_tunnel_dataset_A','axes_convention_id':'VSPAERO_NATIVE',
          'sref':10.,'bref':12.,'cref':2.,'mach':.6,'re_cref':1e7,'alpha_deg':3.,'beta_deg':0.,
          'uncertainty_status':'UNKNOWN','usage_class':'VALIDATION'}
    result=import_reference_csv(db,study_id='VSP3505-VV-001',csv_path=src,metadata=meta)
    qref=db.conn.execute('SELECT quantity_id FROM quantity WHERE dataset_id=?',(result['dataset_id'],)).fetchone()[0]
    ds=insert_dataset(db,study_id='VSP3505-VV-001',kind='DERIVED',source='unit test',method='TEST')
    qpred=insert_quantity(db,ds,'CLi',.404)
    diff=compare_quantities(db,test_id='VAL-01',predicted_id=qpred,reference_id=qref,condition_verified=True)
    assert diff['verdict']=='INCONCLUSIVE' and diff['delta']==pytest.approx(.004)


def test_cd0_derived_dataset_keeps_native_and_dependencies(study, configured_case):
    root,db,_=configured_case
    runid,rd=db.register_run('GRID-01-MEDIUM',root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'test':1})
    (rd/'raw').mkdir()
    shutil.copy2(ROOT/'examples/models/CRM-HS/CRM-HS.polar',rd/'raw'/'test.polar')
    polar=import_polar(db,study_root=root,run_id=runid,path=rd/'raw'/'test.polar')
    cd0=insert_dataset(db,study_id='VSP3505-VV-001',kind='EXTERNAL_CD0',source='test surrogate',
                      metadata={'sref':10.,'mach':0.6,'re_cref':1e7},method='TEST')
    insert_quantity(db,cd0,'CD0_build',.02)
    new=create_drag_dataset(db,study_id='VSP3505-VV-001',polar_dataset=polar['dataset_id'],cd0_dataset=cd0)
    value=db.conn.execute("SELECT canonical_value FROM quantity WHERE dataset_id=? AND qoi='CD_A'",(new['dataset_id'],)).fetchone()[0]
    cdiw=db.conn.execute("SELECT canonical_value FROM quantity WHERE dataset_id=? AND qoi='CDiw'",(polar['dataset_id'],)).fetchone()[0]
    assert value==pytest.approx(cdiw+.02)
    assert db.conn.execute('SELECT COUNT(*) FROM quantity_dependency WHERE derived_quantity_id=?',(new['CD_A_quantity_id'],)).fetchone()[0]==2
    with pytest.raises(ValueError,match='baseline'):
        create_drag_dataset(db,study_id='VSP3505-VV-001',polar_dataset=polar['dataset_id'],cd0_dataset=cd0,cdo_baseline=.008)


def test_evaluate_selects_one_latest_run_not_a_mix(study, capsys, configured_case):
    import datetime
    from openvsppy.vv.__main__ import main
    root,db,_=configured_case
    selected=[]
    for i,value in enumerate((0.111,0.222)):
        rid, rd = db.register_run('GRID-01-MEDIUM',root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'i':i})
        db.transition(rid,'PRECHECKED')
        db.transition(rid,'RUNNING')
        db.transition(rid,'CAPTURED')
        (rd/'raw').mkdir()
        fp=rd/'raw'/'sample.polar'
        shutil.copy2(ROOT/'examples/models/CRM-HS/CRM-HS.polar', fp)
        imp=import_polar(db,study_root=root,run_id=rid,path=fp)
        with db.transaction() as c:
            c.execute("UPDATE quantity SET canonical_value=? WHERE dataset_id=? AND qoi='CLi'",(value,imp['dataset_id']))
        db.transition(rid,'PARSED')
        db.transition(rid,'SUCCEEDED')
        with db.transaction() as c:
            c.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',('GRID-01-MEDIUM',imp['point_id'],datetime.datetime.now(datetime.timezone.utc).isoformat(),'TEST'))
        selected.append((rid,imp['dataset_id']))
    rc=main(['--study',str(root/'study.json'),'--db',str(db.path),'evaluate','--test','GRID-01'])
    assert rc==0
    data=json.loads(capsys.readouterr().out)
    assert data['sources']['GRID-01-MEDIUM']['run_id']==selected[1][0]
    assert data['sources']['GRID-01-MEDIUM']['dataset_id']==selected[1][1]


def test_evaluate_explicit_selection_file(study,tmp_path,capsys, configured_case):
    from openvsppy.vv.__main__ import main
    root,db,_=configured_case
    choices=[]
    for i in range(2):
        rid,runpath=db.register_run('GRID-01-MEDIUM',root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'iteration':i})
        db.transition(rid,'PRECHECKED');db.transition(rid,'RUNNING');db.transition(rid,'CAPTURED')
        raw=runpath/'raw';raw.mkdir()
        file=raw/'sample.polar';shutil.copy2(ROOT/'examples/models/CRM-HS/CRM-HS.polar',file)
        ds=import_polar(db,study_root=root,run_id=rid,path=file)
        db.transition(rid,'PARSED');db.transition(rid,'SUCCEEDED')
        with db.transaction() as conn:
            conn.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',(
                'GRID-01-MEDIUM',ds['point_id'],'2026-10-30T00:00:00+00:00','TEST'))
        choices.append(rid)
    selected=tmp_path/'selected.json'
    selected.write_text(json.dumps({'GRID-01-COARSE':'NO_SUCCESS','GRID-01-MEDIUM':choices[0],
                                   'GRID-01-FINE':'NO_SUCCESS'}))
    rc=main(['--study',str(root/'study.json'),'--db',str(db.path),'evaluate',
             '--test','GRID-01','--selection-file',str(selected)])
    assert rc==0
    assert json.loads(capsys.readouterr().out)['sources']['GRID-01-MEDIUM']['run_id']==choices[0]
    selected.write_text(json.dumps({'GRID-01-MEDIUM':choices[0]}))
    with pytest.raises(ValueError,match='Missing explicit'):
        main(['--study',str(root/'study.json'),'--db',str(db.path),'evaluate',
              '--test','GRID-01','--selection-file',str(selected)])
