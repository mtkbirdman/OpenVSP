"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path
import json
import re
import shutil

from openvsppy.vv.VVAutomation import execute_study
from openvsppy.vv.VVDatabase import VVDatabase
from openvsppy.vv.VVDatabase import digest
from openvsppy.vv.VVPostprocess import import_polar
from openvsppy.vv.VVPostprocess import import_reference_csv
from openvsppy.vv.VVRunner import run_cases

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT
BASE = ROOT / 'vv_study/g103a_baseline.reviewed.json'


def test_empty_study_reports_all_tests_and_blocks_without_binary(tmp_path, study_factory):
    root,db=study_factory(tmp_path / "vv_study")
    result=run_cases(db.path,root/'study.json',['GRID-01-COARSE'],binary_path=None)
    assert result[0]['status']=='BLOCKED'
    assert db.conn.execute('SELECT COUNT(*) FROM run').fetchone()[0]==0
    db.close()
    report=execute_study(root/'study.json',root/'manifest.json',execute=False)
    assert report['case_count']==68
    assert report['figure_count']==0
    content=Path(report['report']).read_text(encoding='utf-8')
    assert len(re.findall(r'^## \d+\. ',content,re.M))==13
    assert len(re.findall(r'^### (VER|GEO|GRID|WAKE|XCD|VAL|STA|CTRL|DRAG|ACC)-\d+',content,re.M))==25
    assert 'グラフ: **生成対象データ不足**' in content
    assert Path(report['output_dir'],'results_snapshot.json').exists()
    assert Path(report['output_dir'],'database_snapshot.sqlite').exists()
    assert 'SQLite frozen snapshot' in content
    blocked=execute_study(root/'study.json',root/'manifest.json',execute=True)
    assert blocked['case_count']==68 and blocked['blocked_or_failed']==68
    assert blocked['successful_or_reused_runs']==0
    progress=Path(blocked['output_dir'],'progress.jsonl').read_text(encoding='utf-8')
    assert progress.count('"case_id"')>=68


def test_frozen_successful_runs_generate_real_review_figures_and_tables(tmp_path, study_factory):
    root,db=study_factory(tmp_path / "vv_study")
    geom=ROOT/'examples/models/CRM-HS/CRM-HS.vsp3'
    polar=ROOT/'examples/models/CRM-HS/CRM-HS.polar'
    db.add_geometry('VSP3505-VV-001',{'geometry_id':'G','stage':'G3','vsp3_source':str(geom)})
    condition={'condition_id':'F','mach':0.6,'re_cref':1e7,'alpha_deg':3.0,'beta_deg':0.0,
        'controls_deg':{},'trim_mode':'FIXED_CONTROL','sref':10.0,'bref':12.0,'cref':2.0,
        'moment_x':0.0,'moment_y':0.0,'moment_z':0.0,'reference_length_unit':'m',
        'axes_convention_id':'VSPAERO_NATIVE'}
    db.add_condition('VSP3505-VV-001',condition)
    saved={}
    for role in ('COARSE','MEDIUM'):
        setting='S-'+role
        db.add_setting('VSP3505-VV-001',{'setting_id':setting,'analysis_mode':'STEADY_SWEEP',
            'mesh_level':role,'requested':{'thin_set_index':4,'thick_set_index':3}})
        cid='GRID-01-'+role
        db.resolve_case(cid,'G','F',setting)
        run,run_dir=db.register_run(cid,root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'test':role})
        db.transition(run,'PRECHECKED');db.transition(run,'RUNNING');db.transition(run,'CAPTURED')
        native=run_dir/'raw'/'sample.polar'
        native.parent.mkdir()
        shutil.copy2(polar,native)
        point=import_polar(db,study_root=root,run_id=run,path=native)
        db.transition(run,'PARSED');db.transition(run,'SUCCEEDED')
        with db.transaction() as conn:
            conn.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',(cid,point['point_id'],'2026-10-31T00:00:00Z','TEST'))
        saved[cid]=run
    raw_hash={r:digest(root/db.row('run','run_id',r)['run_dir']/'raw'/'sample.polar') for r in saved.values()}
    db.close()
    (root/'manifest.json').write_text(json.dumps({'run_selection':saved}),encoding='utf-8')
    result=execute_study(root/'study.json',root/'manifest.json',execute=False)
    assert result['figure_count']>=1
    assert result['run_selection']==saved
    document=Path(result['report']).read_text(encoding='utf-8')
    assert '![FIG-GRID-01-' in document
    from openvsppy.vv.VVDatabase import VVDatabase
    products = Path(result['output_dir'])
    evidence = json.loads((products/'evidence_manifest.json').read_text(encoding='utf8'))
    assert evidence['database_source'] == 'database_snapshot.sqlite'
    for entry in evidence['files']:
        assert digest(products/entry['path']) == entry['sha256']
    frozen = VVDatabase(products/'database_snapshot.sqlite')
    try:
        assert frozen.conn.execute("SELECT COUNT(*) FROM artifact WHERE role='REVIEW_REPORT'").fetchone()[0] == 0
        # Every Artifact in the frozen DB is exported with the same ID set.
        with (products/'tables/artifact.csv').open(encoding='utf-8-sig',newline='') as f:
            rows = list(__import__('csv').DictReader(f))
        assert {r['artifact_id'] for r in rows} == {
            r['artifact_id'] for r in frozen.conn.execute('SELECT artifact_id FROM artifact')}
    finally:
        frozen.close()
    assert result['embedded_figure_count'] <= result['figure_count']
    assert len(list((products/'test_reports').glob('*.md'))) == 25
    assert 'CMiy' in document and 'GRID-01-MEDIUM' in document
    assert all((Path(result['output_dir'])/p['path']).is_file() for p in json.loads(
        (Path(result['output_dir'])/'results_snapshot.json').read_text())['charts'])
    db=VVDatabase(root/'01_database'/'vv.sqlite')
    assert db.conn.execute('SELECT COUNT(*) FROM run').fetchone()[0]==2
    for rid,original in raw_hash.items():
        assert digest(root/db.row('run','run_id',rid)['run_dir']/'raw'/'sample.polar')==original
    db.close()


def test_reference_comparison_and_drag_figures_enter_formal_report(tmp_path, study_factory):
    """Synthetic DB bookkeeping using an existing .polar fixture, not Validation."""
    root,db=study_factory(tmp_path / "vv_study")
    geom=ROOT/'examples/models/CRM-HS/CRM-HS.vsp3'
    db.add_geometry('VSP3505-VV-001',{'geometry_id':'G','stage':'G3','vsp3_source':str(geom)})
    fc={'condition_id':'F','mach':0.6,'re_cref':1e7,'alpha_deg':3.0,'beta_deg':0.0,
        'controls_deg':{},'trim_mode':'FIXED_CONTROL','sref':10.0,'bref':12.0,'cref':2.0,
        'moment_x':0.0,'moment_y':0.0,'moment_z':0.0,'reference_length_unit':'m',
        'axes_convention_id':'VSPAERO_NATIVE'}
    db.add_condition('VSP3505-VV-001',fc)
    db.add_setting('VSP3505-VV-001',{'setting_id':'S','analysis_mode':'STEADY_SWEEP',
        'mesh_level':'MEDIUM','requested':{'thin_set_index':4,'thick_set_index':3}})
    cid='VER-01-COARSE'
    db.resolve_case(cid,'G','F','S')
    run,folder=db.register_run(cid,root,ROOT/'README.md',ROOT/'src/openvsppy/vv/VVRunner.py',{'test':'reference'})
    db.transition(run,'PRECHECKED');db.transition(run,'RUNNING');db.transition(run,'CAPTURED')
    native=folder/'raw'/'example.polar';native.parent.mkdir()
    shutil.copy2(ROOT/'examples/models/CRM-HS/CRM-HS.polar',native)
    point=import_polar(db,study_root=root,run_id=run,path=native)
    db.transition(run,'PARSED');db.transition(run,'SUCCEEDED')
    with db.transaction() as conn:
        conn.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',(cid,point['point_id'],'2026-10-31T00:00:00Z','TEST'))
    db.close()
    for filename,qoi,value in (('theory.csv','CLi',0.39),('cd0.csv','CD0_build',0.005)):
        (root/filename).write_text(f'qoi,value,unit\n{qoi},{value},1\n',encoding='utf8')
        kind='THEORY' if filename=='theory.csv' else 'EXTERNAL_CD0'
        meta={'kind':kind,'source_id':'TH-1' if kind=='THEORY' else 'CD0-1',
            'axes_convention_id':'VSPAERO_NATIVE','sref':10,'bref':12,'cref':2,'mach':.6,
            're_cref':1e7,'alpha_deg':3,'beta_deg':0,'uncertainty_status':'NOT_ASSESSED',
            'usage_class':'REFERENCE','geometry_equivalence_reviewed':True,'reference_reviewed':True}
        (root/(filename+'.json')).write_text(json.dumps(meta),encoding='utf8')
    manifest={'run_selection':{cid:run},'references':[
        {'file':'theory.csv','metadata':'theory.csv.json'},
        {'file':'cd0.csv','metadata':'cd0.csv.json'}],
        'comparisons':[{'test_id':'VER-01','case_id':cid,
                        'reference_source_id':'TH-1','conditions_reviewed':True}],
        'drag_corrections':[{'case_id':cid,'cd0_source_id':'CD0-1',
                             'cdo_baseline':0.005,'baseline_metadata':{'condition_id':'F'}}]}
    (root/'manifest.json').write_text(json.dumps(manifest),encoding='utf8')
    result=execute_study(root/'study.json',root/'manifest.json',execute=False)
    snapshot=json.loads(Path(result['output_dir'],'results_snapshot.json').read_text(encoding='utf8'))
    assert len(snapshot['comparisons'])==1
    assert {'DRAG-01','DRAG-02'} <= {p['test_id'] for p in snapshot['drag_polars']}
    assert len([p for p in snapshot['charts'] if p['kind']=='DRAG_POLAR'])==2
    assert any(p['kind']=='REFERENCE_COMPARISON' for p in snapshot['charts'])
    text=Path(result['report']).read_text(encoding='utf8')
    assert '![FIG-VER-01-' in text
    assert 'ParasiteDrag補正' in (Path(result['output_dir'])/'test_reports/DRAG-01.md').read_text(encoding='utf8')
    # Repeated review/export must not duplicate immutable source datasets,
    # comparisons or drag correction quantities.
    repeated=execute_study(root/'study.json',root/'manifest.json',execute=False)
    reread=json.loads(Path(repeated['output_dir'],'results_snapshot.json').read_text(encoding='utf8'))
    assert all(ref['status']=='REUSED_REFERENCE' for ref in reread['reference_imports'])
    assert reread['derivations'][0]['records'][0]['status']=='REUSED_COMPARISON'
    assert reread['derivations'][1]['status']=='REUSED_CORRECTION'
    db=VVDatabase(root/'01_database'/'vv.sqlite')
    assert db.conn.execute('SELECT COUNT(*) FROM comparison').fetchone()[0]==1
    assert db.conn.execute("SELECT COUNT(*) FROM dataset WHERE method_id='CD0_CORRECTION'").fetchone()[0]==1
    db.close()


def test_approved_manifest_path_and_report_without_solving(tmp_path, study_factory):
    root,db=study_factory(tmp_path / "vv_study")
    db.close()
    # Manifest paths resolve relative to the supplied manifest, never to cwd.
    model=str(REPO/'examples/models/G103A/G103A.vsp3')
    case=json.loads(BASE.read_text())
    case['geometry']['vsp3_source']=model
    (root/'baseline.json').write_text(json.dumps(case))
    (root/'manifest.json').write_text(json.dumps({'baselines':{
        'GEO-04':'baseline.json','GRID-03':'baseline.json','STA-03':'baseline.json'}}))
    output=execute_study(root/'study.json',root/'manifest.json',execute=False)
    assert output['case_count']==69
    db=VVDatabase(root/'01_database/vv.sqlite')
    assert db.row('planned_case','case_id','GRID-03-NREF_1') is not None
    assert db.row('planned_case','case_id','STA-03-P')['geometry_id'] is not None
    assert db.row('planned_case','case_id','GEO-04-THICK')['geometry_id'] is not None
    assert db.conn.execute('SELECT COUNT(*) FROM run').fetchone()[0]==0
    db.close()
    assert Path(output['report']).exists()


def test_automatic_drag_a_b_uses_selected_native_baseline_without_duplication(tmp_path,monkeypatch, study_factory):
    """Integrates existing polar, external-CD0 record and A/B (no solver run)."""
    from openvsppy.vv.VVPostprocess import import_polar, import_reference_csv
    root,db=study_factory(tmp_path / "vv_study")
    study_id='VSP3505-VV-001'
    db.add_geometry(study_id,{'geometry_id':'G','stage':'G3',
        'vsp3_source':str(REPO/'examples/models/CRM-HS/CRM-HS.vsp3')})
    cond={'condition_id':'F','mach':.6,'re_cref':1e7,'alpha_deg':3.,'beta_deg':0.,
        'controls_deg':{},'trim_mode':'FIXED_CONTROL','sref':10.,'bref':12.,'cref':2.,
        'moment_x':0.,'moment_y':0.,'moment_z':0.,
        'reference_length_unit':'m','axes_convention_id':'VSPAERO_NATIVE'}
    db.add_condition(study_id,cond)
    db.add_setting(study_id,{'setting_id':'S','analysis_mode':'STEADY_SWEEP',
        'mesh_level':'MEDIUM','requested':{'thin_set_index':4,'thick_set_index':3}})
    cid='GEO-03-G3';db.resolve_case(cid,'G','F','S')
    run,folder=db.register_run(cid,root,REPO/'README.md',REPO/'src/openvsppy/vv/VVRunner.py',{'test':'drag'})
    for state in ('PRECHECKED','RUNNING','CAPTURED'):
        db.transition(run,state)
    native=folder/'raw'/'example.polar';native.parent.mkdir()
    shutil.copy2(REPO/'examples/models/CRM-HS/CRM-HS.polar',native)
    imported=import_polar(db,study_root=root,run_id=run,path=native)
    for state in ('PARSED','SUCCEEDED'):
        db.transition(run,state)
    with db.transaction() as connection:
        connection.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',
                           (cid,imported['point_id'],'2026-10-31T00:00:00Z','TEST'))
    reference=root/'cd0.csv';reference.write_text('qoi,value,unit\nCD0_build,0.006,1\n')
    metadata={'kind':'EXTERNAL_CD0','source_id':'mock-cd0',
        'axes_convention_id':'VSPAERO_NATIVE','sref':10.,'bref':12.,'cref':2.,
        'mach':.6,'re_cref':1e7,'alpha_deg':3.,'beta_deg':0.,
        'uncertainty_status':'NOT_ASSESSED','usage_class':'REFERENCE'}
    cd0=import_reference_csv(db,study_id=study_id,csv_path=reference,metadata=metadata)
    db.close()
    # Stub solely replaces the unavailable 3.50.5 worker, not the correction formulas.
    monkeypatch.setattr('openvsppy.vv.VVRunner.run_parasite_drag',
                        lambda *args,**kwargs:{'dataset_id':cd0['dataset_id'],'status':'REUSED_PARASITE'})
    manifest={'auto_drag':{
        'source_id':'mock-cd0','model':str(REPO/'examples/models/CRM-HS/CRM-HS.vsp3'),
        'geometry_names':['WingGeom'],'geom_set':5,'velocity_mps':30,
        'baseline_case_id':cid,'source_case_ids':[cid],'methods':['A','B']},
        'run_selection':{cid:run}}
    (root/'manifest.json').write_text(json.dumps(manifest))
    first=execute_study(root/'study.json',root/'manifest.json',execute=True)
    snap=json.loads(Path(first['output_dir'],'results_snapshot.json').read_text())
    assert {x['test_id'] for x in snap['drag_polars']}=={'DRAG-01','DRAG-02'}
    assert {x['qoi'] for x in snap['drag_polars']}=={'CD_A','CD_B'}
    assert first['figure_count']>=2
    db=VVDatabase(root/'01_database/vv.sqlite')
    before=db.conn.execute("SELECT COUNT(*) FROM dataset WHERE method_id='CD0_CORRECTION'").fetchone()[0]
    assert before==2
    db.close()
    second=execute_study(root/'study.json',root/'manifest.json',execute=True)
    db=VVDatabase(root/'01_database/vv.sqlite')
    after=db.conn.execute("SELECT COUNT(*) FROM dataset WHERE method_id='CD0_CORRECTION'").fetchone()[0]
    assert after==before
    assert db.conn.execute('SELECT COUNT(*) FROM run').fetchone()[0]==1
    db.close()
