"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path
import json
import sys

import pytest

from openvsppy.vv.VVAutomation import prepare_cases

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT
BASE = ROOT / 'vv_study/g103a_baseline.reviewed.json'


def baseline_in_tmp(tmp_path):
    case=json.loads(BASE.read_text())
    case['geometry']['vsp3_source']=str(REPO/'examples/models/G103A/G103A.vsp3')
    file=tmp_path/'base.json'
    file.write_text(json.dumps(case))
    return file

def test_geo04_only_wing_changes(tmp_path):
    baseline=baseline_in_tmp(tmp_path)
    prepared=prepare_cases('GEO-04',baseline,tmp_path)
    first,second=[json.loads(Path(row['case_json']).read_text()) for row in prepared]
    a=first['numerical_setting']['requested']['geometry_membership']
    b=second['numerical_setting']['requested']['geometry_membership']
    assert a=={'thin':['WingGeom','HTailGeom','VTailGeom'],'thick':['FuselageGeom']}
    assert b=={'thin':['HTailGeom','VTailGeom'],'thick':['FuselageGeom','WingGeom']}
    assert first['flight_condition']==second['flight_condition']
    assert first['geometry']['vsp3_source']==second['geometry']['vsp3_source']


def test_grid03_one_factor_and_preserve_legacy_plan(tmp_path, study_factory):
    study,db=study_factory(tmp_path / "vv_study")
    baseline=baseline_in_tmp(tmp_path)
    prepared=prepare_cases('GRID-03',baseline,study)
    vals={}
    for record in prepared:
        case=json.loads(Path(record['case_json']).read_text())
        req=case['numerical_setting']['requested']
        vals[case['case_id']]=(req['geom_parms'][0]['value'],req['vspaero']['NRef'])
    assert vals=={'GRID-03-BASE':(0.5,0),'GRID-03-LOCAL_A':(0.25,0),
                  'GRID-03-LOCAL_B':(1.0,0),'GRID-03-NREF_1':(0.5,1)}
    assert db.conn.execute('SELECT COUNT(*) FROM planned_case').fetchone()[0]==68
    assert db.add_extra_planned_case('VSP3505-VV-001','GRID-03','NREF_1')=='GRID-03-NREF_1'
    assert db.conn.execute('SELECT COUNT(*) FROM planned_case').fetchone()[0]==69
    assert list((study/'01_database').glob('vv.sqlite.pre_grid03_*.bak'))
    assert db.add_extra_planned_case('VSP3505-VV-001','GRID-03','NREF_1')=='GRID-03-NREF_1'
    assert db.conn.execute('SELECT COUNT(*) FROM planned_case').fetchone()[0]==69
    db.close()


def test_sta03_three_native_modes(tmp_path):
    src=baseline_in_tmp(tmp_path)
    prepared=prepare_cases('STA-03',src,tmp_path)
    assert len(prepared)==3
    assert [json.loads(Path(r['case_json']).read_text())['numerical_setting']['analysis_mode']
            for r in prepared]==['P','Q','R']
    assert all(json.loads(Path(r['case_json']).read_text())['numerical_setting']
               ['requested']['native_pqr_defaults']['expected_time_steps']==128 for r in prepared)


def test_ver_requires_benchmark_not_substituting_g103a(tmp_path):
    src=baseline_in_tmp(tmp_path)
    with pytest.raises(ValueError,match='benchmark'):
        prepare_cases('VER-01',src,tmp_path,dry_run=True)
    case=json.loads(src.read_text())
    case['benchmark']={'test_id':'VER-01', 'vsp3_source':str(REPO/'examples/models/G103A/G103A.vsp3'),
        'geometry_membership':{'thin':['WingGeom'],'thick':[]},'thin_set_index':4,
        'thick_set_index':-1,'wing_name':'WingGeom'}
    src.write_text(json.dumps(case))
    # A synthetic input demonstrates wiring, not a valid Hershey benchmark.
    out=prepare_cases('VER-01',src,tmp_path,dry_run=True)
    assert len(out)==3 and {x['mesh_level'] for x in out}=={'COARSE','MEDIUM','FINE'}


def test_ver02_geometry_tessellation_avoids_invalid_nref_mesh(tmp_path):
    """NRef=1/2 made VSPAERO identify 2/1 loops, not a valid surface."""
    case = json.loads((ROOT / 'vv_study/ver-02_benchmark.reviewed.json').read_text())
    assert case['flight_condition']['mach'] == 0.001
    # Use an existing fixture as the programmatic model's location; the case
    # generator does not open it, only verifies it exists and fingerprints it.
    case['benchmark']['vsp3_source'] = str(ROOT / 'examples/models/G103A/G103A.vsp3')
    path = tmp_path / 'ver02.json'
    path.write_text(json.dumps(case), encoding='utf-8')
    records = prepare_cases('VER-02', path, tmp_path)
    assert len(records) == 3
    expected = [(20, 21), (30, 31), (40, 41)]
    for record, (u, w) in zip(records, expected):
        item = json.loads(Path(record['case_json']).read_text(encoding='utf-8'))
        requested = item['numerical_setting']['requested']
        assert requested['vspaero']['NRef'] == 0
        assert [(x['name'], x['value']) for x in requested['geom_parms']] == [
            ('Tess_U', u), ('Tess_W', w)]


def test_hershey_benchmark_uses_halfspan_and_verified_symmetry(tmp_path, monkeypatch):
    from openvsppy.vv.worker import api_benchmark

    class Pos:
        def __init__(self, xyz): self.xyz = xyz
        def x(self): return self.xyz[0]
        def y(self): return self.xyz[1]
        def z(self): return self.xyz[2]

    class FakeOpenVSP:
        AR_WSECT_DRIVER = 1
        ROOTC_WSECT_DRIVER = 2
        TIPC_WSECT_DRIVER = 3
        SET_ALL = 0
        def __init__(self): self.parms = {}
        def GetVSPVersion(self): return 'OpenVSP 3.50.5'
        def ClearVSPModel(self): pass
        def AddGeom(self, kind): assert kind == 'WING'; return 'wing'
        def SetGeomName(self, *args): pass
        def SetDriverGroup(self, *args): pass
        def FindParm(self, geom, name, group): return f'{name}:{group}'
        def SetParmVal(self, parm, value): self.parms[parm] = float(value)
        def SetSetFlag(self, *args): pass
        def Update(self): pass
        def GetGeomBBoxMin(self, geom): return Pos((0, -5, 0))
        def GetGeomBBoxMax(self, geom): return Pos((1, 0, 0.12))
        def GetParmVal(self, parm):
            return {'TotalSpan:WingGeom': 10., 'Sym_Planar_Flag:Sym': 2.}[parm]
        def WriteVSPFile(self, file, *args): Path(file).write_text('fake benchmark')

    monkeypatch.setitem(sys.modules, 'openvsp', FakeOpenVSP())
    path = tmp_path / 'hershey.vsp3'
    manifest = tmp_path / 'hershey.generation.json'
    job = tmp_path / 'job.json'
    job.write_text(json.dumps({'kind': 'HERSHEY_BAR', 'output': str(path),
                               'result_manifest': str(manifest)}))
    api_benchmark(job)
    saved = json.loads(manifest.read_text())
    assert saved['bounding_box_xyz'][1] == 5
    assert saved['full_span_m'] == 10
    assert saved['symmetry_flag'] == 2


def test_pqr_and_geometry_recipe_require_review(study,tmp_path):
    root,_=study
    case=json.loads((ROOT/'vv_study/g103a_baseline.example.json').read_text(encoding='utf-8'))
    case['geometry']['vsp3_source']=str(ROOT/'examples/models/G103A/G103A.vsp3')
    path=tmp_path/'base.json';path.write_text(json.dumps(case))
    with pytest.raises(ValueError,match='reviewed'):
        prepare_cases('GEO-03',path,root,dry_run=True)
    case['numerical_setting']['requested']['geometry_membership_reviewed']=True
    case['numerical_setting']['requested']['unsteady_settings_reviewed']=True
    case['numerical_setting']['requested']['vspaero'].update({'Vinf':30,'Rho':1.225})
    path.write_text(json.dumps(case))
    geo=prepare_cases('GEO-03',path,root,dry_run=True)
    assert len(geo)==1 and geo[0]['case_id']=='GEO-03-G3'
    pqr=prepare_cases('STA-03',path,root,dry_run=False)
    assert len(pqr)==3
    for item in pqr:
        saved=json.loads(Path(item['case_json']).read_text())
        assert saved['numerical_setting']['analysis_mode'] in ('P','Q','R')


def test_generated_vspaero_settings_fixture():
    from openvsppy.vv.VVPostprocess import read_vspaero_settings
    path=ROOT/'examples/models/CRM-HS/CRM-HS.vspaero'
    settings=read_vspaero_settings(path)
    assert settings['Sref']==pytest.approx(415.50831178)
    assert settings['WakeIters']==0
    assert settings['NumWakeNodes']==32


def test_prepare_wake_without_solver(study,tmp_path, configured_case):
    from openvsppy.vv.VVAutomation import prepare_cases
    root,db,cond=configured_case
    case=json.loads((ROOT/'vv_study/case.example.json').read_text())
    case['geometry']['vsp3_source']=str(ROOT/'examples/models/CRM-HS/CRM-HS.vsp3')
    case['flight_condition']=cond
    case['numerical_setting']['requested']={'thin_set_index':3,'thick_set_index':4,
                                             'vspaero':{'WakeNumIter':5,'RootWakeNodes':8,'FixedWakeFlag':False}}
    baseline=tmp_path/'base.json';baseline.write_text(json.dumps(case))
    preview=prepare_cases('WAKE-01',baseline,root,dry_run=True)
    assert len(preview)==3
    assert not (root/'generated/WAKE-01').exists()
    prepared=prepare_cases('WAKE-01',baseline,root)
    assert {p['wake']['WakeNumIter'] for p in prepared}=={5,10,20}
    for item in prepared:
        document=json.loads(Path(item['case_json']).read_text())
        db.register_case_document('VSP3505-VV-001',document)
    from openvsppy.vv.VVAutomation import case_difference
    assert case_difference(db,'WAKE-01')==[]
    assert prepare_cases('WAKE-01',baseline,root)==prepared
    grid03=prepare_cases('GRID-03',baseline,root,dry_run=True)
    assert [x['case_id'] for x in grid03]==['GRID-03-BASE','GRID-03-LOCAL_A',
                                           'GRID-03-LOCAL_B','GRID-03-NREF_1']


def test_prepare_ctrl_requires_verified_groups(study,tmp_path, configured_case):
    from openvsppy.vv.VVAutomation import prepare_cases
    root,db,cond=configured_case
    case=json.loads((ROOT/'vv_study/case.example.json').read_text())
    case['geometry']['vsp3_source']=str(ROOT/'examples/models/CRM-HS/CRM-HS.vsp3')
    case['flight_condition']=cond
    case['numerical_setting']['requested']={'control_group_names':{},'control_sign_verified':False}
    baseline=tmp_path/'base.json';baseline.write_text(json.dumps(case))
    with pytest.raises(ValueError,match='approved control'):
        prepare_cases('CTRL-01',baseline,root,dry_run=True)
    case['numerical_setting']['requested']={'control_group_names':{'elevator':'TEST_GROUP'},'control_sign_verified':True}
    baseline.write_text(json.dumps(case))
    items=prepare_cases('CTRL-01',baseline,root)
    assert len(items)==7
    controls=[json.loads(Path(p['case_json']).read_text())['flight_condition']['controls_deg']['elevator'] for p in items]
    assert controls==[-10.,-5.,-2.,0.,2.,5.,10.]


def test_prepare_val_needs_explicit_span(study,tmp_path, configured_case):
    from openvsppy.vv.VVAutomation import prepare_cases
    root,_,cond=configured_case
    case=json.loads((ROOT/'vv_study/case.example.json').read_text())
    case['flight_condition']=cond
    baseline=tmp_path/'base.json';baseline.write_text(json.dumps(case))
    with pytest.raises(ValueError,match='span'):
        prepare_cases('VAL-01',baseline,root,dry_run=True)
    vals=prepare_cases('VAL-01',baseline,root,dry_run=True,span=2.)
    assert len(vals)==3


def test_wake_02_preserves_reviewed_baseline_iterations(study,tmp_path, configured_case):
    from openvsppy.vv.VVAutomation import prepare_cases
    root,db,cond=configured_case
    case=json.loads((ROOT/'vv_study/case.example.json').read_text())
    case['flight_condition']=cond
    case['numerical_setting']['requested']['vspaero']={'WakeNumIter':20,'RootWakeNodes':8}
    baseline=tmp_path/'baseline.json';baseline.write_text(json.dumps(case))
    cases=prepare_cases('WAKE-02',baseline,root,dry_run=True)
    assert [c['wake']['RootWakeNodes'] for c in cases]==[8,12,16]
    assert all(c['wake']['WakeNumIter']==20 for c in cases)
