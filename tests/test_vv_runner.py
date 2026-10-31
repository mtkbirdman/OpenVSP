"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path
import json
import subprocess
import sys

import pytest

from openvsppy.vv.VVRunner import check_case
from openvsppy.vv.worker import api_run
import openvsppy.core.util as util

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('wing_only', [False, True])
def test_control_group_gain_check_uses_active_group_not_generic_parm_label(tmp_path, monkeypatch, wing_only):
    """All 5 native gains report ControlSurfaceGroup, not the indexed group."""
    from openvsppy.vv.worker import api_run
    import openvsppy.core.util as util

    class EndOfScope(Exception): pass
    class FakeOpenVSP:
        MANUAL_REF = 1
        def __init__(self):
            self.root = str(tmp_path)
            self.parms = {}
            self.groups = ['AILERON_GROUP', 'ELEVATOR_GROUP', 'RUDDER_GROUP']
            self.geom_names={'wing':'WingGeom','htail':'HTailGeom',
                             'vtail':'VTailGeom','fuselage':'FuselageGeom'}
        def GetVSPVersion(self): return 'OpenVSP 3.50.5'
        def SetVSPAEROPath(self, value): self.path = value; return True
        def GetVSPAEROPath(self): return self.path
        def ClearVSPModel(self): pass
        def ReadVSPFile(self, value): pass
        def Update(self): pass
        def FindParm(self, container, name, group): return f'{group}:{name}'
        def GetParmVal(self, parm): return self.parms.get(parm, 1. if 'Gain' in parm else 0.)
        def SetParmVal(self, parm, val): self.parms[parm] = float(val)
        def GetNumControlSurfaceGroups(self): return 3
        def GetVSPAEROControlGroupName(self, i): return self.groups[i]
        def GetActiveCSNameVec(self, i): return [
            ['WingGeom_Surf0_AILERON','WingGeom_Surf1_AILERON'],
            ['HTailGeom_Surf0_ELEVATOR','HTailGeom_Surf1_ELEVATOR'],
            ['VTailGeom_Surf0_RUDDER']][i]
        def FindGeoms(self): return list(self.geom_names)
        def GetGeomName(self, gid): return self.geom_names[gid]
        def WriteVSPFile(self, path): Path(path).write_text('effective')
        def GetNumSets(self): return 5
        SET_NONE = -1
        def GetGeomSetAtIndex(self, i): return (['wing'] if wing_only else ['wing','htail','vtail']) if i == 4 else ['fuselage']
        def GetSetName(self, i): return str(i)
        def SetAnalysisInputDefaults(self, value): raise EndOfScope('Gain checks passed')
        def FindContainerParmIDs(self, settings):
            # This falsely enumerates all 5 gains as the *same* group.
            return [f'p{i}' for i in range(5)]
        def GetParmGroupName(self, parm): return 'ControlSurfaceGroup'
        def GetParmName(self, parm): return f'Surf_{parm}_Gain'

    fake = FakeOpenVSP()
    monkeypatch.setitem(sys.modules, 'openvsp', fake)
    monkeypatch.setattr(util, 'find_vspaero_settings_container', lambda vsp: 'settings')
    monkeypatch.setattr(util, 'get_container_parm_value',
                        lambda vsp, container, name: (0, f'core:{name}'))
    binary = tmp_path / 'vspaero.exe'
    binary.write_bytes(b'test executable')
    condition = dict(sref=10., bref=12., cref=1., moment_x=0., moment_y=0., moment_z=0.,
                     controls_deg_json=json.dumps({'aileron': 0., 'elevator': 0., 'rudder': 0.}))
    groups = ['AILERON_GROUP', 'ELEVATOR_GROUP', 'RUDDER_GROUP']
    requested = {'thin_set_index': 4, 'thick_set_index': 3,
                 'control_group_names': dict(zip(('aileron','elevator','rudder'), groups)),
                 'control_gains_reviewed': {
                     'aileron': {'a0_Gain': 1., 'a1_Gain': 1.},
                     'elevator': {'e0_Gain': 1., 'e1_Gain': 1.},
                     'rudder': {'r0_Gain': 1.}}
                 }
    # Simulate the actual default G103A elevator second gain=-1; readback
    # comparison must still hold while all generic Parm labels are identical.
    fake.parms['ControlSurfaceGroup_1:e1_Gain'] = -1.
    requested['control_gains_reviewed']['elevator']['e1_Gain'] = -1.
    job = tmp_path / 'job.json'
    job.write_text(json.dumps({'input_vsp3': 'input.vsp3', 'output_dir': str(tmp_path / 'raw'),
                               'runtime': {'vspaero_binary': str(binary)},
                               'setting': {'requested': requested, 'analysis_mode': 'STEADY_SWEEP'},
                               'condition': condition}))
    with pytest.raises(EndOfScope, match='Gain checks passed'):
        api_run(job)
    assert fake.GetParmVal('ControlSurfaceGroup_0:ActiveFlag')==1.
    for i in (1,2):
        assert fake.GetParmVal(f'ControlSurfaceGroup_{i}:ActiveFlag')==(0 if wing_only else 1.)


@pytest.mark.parametrize("analysis_mode", ["STEADY_SWEEP", "STEADY"])
@pytest.mark.parametrize("ncpu", [1, 4])
def test_reviewed_analysis_modes_and_diagnostics(tmp_path, monkeypatch, analysis_mode, ncpu):
    class StubVSP:
        STABILITY_OFF=0
        STABILITY_DEFAULT=1
        STABILITY_ADJOINT=2
        STABILITY_P_ANALYSIS=3
        STABILITY_Q_ANALYSIS=4
        STABILITY_R_ANALYSIS=5
        MANUAL_REF=0
        SET_NONE=-1
        def __init__(self):
            self.path=''
            self.params={}
            self.inputs={}
            self.analysis_inputs={
                'VSPAEROComputeGeometry': {'GeomSet', 'ThinGeomSet'},
                'VSPAEROSweep': {'GeomSet','ThinGeomSet', 'AlphaStart','AlphaEnd','AlphaNpts',
                    'BetaStart','BetaEnd','BetaNpts','MachStart','MachEnd','MachNpts',
                    'ReCref','ReCrefEnd','ReCrefNpts','WakeNumIter','NumWakeNodes',
                    'FixedWakeFlag','UnsteadyType','NCPU','RedirectFile'} }
        def GetVSPVersion(self): return 'OpenVSP 3.50.5'
        def SetVSPAEROPath(self,path): self.path=path;return True
        def GetVSPAEROPath(self): return self.path
        def ClearVSPModel(self): pass
        def ReadVSPFile(self,path): pass
        def Update(self): pass
        def FindParm(self,container,name,group): return group+':'+name
        def SetParmVal(self,parm,val): self.params[parm]=float(val)
        def GetParmVal(self,parm): return self.params.get(parm,0.)
        def GetNumControlSurfaceGroups(self): return 0
        def WriteVSPFile(self,path): Path(path).write_text('synthetic vsp3')
        def GetNumSets(self):return 8
        def GetGeomSetAtIndex(self,index):return ['fus'] if index==3 else ['wing']
        def GetSetName(self,index):return str(index)
        def SetAnalysisInputDefaults(self,name): self.inputs[name]={}
        def GetAnalysisInputNames(self,name):return list(self.analysis_inputs[name])
        def SetIntAnalysisInput(self,name,key,values,index):self.inputs[name][key]=values
        def SetDoubleAnalysisInput(self,name,key,values,index):self.inputs[name][key]=values
        def SetStringAnalysisInput(self,name,key,values,index):self.inputs[name][key]=values
        def GetStringAnalysisInput(self,name,key):return self.inputs[name][key]
        def GetIntAnalysisInput(self,name,key):return self.inputs[name][key]
        def GetDoubleAnalysisInput(self,name,key):return self.inputs[name][key]
        def ExecAnalysis(self,name):
            if name=='VSPAEROComputeGeometry':
                (tmp_path/'raw/model_effective.vspgeom').write_text('# vspgeom v3\n1\n100 150 4\n')
                (tmp_path/'raw/model_effective.vkey').write_text('# VSPGEOM v3 Tag Key File\nmodel_effective.vspgeom\n2\n')
            return 'STUB_RESULT_NOT_NATIVE'
    monkeypatch.chdir(tmp_path)
    fake=StubVSP()
    monkeypatch.setitem(sys.modules,'openvsp',fake)
    monkeypatch.setattr(util,'find_vspaero_settings_container',lambda vsp:'settings')
    monkeypatch.setattr(util,'get_container_parm_value',
                        lambda vsp,container,name:(None, 'settings:'+name))
    binary=tmp_path/'vspaero.exe';binary.write_bytes(b'STUB-NOT-SOLVER')
    target=tmp_path/'raw';target.mkdir()
    job={
        'input_vsp3':str(tmp_path/'input.vsp3'),'output_dir':str(target),
        'runtime':{'vspaero_binary':str(binary),'ncpu':ncpu},
        'condition':{'sref':17.88185312625759,'bref':17.536552870371967,
            'cref':1.0624862547003904,'moment_x':2.87,'moment_y':0.,'moment_z':.164,
            'controls_deg_json':'{}','alpha_deg':2.,'beta_deg':0.,'mach':.0857102,
            're_cref':2120000.},
        'setting':{'analysis_mode':analysis_mode,
                   'requested':{'thick_set_index':3,'thin_set_index':4,
                                'vspaero':{'WakeNumIter':20,'RootWakeNodes':8}}}
    }
    job_path=tmp_path/'job.json';job_path.write_text(json.dumps(job))
    with pytest.raises(RuntimeError,match='Native results missing after VSPAEROSweep returned'):
        api_run(job_path)
    expected=fake.STABILITY_OFF if analysis_mode=='STEADY_SWEEP' else fake.STABILITY_DEFAULT
    assert fake.inputs['VSPAEROSweep']['UnsteadyType']==[expected]
    assert fake.inputs['VSPAEROSweep']['NCPU']==[ncpu]
    diagnostic=json.loads((target/'api_diagnostics.json').read_text())
    assert diagnostic['state']=='NATIVE_OUTPUT_MISSING'
    assert diagnostic['analysis_mode']==analysis_mode
    assert diagnostic['unsteady_type']==expected
    assert diagnostic['expected_solver_args'][:2]==['-omp',str(ncpu)]
    assert diagnostic['ncpu']==ncpu
    assert diagnostic['state']=='NATIVE_OUTPUT_MISSING'
    assert ('-stab' in diagnostic['expected_solver_args']) == (analysis_mode=='STEADY')
    assert diagnostic['native_outputs_after_sweep']=={}
    assert diagnostic['api_monitor_log'].endswith('vspaero_api_stdout.log')
    assert diagnostic['worker_cwd']==str(target)
    assert 'model_effective.vspgeom' in diagnostic['native_inputs_before_sweep']
    assert not (target/'worker_result.json').exists()


def test_runner_timeout_preserves_log_and_state(study,monkeypatch, configured_case):
    from openvsppy.vv.VVRunner import run_case
    import subprocess
    root,db,_=configured_case
    class FakeProcess:
        pid=1000
        def __init__(self,cmd,**kwargs):
            kwargs['stdout'].write(b'iteration 1\n')
        def wait(self,timeout=None):
            if timeout==15:return 9
            raise subprocess.TimeoutExpired('fake',timeout)
        def kill(self):pass
    monkeypatch.setattr('openvsppy.vv.VVRunner.subprocess.Popen',FakeProcess)
    monkeypatch.setattr('openvsppy.vv.VVRunner.subprocess.run',lambda *args,**kwargs: None)
    result=run_case(db.path,'GRID-01-MEDIUM',root/'study.json',binary_path=ROOT/'README.md',timeout_seconds=1)
    assert result['status']=='INTERRUPTED'
    row=db.row('run','run_id',result['run_id'])
    assert row['status']=='INTERRUPTED'
    assert db.verify_artifacts(root) and all(s=='VERIFIED' for _,s in db.verify_artifacts(root))
    logfile=root/row['run_dir']/'logs/stdout.log'
    assert logfile.read_bytes()==b'iteration 1\n'


def test_runner_worker_failure_preserves_logs(study,monkeypatch, configured_case):
    from openvsppy.vv.VVRunner import run_case
    root,db,_=configured_case
    class FailedProcess:
        pid=3000
        def __init__(self,cmd,**kwargs):kwargs['stderr'].write(b'bad worker\n')
        def wait(self,timeout=None):return 3
    monkeypatch.setattr('openvsppy.vv.VVRunner.subprocess.Popen',FailedProcess)
    result=run_case(db.path,'GRID-01-MEDIUM',root/'study.json',binary_path=ROOT/'README.md')
    assert result['status']=='EXECUTION_FAILED'
    assert 'bad worker' in result['stderr']
    assert all(s=='VERIFIED' for _,s in db.verify_artifacts(root))


def test_explicit_grid_precheck_rejects_both_none(study, configured_case):
    root,db,_=configured_case
    db.add_setting('VSP3505-VV-001',{'setting_id':'SET-NONE','analysis_mode':'STEADY_SWEEP','mesh_level':'G1',
                                     'requested':{'thin_set_index':-1,'thick_set_index':-1}})
    db.resolve_case('GRID-01-MEDIUM','G3','FC-001','SET-NONE')
    assert any('SET_NONE' in error for error in check_case(db,'GRID-01-MEDIUM',root/'study.json',ROOT/'README.md'))


def test_runner_spawn_failure_does_not_leave_running(study,monkeypatch, configured_case):
    from openvsppy.vv.VVRunner import run_case
    root,db,_=configured_case
    def fail_spawn(*args,**kwargs):raise FileNotFoundError('python executable missing')
    monkeypatch.setattr('openvsppy.vv.VVRunner.subprocess.Popen',fail_spawn)
    result=run_case(db.path,'GRID-01-MEDIUM',root/'study.json',binary_path=ROOT/'README.md')
    assert result['status']=='EXECUTION_FAILED'
    assert db.row('run','run_id',result['run_id'])['status']=='EXECUTION_FAILED'
    assert all(s=='VERIFIED' for _,s in db.verify_artifacts(root))


def test_inspect_run_on_interrupted_worker(study,monkeypatch,capsys, configured_case):
    from openvsppy.vv.VVRunner import run_case
    from openvsppy.vv.__main__ import main
    root,db,_=configured_case
    class FailProc:
        def __init__(self,*a,**kwargs):kwargs['stderr'].write(b'failed intentionally')
        def wait(self,timeout=None):return 1
    monkeypatch.setattr('openvsppy.vv.VVRunner.subprocess.Popen',FailProc)
    ret=run_case(db.path,'GRID-01-MEDIUM',root/'study.json',binary_path=ROOT/'README.md')
    assert ret['status']=='EXECUTION_FAILED'
    main(['--study',str(root/'study.json'),'--db',str(db.path),'inspect-run','--run',ret['run_id']])
    report=json.loads(capsys.readouterr().out)
    assert report['status']=='EXECUTION_FAILED'
    assert all(a['status']=='VERIFIED' for a in report['artifacts'])
