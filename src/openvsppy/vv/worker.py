"""Worker executed in separate process with the *target* OpenVSP Python API."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def api_run(job_path):
    job=json.loads(Path(job_path).read_text(encoding='utf-8'))
    import openvsp as vsp  # intentionally absent from offline/unit-test environment
    from openvsppy.core.util import find_vspaero_settings_container, get_container_parm_value
    from openvsppy.vv.VVDatabase import resolve_ncpu
    version=str(vsp.GetVSPVersion())
    if '3.50.5' not in version:
        raise RuntimeError(f'Wrong OpenVSP API: expected 3.50.5, got {version!r}')
    # OpenVSP Python may locate VSPAERO beside the Python executable unless
    # explicitly directed to its solver directory (official API utility).
    binary = Path(job['runtime']['vspaero_binary']).resolve()
    if not binary.is_file():
        raise RuntimeError(f'Missing solver executable: {binary}')
    if not vsp.SetVSPAEROPath(str(binary.parent)):
        raise RuntimeError(f'Cannot select solver directory: {binary.parent}')
    if Path(vsp.GetVSPAEROPath()).resolve() != binary.parent:
        raise RuntimeError('Effective VSPAERO executable directory mismatch')
    vsp.ClearVSPModel()
    vsp.ReadVSPFile(job['input_vsp3'])
    vsp.Update()
    condition=job['condition'];setts=job['setting']['requested']
    ncpu=resolve_ncpu(setts, job['runtime'].get('ncpu'))
    membership=setts.get('geometry_membership')
    if membership is not None:
        if not setts.get('geometry_membership_reviewed',False):
            raise RuntimeError('Explicit geometry selection needs reviewed name-to-set mapping')
        selected={name for group in ('thin','thick') for name in membership[group]}
        if len(selected)!=sum(len(membership[group]) for group in ('thin','thick')):
            raise ValueError('Geometry cannot belong to both Thin and Thick')
        available={vsp.GetGeomName(gid):gid for gid in vsp.FindGeoms()}
        if len(available)!=len(vsp.FindGeoms()) or not selected <= available.keys():
            raise ValueError('Missing or duplicate geometry names in model')
        if 'TailFairing' in selected:raise ValueError('TailFairing cannot be included in aero set')
        for name,gid in available.items():
            for group,index in (('thin',setts['thin_set_index']),('thick',setts['thick_set_index'])):
                if index!=-1:vsp.SetSetFlag(gid,int(index),name in membership[group])
        vsp.Update()

    # Apply reviewed geometric mesh sensitivity inputs before meshing.
    effective_geom=[]
    for item in setts.get('geom_parms',[]):
        ids=vsp.FindGeomsWithName(item['geometry'])
        if len(ids)!=1:raise RuntimeError('Expected unique geometry for mesh setting: '+item['geometry'])
        parm=vsp.FindParm(ids[0],item['name'],item['group'])
        if not parm:raise RuntimeError('Missing mesh parameter: '+str(item))
        vsp.SetParmVal(parm,float(item['value']))
        vsp.Update()
        actual=float(vsp.GetParmVal(parm))
        if abs(actual-float(item['value']))>1e-8:raise RuntimeError('Mesh Parm readback mismatch: '+str(item))
        effective_geom.append({**item,'actual':actual})
    settings_container=find_vspaero_settings_container(vsp)
    if not settings_container:raise RuntimeError('Cannot locate VSPAERO settings ParmContainer')
    def set_checked(name,value):
        _,parm=get_container_parm_value(vsp,settings_container,name)
        if not parm:raise RuntimeError('VSPAERO setting not found: '+name)
        vsp.SetParmVal(parm,float(value))
        actual=float(vsp.GetParmVal(parm))
        if abs(actual-float(value))>1e-7:
            raise RuntimeError(f'{name}: requested {value}, read {actual}')
        return actual
    effective={'geom_parms':effective_geom}
    for name,value in {'GeomSet':setts['thick_set_index'],
                       'ThinGeomSet':setts['thin_set_index'],
                       'Sref':condition['sref'],'bref':condition['bref'],'cref':condition['cref'],
                       'Xcg':condition['moment_x'],'Ycg':condition['moment_y'],'Zcg':condition['moment_z']}.items():
        effective[name]=set_checked(name,value)
    effective['RefFlag']=set_checked('RefFlag',int(vsp.MANUAL_REF))
    # Geometry-dependent control groups must not reference omitted surfaces.
    # Resolve membership using the actual OpenVSP Geometry Sets, not the case label.
    selected_ids = set()
    for index in (int(setts['thin_set_index']), int(setts['thick_set_index'])):
        if index != int(vsp.SET_NONE):
            selected_ids.update(vsp.GetGeomSetAtIndex(index))
    # Minimal test stubs do not implement GetGeomName; production 3.50.5 does.
    active_geometry_names = ({vsp.GetGeomName(gid) for gid in selected_ids}
                             if hasattr(vsp, 'GetGeomName') else None)
    all_geometry_names = ({vsp.GetGeomName(gid) for gid in vsp.FindGeoms()}
                          if active_geometry_names is not None else None)
    controls=json.loads(condition['controls_deg_json'])
    effective['controls']={}

    group_names=[vsp.GetVSPAEROControlGroupName(i)
                 for i in range(vsp.GetNumControlSurfaceGroups())]
    for short,angle in controls.items():
        if short not in setts.get('control_group_names',{}):
            if not float(angle):continue
            raise RuntimeError('No approved Control Surface Group for '+short)
        groupname=setts['control_group_names'][short]
        if group_names.count(groupname)!=1:
            raise RuntimeError('Expected exactly one existing Control Surface Group: '+groupname)
        group_index=group_names.index(groupname)
        # VSPAERO ControlSurfaceGroup is a Parm group of VSPAEROSettings,
        # not a standalone ParmContainer named after the user-visible group.
        parm_group=f'ControlSurfaceGroup_{group_index}'
        active=vsp.FindParm(settings_container,'ActiveFlag',parm_group)
        deflection=vsp.FindParm(settings_container,'DeflectionAngle',parm_group)
        if not active or not deflection:
            raise RuntimeError('Missing Group ActiveFlag/DeflectionAngle: '+groupname)
        surfaces=list(vsp.GetActiveCSNameVec(group_index))
        if not surfaces:raise RuntimeError('Control Group has no active surfaces: '+groupname)
        if active_geometry_names is not None:
            from openvsppy.vv.VVNative import classify_control_surfaces
            present, absent = classify_control_surfaces(
                surfaces, active_geometry_names, all_geometry_names)
            if absent:
                if present:
                    raise RuntimeError('Partial Control Group geometry selection: '+groupname+
                        '; present='+str(present)+'; absent='+str(absent))
                if abs(float(angle))>1e-10:
                    raise RuntimeError('Nonzero deflection for omitted Control Group: '+groupname)
                vsp.SetParmVal(active,0.0)
                vsp.Update()
                if bool(vsp.GetParmVal(active)):
                    raise RuntimeError('Omitted Control Group could not be deactivated: '+groupname)
                effective['controls'][short]={'group':groupname,'active':False,
                    'angle_deg':float(angle),'surfaces':surfaces,'reason':'GEOMETRY_NOT_SELECTED'}
                continue
        approved=setts.get('control_gains_reviewed',{}).get(short)
        if not isinstance(approved,dict) or not approved:
            raise RuntimeError('Missing reviewed per-surface Gain mapping: '+short)
        gains={}
        for name,expected in approved.items():
            parm=vsp.FindParm(settings_container,name,parm_group)
            if not parm:raise RuntimeError('Missing Control Surface Gain Parm: '+name)
            actual=float(vsp.GetParmVal(parm))
            if abs(actual-float(expected))>1e-8:
                raise RuntimeError('Control Gain differs from approved mapping: '+name)
            gains[name]=actual
        # GetParmGroupName returns the generic "ControlSurfaceGroup" for all
        # five G103A gain Parms (including other control groups).  It does NOT
        # return the indexed group name used by FindParm.  Verify the approved
        # per-group gain Parms above, then check cardinality against the API's
        # active surfaces for precisely this control group.
        if len(gains)!=len(surfaces):
            raise RuntimeError('Control Group surface/Gain count differs: '+groupname)
        vsp.SetParmVal(active,1.0)
        vsp.SetParmVal(deflection,float(angle))
        vsp.Update()
        if not bool(vsp.GetParmVal(active)) or abs(vsp.GetParmVal(deflection)-float(angle))>1e-8:
            raise RuntimeError('Control Group ActiveFlag/DeflectionAngle readback mismatch')
        effective['controls'][short]={'group':groupname,'active':True,'angle_deg':float(angle),
                                     'surfaces':surfaces,'gains':gains}

    # Explicit, source-verified optional settings. Values omitted in a case
    # intentionally retain their versioned .vsp3 settings.
    optional={'NRef':'NRef','Symmetry':'Symmetry','ImplicitWake':'ImplicitWake',
              'WakeRelax':'WakeRelax','Vinf':'Vinf','Rho':'Rho'}
    for key,parm_name in optional.items():
        if key in setts.get('vspaero',{}):
            effective[key]=set_checked(parm_name,setts['vspaero'][key])
    vsp.Update()
    # Snapshot the actual modified model in the isolated Run directory.
    raw=Path(job['output_dir']).resolve()
    raw.mkdir(parents=True,exist_ok=True)
    # ProcessUtil::ForkCmd passes a null lpCurrentDirectory to CreateProcess.
    # Thus the nested vspaero.exe inherits the worker cwd.  The previous
    # worker cwd was rundir, unlike the proven direct solver run in raw/.
    # This isolated worker performs no unrelated cases after the chdir.
    os.chdir(raw)
    vsp.WriteVSPFile(str(raw/'model_effective.vsp3'))
    thick=int(setts['thick_set_index']);thin=int(setts['thin_set_index'])
    if thick==thin:raise ValueError('Overlap of Thick/Thin set selection')
    none=int(vsp.SET_NONE)
    if min(thin,thick)<none or max(thin,thick)>=int(vsp.GetNumSets()):
        raise ValueError('VSPAERO Set index outside defined set range')
    thin_members=set() if thin==none else set(vsp.GetGeomSetAtIndex(thin))
    thick_members=set() if thick==none else set(vsp.GetGeomSetAtIndex(thick))
    if not thin_members and not thick_members:
        raise RuntimeError('At least one Thin or Thick set must contain geometry')
    if thin!=none and not thin_members or thick!=none and not thick_members:
        raise RuntimeError('A selected non-NONE set has no geometry')
    if thin_members & thick_members:
        raise RuntimeError(f'Thin/Thick geometry membership overlap: {sorted(thin_members & thick_members)}')
    effective['geometry_sets']={'thin':sorted(thin_members),'thick':sorted(thick_members),
        'thin_name':'SET_NONE' if thin==none else str(vsp.GetSetName(thin)),
        'thick_name':'SET_NONE' if thick==none else str(vsp.GetSetName(thick))}
    analysis='VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(analysis)
    compute=set(vsp.GetAnalysisInputNames(analysis))
    for key,val in {'GeomSet':thick,'ThinGeomSet':thin}.items():
        if key not in compute:raise RuntimeError(f'ComputeGeometry input missing: {key}')
        vsp.SetIntAnalysisInput(analysis,key,[val],0)
    geometry_result_id=vsp.ExecAnalysis(analysis)
    # ComputeGeometry must have generated this Run's actual solver geometry.
    geom_path=raw/'model_effective.vspgeom'
    if not geom_path.is_file() or geom_path.stat().st_size==0:
        raise RuntimeError('VSPAEROComputeGeometry did not generate model_effective.vspgeom')
    # Check the native geometry before starting vspaero.exe; retain diagnostics
    # even when a malformed NRef mesh prevents an executable solver run.
    from openvsppy.vv.VVNative import inspect_native_mesh
    mesh_diagnostics=inspect_native_mesh(geom_path, raw/'model_effective.vkey')
    (raw/'mesh_diagnostics.json').write_text(
        json.dumps(mesh_diagnostics,indent=2),encoding='utf-8')
    if mesh_diagnostics['status']=='INVALID':
        raise RuntimeError('VSPAERO mesh preflight failed: '+
                           '; '.join(mesh_diagnostics['issues']))
    name='VSPAEROSweep';vsp.SetAnalysisInputDefaults(name)
    available=set(vsp.GetAnalysisInputNames(name))
    def integer(key,value):
        if key not in available:raise RuntimeError('Sweep input missing: '+key)
        vsp.SetIntAnalysisInput(name,key,[int(value)],0)
    def real(key,value):
        if key not in available:raise RuntimeError('Sweep input missing: '+key)
        vsp.SetDoubleAnalysisInput(name,key,[float(value)],0)
    integer('GeomSet',thick);integer('ThinGeomSet',thin)
    for pre,value in (('Alpha',condition['alpha_deg']),('Beta',condition['beta_deg']),('Mach',condition['mach'])):
        real(pre+'Start',value);real(pre+'End',value);integer(pre+'Npts',1)
    real('ReCref',condition['re_cref']); real('ReCrefEnd',condition['re_cref']);integer('ReCrefNpts',1)
    wake=setts.get('vspaero',{})
    for req,key in [('WakeNumIter','WakeNumIter'),('RootWakeNodes','NumWakeNodes'),('FixedWakeFlag','FixedWakeFlag')]:
        if req in wake:integer(key,int(wake[req]))
    if job['setting']['analysis_mode']=='STEADY_SWEEP':
        integer('UnsteadyType',vsp.STABILITY_OFF)
    elif job['setting']['analysis_mode']=='STEADY':
        # STA-01's steady stability derivatives intentionally use solver -stab.
        integer('UnsteadyType',vsp.STABILITY_DEFAULT)
    elif job['setting']['analysis_mode']=='ADJOINT':
        integer('UnsteadyType',vsp.STABILITY_ADJOINT)
    elif job['setting']['analysis_mode'] in ('P','Q','R'):
        mode=job['setting']['analysis_mode']
        integer('UnsteadyType',{'P':vsp.STABILITY_P_ANALYSIS,'Q':vsp.STABILITY_Q_ANALYSIS,
                                'R':vsp.STABILITY_R_ANALYSIS}[mode])
    else:
        raise ValueError('Unsupported analysis_mode: '+str(job['setting']['analysis_mode']))
    # Default 1 preserves the G103A baseline; explicit Run/Case NCPU is recorded.
    integer('NCPU',ncpu)
    # Official VSPAEROSweep output redirection: isolate the C++ monitor output
    # from the worker's Python stdout and retain a complete native API log.
    solver_log=raw/'vspaero_api_stdout.log'
    if 'RedirectFile' not in available:
        raise RuntimeError('VSPAEROSweep RedirectFile input missing')
    vsp.SetStringAnalysisInput(name,'RedirectFile',[str(solver_log)],0)
    vsp.Update()
    expected_mode={'STEADY_SWEEP':vsp.STABILITY_OFF,'STEADY':vsp.STABILITY_DEFAULT,
                   'ADJOINT':vsp.STABILITY_ADJOINT,
                   'P':vsp.STABILITY_P_ANALYSIS,'Q':vsp.STABILITY_Q_ANALYSIS,
                   'R':vsp.STABILITY_R_ANALYSIS}[job['setting']['analysis_mode']]
    if vsp.GetStringAnalysisInput(name,'RedirectFile')[0]!=str(solver_log):
        raise RuntimeError('VSPAEROSweep RedirectFile readback mismatch')
    for key,expected in (('UnsteadyType',expected_mode),('NCPU',ncpu)):
        actual=int(vsp.GetIntAnalysisInput(name,key)[0])
        if actual!=int(expected):
            raise RuntimeError(f'VSPAEROSweep {key} readback: expected {expected}, got {actual}')
    for key,value in (('AlphaStart',condition['alpha_deg']),('BetaStart',condition['beta_deg']),('MachStart',condition['mach'])):
        if abs(float(vsp.GetDoubleAnalysisInput(name,key)[0])-float(value))>1e-8:
            raise RuntimeError('Sweep effective input mismatch: '+key)
    effective['sweep']={'analysis_mode':job['setting']['analysis_mode'],
                        'UnsteadyType':int(expected_mode),'NCPU':ncpu,
                        'available_inputs':sorted(available)}

    # Persist evidence BEFORE the native invocation.  If vspaero.exe exits
    # abnormally or the API returns without results, the diagnostic survives.
    from openvsppy.vv.VVDatabase import digest
    from openvsppy.vv.VVNative import required_native_files, read_mode_derivatives
    diagnostic_path=raw/'api_diagnostics.json'
    native_input_suffixes=('.vsp3','.vspgeom','.vspaero','.vkey','.csf',
                           '.ControlSurfaces.taglist','.ALL.taglist')
    def inventory(suffixes):
        return {p.name:{'size_bytes':p.stat().st_size,'sha256':digest(p)}
                for suffix in suffixes
                for p in (raw/('model_effective'+suffix),)
                if p.is_file()}
    diagnostic={
        'state':'BEFORE_SWEEP',
        'openvsp_version':version,
        'geometry_result_id':str(geometry_result_id),
        'analysis_mode':job['setting']['analysis_mode'],
        'unsteady_type':int(expected_mode),
        'ncpu':ncpu,
        'solver_binary':str(binary),
        'solver_sha256':digest(binary),
        'expected_solver_args':['-omp',str(ncpu)] +
            {'STEADY_SWEEP':[],'STEADY':['-stab'],'ADJOINT':['-adjoint'],
             'P':['-pstab'],'Q':['-qstab'],'R':['-rstab']}[job['setting']['analysis_mode']] +
            [str(raw/'model_effective')],
        'api_monitor_log':str(solver_log),
        'expected_args_note':'Derived from source; not an observed child command line',
        'worker_cwd':str(Path.cwd()),
        'worker_pid':os.getpid(),
        'python_executable':sys.executable,
        'environment':{k:os.environ.get(k) for k in
            ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','PATH','VSPAERO_EXE')},
        'native_inputs_before_sweep':{**inventory(native_input_suffixes),
            **{p.name:{'size_bytes':p.stat().st_size,'sha256':digest(p)}
               for p in sorted(raw.glob('model_effective*.tag'))}},
        'native_solver_exit_code':'NOT_EXPOSED_BY_PYTHON_API',
    }
    diagnostic_path.write_text(json.dumps(diagnostic,ensure_ascii=False,indent=2),encoding='utf-8')
    print('VSPAERO_API_DIAGNOSTIC before sweep: mode=%s ncpu=%s input=%s' %
          (job['setting']['analysis_mode'],ncpu,str(raw/'model_effective')),flush=True)
    sweep_result_id=vsp.ExecAnalysis(name)
    diagnostic['state']='RETURNED_FROM_SWEEP'
    diagnostic['sweep_result_id']=str(sweep_result_id)
    diagnostic['native_outputs_after_sweep']=inventory(('.polar','.history','.lod','.adb',
                                                        '.stab','.adjoint.stab','.pstab','.qstab','.rstab'))
    diagnostic_path.write_text(json.dumps(diagnostic,ensure_ascii=False,indent=2),encoding='utf-8')
    # No broad glob across directories: only this Run's basename is accepted.
    # In particular a leftover polar from another Case cannot make this Run succeed.
    mode=job['setting']['analysis_mode']
    required=required_native_files(mode)
    missing=[suffix for suffix in required
             if not (raw/('model_effective'+suffix)).is_file()
             or (raw/('model_effective'+suffix)).stat().st_size==0]
    if missing:
        from openvsppy.vv.VVNative import diagnose_solver_log
        native_failure=diagnose_solver_log(solver_log)
        diagnostic['native_failure']=native_failure
        diagnostic['state']='NATIVE_OUTPUT_MISSING'
        diagnostic_path.write_text(json.dumps(diagnostic,ensure_ascii=False,indent=2),encoding='utf-8')
        raise RuntimeError('Native results missing after VSPAEROSweep returned: '+
            ', '.join(missing)+'; diagnostic='+native_failure['classification']+
            '; stage='+native_failure['last_stage'])
    if mode in ('P','Q','R'):
        # Command-line defaults override these parts of the written .vspaero
        # input (OpenVSP_3.50.5 src/vsp_aero/Solver/vspaero.C, -p/q/rstab).
        # These values are source-derived, NOT read back from the binary.
        diagnostic['source_defined_solver_overrides']={
            'NumberOfTimeSteps':128,'Unsteady_AngleMax_deg':5.0,
            'WakeIterations':25,'NumWakeNodes':32,
            'evidence':'OpenVSP 3.50.5 vspaero.C CLI branches; native readback unverified'}
        effective['native_solver_source_overrides']=diagnostic['source_defined_solver_overrides']
    if mode!='STEADY_SWEEP':
        header, rows, groups=read_mode_derivatives(raw/('model_effective'+required[0]),mode)
        diagnostic['native_derivative_count']=len(rows)
        diagnostic['native_stability_conditions']=header
        diagnostic_path.write_text(json.dumps(diagnostic,ensure_ascii=False,indent=2),encoding='utf-8')
    if mode=='STEADY_SWEEP' and not int(wake.get('FixedWakeFlag',0)) and 'WakeNumIter' in wake:
        from openvsppy.vv.VVPostprocess import read_history
        history_rows=read_history(raw/'model_effective.history')
        last_iteration=int(history_rows[-1]['Iter']) if history_rows else -1
        diagnostic['last_history_iteration']=last_iteration
        diagnostic_path.write_text(json.dumps(diagnostic,ensure_ascii=False,indent=2),encoding='utf-8')
        if last_iteration!=int(wake['WakeNumIter']):
            raise RuntimeError(f'Wake history incomplete: last={last_iteration}, expected={wake["WakeNumIter"]}')
    generated_files=sorted(p.name for p in raw.iterdir() if p.is_file())
    generated={'captured_files':generated_files,
               'vspaero_directory':str(binary.parent),
               'api_diagnostic_file':diagnostic_path.name,
               'native_solver_exit_code':'NOT_EXPOSED_BY_PYTHON_API'}
    settings_files=list(raw.glob('*.vspaero'))
    if len(settings_files)!=1:
        raise RuntimeError(f'Expected one generated .vspaero settings file, got {len(settings_files)}')
    from openvsppy.vv.VVPostprocess import read_vspaero_settings
    generated['solver_settings']=read_vspaero_settings(settings_files[0])
    settings=generated['solver_settings']
    for req,actual in (('Sref','Sref'),('bref','Bref'),('cref','Cref'),
                        ('moment_x','X_cg'),('moment_y','Y_cg'),('moment_z','Z_cg'),
                        ('mach','Mach'),('alpha_deg','AoA'),('beta_deg','Beta'),
                        ('re_cref','ReCref')):
        target=condition[req] if req in condition else effective[req]
        if actual not in settings or abs(float(settings[actual])-float(target))>max(1e-7,abs(float(target))*1e-7):
            raise RuntimeError(f'Generated VSPAERO input mismatch: {req} / {actual}')
    if not int(wake.get('FixedWakeFlag',0)) and 'WakeNumIter' in wake:
        if int(settings.get('WakeIters',-999))!=int(wake['WakeNumIter']):
            raise RuntimeError('Generated WakeIters mismatch')
    if 'RootWakeNodes' in wake and int(settings.get('NumWakeNodes',-999))!=int(wake['RootWakeNodes']):
        raise RuntimeError('Generated NumWakeNodes mismatch')
    for key in ('ImplicitWake','WakeRelax','Symmetry','Vinf','Rho'):
        if key in wake:
            if key not in settings or abs(float(settings[key])-float(wake[key]))>max(1e-8,abs(float(wake[key]))*1e-8):
                raise RuntimeError(f'Generated {key} mismatch')
    result={'openvsp_version':'3.50.5','effective':effective,'generated':generated}
    (raw/'worker_result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')


def api_drag(job_path):
    """Isolated ParasiteDrag execution; no VSPAERO Run success is implied."""
    import openvsp as vsp
    from openvsppy.aerodynamics.AnalysisParasiteDrag import vsp_parasite_drag
    from openvsppy.vv.VVDatabase import digest
    job=json.loads(Path(job_path).read_text(encoding='utf-8'))
    if '3.50.5' not in str(vsp.GetVSPVersion()):
        raise RuntimeError('ParasiteDrag requires OpenVSP Python API 3.50.5')
    vsp.ClearVSPModel();vsp.ReadVSPFile(job['model']);vsp.Update()
    names=job.get('geometry_names')
    if names:
        available={vsp.GetGeomName(gid):gid for gid in vsp.FindGeoms()}
        if set(names)-set(available) or 'TailFairing' in names:
            raise ValueError('ParasiteDrag geometry names are missing or include TailFairing')
        if len(names)!=len(set(names)) or len(available)!=len(vsp.FindGeoms()):
            raise ValueError('Ambiguous geometry names for ParasiteDrag')
        if int(job['geom_set'])<0 or int(job['geom_set'])>=vsp.GetNumSets():
            raise ValueError('ParasiteDrag set index invalid')
        for name,gid in available.items():
            vsp.SetSetFlag(gid,int(job['geom_set']),name in names)
        vsp.Update()
        actual={vsp.GetGeomName(gid) for gid in vsp.GetGeomSetAtIndex(int(job['geom_set']))}
        if actual!=set(names):raise RuntimeError('ParasiteDrag set membership mismatch')
    df=vsp_parasite_drag(vsp,job['velocity_mps'],job['altitude_m'],job['delta_temp_k'],
                         sref=job['sref'],geom_set=job['geom_set'],verbose=0)
    if len(df)!=1:raise RuntimeError('One drag result expected for one condition')
    output=Path(job['output']);output.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(output,index=False)
    row=df.iloc[0].to_dict()
    meta={'openvsp_version':'3.50.5','model_sha256':digest(job['model']),
          'velocity_mps':job['velocity_mps'],'altitude_m':job['altitude_m'],
          'delta_temp_k':job['delta_temp_k'],'geom_set':job['geom_set'],'geometry_names':names,
          'result':{key:float(val) for key,val in row.items()}}
    Path(job['result_manifest']).write_text(json.dumps(meta,indent=2),encoding='utf-8')



def api_benchmark(job_path):
    """Create a versioned VER geometry with OpenVSP 3.50.5, never a G103A proxy.

    The model and its creation metadata are evidence, not a validation result.
    Geometry dimensions are verified by OpenVSP bounding boxes on the target host.
    """
    import openvsp as vsp
    from openvsppy.vv.VVDatabase import digest
    job=json.loads(Path(job_path).read_text(encoding='utf-8'))
    if '3.50.5' not in str(vsp.GetVSPVersion()):
        raise RuntimeError('VER benchmark geometry must be generated by OpenVSP 3.50.5')
    output=Path(job['output']).resolve()
    if output.exists():raise FileExistsError('Benchmark model already exists; do not overwrite: '+str(output))
    vsp.ClearVSPModel()
    kind=job['kind']
    if kind=='HERSHEY_BAR':
        wing=vsp.AddGeom('WING')
        vsp.SetGeomName(wing,'HersheyWing')
        # Follow the OpenVSP official Hershey example: AR driver, 1 m chords,
        # symmetric NACA 0012 and unswept wing. Full AR is 10, half-wing AR=5.
        vsp.SetDriverGroup(wing,1,vsp.AR_WSECT_DRIVER,
                           vsp.ROOTC_WSECT_DRIVER,vsp.TIPC_WSECT_DRIVER)
        for name,group,value in [('Aspect','XSec_1',5),
            ('Root_Chord','XSec_1',1),('Tip_Chord','XSec_1',1),
            ('Sweep','XSec_1',0),('ThickChord','XSecCurve_0',.12),
            ('ThickChord','XSecCurve_1',.12),('LECluster','WingGeom',.2),
            ('TECluster','WingGeom',1)]:
            parm=vsp.FindParm(wing,name,group)
            if not parm:raise RuntimeError('Hershey geometry Parm not found: '+name+'/'+group)
            vsp.SetParmVal(parm,value)
        vsp.SetSetFlag(wing,4,True)
        geom_id=wing
        expected=(1.,10.,.12)
    elif kind=='ELLIPSOID':
        geom_id=vsp.AddGeom('ELLIPSOID')
        vsp.SetGeomName(geom_id,'EllipsoidGeom')
        for name,value in [('A_Radius',1.),('B_Radius',2.),('C_Radius',3.)]:
            parm=vsp.FindParm(geom_id,name,'Design')
            if not parm:raise RuntimeError('Ellipsoid 3.50.5 Parm not found: '+name)
            vsp.SetParmVal(parm,value)
        for name,value in [('Tess_U',40),('Tess_W',41)]:
            parm=vsp.FindParm(geom_id,name,'Shape')
            if not parm:raise RuntimeError('Ellipsoid 3.50.5 tessellation Parm not found: '+name)
            vsp.SetParmVal(parm,value)
        vsp.SetSetFlag(geom_id,3,True)
        expected=(2.,4.,6.)
    else:
        raise ValueError('Unknown benchmark geometry kind: '+str(kind))
    vsp.Update()
    minimum=vsp.GetGeomBBoxMin(geom_id)
    maximum=vsp.GetGeomBBoxMax(geom_id)
    width=[maximum.x()-minimum.x(),maximum.y()-minimum.y(),maximum.z()-minimum.z()]
    # Ellipsoid bounding-box extents must reproduce the exact reference axes.
    if kind=='ELLIPSOID' and any(abs(a-b)>1e-5 for a,b in zip(width,expected)):
        raise RuntimeError('Ellipsoid geometry differs from 1/2/3 semi-axes: '+str(width))
    if kind=='HERSHEY_BAR':
        # GetGeomBBox* reports the unmirrored wing's half span, 5 m on 3.50.5.
        # Do not confuse that with the full, symmetric 10 m reference span.
        span_parm=vsp.FindParm(geom_id,'TotalSpan','WingGeom')
        sym_parm=vsp.FindParm(geom_id,'Sym_Planar_Flag','Sym')
        if abs(width[1]-5.)>1e-3 or not span_parm or not sym_parm or \
                abs(float(vsp.GetParmVal(span_parm))-10.)>1e-3 or \
                not int(round(float(vsp.GetParmVal(sym_parm)))):
            raise RuntimeError('Hershey symmetry/full-span verification failed: '+str(width))
    output.parent.mkdir(parents=True,exist_ok=True)
    vsp.WriteVSPFile(str(output),vsp.SET_ALL)
    if not output.is_file():raise RuntimeError('OpenVSP did not write the benchmark model')
    evidence={'openvsp_version':str(vsp.GetVSPVersion()),'kind':kind,
        'geometry_id':geom_id,'bounding_box_xyz':width,
        'full_span_m':float(vsp.GetParmVal(span_parm)) if kind=='HERSHEY_BAR' else None,
        'symmetry_flag':int(round(vsp.GetParmVal(sym_parm))) if kind=='HERSHEY_BAR' else None,
        'sha256':digest(output),
        'source':'OpenVSP_3.50.5 official benchmark geometry, programmatically constructed',
        'engineering_approval':'NOT_ASSESSED'}
    Path(job['result_manifest']).write_text(json.dumps(evidence,indent=2),encoding='utf-8')


def main():
    p=argparse.ArgumentParser();p.add_argument('--job');p.add_argument('--drag-job');p.add_argument('--benchmark-job');a=p.parse_args()
    if a.benchmark_job:api_benchmark(a.benchmark_job)
    elif a.drag_job:api_drag(a.drag_job)
    elif a.job:api_run(a.job)
    else:p.error('--job or --drag-job required')


if __name__=='__main__':main()
