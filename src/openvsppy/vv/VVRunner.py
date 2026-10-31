"""Isolated OpenVSP 3.50.5 run orchestration with strict preflight checks."""
from __future__ import annotations

import json
import os
import re
import time
import shutil
import subprocess
import sys
from pathlib import Path

from .VVDatabase import VVDatabase, canonical, digest, load_study, REPO_ROOT, resolve_ncpu
from .VVPostprocess import import_polar, import_stab
from .VVNative import required_native_files


def check_case(db, case_id, study_file, binary_path=None):
    root,d,_,_=load_study(study_file)
    case=db.row('planned_case','case_id',case_id)
    errors=[]
    if not case: return ['Unknown Case ID']
    if case['execute_policy']!='EXECUTE': errors.append('Case not approved for execution')
    for k in ('geometry_id','condition_id','setting_id'):
        if case[k] is None: errors.append('Unresolved '+k)
    if errors:return errors
    g=db.row('geometry','geometry_id',case['geometry_id'])
    c=db.row('flight_condition','condition_id',case['condition_id'])
    setting=db.row('numerical_setting','setting_id',case['setting_id'])
    if any(r is None for r in (g,c,setting)):return ['Missing input records']
    geom=json.loads(g['physical_model_json'])
    path=Path(geom['source_file']).resolve()
    if not path.is_file() or digest(path)!=g['source_vsp3_sha256']:
        errors.append('Geometry file absent or hash changed')
    if c['trim_mode']!='FIXED_CONTROL': errors.append('Automatic trim not supported by initial one-point worker')
    if setting['analysis_mode'] not in ('STEADY_SWEEP','STEADY','ADJOINT','P','Q','R'):
        errors.append('This analysis mode is not yet implemented in the worker')
    requested=json.loads(setting['requested_json'])
    try:
        resolve_ncpu(requested)
    except ValueError as exc:
        errors.append(str(exc))
    if any(k not in requested for k in ('thin_set_index','thick_set_index')):
        errors.append('Exact Thin/Thick set indices must be explicitly configured')
    elif requested['thin_set_index']==requested['thick_set_index']:
        errors.append('Thin/Thick sets overlap or both use SET_NONE')
    elif min(requested['thin_set_index'],requested['thick_set_index']) < -1:
        errors.append('Set indices cannot be lower than vsp.SET_NONE (-1)')
    if setting['analysis_mode'] in ('P','Q','R'):
        wake=requested.get('vspaero',{})
        if not requested.get('unsteady_settings_reviewed') or any(float(wake.get(k,0))<=0 for k in ('Vinf','Rho')):
            errors.append('P/Q/R require reviewed unsteady settings and positive Vinf/Rho')
    if requested.get('mesh_overrides'):
        errors.append('Geometry mesh overrides must be pre-applied and saved to a versioned .vsp3')
    controls=json.loads(c['controls_deg_json'])
    if any(abs(float(x))>1e-10 for x in controls.values()):
        if not requested.get('control_group_names') or not requested.get('control_sign_verified',False) or not requested.get('control_gains_reviewed'):
            errors.append('Control Group IDs and physical sign verification are required for nonzero deflection')
    binary=binary_path or d['runtime'].get('vspaero_binary')
    if not binary or not Path(binary).is_file():errors.append('VSPAERO 3.50.5 executable must exist')
    return errors


def notify_progress(callback, event):
    """Best-effort progress transport; solver success never depends on terminal I/O."""
    if callback is not None:
        try:
            callback(event)
        except (OSError, UnicodeError):
            pass


def wait_for_worker(proc, cmd, timeout_seconds, native_log, *,
                    progress_callback=None, verbose=False, heartbeat_seconds=60,
                    poll_seconds=5):
    """Wait for the isolated worker; optionally observe its already-written solver log.

    Never read the child's stdout/stderr pipes: they remain redirected to files.
    Progress reporting is advisory; it cannot alter native success criteria.
    """
    start = time.monotonic()
    next_heartbeat = start + heartbeat_seconds
    position = 0
    pending = b''
    last_iteration = 0
    adjoint_active = False
    last_adjoint_stage = None
    last_adjoint_gmres = -1

    def read_wake():
        nonlocal position, pending, last_iteration, adjoint_active, last_adjoint_stage, last_adjoint_gmres
        try:
            with Path(native_log).open('rb') as stream:
                size = stream.seek(0, 2)
                if size < position:  # Log was replaced/truncated.
                    position, pending = 0, b''
                stream.seek(position)
                block = stream.read()
                position = stream.tell()
        except OSError:  # Solver may not have created its log yet.
            return
        if not block:
            return
        # OpenVSP 3.50.5 writes progress on CR-only lines as well as LF.
        # Reading binary chunks also survives an incomplete final line.
        parts = re.split(rb'[\r\n]+', pending + block)
        pending = parts.pop()[-256:]
        latest = None
        for line in parts:
            if b'Solving adjoint linear systems' in line:
                adjoint_active = True
                notify_progress(progress_callback, {'event':'ADJOINT_STAGE',
                    'stage':'SOLVING_LINEAR_SYSTEMS','elapsed_seconds':round(time.monotonic()-start,1)})
            if adjoint_active:
                stage_match = re.search(rb'(Inviscid (?:Force|Moment)\s+in\s+[xyz])',line,re.I)
                if stage_match:
                    stage = stage_match[1].decode('ascii', errors='replace').strip()
                    if stage != last_adjoint_stage:
                        last_adjoint_stage, last_adjoint_gmres = stage, -1
                        notify_progress(progress_callback,{'event':'ADJOINT_STAGE',
                            'stage':stage,'elapsed_seconds':round(time.monotonic()-start,1)})
                gmres_match = re.search(rb'GMRES Iter:\s*(\d+)',line)
                if gmres_match and last_adjoint_stage is not None:
                    index=int(gmres_match[1])
                    if index>=last_adjoint_gmres+10 or (last_adjoint_gmres<0 and index==0):
                        last_adjoint_gmres=index
                        notify_progress(progress_callback,{'event':'ADJOINT_GMRES',
                            'stage':last_adjoint_stage,'iteration':index,
                            'elapsed_seconds':round(time.monotonic()-start,1)})
            match = re.search(rb'Wake Iter:\s*(\d+)\s*/\s*(\d+)', line)
            if match:
                iteration, total = int(match[1]), int(match[2])
                if iteration > last_iteration:
                    latest = (iteration, total)
                    last_iteration = iteration
        if latest:
            iteration, total = latest
            notify_progress(progress_callback, {'event':'WAKE', 'iteration':iteration, 'total_iterations':total,
                    'elapsed_seconds':round(time.monotonic()-start, 1)})

    while True:
        elapsed = time.monotonic() - start
        remaining = timeout_seconds - elapsed
        if remaining <= 0:
            raise subprocess.TimeoutExpired(cmd, timeout_seconds)
        interval = poll_seconds if verbose else heartbeat_seconds
        try:
            code = proc.wait(timeout=min(interval, remaining))
        except subprocess.TimeoutExpired:
            now = time.monotonic()
            if verbose:
                read_wake()
            if now >= next_heartbeat:
                notify_progress(progress_callback, {'event':'HEARTBEAT','elapsed_seconds':round(now-start, 1)})
                next_heartbeat = now + heartbeat_seconds
            continue
        if verbose:
            read_wake()
        return code


def run_case(db_path, case_id, study_file, *, python_executable=None, binary_path=None, timeout_seconds=3600,
             progress_callback=None, verbose=False, ncpu=None):
    """Execute one approved Case in one Python process. A failed Run is never reused."""
    study_path,d,_,_=load_study(study_file)
    root=study_path.parent
    db=VVDatabase(db_path)
    try:
        errors=check_case(db,case_id,study_file,binary_path)
        if errors:raise ValueError('; '.join(errors))
        case=db.row('planned_case','case_id',case_id)
        geom=db.row('geometry','geometry_id',case['geometry_id'])
        cond=db.row('flight_condition','condition_id',case['condition_id'])
        setting=db.row('numerical_setting','setting_id',case['setting_id'])
        case_requested=json.loads(setting['requested_json'])
        effective_ncpu=resolve_ncpu(case_requested, ncpu)
        requested={'geometry':dict(geom),'condition':dict(cond),'setting':dict(setting)}
        binpath=binary_path or d['runtime']['vspaero_binary']
        runid,rundir=db.register_run(case_id,root,binpath,__file__,requested,ncpu=effective_ncpu)
        inputdir=rundir/'input';raw=rundir/'raw';logs=rundir/'logs'
        for p in (inputdir,raw,logs):p.mkdir()
        geompath=Path(json.loads(geom['physical_model_json'])['source_file'])
        target=inputdir/'model.vsp3';shutil.copy2(geompath,target)
        job={'study_id':d['study_id'],'run_id':runid,'input_vsp3':str(target),
             'output_dir':str(raw),'condition':dict(cond),
             'setting':{**dict(setting),'requested':case_requested},
             'runtime':{'vspaero_binary':str(Path(binpath).resolve()),'ncpu':effective_ncpu,'timeout_seconds':timeout_seconds}}
        job_path=inputdir/'job.json';job_path.write_text(canonical(job),encoding='utf-8')
        if progress_callback:
            notify_progress(progress_callback, {'event':'RUN_ID','run_id':runid,
                'analysis_mode':setting['analysis_mode'],
                'ncpu':effective_ncpu})
        for filename in (target,job_path): db.register_artifact(root,d['study_id'],filename,run_id=runid,role='INPUT')
        db.transition(runid,'PRECHECKED',message='Input and 3.50.5 binary path checked; runtime version checked in worker')
        db.transition(runid,'RUNNING',message='Isolated OpenVSP subprocess started')
        cmd=[python_executable or sys.executable,'-m','openvsppy.vv.worker','--job',str(job_path)]
        env=os.environ.copy()
        repo=REPO_ROOT
        env['PYTHONPATH']=str(repo/'src')+os.pathsep+env.get('PYTHONPATH','')
        env['VSPAERO_EXE']=str(Path(binpath).resolve())
        stdout=logs/'stdout.log';stderr=logs/'stderr.log'
        return_code=None
        failure=None
        # Logs are written during execution, including when the worker times out.
        # On Windows, the worker is one isolated process per Run. Terminate the
        # process tree on interruption so retries cannot overlap a stale solver.
        with stdout.open('wb') as out, stderr.open('wb') as err:
            try:
                proc=subprocess.Popen(cmd,cwd=rundir,env=env,stdout=out,stderr=err,
                                      creationflags=getattr(subprocess,'CREATE_NEW_PROCESS_GROUP',0) if os.name=='nt' else 0)
            except OSError as exc:
                failure='SPAWN_FAILED'
                err.write(str(exc).encode('utf-8',errors='replace'))
            else:
                try:
                    return_code=wait_for_worker(proc, cmd, timeout_seconds,
                        raw/'vspaero_api_stdout.log', progress_callback=progress_callback,
                        verbose=verbose)
                except (subprocess.TimeoutExpired,KeyboardInterrupt) as exc:
                    failure='TIMEOUT' if isinstance(exc,subprocess.TimeoutExpired) else 'OPERATOR_INTERRUPT'
                    if os.name=='nt':
                        subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],
                                       capture_output=True,check=False)
                    else:
                        proc.kill()
                    try:proc.wait(timeout=15)
                    except subprocess.TimeoutExpired:pass
        for log in (stdout,stderr):
            db.register_artifact(root,d['study_id'],log,run_id=runid,role='LOG')
        if failure:
            from .VVNative import diagnose_solver_log
            diagnosis=diagnose_solver_log(raw/'vspaero_api_stdout.log')
            (raw/'runtime_diagnostics.json').write_text(json.dumps({
                'reason':failure,'timeout_seconds':timeout_seconds,
                'native_diagnosis':diagnosis},indent=2),encoding='utf-8')
            if progress_callback:
                notify_progress(progress_callback,{'event':'RUNTIME_FAILURE',
                    'reason':failure,'last_stage':diagnosis['last_stage']})
            # Partial outputs are evidence of an interrupted attempt, not
            # successful numerical results.
            for artifact in raw.iterdir():
                if artifact.is_file():
                    db.register_artifact(root,d['study_id'],artifact,run_id=runid,role='PARTIAL')
            status='EXECUTION_FAILED' if failure=='SPAWN_FAILED' else 'INTERRUPTED'
            db.transition(runid,status,code=failure,message=failure)
            return {'run_id':runid,'status':status,'reason':failure,
                    'native_stage':diagnosis['last_stage']}
        if return_code:
            for artifact in raw.iterdir():
                if artifact.is_file():
                    db.register_artifact(root,d['study_id'],artifact,run_id=runid,role='PARTIAL')
            tail=stderr.read_text(encoding='utf-8',errors='replace')[-2000:]
            db.transition(runid,'EXECUTION_FAILED',code='WORKER_FAILED',message=tail,exit_code=return_code)
            diag_path=raw/'api_diagnostics.json'
            diagnostics=json.loads(diag_path.read_text(encoding='utf-8')) if diag_path.is_file() else {}
            from .VVNative import diagnose_solver_log, classify_worker_error
            failure_diag=diagnostics.get('native_failure') or diagnose_solver_log(raw/'vspaero_api_stdout.log')
            classification=classify_worker_error(tail,failure_diag)
            return {'run_id':runid,'status':'EXECUTION_FAILED',
                    'reason':classification+' ('+failure_diag['last_stage']+')',
                    'stderr':tail[-500:],'api_state':diagnostics.get('state'),
                    'native_failure':failure_diag,
                    'native_outputs':diagnostics.get('native_outputs_after_sweep',{})}
        result_path=raw/'worker_result.json'
        if not result_path.is_file():
            db.transition(runid,'EXECUTION_FAILED',code='NO_RESULT',message='Missing worker_result.json',exit_code=return_code)
            return {'run_id':runid,'status':'EXECUTION_FAILED','reason':'Missing worker_result.json'}
        try:
            result=json.loads(result_path.read_text(encoding='utf-8'))
        except (ValueError,OSError) as exc:
            db.transition(runid,'EXECUTION_FAILED',code='INVALID_MANIFEST',message=str(exc),exit_code=return_code)
            return {'run_id':runid,'status':'EXECUTION_FAILED','reason':'Invalid worker result manifest'}
        if result.get('openvsp_version')!='3.50.5':
            db.transition(runid,'EXECUTION_FAILED',code='VERSION_MISMATCH',message='Worker version did not match 3.50.5',exit_code=return_code)
            return {'run_id':runid,'status':'EXECUTION_FAILED'}
        db.transition(runid,'CAPTURED',effective=result.get('effective',{}),generated=result.get('generated',{}),exit_code=return_code)
        try:
            mode=setting['analysis_mode']
            with db.transaction():
                for p in raw.iterdir():
                    if p.is_file():db.register_artifact(root,d['study_id'],p,run_id=runid,role='NATIVE_OR_MANIFEST')
                suffixes=required_native_files(mode)
                files=[raw/('model_effective'+suffix) for suffix in suffixes]
                if any(not file.is_file() or not file.stat().st_size for file in files):
                    raise ValueError('Missing required native files: '+str([f.name for f in files if not f.is_file() or not f.stat().st_size]))
                if progress_callback:
                    notify_progress(progress_callback, {'event':'OUTPUT','files':
                        {p.name:p.stat().st_size for p in files}})
                if mode=='STEADY_SWEEP':
                    imported=import_polar(db,study_root=root,run_id=runid,path=files[0],
                        tol_deg=d.get('policies',{}).get('reference_angle_tolerance_deg',1e-6))
                    must_have={'CLi','CDi','CLiw','CDiw','CMiy'}
                    native_names={r['qoi'] for r in db.conn.execute(
                        "SELECT qoi FROM quantity WHERE dataset_id=? AND value_status='AVAILABLE'",
                        (imported['dataset_id'],))}
                    if not must_have.issubset(native_names):
                        raise ValueError('Missing required native polar QoIs: '+str(sorted(must_have-native_names)))
                    if not imported['matched']:raise ValueError('Polar condition does not match planned Case')
                else:
                    requested=json.loads(setting['requested_json'])
                    mapping={}
                    if mode=='STEADY':
                        from ..vspaero.VSPAEROStab import read_vspaero_stab
                        native=read_vspaero_stab(files[0])
                        approved=set(requested.get('control_group_names',{}).values())
                        if not approved or set(native.control_groups.values())!=approved:
                            raise ValueError('Native control group labels do not match reviewed groups')
                        mapping=native.control_groups.copy()  # Names observed in file, checked against approved names.
                    imported=import_stab(db,study_root=root,run_id=runid,path=files[0],
                        mode=mode,control_mapping=mapping,
                        method='STAB_NATIVE' if mode=='STEADY' else mode+'_NATIVE')
                db.transition(runid,'PARSED',message='Mode-specific native results parsed and conditions matched')
                db.transition(runid,'SUCCEEDED',message='Parsing and integrity completed; V&V not yet accepted')
                db.conn.execute('INSERT INTO case_run_link VALUES(?,?,?,?)',(
                    case_id,imported['point_id'],
                    __import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),
                    'ONE_POINT_MATCH'))
        except Exception as exc:
            db.transition(runid,'OUTPUT_INVALID',code='OUTPUT_INVALID',message=str(exc))
            return {'run_id':runid,'status':'OUTPUT_INVALID','reason':str(exc)}
        if progress_callback:
            notify_progress(progress_callback, {'event':'DATABASE','run_id':runid,
                'point_id':imported['point_id'],'status':'SUCCEEDED'})
        return {'run_id':runid,'status':'SUCCEEDED','point_id':imported['point_id']}
    finally:
        db.close()


def select_case_timeout(case_id, mode, default_seconds, by_mode=None, by_case=None):
    """Use only explicitly configured per-mode/per-case runtime limits."""
    limits = dict(by_mode or {})
    limits.update(by_case or {})
    value = (by_case or {}).get(case_id, (by_mode or {}).get(mode, default_seconds))
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'Invalid timeout for {case_id} ({mode}): {value!r}')
    return value


def run_cases(db_path, study_file, case_ids, *, binary_path=None, python_executable=None,
              timeout_seconds=3600, resume=True, progress_path=None,
              progress_callback=None, verbose=False, ncpu=None,
              timeouts_by_mode=None, timeouts_by_case=None):
    """Run resolved Cases sequentially; leave every successful raw Run immutable.

    Returns a record for each requested Case. A single failing Case does not
    prevent subsequent independent Cases from running. Never reuses a Run with
    changed solver/execution code or missing/modified raw artifacts.
    """
    from .VVDatabase import execution_fingerprint, postprocess_fingerprint
    db_path=Path(db_path)
    study_file=Path(study_file).resolve()
    study_root=study_file.parent
    _,study_config,_,_=load_study(study_file)
    binary_path=binary_path or study_config.get('runtime',{}).get('vspaero_binary')
    binary=Path(binary_path).resolve() if binary_path else None
    binary_hash=digest(binary) if binary and binary.is_file() else None
    code_hash=execution_fingerprint(REPO_ROOT)
    process_hash=postprocess_fingerprint(REPO_ROOT)
    results=[]
    total=len(case_ids)
    for index,case_id in enumerate(case_ids,1):
        started=time.monotonic()
        if progress_callback:
            notify_progress(progress_callback, {'phase':'RUN','event':'CHECKING','case_id':case_id,
                               'index':index,'total':total})
        db=VVDatabase(db_path)
        try:
            case=db.row('planned_case','case_id',case_id)
            if not case:
                result={'case_id':case_id,'status':'BLOCKED','reason':'Unknown Case ID'}
            elif case['execute_policy']!='EXECUTE' or None in (case['geometry_id'],case['condition_id'],case['setting_id']):
                result={'case_id':case_id,'status':'BLOCKED','reason':'Case is unresolved or is not an executable analysis'}
            elif binary_hash is None:
                result={'case_id':case_id,'status':'BLOCKED','reason':'No VSPAERO binary was supplied or the path is missing'}
            else:
                setting=db.row('numerical_setting','setting_id',case['setting_id'])
                effective_ncpu=resolve_ncpu(json.loads(setting['requested_json']),ncpu)
                case_timeout=select_case_timeout(case_id,setting['analysis_mode'],
                    timeout_seconds,timeouts_by_mode,timeouts_by_case)
                previous=db.conn.execute('''SELECT r.run_id,r.binary_sha256,r.environment_json
                    FROM run r JOIN case_point pt ON pt.run_id=r.run_id
                    JOIN case_run_link l ON l.point_id=pt.point_id
                    WHERE l.case_id=? AND r.status='SUCCEEDED' AND pt.point_status='MATCHED'
                      AND r.geometry_id=? AND r.condition_id=? AND r.setting_id=?
                    ORDER BY r.rowid DESC''',
                    (case_id,case['geometry_id'],case['condition_id'],case['setting_id'])).fetchall()
                result=None
                if resume and not check_case(db,case_id,study_file,binary):
                    for old in previous:
                        env=json.loads(old['environment_json'])
                        if (old['binary_sha256']!=binary_hash or
                            env.get('execution_source_sha256')!=code_hash or
                            env.get('ncpu')!=effective_ncpu):
                            continue
                        artifacts=list(db.conn.execute(
                            'SELECT relative_path,sha256 FROM artifact WHERE run_id=?',
                            (old['run_id'],)))
                        if not artifacts or not all((study_root/a['relative_path']).is_file()
                            and digest(study_root/a['relative_path'])==a['sha256'] for a in artifacts):
                            continue
                        if env.get('postprocess_source_sha256')!=process_hash:
                            update=reprocess_run(db,study_root=study_root,run_id=old['run_id'])
                            result={'case_id':case_id,'run_id':old['run_id'],**update}
                        else:
                            result={'case_id':case_id,'run_id':old['run_id'],'status':'REUSED','ncpu':effective_ncpu}
                        break
                if result is None:
                    # Release SQLite before the worker creates its own connection.
                    db.close();db=None
                    if progress_callback:
                        notify_progress(progress_callback, {'phase':'RUN','event':'STARTED','case_id':case_id,
                                           'index':index,'total':total,
                                           'timeout_seconds':case_timeout})
                    def report_native(entry):
                        if progress_callback:
                            notify_progress(progress_callback, {'phase':'RUN','case_id':case_id,
                                'index':index,'total':total, **entry})
                    result=run_case(db_path,case_id,study_file,
                                    python_executable=python_executable,
                                    binary_path=str(binary),timeout_seconds=case_timeout,
                                    progress_callback=report_native, verbose=verbose, ncpu=effective_ncpu)
                    result['case_id']=case_id
                    result['ncpu']=effective_ncpu
        except (Exception) as exc:
            result={'case_id':case_id,'status':'BLOCKED','reason':f'{type(exc).__name__}: {exc}'}
        finally:
            if db is not None: db.close()
        results.append(result)
        if progress_callback:
            notify_progress(progress_callback, {'phase':'RUN','event':'FINISHED','index':index,'total':total,
                'elapsed_seconds':round(time.monotonic()-started, 1),**result})
        if progress_path:
            progress=Path(progress_path);progress.parent.mkdir(parents=True,exist_ok=True)
            with progress.open('a',encoding='utf-8') as file:
                file.write(json.dumps(result,ensure_ascii=False)+'\n')
    return results


def run_parasite_drag(db, study_root, *, model, velocity_mps, altitude_m,
                      delta_temp_k, sref, geom_set, metadata,
                      python_executable=None, timeout_seconds=3600, geometry_names=None):
    """Isolated reviewed 3.50.5 ParasiteDrag, with immutable native evidence."""
    import csv
    import math
    import uuid
    from .VVPostprocess import import_reference_csv
    root=Path(study_root).resolve()
    model=Path(model).resolve()
    if not model.is_file():raise FileNotFoundError(model)
    if not math.isfinite(float(velocity_mps)) or float(velocity_mps)<=0 or float(sref)<=0:
        raise ValueError('Positive velocity and reference area required')
    if metadata.get('kind')!='EXTERNAL_CD0' or not metadata.get('reference_reviewed'):
        raise ValueError('Reviewed EXTERNAL_CD0 metadata is required')
    existing=[row for row in db.conn.execute(
        "SELECT dataset_id,metadata_json FROM dataset WHERE kind='EXTERNAL_CD0'")
        if json.loads(row['metadata_json']).get('source_id')==metadata.get('source_id')]
    if existing:
        if len(existing)!=1:raise ValueError('Ambiguous existing ParasiteDrag source_id')
        dsid=existing[0]['dataset_id']
        candidates=[root/a['relative_path'] for a in db.conn.execute(
            "SELECT relative_path FROM artifact WHERE dataset_id=? AND role='PARASITE_RAW'",(dsid,))
            if Path(a['relative_path']).name=='manifest.json']
        if len(candidates)!=1 or not candidates[0].is_file():
            raise ValueError('ParasiteDrag source_id already exists but no reusable native manifest')
        old=json.loads(candidates[0].read_text(encoding='utf-8'))
        required={'velocity_mps':velocity_mps,'altitude_m':altitude_m,
                  'delta_temp_k':delta_temp_k,'geom_set':geom_set,'geometry_names':geometry_names}
        if (old.get('model_sha256')!=digest(model) or
            any((old.get(k)!=v if k=='geometry_names' else
                 not math.isclose(float(old[k]),float(v),rel_tol=1e-9,abs_tol=1e-9))
                for k,v in required.items()) or
            not math.isclose(float(old['result']['Sref']),float(sref),rel_tol=1e-9)):
            raise ValueError('ParasiteDrag source_id was reused with changed physical inputs')
        return {'status':'REUSED_PARASITE','dataset_id':dsid,
                'result':old['result'],'native':str(candidates[0].parent)}
    folder=root/'02_reference'/'external_cd0'/('parasite_'+uuid.uuid4().hex[:16])
    folder.mkdir(parents=True,exist_ok=False)
    job={'model':str(model),'velocity_mps':float(velocity_mps),
         'altitude_m':float(altitude_m),'delta_temp_k':float(delta_temp_k),
         'sref':float(sref),'geom_set':int(geom_set),
         'geometry_names':geometry_names,
         'output':str(folder/'parasite_native.csv'),
         'result_manifest':str(folder/'manifest.json')}
    (folder/'job.json').write_text(canonical(job),encoding='utf-8')
    repo=REPO_ROOT
    env=os.environ.copy();env['PYTHONPATH']=str(repo/'src')+os.pathsep+env.get('PYTHONPATH','')
    with (folder/'stdout.log').open('wb') as out,(folder/'stderr.log').open('wb') as err:
        proc=subprocess.run([python_executable or sys.executable,'-m','openvsppy.vv.worker',
                             '--drag-job',str(folder/'job.json')],cwd=folder,
                             env=env,stdout=out,stderr=err,timeout=timeout_seconds)
    if proc.returncode:
        raise RuntimeError(f'ParasiteDrag worker failed; native evidence preserved: {folder}')
    native=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
    result=native['result']
    if not math.isclose(result['Mach'],float(metadata['mach']),rel_tol=1e-4,abs_tol=1e-6):
        raise ValueError('ParasiteDrag Mach differs from reviewed reference metadata: '+str(folder))
    if not math.isclose(result['Sref'],float(metadata['sref']),rel_tol=1e-7):
        raise ValueError('ParasiteDrag Sref differs from reviewed reference metadata: '+str(folder))
    ref_csv=folder/'cd0_import.csv'
    with ref_csv.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['qoi','value','unit'])
        writer.writerow(['CD0_build',result['CD0'],'1'])
    metadata={**metadata,'source_sha256':digest(ref_csv),
             'parasite_native_sha256':digest(folder/'parasite_native.csv'),
             'parasite_manifest_sha256':digest(folder/'manifest.json'),
             'model_sha256':native['model_sha256'],
             're_cref_source':'REVIEWED_EXTERNAL_METADATA'}
    study=db.conn.execute('SELECT study_id FROM study LIMIT 1').fetchone()
    if not study:raise ValueError('Study not initialized')
    imported=import_reference_csv(db,study_id=study['study_id'],csv_path=ref_csv,metadata=metadata)
    for artifact in folder.iterdir():
        if artifact.is_file():
            db.register_artifact(root,study['study_id'],artifact,
                                 dataset_id=imported['dataset_id'],role='PARASITE_RAW')
    return {'dataset_id':imported['dataset_id'],'result':result,'native':str(folder)}
