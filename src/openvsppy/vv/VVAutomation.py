"""Prepare reviewed V&V Cases and execute the study through report generation.

The 25 Test IDs are defined by the plan and SQLite catalog. Only experiments
with explicitly reviewed inputs can run; unknown inputs remain blockers.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

from .VVDatabase import digest, CASE_ROLES, REPO_ROOT


# Initial candidate levels, not accepted grid convergence settings.
GRID_LEVELS = {'COARSE': 21, 'MEDIUM': 31, 'FINE': 49}
WAKE_STUDIES = {
    'WAKE-01': [('ITER_5', 5, 8), ('ITER_10', 10, 8), ('ITER_20', 20, 8)],
    'WAKE-02': [('NODES_8', None, 8), ('NODES_12', None, 12), ('NODES_16', None, 16)],
}


def prepare_cases(test_id, base_file, study_root, *, dry_run=False, levels=None, span=None):
    """Generate explicit Case JSON and optionally .vsp3 grid levels.

    The output directory is dedicated to the test. Existing differing files
    are never silently overwritten. A repeated identical prepare is safe.
    """
    source=Path(base_file).resolve()
    base=json.loads(source.read_text(encoding='utf-8-sig'))
    root=Path(study_root).resolve()
    output=root/'generated'/test_id
    case_specs=[]
    grid=test_id in ('GRID-01','GRID-02')
    if grid:
        levels=levels or GRID_LEVELS
        case_specs=[(level,None,None) for level in ('COARSE','MEDIUM','FINE')]
    elif test_id=='GRID-03':
        # BASE, LOCAL_A, LOCAL_B isolate LECluster; NREF_1 changes only NRef.
        case_specs=[(role,None,None) for role in (*CASE_ROLES[test_id],'NREF_1')]
    elif test_id in ('VER-01','VER-02'):
        case_specs=[(role,None,None) for role in CASE_ROLES[test_id]]
    elif test_id in WAKE_STUDIES:
        case_specs=WAKE_STUDIES[test_id]
    elif test_id in ('GEO-01','GEO-02','GEO-03','GEO-04'):
        case_specs=[(role,None,None) for role in CASE_ROLES[test_id]]
    elif test_id=='WAKE-03':
        case_specs=[(role,None,None) for role in CASE_ROLES[test_id]]
    elif test_id=='STA-03':
        case_specs=[(role,None,None) for role in CASE_ROLES[test_id]]
    elif test_id=='STA-02':
        case_specs=[(role,None,None) for role in CASE_ROLES[test_id]]
    elif test_id=='STA-01':
        case_specs=[('STEADY',None,None),('ADJOINT',None,None)]
    elif test_id in ('CTRL-01','CTRL-02','CTRL-03'):
        group={'CTRL-01':'elevator','CTRL-02':'aileron','CTRL-03':'rudder'}[test_id]
        req=base['numerical_setting']['requested']
        if not req.get('control_sign_verified') or group not in req.get('control_group_names',{}):
            raise ValueError(f'{test_id}: approved control group and physical sign mapping required in baseline')
        case_specs=[(f'DELTA_{angle:+d}',None,None) for angle in (-10,-5,-2,0,2,5,10)]
    elif test_id in ('VAL-01','VAL-02'):
        if span is None or not 0<float(span)<=20:
            raise ValueError('VAL requires --span from an approved flight condition study (0 < span <= 20 deg)')
        label='ALPHA' if test_id=='VAL-01' else 'BETA'
        case_specs=[(f'{label}_NEG',None,None),(f'{label}_0',None,None),(f'{label}_POS',None,None)]
    else:
        raise ValueError(f'{test_id}: no approved automatic input-generation recipe; use reviewed Case JSON.')
    prepared=[]
    for role,iterations,nodes in case_specs:
        case=copy.deepcopy(base)
        case['case_id']=f'{test_id}-{role}'
        geom=case['geometry'];setting=case['numerical_setting']
        if test_id.startswith('VER-'):
            benchmark=base.get('benchmark')
            if not isinstance(benchmark,dict) or benchmark.get('test_id')!=test_id:
                raise ValueError(f'{test_id}: reviewed benchmark model and method must be supplied')
            model_file=Path(benchmark['vsp3_source'])
            model_file=(source.parent/model_file).resolve() if not model_file.is_absolute() else model_file.resolve()
            if not model_file.is_file():raise FileNotFoundError(model_file)
            geom['vsp3_source']=str(model_file)
            geom['geometry_id']=f'{test_id}-{role}-G-{digest(model_file)[:12]}'
            geom['physical_model']['purpose']=f'{test_id}_BENCHMARK_UNVERIFIED'
            geom['representation']['case_role']=role
            setting['requested']['geometry_membership']=benchmark['geometry_membership']
            setting['requested']['geometry_membership_reviewed']=True
            setting['requested']['thin_set_index']=int(benchmark['thin_set_index'])
            setting['requested']['thick_set_index']=int(benchmark['thick_set_index'])
            setting['mesh_level']=role
            if test_id=='VER-01':
                # The baseline has already supplied benchmark-specific reference dimensions.
                setting['requested']['geom_parms']=[{'geometry':benchmark['wing_name'],
                    'name':'Tess_W','group':'Shape',
                    'value':{'COARSE':21,'MEDIUM':31,'FINE':49}[role]}]
            else:
                # Windows 3.50.5 native logs: NRef=1/2 on the reference
                # ellipsoid produced 2/1 unassociated loops and no .polar.
                # Resolve geometry convergence with explicit surface tessellation
                # instead.  40/41 is the original, successfully solved level.
                tess={'COARSE':(20,21),'MEDIUM':(30,31),'FINE':(40,41)}[role]
                setting['requested']['geom_parms']=[
                    {'geometry':'EllipsoidGeom','name':'Tess_U','group':'Shape','value':tess[0]},
                    {'geometry':'EllipsoidGeom','name':'Tess_W','group':'Shape','value':tess[1]}]
                setting['requested']['vspaero']['NRef']=0
        elif test_id=='GRID-03':
            wing={'BASE':0.5,'LOCAL_A':0.25,'LOCAL_B':1.0,'NREF_1':0.5}[role]
            setting['requested']['geom_parms']=[{'geometry':'WingGeom','name':'LECluster',
                                               'group':'WingGeom','value':wing}]
            setting['requested']['vspaero']['NRef']=1 if role=='NREF_1' else 0
            geom['geometry_id']=f'GRID-03-{role}-G-{digest(geom["vsp3_source"])[:12]}'
            geom['physical_model']['purpose']='GRID03_LE_CLUSTER_NREF_SENSITIVITY'
            geom['representation']['case_role']=role
        elif grid:
            model=output/f'{test_id}_{role}.vsp3'
            origin=Path(geom['vsp3_source']).resolve()
            if not origin.is_file():raise FileNotFoundError(origin)
            geom['geometry_id']=f'{test_id}-{role}-G-{digest(origin)[:12]}'
            geom['vsp3_source']=str(model)
            geom['remarks']=f'Grid study: baseline_sha256={digest(origin)}; reference_tess_w={levels[role]}'
            geom.setdefault('physical_model',{})['purpose']=f'{test_id}_MESH_SENSITIVITY'
            geom['physical_model']['aircraft_variant']='G103A_TWIN_II_ACRO'
            setting['mesh_level']=role
        elif test_id.startswith('GEO-'):
            configurations={
              'G1':{'thin':['WingGeom'],'thick':[]},
              'G2':{'thin':['WingGeom'],'thick':['FuselageGeom']},
              'G3':{'thin':['WingGeom','HTailGeom','VTailGeom'],'thick':['FuselageGeom']},
              # GEO-04: reference M0 vs wing-only Thick M1, never all-geometry Thick.
              'THIN':{'thin':['WingGeom','HTailGeom','VTailGeom'],'thick':['FuselageGeom']},
              'THICK':{'thin':['HTailGeom','VTailGeom'],'thick':['FuselageGeom','WingGeom']},
            }
            if not setting['requested'].get('geometry_membership_reviewed'):
                raise ValueError('GEO requires explicit reviewed geometry name/set mapping in baseline')
            setting['requested']['geometry_membership']=configurations[role]
            setting['requested']['thin_set_index']=(-1 if not configurations[role]['thin'] else base['numerical_setting']['requested']['thin_set_index'])
            setting['requested']['thick_set_index']=(-1 if not configurations[role]['thick'] else base['numerical_setting']['requested']['thick_set_index'])
            geom['geometry_id']=f'{test_id}-{role}-G'
            geom['physical_model']['purpose']=f'GEO_{role}_WING_THICK_COMPARISON' if test_id=='GEO-04' else f'GEO_{role}_REVIEWED'
            geom['representation']['case_role']=role
        elif test_id=='WAKE-03':
            setting['requested'].setdefault('vspaero',{})['ImplicitWake']=role=='IMPLICIT'
            if 'WakeNumIter' not in setting['requested']['vspaero']:
                raise ValueError('WAKE-03 needs a reviewed common wake iteration count')
        elif test_id=='STA-03':
            if not setting['requested'].get('unsteady_settings_reviewed'):
                raise ValueError('STA-03 needs reviewed Vinf/Rho and native P/Q/R solver settings')
            setting['analysis_mode']=role
            # Native 3.50.5 solver CLI uses its own time defaults for P/Q/R.
            # Do not claim a configurable time step unless its input is verified.
        elif test_id=='STA-02':
            if not setting['requested'].get('stability_sensitivity_reviewed'):
                raise ValueError('STA-02 requires reviewed Mesh/Wake sensitivity levels')
            setting['analysis_mode']='ADJOINT'
            variants=setting['requested']['stability_sensitivity_reviewed']
            if role not in variants:raise ValueError('Missing STA-02 settings for '+role)
            if role=='MESH_FINE':
                alternate=variants[role].get('vsp3_source')
                if not alternate or not Path(alternate).is_file():
                    raise ValueError('STA-02 MESH_FINE requires an existing reviewed refined .vsp3')
                geom['vsp3_source']=str(Path(alternate).resolve())
                geom['geometry_id']=f'STA-02-{role}-G-{digest(alternate)[:12]}'
            setting['requested'].setdefault('vspaero',{}).update({k:v for k,v in variants[role].items() if k!='vsp3_source'})
        elif test_id in WAKE_STUDIES:
            if iterations is None:
                iterations=setting['requested'].get('vspaero',{}).get('WakeNumIter')
                if iterations is None:raise ValueError('WAKE-02 requires reviewed baseline WakeNumIter')
            setting['requested'].setdefault('vspaero',{})['WakeNumIter']=iterations
            setting['requested']['vspaero']['RootWakeNodes']=nodes
            setting['requested']['vspaero']['FixedWakeFlag']=False
        elif test_id=='STA-01':
            setting['analysis_mode']=role
        elif test_id.startswith('CTRL-'):
            control={'CTRL-01':'elevator','CTRL-02':'aileron','CTRL-03':'rudder'}[test_id]
            angle=int(role.split('_')[1])
            case['flight_condition']['controls_deg'][control]=float(angle)
            case['flight_condition']['condition_id']=f'{test_id}-{role}-FC'
        elif test_id.startswith('VAL-'):
            axis='alpha_deg' if test_id=='VAL-01' else 'beta_deg'
            increment=-float(span) if role.endswith('NEG') else (float(span) if role.endswith('POS') else 0.0)
            case['flight_condition'][axis]=float(base['flight_condition'][axis])+increment
            case['flight_condition']['condition_id']=f'{test_id}-{role}-FC'
        # Unique IDs prevent previously registered ID mutation.
        if test_id=='GRID-03' and role=='NREF_1':
            case['case_id']='GRID-03-NREF_1'
        setting['setting_id']=f'{test_id}-{role}-S'
        path=output/f'{case["case_id"]}.json'
        prepared.append({'case_id':case['case_id'],'case_json':str(path),
                         'vsp3_source':geom['vsp3_source'],
                         'wake':setting['requested'].get('vspaero',{}),
                         'mesh_level':setting['mesh_level']})
        if dry_run:continue
        output.mkdir(parents=True,exist_ok=True)
        if grid:
            if not model.is_file():
                # Existing surface tessellation algorithm is the source of
                # truth. It is not a proof of volume/grid quality.
                from ..vspaero.VSPAEROMesh import equalize_vspaero_tessellation
                metrics=equalize_vspaero_tessellation(
                    input_vsp3_path=origin,output_vsp3_path=model,
                    reference_wing_name='WingGeom',reference_tess_w=levels[role],
                    geom_names=['WingGeom','HTailGeom','VTailGeom','FuselageGeom'],
                    tolerance=0.10,max_iterations=4)
                metrics.to_csv(output/f'{test_id}_{role}_mesh_metrics.csv',index=False)
            geom['remarks']+=f'; generated_sha256={digest(model)}'
        content=json.dumps(case,indent=2,ensure_ascii=False)+'\n'
        if path.is_file() and path.read_text(encoding='utf-8')!=content:
            raise ValueError(f'Prepared Case changed; choose a new output directory or version: {path}')
        if not path.is_file():path.write_text(content,encoding='utf-8')
    return prepared


def case_difference(db, test_id):
    """Check immutable comparison inputs for supported series.

    GRID varies geometry/setting but holds flight condition and wake settings.
    WAKE varies only requested wake parameters, holding all other input IDs.
    """
    rows=list(db.conn.execute('SELECT * FROM planned_case WHERE test_id=? ORDER BY case_id',(test_id,)))
    if test_id not in ('GRID-01','GRID-02','WAKE-01','WAKE-02','WAKE-03'):
        return ['No automated comparison equivalence rule for this Test ID']
    if any(None in (r['geometry_id'],r['condition_id'],r['setting_id']) for r in rows):
        return ['One or more Cases remain unresolved']
    conditions={r['condition_id'] for r in rows}
    if len(conditions)!=1:return ['Flight conditions differ across Cases']
    requested=[json.loads(db.row('numerical_setting','setting_id',r['setting_id'])['requested_json']) for r in rows]
    errors=[]
    if test_id.startswith('WAKE'):
        if len({r['geometry_id'] for r in rows})!=1:errors.append('WAKE geometry differs')
        variation={'WAKE-01':'WakeNumIter','WAKE-02':'RootWakeNodes','WAKE-03':'ImplicitWake'}[test_id]
        norm=[]
        for s in requested:
            s.get('vspaero',{}).pop(variation,None)
            norm.append(s)
        if any(s!=norm[0] for s in norm[1:]):errors.append('WAKE has unrelated setting differences')
    else:
        # Representation changes, not flight/solver settings, are intended.
        for s in requested:
            s.pop('mesh_overrides',None)
        if any(s!=requested[0] for s in requested[1:]):errors.append('GRID has unrelated solver setting differences')
    return errors


def execute_study(study_file, manifest_file, *, db_path=None, binary_path=None,
                  python_executable=None, timeout_seconds=3600, resume=True,
                  execute=True, verbose=False, ncpu=None):
    """One readable entry point: prepare, run, collect, plot, and report.

    The manifest supplies *approved* recipes and/or reviewed Case documents.
    Unspecified inputs stay unresolved. Successful raw Runs are never edited.
    Every invocation writes a new, immutable-by-convention reviewer directory.
    """
    import csv
    import shutil
    import uuid
    from datetime import datetime, timezone
    from .VVDatabase import VVDatabase, load_study
    from .VVRunner import run_cases, run_parasite_drag
    from .VVPostprocess import (collect_study_results, evaluate_control_sweep,
                                richardson_gci, import_reference_csv,
                                compare_reference_dataset, create_drag_dataset)
    from .VVPlots import create_study_plots
    from .VVReport import write_evaluation_report, export_csv

    study_path,study,_,_=load_study(study_file)
    root=study_path.parent
    manifest_path=Path(manifest_file).resolve()
    config=json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    db_path=Path(db_path) if db_path else root/'01_database'/'vv.sqlite'
    # Directories from prior executions remain untouched, even on rerun.
    name=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
    destination=root/'04_reports'/('execute_'+name)
    destination.mkdir(parents=True,exist_ok=False)
    progress=destination/'progress.jsonl'
    def display(entry):
        phase, event = entry.get('phase'), entry.get('event')
        message = None
        if phase == 'RUN' and entry.get('case_id'):
            case = entry.get('case_id','')
            label = f"{entry.get('index', '?')}/{entry.get('total', '?')} {case}"
            if event == 'CHECKING':
                message = f'[CASE] {label} CHECKING'
            elif event == 'STARTED':
                message = f'[RUN] {label} STARTED (timeout={entry.get("timeout_seconds", "default")}s)'
            elif event == 'FINISHED':
                message = f'[RUN] {label} {entry["status"]} ({entry["elapsed_seconds"]:.1f}s)'
                if entry['status'] not in ('SUCCEEDED','REUSED','REPROCESSED','REUSED_POSTPROCESS'):
                    message += f" - {entry.get('reason',entry.get('stderr',''))[:200]}"
            elif verbose and event == 'ADJOINT_STAGE':
                message = f'[ADJOINT] {case} stage: {entry["stage"]}'
            elif verbose and event == 'ADJOINT_GMRES':
                message = f'[ADJOINT] {case} {entry["stage"]} GMRES {entry["iteration"]}'
            elif verbose and event == 'RUNTIME_FAILURE':
                message = f'[DIAGNOSTIC] {case} {entry["reason"]}: {entry["last_stage"]}'
            elif event == 'HEARTBEAT':
                message = f'[RUN] {label} running ({entry["elapsed_seconds"]:.0f}s elapsed)'
            elif verbose and event == 'RUN_ID':
                message = f'[RUN] {case} Run ID: {entry["run_id"]}, mode: {entry["analysis_mode"]}'
                if entry.get('ncpu') is not None:
                    message += f', NCPU: {entry["ncpu"]}'
            elif verbose and event == 'WAKE':
                message = f'[VSPAERO] {case} Wake {entry["iteration"]}/{entry["total_iterations"]} started'
            elif verbose and event == 'OUTPUT':
                message = f'[OUTPUT] {case} ' + ', '.join(f'{name} ({size} B)' for name,size in entry['files'].items())
            elif verbose and event == 'DATABASE':
                message = f'[DATABASE] {case} {entry["status"]} (Point {entry["point_id"]})'
        elif event in ('STARTED','FINISHED'):
            message = f'[{phase}] {entry.get("message",event)}'
        if message:
            try:
                print(message,file=sys.stderr,flush=True)
            except (OSError,UnicodeError):
                pass
    def path_of(value):
        p=Path(value)
        return (manifest_path.parent/p).resolve() if not p.is_absolute() else p.resolve()
    def log(entry):
        with progress.open('a',encoding='utf-8') as stream:
            stream.write(json.dumps(entry,ensure_ascii=False)+'\n')
        display(entry)

    log({'phase':'WORKFLOW','event':'STARTED',
         'message':f'execute-all started; progress: {progress}'})
    log({'phase':'PREPARE','event':'STARTED','message':'Preparing reviewed Cases'})
    db=VVDatabase(db_path)
    prepared=[]
    try:
        if db.row('study','study_id',study['study_id']) is None:
            raise ValueError('Study DB is not initialized. Run: python -m openvsppy.vv init')
        # 1. Register explicitly reviewed Case JSON files; never modify resolved Case IDs.
        for source in config.get('case_files',[]):
            try:
                definition=json.loads(path_of(source).read_text(encoding='utf-8-sig'))
                cid=definition['case_id']
                existing=db.row('planned_case','case_id',cid)
                if existing is None:raise ValueError('Not in the plan: '+cid)
                if existing['geometry_id'] is not None:
                    # Validate immutable content too, not only the three IDs.
                    db.add_geometry(study['study_id'],definition['geometry'])
                    db.add_condition(study['study_id'],definition['flight_condition'])
                    db.add_setting(study['study_id'],definition['numerical_setting'])
                    expected=(definition['geometry']['geometry_id'],
                              definition['flight_condition']['condition_id'],
                              definition['numerical_setting']['setting_id'])
                    actual=(existing['geometry_id'],existing['condition_id'],existing['setting_id'])
                    if expected!=actual:raise ValueError('Resolved Case ID differs from provided Case: '+cid)
                    status='ALREADY_RESOLVED'
                else:
                    db.register_case_document(study['study_id'],definition)
                    status='REGISTERED'
                prepared.append({'case_id':cid,'status':status})
            except Exception as exc:
                prepared.append({'file':source,'status':'BLOCKED','reason':str(exc)})
            log({'phase':'PREPARE','result':prepared[-1]})
        # 1b. Create the two official benchmark geometries with the same 3.50.5
        # Python API as the analysis worker. Only new versioned files are written.
        for spec in config.get('benchmark_models',[]):
            target=path_of(spec['output'])
            manifest=target.with_suffix('.generation.json')
            record={'kind':spec['kind'],'output':str(target),'result_manifest':str(manifest)}
            try:
                if target.exists():
                    if not manifest.is_file() or json.loads(manifest.read_text(
                        encoding='utf-8'))['sha256']!=digest(target):
                        raise ValueError('Benchmark exists without matching 3.50.5 creation evidence')
                    status='REUSED_BENCHMARK'
                elif not execute:
                    raise RuntimeError('Benchmark model absent; no-run mode cannot generate it')
                else:
                    import os,subprocess
                    target.parent.mkdir(parents=True,exist_ok=True)
                    job_file=target.with_suffix('.job.json')
                    job_file.write_text(json.dumps(record,indent=2),encoding='utf-8')
                    env=os.environ.copy()
                    repo=REPO_ROOT
                    env['PYTHONPATH']=str(repo/'src')+os.pathsep+env.get('PYTHONPATH','')
                    with target.with_suffix('.stdout.log').open('wb') as stdout, \
                         target.with_suffix('.stderr.log').open('wb') as stderr:
                        completed=subprocess.run([python_executable or sys.executable,
                            '-m','openvsppy.vv.worker','--benchmark-job',str(job_file)],
                            cwd=repo,env=env,stdout=stdout,stderr=stderr,
                            timeout=timeout_seconds)
                    if completed.returncode or not target.is_file():
                        raise RuntimeError('Benchmark generation failed; inspect '+str(stderr.name))
                    status='GENERATED_BENCHMARK'
                prepared.append({'test_id':spec['test_id'],'status':status,
                                 'benchmark_vsp3':str(target)})
            except Exception as exc:
                prepared.append({'test_id':spec['test_id'],'status':'BLOCKED',
                                 'reason':str(exc)})
            log({'phase':'BENCHMARK','result':prepared[-1]})
        # 2. Generate only the recipes the user explicitly configured.
        for test_id,source in config.get('baselines',{}).items():
            if test_id not in CASE_ROLES:
                prepared.append({'test_id':test_id,'status':'BLOCKED','reason':'Unknown Test ID'})
                continue
            try:
                cases=prepare_cases(test_id,path_of(source),root,
                                    span=config.get('spans_deg',{}).get(test_id),
                                    levels=config.get('grid_levels',{}).get(test_id))
                for item in cases:
                    if item['case_id']=='GRID-03-NREF_1' and db.row('planned_case','case_id',item['case_id']) is None:
                        db.add_extra_planned_case(study['study_id'],'GRID-03','NREF_1')
                    definition=json.loads(Path(item['case_json']).read_text(encoding='utf-8'))
                    old=db.row('planned_case','case_id',item['case_id'])
                    if old['geometry_id'] is None:
                        db.register_case_document(study['study_id'],definition)
                        status='REGISTERED'
                    else:
                        db.add_geometry(study['study_id'],definition['geometry'])
                        db.add_condition(study['study_id'],definition['flight_condition'])
                        db.add_setting(study['study_id'],definition['numerical_setting'])
                        expected=(definition['geometry']['geometry_id'],
                                  definition['flight_condition']['condition_id'],
                                  definition['numerical_setting']['setting_id'])
                        current=(old['geometry_id'],old['condition_id'],old['setting_id'])
                        if current!=expected:raise ValueError('Resolved Case differs: '+item['case_id'])
                        status='ALREADY_RESOLVED'
                    prepared.append({'case_id':item['case_id'],'status':status})
            except Exception as exc:
                prepared.append({'test_id':test_id,'status':'BLOCKED','reason':str(exc)})
            log({'phase':'PREPARE','result':prepared[-1]})
        log({'phase':'PREPARE','event':'FINISHED',
             'message':f'Case preparation completed ({len(prepared)} records)'})
        # 3. Import referenced data when provided. Unprovided data remains missing.
        imports=[]
        for item in config.get('references',[]):
            try:
                source=path_of(item['file']);metadata=json.loads(path_of(item['metadata']).read_text(encoding='utf-8'))
                if not source.is_file():raise FileNotFoundError(source)
                folder=root/'02_reference'/'imported'
                folder.mkdir(parents=True,exist_ok=True)
                target=folder/(digest(source)[:12]+'_'+source.name)
                if not target.is_file():shutil.copy2(source,target)
                if digest(source)!=digest(target):raise ValueError('Reference content hash mismatch')
                metadata['source_sha256']=digest(target)
                prior=[row for row in db.conn.execute(
                    'SELECT dataset_id,metadata_json FROM dataset WHERE kind=?',(metadata['kind'],))
                    if json.loads(row['metadata_json']).get('source_id')==metadata['source_id']]
                if prior:
                    if len(prior)!=1 or json.loads(prior[0]['metadata_json'])!=metadata:
                        raise ValueError('Reference source_id was reused with different metadata or evidence: '+metadata['source_id'])
                    record={'dataset_id':prior[0]['dataset_id'],
                            'count':db.conn.execute('SELECT COUNT(*) FROM quantity WHERE dataset_id=?',
                                (prior[0]['dataset_id'],)).fetchone()[0]}
                    status='REUSED_REFERENCE'
                else:
                    record=import_reference_csv(db,study_id=study['study_id'],csv_path=target,metadata=metadata)
                    status='IMPORTED'
                db.register_artifact(root,study['study_id'],target,
                                     dataset_id=record['dataset_id'],role='REFERENCE_RAW')
                imports.append({'status':status,'file':item['file'],**record})
            except Exception as exc:
                imports.append({'status':'BLOCKED','file':item.get('file'), 'reason':str(exc)})
            log({'phase':'REFERENCE','result':imports[-1]})
        log({'phase':'REFERENCE','event':'FINISHED','message':f'{len(imports)} reference imports checked'})
        # 4. Execute all resolved, executable planned Cases. No binary: capture BLOCKED.
        ids=[row['case_id'] for row in db.conn.execute('SELECT case_id FROM planned_case ORDER BY test_id,role')]
    finally:
        db.close()
    log({'phase':'RUN','event':'STARTED','message':f'{len(ids)} Cases scheduled'})
    executions=(run_cases(db_path,study_path,ids,binary_path=binary_path,
                          python_executable=python_executable,timeout_seconds=timeout_seconds,
                          resume=resume,progress_callback=log,verbose=verbose,ncpu=ncpu,
                          timeouts_by_mode=config.get('timeouts_seconds_by_mode'),
                          timeouts_by_case=config.get('timeouts_seconds_by_case')) if execute else
                [{'case_id':cid,'status':'NOT_EXECUTED','reason':'--no-run requested'} for cid in ids])
    log({'phase':'RUN','event':'FINISHED',
         'message':f'{len(executions)} Case results collected' + (' (no-run)' if not execute else '')})

    db=VVDatabase(db_path)
    try:
        drag_runs=[]
        for item in config.get('parasite_drag_runs',[]):
            if not execute:
                result={'status':'NOT_EXECUTED','reason':'--no-run requested'}
            else:
                try:
                    result=run_parasite_drag(db,root,model=path_of(item['model']),
                        velocity_mps=item['velocity_mps'],altitude_m=item.get('altitude_m',0),
                        delta_temp_k=item.get('delta_temp_k',0),sref=item['sref'],
                        geom_set=item['geom_set'],
                        metadata=json.loads(path_of(item['metadata']).read_text(encoding='utf-8')),
                        python_executable=python_executable,timeout_seconds=timeout_seconds)
                    result={'status':'IMPORTED',**result}
                except Exception as exc:
                    result={'status':'BLOCKED','reason':str(exc)}
            drag_runs.append(result);log({'phase':'PARASITE_DRAG','result':result})

        def resolve_dataset(item, dataset_key, case_key, source_key, kind):
            # Prefer explicitly identified datasets; otherwise use the selected
            # successful Case Run or exactly one reviewed reference source ID.
            if item.get(dataset_key):return item[dataset_key]
            if case_key and case_key in item:
                cid=item[case_key]
                selection=config.get('run_selection',{}).get(cid)
                rows=db.conn.execute("""SELECT ds.dataset_id FROM planned_case pc
                    JOIN case_run_link l ON l.case_id=pc.case_id
                    JOIN case_point pt ON pt.point_id=l.point_id
                    JOIN run r ON r.run_id=pt.run_id
                    JOIN dataset ds ON ds.point_id=pt.point_id
                    WHERE pc.case_id=? AND r.status='SUCCEEDED' AND pt.point_status='MATCHED'
                      AND ds.method_id IN ('POLAR_NATIVE','POLAR_REPROCESS')
                      AND (? IS NULL OR r.run_id=?)
                    ORDER BY r.rowid DESC, ds.rowid DESC LIMIT 1""",
                    (cid,selection,selection)).fetchall()
                if len(rows)!=1:raise ValueError('No successful selected polar dataset: '+cid)
                return rows[0]['dataset_id']
            if source_key and source_key in item:
                matches=[ds['dataset_id'] for ds in db.conn.execute(
                    'SELECT dataset_id,metadata_json FROM dataset WHERE kind=?',(kind,))
                    if json.loads(ds['metadata_json']).get('source_id')==item[source_key]]
                if len(matches)!=1:
                    raise ValueError(f'Expected one {kind} dataset for {source_key}: {item[source_key]} (found {len(matches)})')
                return matches[0]
            raise ValueError(f'Supply {dataset_key}, {case_key}, or {source_key}')

        # 5. Optional reference/drag calculations, never infer missing pairings.
        processed=[]
        for item in config.get('comparisons',[]):
            try:
                records=compare_reference_dataset(db,test_id=item['test_id'],
                    predicted_dataset=resolve_dataset(item,'predicted_dataset','case_id','', 'VSPAERO'),
                    reference_dataset=resolve_dataset(item,'reference_dataset','', 'reference_source_id',
                        {'VER-01':'THEORY','VER-02':'THEORY','XCD-01':'PANEL','XCD-02':'NS',
                         'VAL-01':'WIND_TUNNEL','VAL-02':'WIND_TUNNEL'}[item['test_id']]),
                    conditions_reviewed=item.get('conditions_reviewed',False))
                result={'status':'COMPUTED','test_id':item['test_id'],'records':records}
            except Exception as exc:
                result={'status':'BLOCKED','test_id':item.get('test_id'),'reason':str(exc)}
            processed.append(result);log({'phase':'COMPARE','result':result})
        for item in config.get('drag_corrections',[]):
            try:
                polar=resolve_dataset(item,'polar_dataset','case_id','', 'VSPAERO')
                cd0=resolve_dataset(item,'cd0_dataset','', 'cd0_source_id','EXTERNAL_CD0')
                candidate={'polar_dataset':polar,'cd0_dataset':cd0}
                if item.get('cdo_baseline') is not None:
                    candidate.update({'cdo_baseline':item['cdo_baseline'],
                                      'baseline_metadata':item.get('baseline_metadata')})
                prior=[d for d in db.conn.execute(
                    "SELECT dataset_id,metadata_json FROM dataset WHERE method_id='CD0_CORRECTION'")
                    if json.loads(d['metadata_json'])==candidate]
                if len(prior)>1:raise ValueError('Multiple existing identical drag corrections; resolve before reuse')
                if prior:
                    result={'status':'REUSED_CORRECTION','dataset_id':prior[0]['dataset_id']}
                else:
                    result={'status':'COMPUTED',**create_drag_dataset(db,study_id=study['study_id'],
                        polar_dataset=polar,cd0_dataset=cd0,
                        cdo_baseline=item.get('cdo_baseline'),
                        baseline_metadata=item.get('baseline_metadata'))}
            except Exception as exc:
                result={'status':'BLOCKED','reason':str(exc)}
            processed.append(result);log({'phase':'DRAG','result':result})
        # 5b. Reviewed ParasiteDrag from the same G103A baseline; derive A/B from
        # explicitly identified solver datasets. Never infer skin-friction inputs.
        automatic=config.get('auto_drag')
        if automatic:
            try:
                if not execute:
                    raise RuntimeError('ParasiteDrag calculation requires execution (--no-run)')
                baseline_case=automatic['baseline_case_id']
                baseline_polar=resolve_dataset({'case_id':baseline_case},'polar_dataset','case_id','', 'VSPAERO')
                base_ds=db.row('dataset','dataset_id',baseline_polar)
                flight=db.row('flight_condition','condition_id',base_ds['condition_id'])
                model=path_of(automatic['model'])
                cd_metadata={'kind':'EXTERNAL_CD0','source_id':automatic['source_id'],
                    'axes_convention_id':flight['axes_convention_id'],
                    'sref':flight['sref'],'bref':flight['bref'],'cref':flight['cref'],
                    'mach':flight['mach'],'re_cref':flight['re_cref'],
                    'alpha_deg':flight['alpha_deg'],'beta_deg':flight['beta_deg'],
                    'uncertainty_status':'NOT_ASSESSED','usage_class':'REFERENCE',
                    'reference_reviewed':True,'geometry_equivalence_reviewed':False,
                    'physical_approval':'NOT_ASSESSED',
                    'drag_surface_selection':automatic['geometry_names']}
                drag_result=run_parasite_drag(db,root,model=model,
                    velocity_mps=automatic['velocity_mps'],
                    altitude_m=automatic.get('altitude_m',0),
                    delta_temp_k=automatic.get('delta_temp_k',0),
                    sref=flight['sref'],geom_set=automatic['geom_set'],
                    geometry_names=automatic['geometry_names'],metadata=cd_metadata,
                    python_executable=python_executable,timeout_seconds=timeout_seconds)
                cd0_id=drag_result['dataset_id']
                drag_runs.append({'status':'IMPORTED' if drag_result.get('status')!='REUSED_PARASITE'
                    else 'REUSED_PARASITE',**drag_result})
                log({'phase':'PARASITE_DRAG_AUTO','result':drag_runs[-1]})
                cdo_record=db.conn.execute("SELECT canonical_value FROM quantity WHERE "
                    "dataset_id=? AND qoi='CDo' AND variant='DEFAULT' AND value_status='AVAILABLE'",
                    (baseline_polar,)).fetchone()
                if cdo_record is None:
                    raise ValueError('Method B requires native baseline CDo in selected polar')
                baseline_value=float(cdo_record['canonical_value'])
                for cid in automatic['source_case_ids']:
                    try:
                        polar_id=resolve_dataset({'case_id':cid},'polar_dataset','case_id','', 'VSPAERO')
                        polar_ds=db.row('dataset','dataset_id',polar_id)
                        cond=db.row('flight_condition','condition_id',polar_ds['condition_id'])
                        # ParasiteDrag CD0 is only reusable at the identical freestream
                        # and reference dimensions. Different AoA is allowed for polars.
                        import math
                        if (not math.isclose(cond['mach'],flight['mach'],abs_tol=1e-6)
                            or not math.isclose(cond['re_cref'],flight['re_cref'],rel_tol=1e-7)
                            or not math.isclose(cond['sref'],flight['sref'],rel_tol=1e-8)):
                            raise ValueError('CD0 source conditions differ from solver Case')
                        for method in automatic['methods']:
                            baseline_meta=({'condition_id':base_ds['condition_id'],
                                'dataset_id':baseline_polar,
                                'approval_status':'EXPLORATORY_NOT_ENGINEERING_APPROVED'}
                                if method=='B' else None)
                            expected={'polar_dataset':polar_id,'cd0_dataset':cd0_id}
                            if method=='B':expected.update({'cdo_baseline':baseline_value,
                                'baseline_metadata':baseline_meta})
                            existing=[row['dataset_id'] for row in db.conn.execute(
                                "SELECT dataset_id,metadata_json FROM dataset WHERE method_id='CD0_CORRECTION'")
                                if json.loads(row['metadata_json'])==expected]
                            if len(existing)>1:raise ValueError('Ambiguous repeated drag correction')
                            if existing:
                                result={'status':'REUSED_CORRECTION','dataset_id':existing[0]}
                            else:
                                # Method A-only dataset and method A+B dataset are kept
                                # independent so changing the baseline cannot mutate A.
                                result={'status':'COMPUTED',**create_drag_dataset(db,
                                    study_id=study['study_id'],polar_dataset=polar_id,
                                    cd0_dataset=cd0_id,
                                    cdo_baseline=baseline_value if method=='B' else None,
                                    baseline_metadata=baseline_meta)}
                            result.update({'test_id':'DRAG-01' if method=='A' else 'DRAG-02',
                                           'case_id':cid,'method':method})
                            processed.append(result);log({'phase':'DRAG_AUTO','result':result})
                    except Exception as exc:
                        result={'status':'BLOCKED','case_id':cid,'reason':str(exc)}
                        processed.append(result);log({'phase':'DRAG_AUTO','result':result})
            except Exception as exc:
                result={'status':'BLOCKED','reason':str(exc)}
                drag_runs.append(result);log({'phase':'DRAG_AUTO','result':result})
        log({'phase':'AGGREGATE','event':'STARTED','message':'Collecting results and numerical comparisons'})
        # 6. Freeze which successful Run and Dataset was used for *this* report.
        snapshot=collect_study_results(db,root,selected_runs=config.get('run_selection'))
        from .VVDatabase import source_fingerprint, VISUALIZATION_SOURCES
        snapshot['execution_manifest_sha256']=digest(manifest_path)
        snapshot['visualization_source_sha256']=source_fingerprint(
            REPO_ROOT,VISUALIZATION_SOURCES)
        snapshot['preparation']=prepared
        snapshot['execution']=executions
        snapshot['reference_imports']=imports
        snapshot['parasite_drag_runs']=drag_runs
        snapshot['derivations']=processed
        numerical={}
        for test in CASE_ROLES:
            cases=[c for c in snapshot['cases'] if c['test_id']==test and c['run_id']]
            if test.startswith('CTRL-'):
                results={}
                for qoi in sorted({q['qoi'] for c in cases for q in c['quantities']
                    if q['variant']=='DEFAULT' and q['status']=='AVAILABLE'}):
                    samples={int(c['role'].split('_')[-1]):next(
                        (q['value'] for q in c['quantities']
                         if q['qoi']==qoi and q['variant']=='DEFAULT' and q['status']=='AVAILABLE'),None)
                        for c in cases if c['role'].startswith('DELTA_')}
                    samples={k:v for k,v in samples.items() if v is not None}
                    if len(samples)>=2:
                        results[qoi]=evaluate_control_sweep(samples,qoi)
                if results:numerical[test]=results
            elif test.startswith(('GRID-','WAKE-','GEO-','STA-')):
                results={}
                for qoi in sorted({q['qoi'] for c in cases for q in c['quantities']
                    if q['variant']=='DEFAULT' and q['status']=='AVAILABLE'}):
                    pairs=[]
                    for c in cases:
                        q=next((q for q in c['quantities'] if q['qoi']==qoi
                            and q['variant']=='DEFAULT' and q['status']=='AVAILABLE'),None)
                        if q:pairs.append({'case_id':c['case_id'],'value':q['value'],
                                           'quantity_id':q['quantity_id']})
                    if len(pairs)>=2:
                        results[qoi]={'samples':pairs,
                            'successive_differences':[pairs[i+1]['value']-pairs[i]['value']
                                                     for i in range(len(pairs)-1)]}
                        widths=config.get('grid_widths',{}).get(test)
                        if test.startswith('GRID-') and len(pairs)==3 and widths:
                            results[qoi]['gci_candidate']=richardson_gci(
                                [p['value'] for p in pairs],widths)
                if results:numerical[test]=results
        snapshot['numerical']=numerical
        snapshot_path=destination/'results_snapshot.json'
        snapshot_path.write_text(json.dumps(snapshot,ensure_ascii=False,indent=2),encoding='utf-8')
        with (destination/'quantities.csv').open('w',newline='',encoding='utf-8-sig') as stream:
            writer=csv.writer(stream)
            writer.writerow(['test_id','case_id','run_id','dataset_id','quantity_id',
                             'qoi','variant','component','value','unit','status'])
            for case in snapshot['cases']:
                for q in case['quantities']:
                    writer.writerow([case['test_id'],case['case_id'],case['run_id'],
                        q['dataset_id'],q['quantity_id'],q['qoi'],q['variant'],
                        q['component'],q['value'],q['unit'],q['status']])
        with (destination/'history.csv').open('w',newline='',encoding='utf-8-sig') as stream:
            writer=csv.writer(stream)
            writer.writerow(['test_id','case_id','run_id','solver_case','iteration','qoi','value'])
            for case in snapshot['cases']:
                for point in case['history']:
                    for key,value in point.items():
                        if key not in ('case','iteration') and isinstance(value,(int,float)):
                            writer.writerow([case['test_id'],case['case_id'],case['run_id'],
                                             point['case'],point['iteration'],key,value])
        log({'phase':'AGGREGATE','event':'FINISHED','message':'CSV and results snapshot written'})
        log({'phase':'REPORT','event':'STARTED','message':'Generating figures, tables and review report'})
        # 7. Generate figures and the populated, but not technically approved, report.
        charts=create_study_plots(snapshot,destination/'figures')
        snapshot['charts']=charts
        snapshot_path.write_text(json.dumps(snapshot,ensure_ascii=False,indent=2),encoding='utf-8')
        for file in (snapshot_path,destination/'quantities.csv',destination/'history.csv'):
            db.register_artifact(root,study['study_id'],file,role='REVIEW_DATA')
        for chart in charts:
            db.register_artifact(root,study['study_id'],chart['path'],role='REVIEW_FIGURE')
        # Copy SQLite with the official backup API, not a raw file copy of
        # a potentially WAL-enabled live database. This snapshot is immutable.
        import sqlite3
        db_snapshot=destination/'database_snapshot.sqlite'
        with sqlite3.connect(db_snapshot) as backup:
            db.conn.backup(backup)
        db.register_artifact(root,study['study_id'],db_snapshot,role='REVIEW_DB_SNAPSHOT')
        report=write_evaluation_report(db,destination/'evaluation_report.md',
                                       selected_runs=config.get('run_selection'),
                                       results=snapshot,database_snapshot=db_snapshot)
        export_csv(db,destination/'tables')
        report_id=db.register_artifact(root,study['study_id'],report,role='REVIEW_REPORT')
        log({'phase':'REPORT','event':'FINISHED','message':'REVIEW_REPORT_GENERATED'})
        log({'phase':'SUMMARY','event':'FINISHED',
             'message':f"{len(ids)} Cases; {sum(x['status'] in ('SUCCEEDED','REUSED','REPROCESSED','REUSED_POSTPROCESS') for x in executions)} succeeded/reused; {sum(x['status'] not in ('SUCCEEDED','REUSED','REPROCESSED','REUSED_POSTPROCESS') for x in executions)} blocked/failed"})
        return {'status':'REVIEW_REPORT_GENERATED','output_dir':str(destination),
                'report':str(report),'report_artifact_id':report_id,
                'case_count':len(ids),'successful_or_reused_runs':sum(
                    1 for x in executions if x['status'] in ('SUCCEEDED','REUSED','REPROCESSED','REUSED_POSTPROCESS')),
                'blocked_or_failed':sum(1 for x in executions if x['status'] not in
                    ('SUCCEEDED','REUSED','REPROCESSED','REUSED_POSTPROCESS')),
                'figure_count':len(charts),'run_selection':{
                    row['case_id']:row['run_id'] for row in snapshot['cases'] if row['run_id']}}
    finally:
        db.close()
