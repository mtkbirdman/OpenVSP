"""CLI: python -m openvsppy.vv --help. All metadata operations work without OpenVSP."""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from openvsppy.vv.VVDatabase import VVDatabase, load_study, digest, execution_fingerprint, postprocess_fingerprint, CASE_ROLES
from openvsppy.vv.VVRunner import check_case, run_case, run_cases, run_parasite_drag
from openvsppy.vv.VVAutomation import prepare_cases, case_difference, execute_study
from openvsppy.vv.VVPostprocess import import_reference_csv, compare_quantities, screen_convergence, create_drag_dataset, read_history, read_mesh_stats, compare_reference_dataset, richardson_gci, evaluate_control_sweep, store_assessment, reprocess_run, evaluate_test
from openvsppy.vv.VVReport import export_csv, write_status, write_evaluation_report


def positive_ncpu(value):
    try:
        number=int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError('NCPU must be a positive integer') from exc
    if number < 1:
        raise argparse.ArgumentTypeError('NCPU must be >= 1')
    return number


def main(argv=None):
    p=argparse.ArgumentParser(description='OpenVSP 3.50.5 V&V Study DB / Runner')
    p.add_argument('--study',default=str(Path(__file__).resolve().parents[3] / 'vv_study' / 'study.json'),help='Study configuration path')
    p.add_argument('--db',default=None,help='SQLite path (default vv_study/01_database/vv.sqlite)')
    sub=p.add_subparsers(dest='command',required=True)
    sub.add_parser('init',help='Initialize SQLite for the 25 approved tests')
    sub.add_parser('plan',help='List approved Test IDs and planned Case roles')
    reg=sub.add_parser('register-case',help='Register explicit full input JSON then resolve one Case ID')
    reg.add_argument('--file',required=True)
    check=sub.add_parser('check');check.add_argument('--case',required=True)
    check.add_argument('--binary')
    run=sub.add_parser('run');run.add_argument('--case',required=True)
    run.add_argument('--binary');run.add_argument('--python');run.add_argument('--timeout',type=int,default=3600)
    run.add_argument('--ncpu',type=positive_ncpu,help='OpenMP threads (overrides Case NCPU; default 1)')
    imp=sub.add_parser('import-reference');imp.add_argument('--file',required=True);imp.add_argument('--metadata',required=True)
    evaluate=sub.add_parser('evaluate');evaluate.add_argument('--test',required=True)
    evaluate.add_argument('--selection-file',help='JSON: case_id -> run_id')
    evaluate.add_argument('--label',help='Freeze explicit Run selection under this label')
    gci=sub.add_parser('gci',help='Calculate provisional GCI only when representative grid h is supplied')
    gci.add_argument('--values',nargs=3,required=True,type=float)
    gci.add_argument('--widths',nargs=3,required=True,type=float)
    prep=sub.add_parser('prepare',help='Generate versioned GRID/WAKE case files from baseline')
    prep.add_argument('--test',required=True);prep.add_argument('--base',required=True)
    prep.add_argument('--dry-run',action='store_true');prep.add_argument('--register',action='store_true')
    prep.add_argument('--span',type=float,help='VAL +/- alpha/beta degrees (approved value required)')
    batch=sub.add_parser('run-test',help='Sequentially execute resolved Cases for one Test ID')
    batch.add_argument('--test',required=True);batch.add_argument('--binary',required=True)
    batch.add_argument('--python');batch.add_argument('--timeout',type=int,default=3600)
    batch.add_argument('--resume',action='store_true');batch.add_argument('--jobs',type=int,default=1)
    batch.add_argument('--ncpu',type=positive_ncpu,help='OpenMP threads (overrides Case NCPU; default 1)')
    full=sub.add_parser('execute-all',help='Prepare approved inputs, run all Cases, aggregate, plot and write a 13-chapter reviewer report')
    full.add_argument('--config',required=True,help='User-reviewed execution manifest JSON')
    full.add_argument('--binary',help='VSPAERO 3.50.5 executable')
    full.add_argument('--python',help='Python with OpenVSP 3.50.5 installed')
    full.add_argument('--timeout',type=int,default=3600)
    full.add_argument('--no-resume',action='store_true')
    full.add_argument('--no-run',action='store_true',help='Aggregate/report existing DB only')
    full.add_argument('--ncpu',type=positive_ncpu,help='OpenMP threads (overrides each Case NCPU; default 1)')
    full.add_argument('--verbose',action='store_true',help='Show native Wake progress and detailed Run information')
    inventory=sub.add_parser('assess-all',help='Inventory evidence for all 25 Test IDs; no automatic approval')
    inventory.add_argument('--output',help='Optional JSON file for evidence inventory')
    inventory.add_argument('--persist',action='store_true',help='Persist unapproved evidence inventory as assessment records')
    reproc=sub.add_parser('reprocess-run',help='Recalculate postprocessed datasets from a successful immutable Run')
    reproc.add_argument('--run',required=True)
    inspect=sub.add_parser('inspect-run',help='Audit a Run and parse native history/mesh diagnostics')
    inspect.add_argument('--run',required=True)
    report=sub.add_parser('report',help='Generate 25-Test evidence report (not an approval)')
    report.add_argument('--output',help='Markdown output file')
    cmp=sub.add_parser('compare');cmp.add_argument('--test',required=True);cmp.add_argument('--predicted',required=True);cmp.add_argument('--reference',required=True);cmp.add_argument('--rule');cmp.add_argument('--conditions-verified',action='store_true')
    cmp_batch=sub.add_parser('compare-test',help='Match all QoIs across explicit, reviewed VSPAERO/reference dataset pair')
    cmp_batch.add_argument('--test',required=True)
    cmp_batch.add_argument('--polar-dataset',required=True);cmp_batch.add_argument('--reference-dataset',required=True)
    cmp_batch.add_argument('--conditions-reviewed',action='store_true')
    cmp_batch.add_argument('--rule')
    drag_exec=sub.add_parser('drag-run',help='Isolated 3.50.5 ParasiteDrag; import CD0 with explicit metadata')
    drag_exec.add_argument('--model',required=True);drag_exec.add_argument('--velocity-mps',required=True,type=float)
    drag_exec.add_argument('--altitude-m',type=float,default=0);drag_exec.add_argument('--delta-temp-k',type=float,default=0)
    drag_exec.add_argument('--sref',required=True,type=float);drag_exec.add_argument('--geom-set',required=True,type=int)
    drag_exec.add_argument('--metadata',required=True,help='Reviewed EXTERNAL_CD0 metadata JSON')
    drag_exec.add_argument('--python',help='3.50.5 Conda Python interpreter')
    drag=sub.add_parser('drag');drag.add_argument('--polar-dataset',required=True);drag.add_argument('--cd0-dataset',required=True)
    drag.add_argument('--cdo-baseline',type=float);drag.add_argument('--baseline-condition-id')
    exp=sub.add_parser('export',help='Generate reviewer CSV sheets and status Markdown')
    exp.add_argument('--output',help='Separate output directory for each DB/study')
    sub.add_parser('verify-artifacts',help='Rehash registered raw assets')
    args=p.parse_args(argv)
    study_path,study,_,_=load_study(args.study)
    root=study_path.parent
    db_path=Path(args.db) if args.db else root/'01_database'/'vv.sqlite'
    if args.command=='execute-all':
        # The workflow owns its DB connections, including those used by workers.
        result=execute_study(study_path,args.config,db_path=db_path,binary_path=args.binary,
            python_executable=args.python,timeout_seconds=args.timeout,
            resume=not args.no_resume,execute=not args.no_run,verbose=args.verbose,ncpu=args.ncpu)
        print(json.dumps(result,ensure_ascii=False,indent=2))
        return 0  # The report preserves every blocked/failed Case without declaring success
    if args.command=='init':
        db=VVDatabase(db_path)
        try:
            n=db.initialize(study_path,root/'01_database'/'schema.sql')
        finally:
            db.close()
        print(json.dumps({'database':str(db_path),'test_ids':n,
                          'cases':sum(map(len,CASE_ROLES.values()))},ensure_ascii=False))
        return 0
    db=VVDatabase(db_path)
    try:
        study_record=db.row('study','study_id',study['study_id'])
        if study_record is None:
            raise ValueError('Study is not initialized. Use: python -m openvsppy.vv init')
        if study_record['plan_sha256']!=study['plan']['sha256']:
            raise ValueError('Plan SHA-256 differs from initialized study.')
        if args.command=='plan':
            for r in db.conn.execute('SELECT case_id,execute_policy FROM planned_case ORDER BY test_id,role'):
                print(f'{r["case_id"]:28} {r["execute_policy"]}')
        elif args.command=='register-case':
            definition=json.loads(Path(args.file).read_text(encoding='utf-8'))
            db.register_case_document(study['study_id'],definition)
            print('Resolved',definition['case_id'])
        elif args.command=='check':
            errs=check_case(db,args.case,study_path,args.binary)
            print('PASS: static precheck (not an OpenVSP execution)' if not errs else 'PRECHECK FAILED:\n'+'\n'.join(errs))
            return int(bool(errs))
        elif args.command=='run':
            db.close();db=None
            output=run_case(db_path,args.case,study_path,python_executable=args.python,
                            binary_path=args.binary,timeout_seconds=args.timeout,ncpu=args.ncpu)
            print(json.dumps(output,ensure_ascii=False,indent=2))
            return 0 if output['status']=='SUCCEEDED' else 2
        elif args.command=='import-reference':
            metadata=json.loads(Path(args.metadata).read_text(encoding='utf-8'))
            original=Path(args.file).resolve()
            if not original.is_file():raise FileNotFoundError(original)
            path=root/'02_reference'/('external_cd0' if metadata.get('kind')=='EXTERNAL_CD0' else {
                'WIND_TUNNEL':'wind_tunnel','PANEL':'panel','NS':'ns','THEORY':'research'}.get(metadata.get('kind'),'research'))
            path.mkdir(parents=True,exist_ok=True)
            dst=path/(digest(original)[:12]+'_'+original.name)
            if not dst.is_file():shutil.copy2(original,dst)
            metadata['source_sha256']=digest(dst)
            output=import_reference_csv(db,study_id=study['study_id'],csv_path=dst,metadata=metadata)
            db.register_artifact(root,study['study_id'],dst,dataset_id=output['dataset_id'],role='REFERENCE_RAW')
            print(json.dumps(output,indent=2))
        elif args.command=='compare':
            print(json.dumps(compare_quantities(db,test_id=args.test,predicted_id=args.predicted,
               reference_id=args.reference,rule_id=args.rule,condition_verified=args.conditions_verified),indent=2))
        elif args.command=='compare-test':
            results=compare_reference_dataset(db,test_id=args.test,predicted_dataset=args.polar_dataset,
                reference_dataset=args.reference_dataset,conditions_reviewed=args.conditions_reviewed,
                rule_id=args.rule)
            print(json.dumps({'test_id':args.test,'comparisons':results},indent=2))
        elif args.command=='drag-run':
            metadata=json.loads(Path(args.metadata).read_text(encoding='utf-8'))
            result=run_parasite_drag(db,root,model=args.model,
                velocity_mps=args.velocity_mps,altitude_m=args.altitude_m,
                delta_temp_k=args.delta_temp_k,sref=args.sref,geom_set=args.geom_set,
                metadata=metadata,python_executable=args.python)
            print(json.dumps(result,ensure_ascii=False,indent=2))
        elif args.command=='drag':
            baseline_meta={'condition_id':args.baseline_condition_id} if args.baseline_condition_id else None
            print(json.dumps(create_drag_dataset(db,study_id=study['study_id'],
                  polar_dataset=args.polar_dataset,cd0_dataset=args.cd0_dataset,
                  cdo_baseline=args.cdo_baseline,baseline_metadata=baseline_meta),indent=2))
        elif args.command=='gci':
            print(json.dumps(richardson_gci(args.values,args.widths),indent=2))
        elif args.command=='prepare':
            if args.register and args.dry_run:
                raise ValueError('--register and --dry-run cannot be combined')
            items=prepare_cases(args.test,args.base,root,dry_run=args.dry_run,span=args.span)
            if args.register:
                for item in items:
                    if item['case_id']=='GRID-03-NREF_1':
                        db.add_extra_planned_case(study['study_id'],'GRID-03','NREF_1')
                    definition=json.loads(Path(item['case_json']).read_text(encoding='utf-8'))
                    prior=db.row('planned_case','case_id',item['case_id'])
                    if prior and prior['geometry_id']:
                        if prior['geometry_id']!=definition['geometry']['geometry_id'] or prior['setting_id']!=definition['numerical_setting']['setting_id']:
                            raise ValueError('Existing Case differs; do not mutate a previous execution: '+item['case_id'])
                        item['registration']='ALREADY_RESOLVED'
                    else:
                        db.register_case_document(study['study_id'],definition)
                        item['registration']='REGISTERED'
            print(json.dumps({'test_id':args.test,'dry_run':args.dry_run,'cases':items},ensure_ascii=False,indent=2))
        elif args.command=='run-test':
            if args.jobs!=1:raise ValueError('Only --jobs 1 is implemented; OpenVSP uses global state')
            cases=list(db.conn.execute('SELECT case_id FROM planned_case WHERE test_id=? ORDER BY role',(args.test,)))
            if not cases:raise ValueError('Unknown Test ID: '+args.test)
            db.close();db=None
            results=run_cases(db_path,study_path,[r['case_id'] for r in cases],
                binary_path=args.binary,python_executable=args.python,
                timeout_seconds=args.timeout,resume=args.resume,ncpu=args.ncpu,
                progress_path=root/'03_run_raw'/'batch_progress.jsonl')
            print(json.dumps({'test_id':args.test,'results':results},indent=2,ensure_ascii=False))
            return 2 if any(r['status'] not in ('SUCCEEDED','REUSED','REPROCESSED','REUSED_POSTPROCESS') for r in results) else 0
        elif args.command=='assess-all':
            records=[]
            for t in db.conn.execute('SELECT * FROM test_definition ORDER BY test_id'):
                total,resolved=db.conn.execute('''SELECT COUNT(*),SUM(CASE WHEN geometry_id IS NOT NULL AND condition_id IS NOT NULL AND setting_id IS NOT NULL THEN 1 ELSE 0 END) FROM planned_case WHERE test_id=?''',(t['test_id'],)).fetchone()
                runs=db.conn.execute('''SELECT COUNT(DISTINCT r.run_id) FROM planned_case pc JOIN case_run_link l ON l.case_id=pc.case_id JOIN case_point pt ON pt.point_id=l.point_id JOIN run r ON r.run_id=pt.run_id WHERE pc.test_id=? AND r.status='SUCCEEDED' ''',(t['test_id'],)).fetchone()[0]
                comparisons=db.conn.execute('SELECT COUNT(*) FROM comparison WHERE test_id=?',(t['test_id'],)).fetchone()[0]
                assessments=db.conn.execute('SELECT COUNT(*) FROM assessment WHERE test_id=?',(t['test_id'],)).fetchone()[0]
                records.append({'test_id':t['test_id'],'requirement':t['requirement'],
                    'planned_cases':total,'resolved_cases':resolved or 0,'successful_runs':runs,
                    'comparisons':comparisons,'assessments':assessments,
                    'engineering_verdict':'REVIEW_REQUIRED' if assessments else 'INSUFFICIENT_EVIDENCE'})
            if args.persist:
                from openvsppy.vv.VVPostprocess import store_assessment
                for item in records:
                    store_assessment(db,study_id=study['study_id'],test_id=item['test_id'],qoi='EVIDENCE_COMPLETENESS',
                                     evidence=[],reason='Automated inventory; no aerodynamic certification. '+str(item))
            payload={'status':'EVIDENCE_INVENTORY_ONLY','test_ids':len(records),'tests':records}
            if args.output:
                Path(args.output).parent.mkdir(parents=True,exist_ok=True)
                Path(args.output).write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
            print(json.dumps(payload,indent=2,ensure_ascii=False))
        elif args.command=='reprocess-run':
            print(json.dumps(reprocess_run(db,study_root=root,run_id=args.run),indent=2))
        elif args.command=='inspect-run':
            run_record=db.row('run','run_id',args.run)
            if run_record is None:raise ValueError('Unknown Run ID: '+args.run)
            rd=root/run_record['run_dir'];audit=[]
            for a in db.conn.execute('SELECT * FROM artifact WHERE run_id=? ORDER BY relative_path',(args.run,)):
                path=root/a['relative_path']
                status='VERIFIED' if path.is_file() and digest(path)==a['sha256'] else 'CHANGED_OR_MISSING'
                audit.append({'path':a['relative_path'],'status':status})
            history=list((rd/'raw').glob('*.history')) if (rd/'raw').is_dir() else []
            log=rd/'raw'/'vspaero_api_stdout.log'
            if not log.is_file():log=rd/'logs'/'stdout.log'
            iters=read_history(history[0]) if len(history)==1 else []
            print(json.dumps({'run_id':args.run,'status':run_record['status'],
                'artifacts':audit,'history':{'iterations':len(iters),'final':iters[-1] if iters else None},
                'mesh':read_mesh_stats(log) if log.is_file() else {}},ensure_ascii=False,indent=2))
        elif args.command=='report':
            outfile=Path(args.output) if args.output else root/'04_reports'/'evaluation_report.md'
            print(json.dumps({'report':str(write_evaluation_report(db,outfile))},ensure_ascii=False))
        elif args.command=='evaluate':
            selection=json.loads(Path(args.selection_file).read_text(encoding='utf-8')) if args.selection_file else None
            print(json.dumps(evaluate_test(db,study['study_id'],args.test,
                selection=selection,label=args.label),ensure_ascii=False,indent=2))
        elif args.command=='export':
            output=Path(args.output) if args.output else root/'04_reports';counts=export_csv(db,output/'tables')
            summary=write_status(db,output/'status.md')
            print(json.dumps({'csv_tables':counts,'summary':str(summary)},ensure_ascii=False))
        elif args.command=='verify-artifacts':
            results=db.verify_artifacts(root)
            for path,status in results:print(status,path)
            return 2 if any(status!='VERIFIED' for _,status in results) else 0
        return 0
    finally:
        if db is not None:db.close()



if __name__ == "__main__":
    raise SystemExit(main())
