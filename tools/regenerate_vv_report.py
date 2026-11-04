#!/usr/bin/env python3
"""Generate a v20261104.1 reviewer report from historical V&V evidence.

Read the original SQLite in SQLite read-only mode, freeze it via backup, then
collect results/plots/tables/reports only on the frozen DB. This does NOT run
OpenVSP, migrate DBs, register artifacts, or update the historical study.
"""
from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path


def run():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', required=True, type=Path, help='Git working tree containing src/openvsppy')
    ap.add_argument('--study-root', required=True, type=Path, help='Historical vv_study root containing 03_run_raw/')
    ap.add_argument('--db', required=True, type=Path, help='Historical SQLite; never modified')
    ap.add_argument('--output', type=Path, help='Optional new output directory (must not exist)')
    ap.add_argument('--selection', type=Path, help='Optional JSON mapping of case IDs to successful run IDs')
    ap.add_argument('--report-profile', type=Path, help='Optional explicit review metadata and display baseline (not approval)')
    args = ap.parse_args()

    repo, root, source_db = (p.resolve() for p in (args.repo, args.study_root, args.db))
    if not (repo/'src/openvsppy/vv/VVReport.py').is_file():
        ap.error('Not a compatible OpenVSP repository: ' + str(repo))
    if not (repo/'docs/vv/templates/evaluation_report.md').is_file():
        ap.error('Report template missing in repository: ' + str(repo))
    if not source_db.is_file():
        ap.error('Historical DB not found: ' + str(source_db))
    if not (root/'03_run_raw').is_dir():
        ap.error('Run archive missing: ' + str(root/'03_run_raw'))
    if not source_db.is_relative_to(root):
        ap.error('DB must be located under the supplied historical study root')
    if source_db.name != 'vv_20261031_10.sqlite':
        print('WARNING: DB name differs from expected 20261031_10; verify selection.', file=sys.stderr)

    sys.path.insert(0, str(repo/'src'))
    from openvsppy.vv.VVDatabase import (VVDatabase, CASE_ROLES, digest,
        source_fingerprint, VISUALIZATION_SOURCES)
    from openvsppy.vv.VVPostprocess import collect_study_results, evaluate_control_sweep
    from openvsppy.vv.VVPlots import create_study_plots
    from openvsppy.vv.VVReport import export_csv, write_evaluation_report, write_status

    report_profile = None
    if args.report_profile:
        report_profile = json.loads(args.report_profile.read_text(encoding='utf-8-sig'))
    selected = None
    if args.selection:
        selected = json.loads(args.selection.read_text(encoding='utf-8-sig'))
        if not isinstance(selected, dict) or not all(isinstance(k, str) and isinstance(v, str)
                                                      for k, v in selected.items()):
            raise ValueError('--selection must contain {"CASE_ID":"RUN_ID",...}')

    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    dest = args.output.resolve() if args.output else root/'04_reports'/f'offline_regenerated_{stamp}_{uuid.uuid4().hex[:8]}'
    if dest.exists():
        raise FileExistsError(f'Refusing to overwrite an existing output directory: {dest}')
    if dest == root or dest == repo:
        raise ValueError('Refusing to use a source directory as output')
    if dest.is_relative_to(root/'03_run_raw') or dest.is_relative_to(root/'01_database'):
        raise ValueError('Output cannot be under original run or DB directories')
    dest.mkdir(parents=True, exist_ok=False)
    print('Output:', dest, flush=True)

    # Do not load study.json: older Study metadata references an older Plan
    # with a different hash/name. Preserve the original Plan hash from SQLite.
    # Source is opened with mode=ro; backup() produces a separate frozen DB.
    old_db_sha256_before = digest(source_db)
    frozen_path = dest/'database_snapshot.sqlite'
    source = sqlite3.connect(source_db.as_uri() + '?mode=ro', uri=True)
    try:
        source.execute('PRAGMA query_only=ON')
        integrity = source.execute('PRAGMA quick_check').fetchone()[0]
        if integrity != 'ok':
            raise ValueError(f'Original SQLite quick_check failed: {integrity}')
        with sqlite3.connect(frozen_path) as backup:
            source.backup(backup)
    finally:
        source.close()
    print('Frozen SQLite:', frozen_path, flush=True)

    db = VVDatabase(frozen_path)
    db.conn.execute('PRAGMA query_only=ON')
    try:
        if db.conn.execute('PRAGMA quick_check').fetchone()[0] != 'ok':
            raise ValueError('Frozen DB integrity check failed')
        if db.conn.execute('PRAGMA foreign_key_check').fetchone():
            raise ValueError('Frozen DB foreign_key_check failed')
        study = db.conn.execute('SELECT study_id,plan_sha256 FROM study').fetchone()
        if study is None:
            raise ValueError('The source DB contains no registered Study')
        print('Study:', study['study_id'], '/ historical plan:', study['plan_sha256'][:12], flush=True)
        snapshot = collect_study_results(db, root, selected_runs=selected)
        snapshot['schema_version'] = 'REPORT_SNAPSHOT_V2_RELATIVE_PATHS'
        snapshot['execution_manifest_sha256'] = 'NOT_USED_READ_ONLY_OFFLINE_REGENERATION'
        snapshot['visualization_source_sha256'] = source_fingerprint(repo, VISUALIZATION_SOURCES)
        snapshot['preparation'] = []
        snapshot['execution'] = []
        snapshot['reference_imports'] = []
        snapshot['parasite_drag_runs'] = []
        snapshot['derivations'] = []
        snapshot['numerical'] = _numerical_screen(snapshot, CASE_ROLES, evaluate_control_sweep)
        snapshot['offline_regeneration'] = {
            'solver_executed': False,
            'historical_db_path': str(source_db),
            'historical_db_sha256_at_start': old_db_sha256_before,
            'source_plan_sha256_from_db': study['plan_sha256'],
            'warning': 'No Study config/Plan replacement; no Native Run or derivative reprocessing.',
        }
        print('Aggregated selected Runs and Quantities.', flush=True)
        charts = create_study_plots(snapshot, dest/'figures')
        snapshot['charts'] = charts
        (dest/'results_snapshot.json').write_text(json.dumps(snapshot, ensure_ascii=False, indent=2), encoding='utf-8')
        print('Generated figures:', len(charts), flush=True)

        with (dest/'quantities.csv').open('w', newline='', encoding='utf-8-sig') as file:
            w = csv.writer(file)
            w.writerow(['test_id','case_id','run_id','dataset_id','quantity_id',
                        'qoi','variant','component','value','unit','status'])
            for case in snapshot['cases']:
                for q in case['quantities']:
                    w.writerow([case['test_id'],case['case_id'],case['run_id'],q['dataset_id'],
                        q['quantity_id'],q['qoi'],q['variant'],q['component'],q['value'],q['unit'],q['status']])
        with (dest/'history.csv').open('w', newline='', encoding='utf-8-sig') as file:
            w = csv.writer(file)
            w.writerow(['test_id','case_id','run_id','solver_case','iteration','qoi','value'])
            for case in snapshot['cases']:
                for point in case['history']:
                    for name, value in point.items():
                        if name not in ('case','iteration') and isinstance(value, (int, float)):
                            w.writerow([case['test_id'],case['case_id'],case['run_id'],
                                        point.get('case'),point.get('iteration'),name,value])
        with (dest/'figure_manifest.csv').open('w', newline='', encoding='utf-8-sig') as file:
            w = csv.writer(file)
            w.writerow(['figure_id','test_id','plot_kind','qoi','x_label','x_unit',
                'y_label','y_unit','source_quantity_ids','source_dataset_ids','path',
                'sha256','display_priority','generation_source_sha256','note'])
            for chart in charts:
                image = dest/chart['path']
                if not image.is_file():
                    raise FileNotFoundError(image)
                sources = chart['sources']
                w.writerow([chart['figure_id'],chart['test_id'],chart['kind'],chart['qoi'],
                    chart.get('x_label',''),chart.get('x_unit',''),chart.get('y_label',''),
                    chart.get('y_unit',''),
                    json.dumps(sorted({x['quantity_id'] for x in sources if x.get('quantity_id')})),
                    json.dumps(sorted({x['dataset_id'] for x in sources if x.get('dataset_id')})),
                    chart['path'],digest(image),chart.get('display_priority','DIAGNOSTIC'),
                    snapshot['visualization_source_sha256'],chart.get('note','')])
        export_csv(db, dest/'tables')
        write_status(db, dest/'study_status.md')
        if args.report_profile:
            import shutil
            shutil.copyfile(args.report_profile, dest/'report_profile_input.json')
        write_evaluation_report(db, dest/'evaluation_report.md', results=snapshot,
                                database_snapshot=frozen_path, report_profile=report_profile)
        all_reports = list((dest/'test_reports').glob('*.md'))
        if len(all_reports) != 25:
            raise ValueError(f'Expected 25 per-test reports, got {len(all_reports)}')
        all_quantity_ids = {q['quantity_id'] for case in snapshot['cases'] for q in case['quantities']}
        report_text = '\n'.join(path.read_text(encoding='utf-8') for path in all_reports)
        absent = [qid for qid in all_quantity_ids if qid not in report_text]
        if absent:
            raise ValueError(f'Quantity IDs missing from per-test Markdown: {absent[:5]}')
        for fig in charts:
            if not (dest/fig['path']).is_file():
                raise FileNotFoundError('Missing linked chart ' + fig['path'])
    finally:
        db.close()

    if digest(source_db) != old_db_sha256_before:
        raise RuntimeError('Source DB changed during generation; review concurrent updates.')
    evidence = []
    for path in sorted(dest.rglob('*')):
        if path.is_file() and path.name != 'evidence_manifest.json':
            evidence.append({'path': path.relative_to(dest).as_posix(),
                             'sha256': digest(path), 'size_bytes':path.stat().st_size})
    (dest/'evidence_manifest.json').write_text(json.dumps({
        'schema_version': 'REVIEW_EVIDENCE_MANIFEST_V1',
        'database_source':'database_snapshot.sqlite',
        'historical_db_source_sha256': old_db_sha256_before,
        'visualization_source_sha256': snapshot['visualization_source_sha256'],
        'files': evidence,
    }, ensure_ascii=False, indent=2), encoding='utf-8')

    print('PASS: Offline reviewer report generated.', flush=True)
    print('Case count:', len(snapshot['cases']), flush=True)
    print('Selected successful Runs:', sum(bool(c['run_id']) for c in snapshot['cases']), flush=True)
    print('Quantity ID count:', len(all_quantity_ids), flush=True)
    print('Figures:', len(charts), flush=True)
    print('Per-test Markdown:', len(all_reports), flush=True)
    print('Main report:', dest/'evaluation_report.md', flush=True)
    print('Original DB SHA-256 unchanged:', old_db_sha256_before, flush=True)


def _numerical_screen(snapshot, case_roles, evaluate_control_sweep):
    """Replicate existing report-only screening. No engineering PASS/FAIL."""
    numerical = {}
    for test in case_roles:
        cases = [c for c in snapshot['cases'] if c['test_id'] == test and c['run_id']]
        if test.startswith('CTRL-'):
            results = {}
            for qoi in sorted({q['qoi'] for c in cases for q in c['quantities']
                               if q['variant'] == 'DEFAULT' and q['status'] == 'AVAILABLE'}):
                samples = {}
                for c in cases:
                    if not c['role'].startswith('DELTA_'):
                        continue
                    value = next((q['value'] for q in c['quantities']
                        if q['qoi'] == qoi and q['variant'] == 'DEFAULT'
                        and q['status'] == 'AVAILABLE'), None)
                    if value is not None:
                        samples[int(c['role'].split('_')[-1])] = value
                if len(samples) >= 2:
                    results[qoi] = evaluate_control_sweep(samples, qoi)
            if results:
                numerical[test] = results
        elif test.startswith(('GRID-', 'WAKE-', 'GEO-', 'STA-')):
            results = {}
            for qoi in sorted({q['qoi'] for c in cases for q in c['quantities']
                               if q['variant'] == 'DEFAULT' and q['status'] == 'AVAILABLE'}):
                pairs = []
                for c in cases:
                    q = next((q for q in c['quantities']
                        if q['qoi'] == qoi and q['variant'] == 'DEFAULT'
                        and q['status'] == 'AVAILABLE'), None)
                    if q:
                        pairs.append({'case_id':c['case_id'],'value':q['value'],
                                      'quantity_id':q['quantity_id']})
                if len(pairs) >= 2:
                    results[qoi] = {
                        'samples':pairs,
                        'successive_differences':[
                            pairs[i+1]['value'] - pairs[i]['value']
                            for i in range(len(pairs)-1)],
                    }
            if results:
                numerical[test] = results
    return numerical


if __name__ == '__main__':
    try:
        run()
    except Exception as exc:
        print('FAILED: '+str(exc), file=sys.stderr)
        raise
