"""Deterministic read-only exports for review (CSV/Markdown)."""
from __future__ import annotations
import csv
from pathlib import Path

TABLES=('test_definition','planned_case','geometry','flight_condition','numerical_setting',
        'run','case_point','case_run_link','dataset','quantity','acceptance_rule',
        'comparison','assessment','artifact','run_event')


def export_csv(db, destination):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    reports={}
    for table in TABLES:
        rows=db.conn.execute(f'SELECT * FROM {table}').fetchall()
        headers=[col[1] for col in db.conn.execute(f'PRAGMA table_info({table})')]
        output=destination/(table+'.csv')
        with output.open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.writer(f);w.writerow(headers);w.writerows([list(row) for row in rows])
        reports[table]=len(rows)
    return reports


def write_status(db,destination):
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    runs=list(db.conn.execute('SELECT status,COUNT(*) FROM run GROUP BY status'))
    tests=list(db.conn.execute('SELECT t.test_id,t.title,COUNT(p.case_id) n, SUM(CASE WHEN p.execute_policy=\'UNRESOLVED\' THEN 1 ELSE 0 END) missing FROM test_definition t LEFT JOIN planned_case p ON p.test_id=t.test_id GROUP BY t.test_id ORDER BY t.test_id'))
    with destination.open('w',encoding='utf-8') as f:
        f.write('# VSPAERO V&V Study status\n\n')
        f.write('This report tracks completeness. It does not certify aerodynamic accuracy.\n\n')
        f.write('## Runs\n\n')
        for row in runs:f.write(f'- {row[0]}: {row[1]}\n')
        if not runs:f.write('- No VSPAERO Run executed or imported.\n')
        f.write('\n## Test coverage\n\n| Test ID | Title | Planned Cases | Unresolved |\n|---|---|---:|---:|\n')
        for r in tests:f.write(f'| {r[0]} | {r[1]} | {r[2]} | {r[3] or 0} |\n')
        f.write('\nNo engineering acceptance is inferred from unresolved or merely succeeded Runs.\n')
    return destination


def write_evaluation_report(db, destination, *, study_path=None, selected_runs=None, results=None, database_snapshot=None):
    """Use official 13-chapter v0.1 template; insert DB-grounded evidence.

    All placeholders not supported by evidence remain unfilled on purpose.
    The output is a reviewer draft, never an engineering certification.
    """
    from .VVDatabase import REPO_ROOT, digest, postprocess_fingerprint
    import json
    import re
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    template=(REPO_ROOT/'docs'/'vv'/
             'evaluation_report_template.md')
    document=template.read_text(encoding='utf-8')
    study=db.conn.execute('SELECT * FROM study LIMIT 1').fetchone()
    source_info=['**自動集計証拠（未承認）**','',
        f'- Study ID: `{study["study_id"]}`',
        f'- Plan v{study["plan_version"]}: `{study["plan_sha256"]}`',
        f'- SQLite SHA-256 (report-time source): `{digest(database_snapshot or db.path)}`',
        f'- Postprocessor source: `{postprocess_fingerprint(REPO_ROOT)}`',
        '- 全数値はSQLiteのRun/Point/Dataset/Quantity/Comparison/Assessment原記録を参照。',
        '- **使用認定：INSUFFICIENT_EVIDENCE**。実行成功、収束、Validationは別の判定。','']
    if results is not None:
        source_info.extend([f'- Execution manifest SHA-256: `{results.get("execution_manifest_sha256","—")}`',
                            f'- Figure/report code SHA-256: `{results.get("visualization_source_sha256","—")}`',
                            '- Run choice is fixed in `results_snapshot.json`; native values retain their original convention.', ''])
    if database_snapshot is not None:
        source_info.extend([f'- SQLite frozen snapshot: `{Path(database_snapshot).name}`',
                            '- The live DB may change after report generation; use the frozen snapshot for exact replication.', ''])
    chosen=selected_runs or {}
    rows=['| Test ID | Required | Cases | Resolved | SUCCEEDED Run | Comparison | Assessment |',
          '|---|---|---:|---:|---:|---:|---:|']
    for test in db.conn.execute('SELECT test_id,requirement FROM test_definition ORDER BY test_id'):
        t=test['test_id']
        total,resolved=db.conn.execute("""SELECT COUNT(*),SUM(CASE WHEN geometry_id IS NOT NULL AND condition_id IS NOT NULL AND setting_id IS NOT NULL THEN 1 ELSE 0 END)
            FROM planned_case WHERE test_id=?""",(t,)).fetchone()
        ok=db.conn.execute("""SELECT COUNT(DISTINCT run.run_id) FROM run
            JOIN case_point pt ON pt.run_id=run.run_id
            JOIN case_run_link l ON l.point_id=pt.point_id
            JOIN planned_case c ON c.case_id=l.case_id
            WHERE c.test_id=? AND run.status='SUCCEEDED' AND pt.point_status='MATCHED' """,(t,)).fetchone()[0]
        cmp=db.conn.execute('SELECT COUNT(*) FROM comparison WHERE test_id=?',(t,)).fetchone()[0]
        ass=db.conn.execute('SELECT COUNT(*) FROM assessment WHERE test_id=?',(t,)).fetchone()[0]
        rows.append(f'| {t} | {test["requirement"]} | {total} | {resolved or 0} | {ok} | {cmp} | {ass} |')
    sections={'1':source_info,'5':['**自動生成：25試験の証拠充足表**','']+rows}
    chapter={'6':('VER',),'7':('GEO','GRID','WAKE'),
             '8':('XCD',),'9':('VAL',),'10':('STA','CTRL'),
             '11':('DRAG',),'12':('ACC',)}
    for sec,families in chapter.items():
        records=['**自動生成：比較・Assessment証拠一覧（工学的受容を意味しません）**','',
                 '| Kind | ID | Test ID | Value / Status | Source |',
                 '|---|---|---|---|---|']
        for c in db.conn.execute("""SELECT c.*,p.dataset_id AS predicted_dataset,r.dataset_id AS reference_dataset
            FROM comparison c JOIN quantity p ON p.quantity_id=c.predicted_quantity_id
            JOIN quantity r ON r.quantity_id=c.reference_quantity_id ORDER BY c.test_id,c.created_utc"""):
            if c['test_id'].split('-')[0] in families:
                records.append(f'| Comparison | {c["comparison_id"]} | {c["test_id"]} | '
                               f'delta={c["delta_value"]}; {c["verdict"]} | '
                               f'{c["predicted_quantity_id"]}, {c["reference_quantity_id"]} / '
                               f'{c["predicted_dataset"]}, {c["reference_dataset"]} |')
        for a in db.conn.execute('SELECT * FROM assessment ORDER BY test_id,created_utc'):
            if a['test_id'] and a['test_id'].split('-')[0] in families:
                records.append(f'| Assessment | {a["assessment_id"]} | {a["test_id"]} | '
                               f'{a["disposition"]} | {a["evidence_ids_json"]} |')
        if len(records)==4:records.append('| 未取得 | — | — | INSUFFICIENT_EVIDENCE | 外部資料・承認待ち |')
        sections[sec]=records
    if results is not None:
        # Results and figures are a frozen snapshot of explicitly identified Runs.
        # These tables are numerical evidence, not engineering conclusions.
        for sec,families in chapter.items():
            blocks=sections[sec]
            blocks.extend(['', '**自動生成：解析条件・数値結果・グラフ（技術判定なし）**', ''])
            for test in (t['test_id'] for t in db.conn.execute('SELECT test_id FROM test_definition ORDER BY test_id')
                         if t['test_id'].split('-')[0] in families):
                cases=[c for c in results['cases'] if c['test_id']==test]
                blocks.extend([f'### {test} — 解析・結果一覧', '',
                    '| Case ID | 状態 | Run ID | 使用したDataset |',
                    '|---|---|---|---|'])
                for case in cases:
                    ds=', '.join(d['dataset_id'] for d in case['datasets']) or '—'
                    blocks.append(f'| {case["case_id"]} | {case["status"]} | {case["run_id"] or "—"} | {ds} |')
                blocks.extend(['', '| Case ID | QoI | Variant | Value | Unit | Quantity ID |',
                               '|---|---|---|---:|---|---|'])
                blocks.extend(['', '| Case ID | Mesh/History Metric | Value | Run ID |',
                               '|---|---|---:|---|'])
                for case in cases:
                    for metric,value in case['mesh'].items():
                        blocks.append(f'| {case["case_id"]} | {metric} | {value} | {case["run_id"]} |')
                    if case['history']:
                        blocks.append(f'| {case["case_id"]} | history_records | '
                                      f'{len(case["history"])} | {case["run_id"]} |')
                blocks.append('')
                has_quantities=False
                for case in cases:
                    for q in case['quantities']:
                        has_quantities=True
                        label=q['qoi']+((' [native; unverified]') if q['status']!='AVAILABLE' else '')
                        blocks.append(f'| {case["case_id"]} | {label} | {q["variant"]} | '
                                      f'{q["value"]:.9g} | {q["unit"]} | {q["quantity_id"]} |')
                if not has_quantities:
                    blocks.append('| — | — | — | — | — | — |')
                drag_data=[r for r in results.get('drag_polars',[]) if r['test_id']==test]
                if drag_data:
                    blocks.extend(['', '**ParasiteDrag由来の補正結果（係数の物理的採否は未判定）**', '',
                        '| Lift QoI | CL | CDiw | Derived CD | Dataset ID | Quantity ID |',
                        '|---|---:|---:|---:|---|---|'])
                    for r in drag_data:
                        blocks.append(f'| {r["cl_qoi"] or "—"} | {r["cl"] if r["cl"] is not None else "—"} | '
                                      f'{r["cd_native_wake"] if r["cd_native_wake"] is not None else "—"} | '
                                      f'{r["cd"]} | {r["dataset_id"]} | {r["quantity_id"]} |')
                numerical=results.get('numerical',{}).get(test)
                if numerical:
                    blocks.extend(['', '**数値集計（物理的な採否判定ではありません）**', '',
                                   '```json',json.dumps(numerical,ensure_ascii=False,indent=2),'```'])
                figures=[chart for chart in results.get('charts',[]) if chart['test_id']==test]
                for index,chart in enumerate(figures,1):
                    relative=Path(chart['path']).resolve().relative_to(destination.parent.resolve()).as_posix()
                    caption=f'{test} 図{index}: {chart["kind"]} / {chart["qoi"]}'
                    blocks.extend(['',f'![{caption}]({relative})', '',
                        f'**{caption}** — 出典: `{json.dumps(chart["sources"],ensure_ascii=False)}`', ''])
                if not figures:blocks.extend(['', 'グラフ: **生成対象データ不足**', ''])
                failures=[r for r in results.get('execution',[]) if r['case_id'] in
                          {c['case_id'] for c in cases} and r['status'] not in
                          ('SUCCEEDED','REUSED','REPROCESSED','REUSED_POSTPROCESS')]
                for record in failures:
                    blocks.append(f'- 未実行・失敗: `{record["case_id"]}` / {record["status"]}` '
                                  f'— {record.get("reason","実行結果を確認")}' )
                blocks.append('')
    sources=['**自動生成：解析原本ID・SHA-256**','',
             '| Artifact ID | Run ID | Dataset ID | File | SHA-256 |',
             '|---|---|---|---|---|']
    for a in db.conn.execute('SELECT * FROM artifact ORDER BY artifact_id'):
        sources.append(f'| {a["artifact_id"]} | {a["run_id"] or "—"} | '
                       f'{a["dataset_id"] or "—"} | `{a["relative_path"]}` | `{a["sha256"]}` |')
    sections['13']=['**自動出典追跡：本書はレビュー草稿です。**','']+sources
    for number,addition in sections.items():
        m=re.search(r'^## '+number+r'\. .+$',document,re.MULTILINE)
        if not m:raise ValueError('Missing formal report chapter '+number)
        document=document[:m.end()]+ '\n\n'+'\n'.join(addition)+'\n'+document[m.end():]
    destination.write_text(document,encoding='utf-8')
    return destination
