"""Deterministic reviewer exports; the output is not engineering approval.

The immutable database snapshot and selected Run result model are the sole
sources of report data. Review documents are presentation, not source records.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

TABLES = ('test_definition', 'planned_case', 'geometry', 'flight_condition',
          'numerical_setting', 'run', 'case_point', 'case_run_link', 'dataset',
          'quantity', 'acceptance_rule', 'comparison', 'assessment', 'artifact',
          'run_event')
CHAPTER_FAMILIES = {
    '6': ('VER',), '7': ('GEO', 'GRID', 'WAKE'),
    '8': ('XCD',), '9': ('VAL',), '10': ('STA', 'CTRL'),
    '11': ('DRAG',), '12': ('ACC',),
}


def _cell(value):
    """Escape GFM table cells; never permit injected table rows."""
    if value is None:
        return '—'
    s = str(value).replace('\\', '\\\\').replace('|', '\\|')
    return ' '.join(s.splitlines()) or '—'


def _value(value):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return '—'
    if isinstance(value, (float, int)) and not isinstance(value, bool):
        return format(value, '.7g')
    return _cell(value)


def _table(headers, rows, numeric=()):
    numeric = set(numeric)
    output = ['| ' + ' | '.join(_cell(v) for v in headers) + ' |',
              '| ' + ' | '.join('---:' if i in numeric else '---'
                               for i in range(len(headers))) + ' |']
    for row in rows:
        if len(row) != len(headers):
            raise ValueError(f'GFM table width mismatch: {len(row)} != {len(headers)}')
        output.append('| ' + ' | '.join(_cell(v) for v in row) + ' |')
    return '\n'.join(output) + '\n'


def export_csv(db, destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    reports = {}
    for name in TABLES:
        rows = db.conn.execute(f'SELECT * FROM {name} ORDER BY rowid').fetchall()
        headers = [col[1] for col in db.conn.execute(f'PRAGMA table_info({name})')]
        path = destination / (name + '.csv')
        with path.open('w', newline='', encoding='utf-8-sig') as handle:
            writer = csv.writer(handle)
            writer.writerow(headers)
            writer.writerows([list(row) for row in rows])
        reports[name] = len(rows)
    return reports


def write_status(db, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    runs = list(db.conn.execute('SELECT status,COUNT(*) FROM run GROUP BY status'))
    tests = list(db.conn.execute('''SELECT t.test_id,t.title,COUNT(p.case_id) n,
        SUM(CASE WHEN p.execute_policy='UNRESOLVED' THEN 1 ELSE 0 END) missing
        FROM test_definition t LEFT JOIN planned_case p ON p.test_id=t.test_id
        GROUP BY t.test_id ORDER BY t.test_id'''))
    parts = ['# VSPAERO V&V Study status', '',
             'Run success is not numerical or engineering approval.', '', '## Runs', '']
    parts += [f'- {row[0]}: {row[1]}' for row in runs] or ['- No Runs recorded.']
    parts += ['', '## Test coverage', '',
              _table(('Test ID', 'Title', 'Planned Cases', 'Unresolved'),
                     ((r['test_id'], r['title'], r['n'], r['missing'] or 0)
                      for r in tests), (2, 3))]
    destination.write_text('\n'.join(parts), encoding='utf-8')
    return destination


def _test_model(db, snapshot):
    cases_by_test = {}
    if snapshot is not None:
        for case in snapshot['cases']:
            cases_by_test.setdefault(case['test_id'], []).append(case)
    tests = []
    for row in db.conn.execute('SELECT test_id,title,requirement FROM test_definition ORDER BY test_id'):
        tid = row['test_id']
        if snapshot is not None:
            cases = cases_by_test.get(tid, [])
            # A planned Case must never silently disappear from a frozen result set.
            registered = db.conn.execute('SELECT COUNT(*) FROM planned_case WHERE test_id=?',
                                         (tid,)).fetchone()[0]
            if len(cases) != registered:
                raise ValueError(f'Case mismatch between snapshot and frozen DB: {tid}')
            for c in cases:
                run_id = c.get('run_id')
                if (c['status'] == 'SUCCEEDED') != bool(run_id):
                    raise ValueError(f'Inconsistent selected Run status: {c["case_id"]}')
                if run_id:
                    verified = db.conn.execute('''SELECT 1 FROM run r
                        JOIN case_point pt ON pt.run_id=r.run_id
                        JOIN case_run_link link ON link.point_id=pt.point_id
                        WHERE r.run_id=? AND link.case_id=?
                        AND r.status='SUCCEEDED' AND pt.point_status='MATCHED'
                        LIMIT 1''', (run_id, c['case_id'])).fetchone()
                    if not verified:
                        raise ValueError(f'Snapshot Run not in frozen database: {c["case_id"]} {run_id}')
            # This audit checks that displayed quantities exist in the same
            # captured database. Dataset and Case linkage is established when
            # collect_study_results() selects linked matched Run points.
            qids = {q['quantity_id'] for c in cases for q in c['quantities']}
            if qids:
                found = {r[0] for r in db.conn.execute(
                    'SELECT q.quantity_id FROM quantity q JOIN dataset d '
                    'ON d.dataset_id=q.dataset_id WHERE d.study_id=?',
                    (db.conn.execute('SELECT study_id FROM study LIMIT 1').fetchone()[0],))}
                if not qids <= found:
                    raise ValueError('Snapshot Quantity missing from frozen database: ' +
                                     repr(sorted(qids - found)[:5]))
            success = sum(bool(c['status'] == 'SUCCEEDED' and c.get('run_id')) for c in cases)
            unresolved = sum(c['status'] == 'UNRESOLVED' for c in cases)
            total = len(cases)
        else:
            total, unresolved = db.conn.execute('''SELECT COUNT(*),SUM(CASE WHEN
                execute_policy='UNRESOLVED' THEN 1 ELSE 0 END) FROM planned_case
                WHERE test_id=?''', (tid,)).fetchone()
            unresolved = unresolved or 0
            success = db.conn.execute('''SELECT COUNT(DISTINCT l.case_id) FROM case_run_link l
                JOIN case_point pt ON pt.point_id=l.point_id
                JOIN run r ON r.run_id=pt.run_id
                JOIN planned_case pc ON pc.case_id=l.case_id
                WHERE pc.test_id=? AND pt.point_status='MATCHED'
                AND r.status='SUCCEEDED' ''', (tid,)).fetchone()[0]
            cases = []
        failed = total - success - unresolved
        if failed < 0:
            raise ValueError(f'Invalid Case counts for {tid}')
        state = ('未解決' if unresolved == total else
                 '予定Caseの解析成功' if success == total and total else
                 '一部実行済み' if success else
                 '未実行または実行未成功')
        tests.append({'test_id': tid, 'title': row['title'], 'requirement': row['requirement'],
                      'total': total, 'success': success, 'failed': failed,
                      'unresolved': unresolved, 'state': state, 'cases': cases,
                      'engineering': 'NOT_EVALUATED'})
    return tests


def _render_test_report(test, charts, snapshot, destination):
    """Write all selected Quantity values, preserving Method and Dataset identity."""
    path = destination / 'test_reports' / (test['test_id'] + '.md')
    path.parent.mkdir(parents=True, exist_ok=True)
    tid = test['test_id']
    lines = [f'# {tid} — {_cell(test["title"])}', '',
             '> 自動生成レビュー草稿。解析成功はV&V合格を意味しません。', '',
             '## 実施状況', '',
             _table(('計画Case', '採用成功Case', '未成功Case', '未解決Case', '工学評価'),
                    [(test['total'], test['success'], test['failed'],
                      test['unresolved'], test['engineering'])], (0, 1, 2, 3)),
             '', '## Case別状況', '']
    cases = test['cases']
    if cases:
        exec_by_case = {}
        for record in snapshot.get('execution', []):
            exec_by_case.setdefault(record.get('case_id'), []).append(record)
        rows = []
        for c in cases:
            attempts = exec_by_case.get(c['case_id'], [])
            latest = attempts[-1] if attempts else {}
            reason = '; '.join(c.get('issues', []))
            if not c.get('run_id'):
                reason = reason or str(latest.get('reason') or c['status'])
            rows.append((c['case_id'], c['role'], c['status'], c.get('run_id') or '—',
                         len(c['datasets']), latest.get('status', 'NOT_RECORDED'), reason or '—'))
        lines.append(_table(('Case ID', 'Role', '状態', '採用Run ID', 'Dataset数',
                             'ワークフロー状態', '未実行・不成功理由'), rows, (4,)))
    else:
        lines.append('Caseの固定スナップショットはありません。\n')
    lines += ['', '## 全Quantity（元数値は quantities.csv / SQLiteを参照）', '']
    quantity_rows = []
    for case in cases:
        methods = {d['dataset_id']: d['method_id'] for d in case['datasets']}
        for q in case['quantities']:
            quantity_rows.append((case['case_id'], q['qoi'], q['variant'], q['component'],
                                  _value(q['value']), q['unit'], q['status'],
                                  methods.get(q['dataset_id'], 'UNKNOWN'),
                                  q['dataset_id'], q['quantity_id']))
    lines.append(_table(('Case ID', 'QoI', 'Variant', 'Component', 'Value', 'Unit', 'Status',
                         'Method', 'Dataset ID', 'Quantity ID'), quantity_rows, (4,))
                 if quantity_rows else '取得済みQuantityなし。\n')
    lines += ['', '## メッシュ・履歴', '']
    diagnostics = []
    for c in cases:
        for metric, value in c.get('mesh', {}).items():
            diagnostics.append((c['case_id'], metric, _value(value), c.get('run_id')))
        if c.get('history'):
            diagnostics.append((c['case_id'], 'history_records', len(c['history']), c.get('run_id')))
    lines.append(_table(('Case ID', 'Metric', 'Value', 'Run ID'), diagnostics, (2,))
                 if diagnostics else 'メッシュ・履歴情報なし。\n')
    drag = [r for r in snapshot.get('drag_polars', []) if r['test_id'] == tid]
    if drag:
        lines += ['', '## ParasiteDrag補正（未承認）', '',
                  _table(('Lift QoI', 'CL', 'Native CDiw', 'Derived CD', 'Dataset ID', 'Quantity ID'),
                         [(r['cl_qoi'], _value(r['cl']), _value(r['cd_native_wake']),
                           _value(r['cd']), r['dataset_id'], r['quantity_id']) for r in drag], (1, 2, 3))]
    numerical = snapshot.get('numerical', {}).get(tid)
    if numerical:
        lines += ['', '## 数値スクリーニング（正式な工学判定ではない）', '',
                  '```json', json.dumps(numerical, ensure_ascii=False, indent=2), '```', '']
    related = [x for x in charts if x['test_id'] == tid]
    lines += ['', f'## 図表（{len(related)}点）', '']
    if related:
        for fig in related:
            # This report is one directory below the main document.
            image = '../' + fig['path']
            lines += [f'### {fig["figure_id"]} — {_cell(fig["qoi"])}', '',
                      f'![{fig["figure_id"]}]({image})', '',
                      f'{fig["kind"]}。横軸：{_cell(fig["x_label"])}。'
                      f'縦軸：{_cell(fig["y_label"])} [{_cell(fig["y_unit"])}]。'
                      f'{_cell(fig.get("note") or "")}\n',
                      '詳細出典：[`figure_manifest.csv`](../figure_manifest.csv)。', '']
    else:
        lines += ['図：**生成対象データ不足**（適切な同一QoI／方法の比較系列なし）。', '']
    related_cmp = [x for x in snapshot.get('comparisons', []) if x['test_id'] == tid]
    lines += ['', '## 参照比較（採用Runとの対応が確認されたもの）', '']
    lines.append(_table(('Comparison ID', 'QoI', 'Predicted', 'Reference', 'Delta'),
                        [(x['comparison_id'], x['qoi'], _value(x['predicted']),
                          _value(x['reference']), _value(x['delta'])) for x in related_cmp], (2, 3, 4))
                 if related_cmp else '該当する参照比較なし。正式なPASS/FAILは未評価。\n')
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return path


def write_evaluation_report(db, destination, *, study_path=None, selected_runs=None,
                            results=None, database_snapshot=None):
    """Populate explicit markers in the official 13-chapter template.

    The caller must supply a *frozen database* when database_snapshot is set.
    No SQL query accesses the mutable live database in the integrated workflow.
    """
    from .VVDatabase import REPO_ROOT, digest, postprocess_fingerprint

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    template = REPO_ROOT / 'docs/vv/templates/evaluation_report.md'
    document = template.read_text(encoding='utf-8')
    study = db.conn.execute('SELECT * FROM study LIMIT 1').fetchone()
    if study is None:
        raise ValueError('No registered Study in report database')
    tests = _test_model(db, results)
    total = sum(t['total'] for t in tests)
    success = sum(t['success'] for t in tests)
    unresolved = sum(t['unresolved'] for t in tests)
    sources = ['**自動集計証拠（未承認・技術評価なし）**', '',
               _table(('項目', '固定された値'), [
                   ('Study ID', study['study_id']),
                   ('対象OpenVSP版', study['target_openvsp_version']),
                   ('計画書版', study['plan_version']),
                   ('計画書SHA-256', study['plan_sha256']),
                   ('計画Case', total),
                   ('採用成功Case', success),
                   ('成功Runなし・未解決を除く', total - success - unresolved),
                   ('未解決Case', unresolved),
                   ('技術評価', 'NOT_EVALUATED'),
                   ('使用認定', 'INSUFFICIENT_EVIDENCE'),
                   ('選択方針', results.get('selection_policy', 'DBに記録された成功Run')
                    if results else 'DBに記録された成功Run'),
               ])]
    if results is not None:
        sources += ['採用Runに紐付かない過去のComparison：' +
                    str(len(results.get('excluded_historical_comparison_ids', []))) + '件（本文の評価から除外）。',
                    '', '全数値は固定された採用RunのQuantityを参照し、'
                    '過去のRun総数および解析成功回数とは区別しています。', '',
                    f'- Execution manifest SHA-256: `{results.get("execution_manifest_sha256", "未記録")}`',
                    f'- Visualization source SHA-256: `{results.get("visualization_source_sha256", "未記録")}`',
                    '- Run choice is frozen in `results_snapshot.json`.', '']
    if database_snapshot is not None:
        sources += [f'- SQLite frozen snapshot: `{Path(database_snapshot).name}`',
                    f'- SQLite frozen SHA-256: `{digest(database_snapshot)}`',
                    '- The live DB can change after generation; the frozen database is the report source.', '']

    state_rows = [(f'[{t["test_id"]}](test_reports/{t["test_id"]}.md)' if results is not None
                   else t['test_id'], t['title'], t['requirement'], t['total'],
                   t['success'], t['failed'], t['unresolved'], t['state'], t['engineering'])
                  for t in tests]
    states = ['**計画・解析実行・工学判定は独立した状態です。**', '',
              _table(('Test ID', '試験名', '区分', '計画', '成功', '未成功',
                      '未解決', '実行状況', '工学評価'), state_rows, (3, 4, 5, 6)),
              '', '成功はSUCCEEDEDの採用Runが存在するCase数であり、V&VのPASSではありません。', '']
    by_test = {t['test_id']: t for t in tests}
    replacements = {'1': '\n'.join(sources), '5.1': '\n'.join(states)}
    chart_list = results.get('charts', []) if results is not None else []
    for number, families in CHAPTER_FAMILIES.items():
        content = ['**自動生成：採用Runの解析証拠（工学的受容を意味しません）**', '']
        for tid, t in by_test.items():
            if tid.split('-')[0] not in families:
                continue
            content += [f'### {tid} — {_cell(t["title"])}', '',
                        f'計画 {t["total"]} Case、採用成功 {t["success"]} Case、'
                        f'未成功 {t["failed"]} Case、未解決 {t["unresolved"]} Case。'
                        f'実行状態：**{t["state"]}**。工学評価：**NOT_EVALUATED**。', '']
            if results is not None:
                used_ids = [c['case_id'] for c in t['cases'] if c.get('run_id')]
                if used_ids:
                    content += ['採用成功Case：' + ', '.join('`' + _cell(cid) + '`' for cid in used_ids), '']
                content += [f'全数値・メッシュ・履歴・出典：'
                            f'[{tid}の詳細レポート](test_reports/{tid}.md)。', '']
                figures = [c for c in chart_list if c['test_id'] == tid]
                important = [c for c in figures if c.get('display_priority') == 'PRIMARY'][:3]
                for chart in important:
                    content += [f'![{_cell(chart["figure_id"])}]({chart["path"]})', '',
                                f'**{chart["figure_id"]}**：{_cell(chart["qoi"])}。'
                                f'横軸：{_cell(chart["x_label"])}。'
                                f'縦軸：{_cell(chart["y_label"])} [{_cell(chart["y_unit"])}]。'
                                '出典は[`figure_manifest.csv`](figure_manifest.csv)を参照。', '']
                if not important:
                    content += ['グラフ: **生成対象データ不足**（または判定用途に適した図なし）。', '']
                if len(figures) > len(important):
                    content += [f'その他の図 {len(figures) - len(important)}点は試験別詳細レポートを参照。', '']
        replacements[number] = '\n'.join(content)

    arts = db.conn.execute('SELECT COUNT(*) FROM artifact').fetchone()[0]
    replacements['13'] = '\n'.join([
        '**自動生成：固定SQLite内の原本証拠一覧（未承認）**', '',
        f'- 固定SQLite内に登録されたArtifact：**{arts}件**。',
        '- 完全なArtifact ID、Dataset ID、相対パス、SHA-256は'
        ' [`tables/artifact.csv`](tables/artifact.csv)から確認してください。',
        '- 固定SQLite自身、Markdown、図表のハッシュは'
        ' [`evidence_manifest.json`](evidence_manifest.json)に別途記載します。',
        '- 生成報告書からの自己参照を防ぐため、固定DB自身はその内部Artifactに含めません。',
    ])

    # Known document-control metadata must not remain an unfilled template.
    db_hash = digest(database_snapshot) if database_snapshot is not None else None
    meta = {
        'NOTICE': '**自動生成レビュー草稿／工学評価未承認（テンプレート v0.1）**',
        'STUDY': _cell(study['study_id']),
        'DB': (f'`{Path(database_snapshot).name}` / SHA-256 `{db_hash}`'
               if db_hash else '未作成（単独レポート出力）'),
        'POST': '後処理コード指紋 `' + postprocess_fingerprint(REPO_ROOT) + '`',
        'PLAN_HASHES': ('計画 `' + study['plan_sha256'] +
                        '` / Study設定 `' + study['config_sha256'] + '`'),
        'DB_SHA': ('`' + db_hash + '`' if db_hash else '未作成'),
    }
    for key, value in meta.items():
        marker = '<!-- VV:META:' + key + ' -->'
        if document.count(marker) != 1:
            raise ValueError('Missing/duplicate metadata marker: ' + marker)
        document = document.replace(marker, value, 1)
    for section, content in replacements.items():
        marker = f'<!-- VV:CHAPTER:{section} -->'
        if document.count(marker) != 1:
            raise ValueError(f'Missing/duplicate report template marker: {marker}')
        document = document.replace(marker, content, 1)
    if '<!-- VV:CHAPTER:' in document or '<!-- VV:META:' in document:
        raise ValueError('Unpopulated report template marker')
    destination.write_text(document, encoding='utf-8')
    if results is not None:
        for test in tests:
            _render_test_report(test, chart_list, results, destination.parent)
    return destination
