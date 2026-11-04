"""Offline regression for portable, readable, provenance-preserving V&V reports."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from pathlib import Path

import pytest

from openvsppy.vv.VVDatabase import VVDatabase, digest
from openvsppy.vv.VVPlots import _figure_name, create_study_plots
from openvsppy.vv.VVPostprocess import collect_study_results
from openvsppy.vv.VVReport import _cell, _table, write_evaluation_report


def sample(test, role, n, *, method='POLAR_NATIVE', extra=False, value=None):
    cid = test + '-' + role
    ds = 'DS-' + cid + ('-extra' if extra else '')
    return {'test_id': test, 'case_id': cid, 'role': role, 'status': 'SUCCEEDED',
            'run_id': 'RUN-' + cid, 'datasets': [{'dataset_id': ds, 'method_id': method}],
            'quantities': [{'dataset_id': ds, 'quantity_id': 'Q-' + ds, 'qoi': 'CMiy',
                            'variant': 'DEFAULT', 'component': 'TOTAL', 'unit': '1',
                            'status': 'AVAILABLE', 'value': n if value is None else value}],
            'mesh': {}, 'history': [], 'issues': []}


def _links(md):
    return re.findall(r'(?<!!)\[[^\]]+\]\(([^)]+)\)|!\[[^\]]+\]\(([^)]+)\)', md)


def test_gfm_cell_escaping_and_width():
    text = _table(('name', 'comment'), [('a|b', 'line1\nline2\\name')])
    assert 'a\\|b' in text
    assert 'line1 line2\\\\name' in text
    parser = pytest.importorskip('markdown_it').MarkdownIt('default').enable('table')
    tokens = parser.parse(text)
    assert len([t for t in tokens if t.type == 'td_open']) == 2
    assert _cell('') == '—'


def test_marker_missing_fails_clearly(study, tmp_path, monkeypatch):
    from openvsppy.vv.VVDatabase import REPO_ROOT
    template = REPO_ROOT / 'docs/vv/templates/evaluation_report.md'
    original_read = Path.read_text

    def read_with_removed_marker(path, *args, **kwargs):
        text = original_read(path, *args, **kwargs)
        if Path(path) == template:
            return text.replace('<!-- VV:CHAPTER:5.1 -->', '', 1)
        return text

    monkeypatch.setattr(Path, 'read_text', read_with_removed_marker)
    _, db = study
    with pytest.raises(ValueError, match='Missing/duplicate report template marker'):
        write_evaluation_report(db, tmp_path/'report.md')


def test_25_case_reports_and_unapproved_study(study, tmp_path):
    root, db = study
    snapshot = collect_study_results(db, root)
    snapshot['charts'] = []
    output = tmp_path/'evaluation_report.md'
    write_evaluation_report(db, output, results=snapshot)
    main = output.read_text(encoding='utf-8')
    assert len(re.findall(r'^## \d+\. ', main, re.M)) == 13
    assert len(list((tmp_path/'test_reports').glob('*.md'))) == 25
    assert 'GRID-01' in main and '| 未実施 | NOT_EVALUATED |' not in main
    assert '工学評価：**PASS**' not in main
    assert '<!-- VV:' not in main
    assert '<!-- VV:' not in (tmp_path/'test_reports/GRID-01.md').read_text(encoding='utf-8')
    for match in _links(main):
        url = next(v for v in match if v)
        if url.startswith(('test_reports/', 'figures/')):
            assert (tmp_path/url).exists(), url


def test_numeric_ctrl_positions_are_physical(monkeypatch, tmp_path):
    import matplotlib.axes
    positions = []
    old = matplotlib.axes.Axes.plot

    def track(ax, *args, **kwargs):
        if len(args) >= 2 and len(args[0]) == 7:
            positions.append(tuple(args[0]))
        return old(ax, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, 'plot', track)
    angles = (-10, -5, -2, 0, 2, 5, 10)
    cases = [sample('CTRL-01', f'DELTA_{a:+d}', i) for i, a in enumerate(angles)]
    charts = create_study_plots({'cases': cases, 'comparisons': []}, tmp_path/'figures')
    assert positions == [tuple(float(a) for a in angles)]
    assert len(charts) == 1 and charts[0]['x_label'] == 'Control command [deg]'
    assert charts[0]['path'].startswith('figures/')


def test_incomplete_numeric_series_not_interpolated(monkeypatch, tmp_path):
    import matplotlib.axes
    seen = []
    old = matplotlib.axes.Axes.plot

    def track(ax, *args, **kwargs):
        seen.append(args)
        return old(ax, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, 'plot', track)
    angles = (-10, -5, -2, 0, 2, 5, 10)
    cases = [sample('CTRL-02', f'DELTA_{a:+d}', i) for i, a in enumerate(angles)]
    # One succeeded Case has no CMiy; it must not become an interpolated value.
    cases[3]['quantities'] = []
    charts = create_study_plots({'cases': cases, 'comparisons': []}, tmp_path/'figures')
    assert len(charts) == 1 and '1 Case(s) omitted' in charts[0]['note']
    assert not seen


def test_missing_case_run_does_not_make_false_continuity(monkeypatch, tmp_path):
    import matplotlib.axes
    calls = []
    previous = matplotlib.axes.Axes.plot

    def spy(ax, *args, **kwargs):
        calls.append(args)
        return previous(ax, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, 'plot', spy)
    angles = (-10, -5, -2, 0, 2, 5, 10)
    cases = [sample('CTRL-03', f'DELTA_{a:+d}', i) for i, a in enumerate(angles)]
    cases[3]['run_id'] = None
    cases[3]['status'] = 'NO_SUCCESSFUL_RUN'
    plots = create_study_plots({'cases': cases, 'comparisons': []}, tmp_path/'figures')
    assert len(plots) == 1
    assert '1 Case(s) omitted' in plots[0]['note']
    assert not calls


def test_categorical_geo_does_not_connect_points(monkeypatch, tmp_path):
    import matplotlib.axes
    old = matplotlib.axes.Axes.plot

    def forbidden(ax, *args, **kwargs):
        pytest.fail('GEO categorical comparison must not use a connected line')

    monkeypatch.setattr(matplotlib.axes.Axes, 'plot', forbidden)
    cases = [sample('GEO-04', role, i) for i, role in enumerate(('THIN', 'THICK'))]
    charts = create_study_plots({'cases': cases, 'comparisons': []}, tmp_path/'figures')
    assert len(charts) == 1 and charts[0]['x_label'] == 'Case role (categorical)'
    monkeypatch.setattr(matplotlib.axes.Axes, 'plot', old)


def test_identical_qoi_multiple_methods_and_datasets_are_preserved(tmp_path):
    cases = [sample('GEO-04', role, i) for i, role in enumerate(('THIN', 'THICK'))]
    for i, case in enumerate(cases):
        supplement = sample('GEO-04', case['role'], i + 10, method='POLAR_REPROCESS', extra=True)
        case['datasets'].extend(supplement['datasets'])
        case['quantities'].extend(supplement['quantities'])
    charts = create_study_plots({'cases': cases, 'comparisons': []}, tmp_path/'figures')
    assert len(charts) == 2
    assert {x['method_id'] for c in charts for x in c['sources']} == {'POLAR_NATIVE', 'POLAR_REPROCESS'}
    assert len({c['path'] for c in charts}) == 2


def test_same_method_duplicate_dataset_gets_separate_series(tmp_path):
    cases = [sample('GEO-04', role, i) for i, role in enumerate(('THIN', 'THICK'))]
    for i, case in enumerate(cases):
        supplement = sample('GEO-04', case['role'], i + 5, extra=True)
        case['datasets'].extend(supplement['datasets'])
        case['quantities'].extend(supplement['quantities'])
    charts = create_study_plots({'cases': cases, 'comparisons': []}, tmp_path/'figures')
    assert len(charts) == 2
    assert len({q['quantity_id'] for c in charts for q in c['sources']}) == 4


def test_figure_names_stable_and_collision_resistant():
    prefix = 'a' * 75
    a = _figure_name('CTRL-01', 'CASE_SWEEP', (prefix + 'X', 'DEFAULT'))
    b = _figure_name('CTRL-01', 'CASE_SWEEP', (prefix + 'Y', 'DEFAULT'))
    assert a != b and a == _figure_name('CTRL-01', 'CASE_SWEEP', (prefix + 'X', 'DEFAULT'))


def test_nonfinite_qoi_not_plotted(tmp_path):
    cases = [sample('GEO-04', role, i, value=float('nan')) for i, role in enumerate(('THIN', 'THICK'))]
    charts = create_study_plots({'cases': cases, 'comparisons': []}, tmp_path/'figures')
    assert not charts


def test_report_data_numeric_values_round_for_display_only():
    from openvsppy.vv.VVReport import _value
    x = 0.00101790435
    assert _value(x) != str(x)
    assert x == 0.00101790435


def test_template_metadata_is_not_falsely_approved(study, tmp_path):
    _, db = study
    file = tmp_path/'report.md'
    write_evaluation_report(db, file)
    text = file.read_text(encoding='utf-8')
    assert 'VSP3505-VV-001' in text
    assert '| 報告対象Study ID | 《未設定' not in text
    assert '| 対象DBスナップショット | 未作成' in text
    assert '未発行' in text


def test_portable_links_after_move(tmp_path):
    # Source folder is immaterial: only POSIX relative paths are serialized.
    original = tmp_path/'a'; cloned = tmp_path/'b'
    original.mkdir(); (original/'figures').mkdir()
    (original/'figures/demo.png').write_bytes(b'test')
    cases = [sample('GEO-04', role, i) for i, role in enumerate(('THIN', 'THICK'))]
    charts = create_study_plots({'cases': cases, 'comparisons': []}, original/'figures')
    shutil.copytree(original, cloned)
    assert all((cloned/chart['path']).exists() for chart in charts)
    assert all(not Path(chart['path']).is_absolute() for chart in charts)
