"""Offline V&V figures. Plot positions follow measured inputs, never Case order.

Figures convey evidence; they do not imply numerical or engineering approval.
The complete Quantity collection remains available in CSV/SQLite even when a
particular Quantity is not scientifically useful as a standalone figure.
"""
from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path

# Deliberately small, reviewed *display* selection. This does not filter data.
PRIMARY_QOI = {
    'VER': ('CLi', 'CDi', 'CMiy'),
    'GEO': ('CLi', 'CDiw', 'CMiy'),
    'GRID': ('CLi', 'CDiw', 'CMiy'),
    'WAKE': ('CLi', 'CDiw', 'CMiy'),
    'XCD': ('CLi', 'CDiw', 'CMiy'),
    'VAL': ('CLi', 'CDiw', 'CMiy'),
    'STA': ('CMm', 'CMl', 'CMn', 'CMiy', 'CLi', 'CDiw'),
    'CTRL': ('CMiy', 'CMix', 'CMiz'),
    'DRAG': ('CD_A', 'CD_B', 'CDiw'),
}


def _numeric_x(test, case):
    role = case['role']
    if test.startswith('CTRL-') and re.fullmatch(r'DELTA_[+-]?\d+(?:\.\d+)?', role):
        return float(role.split('_', 1)[1]), 'Control command [deg]'
    if test == 'WAKE-01' and re.fullmatch(r'ITER_\d+', role):
        return float(role[5:]), 'Requested wake iterations [1]'
    if test == 'WAKE-02' and re.fullmatch(r'NODES_\d+', role):
        return float(role[6:]), 'Requested wake nodes [1]'
    if test.startswith(('GRID-', 'VER-')):
        n = case.get('mesh', {}).get('triangles_total')
        if isinstance(n, (int, float)) and math.isfinite(n) and n > 0:
            return float(n), 'Observed mesh triangles [1]'
    if test.startswith(('XCD-', 'VAL-')):
        condition = case.get('flight_condition', {})
        name = 'beta_deg' if test == 'VAL-02' else 'alpha_deg'
        x = condition.get(name)
        if isinstance(x, (int, float)) and math.isfinite(x):
            return float(x), name + ' [deg]'
    return None, 'Case role (categorical)'


def _figure_name(test, kind, key):
    token = hashlib.sha256(repr((test, kind, key)).encode('utf-8')).hexdigest()[:12]
    slug = re.sub('[^A-Za-z0-9_-]', '_', str(key[0]))[:40]
    return f'{test}_{kind}_{slug}_{token}.png'


def _is_number(v):
    return isinstance(v, (float, int)) and not isinstance(v, bool) and math.isfinite(v)


def create_study_plots(snapshot, destination):
    """Return portable figure metadata, including all source Quantity/Dataset IDs.

    Multi-method and multi-dataset observations are distinct series. Never use
    dictionary assignment to collapse equal QoI names in different datasets.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    charts = []
    cases_by_test = {}
    planned_counts = {}
    for case in snapshot['cases']:
        planned_counts[case['test_id']] = planned_counts.get(case['test_id'], 0) + 1
        if case.get('run_id'):
            cases_by_test.setdefault(case['test_id'], []).append(case)

    def add(test, kind, qoi, fig, sources, key, *, x_label, x_unit='', y_label='', y_unit='', note=''):
        name = _figure_name(test, kind, key)
        image = destination / name
        fig.savefig(image, dpi=150)
        plt.close(fig)
        charts.append({
            'figure_id': 'FIG-' + test + '-' + hashlib.sha256(name.encode()).hexdigest()[:10],
            'test_id': test, 'kind': kind, 'qoi': qoi,
            'path': (Path(destination.name) / name).as_posix(),
            'sources': sources, 'x_label': x_label, 'x_unit': x_unit,
            'y_label': y_label, 'y_unit': y_unit, 'note': note,
            'display_priority': 'PRIMARY' if kind != 'HISTORY' else 'DIAGNOSTIC',
        })

    for test, cases in sorted(cases_by_test.items()):
        family = test.split('-')[0]
        selected = PRIMARY_QOI.get(family, ())
        groups = {}
        for case in cases:
            methods = {d['dataset_id']: d['method_id'] for d in case.get('datasets', [])}
            for q in case['quantities']:
                if q['qoi'] not in selected or not _is_number(q.get('value')):
                    continue
                if q.get('status') not in ('AVAILABLE', 'NATIVE_ONLY'):
                    continue
                key = (q['qoi'], q['variant'], q['component'], q['unit'],
                       q['status'], methods.get(q['dataset_id'], 'UNKNOWN'))
                groups.setdefault(key, {}).setdefault(case['case_id'], []).append((case, q))

        for key, observations in sorted(groups.items(), key=lambda x: repr(x[0])):
            # Datasets of the same method can also be multiple per Case: keep
            # every occurrence in a separately labelled replicate series.
            copies = max(len(items) for items in observations.values())
            for replicate in range(copies):
                ordered = [observations[c['case_id']][replicate]
                           for c in cases if len(observations.get(c['case_id'], [])) > replicate]
                if len(ordered) < 2:
                    continue
                xs_and_labels = [_numeric_x(test, case) for case, _ in ordered]
                numeric = all(x is not None for x, _ in xs_and_labels)
                # Repeated physical X values do not define a single-valued curve.
                if numeric and len({x for x, _ in xs_and_labels}) != len(ordered):
                    numeric = False
                fig, ax = plt.subplots(figsize=(8.5, 4.5))
                values = [float(q['value']) for _, q in ordered]
                if numeric:
                    xs = [item[0] for item in xs_and_labels]
                    ordered_plot = sorted(zip(xs, values), key=lambda x: x[0])
                    xplot, yplot = zip(*ordered_plot)
                    # Do not interpolate through absent/failed planned Cases.
                    if len(ordered) == planned_counts[test]:
                        ax.plot(xplot, yplot, '-o', linewidth=1.3)
                    else:
                        ax.scatter(xplot, yplot)
                    ax.set_xlabel(xs_and_labels[0][1])
                else:
                    ax.scatter(range(len(ordered)), values)
                    ax.set_xticks(range(len(ordered)), [c['role'] for c, _ in ordered],
                                  rotation=30, ha='right')
                    ax.set_xlabel('Case role (categorical)')
                ax.set_ylabel(f'{key[0]} [{key[3]}]' +
                              (' (native, unverified)' if key[4] == 'NATIVE_ONLY' else ''))
                ax.set_title(f'{test}: {key[0]} / {key[1]} / {key[2]}\n'
                             f'Method: {key[5]} / {key[4]}' +
                             (f' / series {replicate + 1}' if copies > 1 else ''))
                ax.grid(True, alpha=0.25)
                fig.tight_layout()
                sources = [{'case_id': c['case_id'], 'run_id': c['run_id'],
                            'quantity_id': q['quantity_id'], 'dataset_id': q['dataset_id'],
                            'method_id': key[5]} for c, q in ordered]
                add(test, 'CASE_SWEEP', key[0], fig, sources, key + (replicate,),
                    x_label=xs_and_labels[0][1] if numeric else 'Case role (categorical)',
                    x_unit='deg' if numeric and 'deg' in xs_and_labels[0][1] else '',
                    y_label=key[0], y_unit=key[3],
                    note=f'{planned_counts[test] - len(ordered)} Case(s) omitted: no adopted Run or matching Quantity.'
                         if len(ordered) < planned_counts[test] else '')

        if test.startswith('WAKE'):
            for case in cases:
                history = case.get('history', [])
                if not history:
                    continue
                fields = [name for name in ('CLi', 'CDiw', 'CMiy', 'CLiw', 'CDi',
                                             'l2_residual', 'max_residual')
                          if name in history[0] and any(_is_number(p.get(name)) for p in history)]
                if not fields:
                    continue
                fig, axes = plt.subplots(len(fields), 1, figsize=(8, 2.3 * len(fields)),
                                         sharex=True)
                if len(fields) == 1:
                    axes = [axes]
                for ax, name in zip(axes, fields):
                    points = [(p['iteration'], p[name]) for p in history
                              if _is_number(p.get('iteration')) and _is_number(p.get(name))]
                    if points:
                        ax.plot(*zip(*points), '-o', ms=2)
                    ax.set_ylabel(name)
                    ax.grid(True, alpha=.25)
                axes[-1].set_xlabel('Native iteration [1]')
                fig.suptitle(f'{case["case_id"]}: native history (no convergence assertion)')
                fig.tight_layout()
                add(test, 'HISTORY', ', '.join(fields), fig,
                    [{'case_id': case['case_id'], 'run_id': case['run_id']}],
                    (case['case_id'],), x_label='Native iteration', y_label=', '.join(fields),
                    note='Native histories; stopping criteria not independently assessed.')

    for test in ('DRAG-01', 'DRAG-02'):
        items = [dict(r) for r in snapshot.get('drag_polars', [])
                 if r['test_id'] == test and _is_number(r.get('cl')) and _is_number(r.get('cd'))]
        if not items or len({r['cl_qoi'] for r in items}) != 1:
            continue
        items.sort(key=lambda x: x['cl'])
        fig, ax = plt.subplots(figsize=(8, 4.5))
        if all(_is_number(r.get('cd_native_wake')) for r in items):
            ax.plot([r['cl'] for r in items], [r['cd_native_wake'] for r in items],
                    '-o', label='Native CDiw')
        ax.plot([r['cl'] for r in items], [r['cd'] for r in items],
                '-s', label='Corrected ' + items[0]['qoi'])
        ax.set_xlabel(items[0]['cl_qoi'] + ' [1]')
        ax.set_ylabel('CD [1]')
        ax.set_title(f'{test}: derived drag polar (not approved)')
        ax.grid(True, alpha=.25)
        ax.legend()
        fig.tight_layout()
        add(test, 'DRAG_POLAR', items[0]['qoi'], fig,
            [{'dataset_id': r['dataset_id'], 'quantity_id': r['quantity_id'],
              'polar_dataset_id': r['polar_dataset_id'],
              'cd0_dataset_id': r['cd0_dataset_id']} for r in items],
            (items[0]['qoi'], items[0]['cl_qoi']), x_label=items[0]['cl_qoi'],
            x_unit='1', y_label=items[0]['qoi'], y_unit='1')

    for test in sorted({x['test_id'] for x in snapshot.get('comparisons', [])}):
        for qoi in sorted({x['qoi'] for x in snapshot['comparisons'] if x['test_id'] == test}):
            pairs = [x for x in snapshot['comparisons'] if x['test_id'] == test
                     and x['qoi'] == qoi and _is_number(x.get('predicted'))
                     and _is_number(x.get('reference'))]
            if not pairs or len({x['unit'] for x in pairs}) != 1:
                continue
            fig, ax = plt.subplots(figsize=(8, 4))
            indexes = list(range(len(pairs)))
            ax.scatter(indexes, [x['predicted'] for x in pairs], marker='o', label='VSPAERO')
            ax.scatter(indexes, [x['reference'] for x in pairs], marker='s', label='Reference')
            ax.set_xticks(indexes, [x['comparison_id'][-8:] for x in pairs], rotation=35)
            ax.set_xlabel('Comparison ID (categorical)')
            ax.set_ylabel(f'{qoi} [{pairs[0]["unit"]}]')
            ax.set_title(f'{test}: reference comparison (condition review required)')
            ax.legend()
            ax.grid(True, alpha=.25)
            fig.tight_layout()
            add(test, 'REFERENCE_COMPARISON', qoi, fig,
                [{'comparison_id': x['comparison_id'],
                  'predicted_dataset': x['predicted_dataset'],
                  'reference_dataset': x['reference_dataset']} for x in pairs],
                (qoi, tuple(x['comparison_id'] for x in pairs)),
                x_label='Comparison ID (categorical)', y_label=qoi, y_unit=pairs[0]['unit'])
    return charts
