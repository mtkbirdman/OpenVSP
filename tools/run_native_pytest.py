"""Explicit, subprocess-isolated runner for expensive OpenVSP native pytest cases.

Each pytest node runs in its own subprocess and writes a persistent log. This
prevents one native crash or timeout from stopping the remaining test cases.
Use pytest -m directly for ordinary interactive debugging (without deadlines).
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import io
import json
import os
import platform
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
MARKERS = {
    'native_openvsp': 'native_openvsp',
    'native_vspaero': 'native_vspaero',
    'native_solver': 'native_solver',
    'gui': 'gui',
    'all-native': 'native_openvsp or native_vspaero or native_solver',
}


def collect_nodeids(marker):
    """Collect only; never import the OpenVSP native extension at collection time."""
    import pytest

    class Collector:
        def __init__(self):
            self.nodeids = []

        def pytest_collection_finish(self, session):
            self.nodeids = [item.nodeid for item in session.items]

    collector = Collector()
    # Pytest writes a collection summary; it is intentionally not a test log.
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        result = pytest.main(['--collect-only', '-q', '-m', marker, 'tests'], plugins=[collector])
    if int(result) not in (0, 5):
        raise RuntimeError(f'Pytest native collection failed: exit code {result}')
    return collector.nodeids


def stop_process_tree(proc):
    """Terminate solver descendants as well as the supervising pytest process."""
    if proc.poll() is not None:
        return
    if os.name == 'nt':
        # taskkill /T also stops the external VSPAERO solver started by OpenVSP.
        subprocess.run(['taskkill', '/PID', str(proc.pid), '/T', '/F'],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        if proc.poll() is None:
            proc.kill()
    else:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except OSError:
            proc.kill()
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()


def execute_command(command, log_path, timeout, heartbeat=60, env=None, cwd=ROOT):
    """Run one test command in an isolated process with a hard time budget."""
    start = time.monotonic()
    creationflags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == 'nt' else 0
    with open(log_path, 'wb') as log:
        proc = subprocess.Popen(
            command, cwd=cwd, env=env, stdin=subprocess.DEVNULL,
            stdout=log, stderr=subprocess.STDOUT,
            start_new_session=(os.name != 'nt'), creationflags=creationflags,
        )
        while True:
            elapsed = time.monotonic() - start
            remaining = timeout - elapsed
            if remaining <= 0:
                stop_process_tree(proc)
                return 'TIMEOUT', time.monotonic() - start, proc.returncode
            try:
                code = proc.wait(timeout=min(heartbeat, remaining))
                return ('PASSED' if code == 0 else 'FAILED'), time.monotonic() - start, code
            except subprocess.TimeoutExpired:
                print(f'  [ACTIVE] {time.monotonic()-start:.0f}s elapsed', flush=True)
            except BaseException:
                stop_process_tree(proc)
                raise


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--marker', choices=sorted(MARKERS), default='native_solver')
    ap.add_argument('--timeout', type=float, default=3600,
                    help='seconds allowed for each pytest case (default: 3600)')
    ap.add_argument('--heartbeat', type=float, default=60)
    ap.add_argument('--report-dir', type=Path, help='directory to preserve logs and summary')
    ap.add_argument('--nodeid', help='run only a single exact collected test nodeid')
    ap.add_argument('--list-only', action='store_true')
    args = ap.parse_args(argv)
    if args.timeout <= 0 or args.heartbeat <= 0:
        ap.error('--timeout and --heartbeat must be positive')

    try:
        nodeids = collect_nodeids(MARKERS[args.marker])
    except RuntimeError as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2
    if args.nodeid:
        nodeids = [n for n in nodeids if n == args.nodeid]
    if not nodeids:
        print('No native cases selected; check --marker or --nodeid', file=sys.stderr)
        return 2
    if args.list_only:
        print('\n'.join(nodeids))
        return 0
    stamp = dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    report = (args.report_dir or ROOT / 'native_test_results' / stamp).resolve()
    report.mkdir(parents=True, exist_ok=False)
    (report / 'scratch').mkdir()
    records = []
    try:
        for i, nodeid in enumerate(nodeids, start=1):
            log_path = report / f'{i:03d}.log'
            print(f'[START] {i}/{len(nodeids)} {nodeid}', flush=True)
            command = [sys.executable, '-m', 'pytest', '-o', 'addopts=',
                       '-q', '-rs', f'--basetemp={report / 'scratch' / f"{i:03d}"}', nodeid]
            try:
                status, elapsed, returncode = execute_command(
                    command, log_path, timeout=args.timeout, heartbeat=args.heartbeat,
                    env=os.environ.copy(),
                )
            except OSError as exc:
                status, elapsed, returncode = 'LAUNCH_FAILED', 0, None
                log_path.write_text(f'{type(exc).__name__}: {exc}\n', encoding='utf-8')
            records.append(dict(test=nodeid, status=status, elapsed_seconds=round(elapsed, 3),
                                returncode=returncode, log=str(log_path)))
            print(f'[END] {status} ({elapsed:.1f}s) {nodeid}  log={log_path}', flush=True)
    except KeyboardInterrupt:
        print('Interrupted; completed test logs retained', file=sys.stderr)
        return 130
    finally:
        (report / 'summary.json').write_text(json.dumps({
            'version': '2026.10.31.9', 'marker': args.marker,
            'python': sys.executable, 'platform': platform.platform(),
            'vspaero_binary': os.environ.get('VSPAERO_BINARY'),
            'timeout_seconds_per_case': args.timeout, 'results': records,
        }, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'[SUMMARY] {sum(r["status"] == "PASSED" for r in records)}/{len(records)} PASSED; '
          f'report={report}', flush=True)
    return 0 if all(r['status'] == 'PASSED' for r in records) else 1


if __name__ == '__main__':
    sys.exit(main())
