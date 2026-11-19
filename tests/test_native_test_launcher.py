"""Offline coverage for collection and subprocess timeout of the native test runner."""
import importlib.util
from pathlib import Path
import sys
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('native_pytest_runner', ROOT / 'tools/run_native_pytest.py')
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


def test_native_marker_collects_without_importing_openvsp():
    counts = {'native_vspaero': 7, 'native_openvsp': 4, 'native_solver': 1}
    for mark, expected in counts.items():
        run = subprocess.run(
            [sys.executable, str(ROOT / 'tools/run_native_pytest.py'),
             '--marker', mark, '--list-only'],
            cwd=ROOT, capture_output=True, text=True, timeout=30,
        )
        assert run.returncode == 0, run.stderr
        nodeids = [line for line in run.stdout.splitlines() if line.strip()]
        assert len(nodeids) == expected
        assert all('::test_' in name for name in nodeids)
        if mark == 'native_vspaero':
            assert any('adjoint' in name for name in nodeids)


def test_subprocess_exit_status_and_timeout(tmp_path):
    status, elapsed, code = runner.execute_command(
        [sys.executable, '-c', 'print("ok")'], tmp_path / 'pass.log', timeout=5,
    )
    assert status == 'PASSED' and code == 0
    assert 'ok' in (tmp_path / 'pass.log').read_text()
    status, elapsed, code = runner.execute_command(
        [sys.executable, '-c', 'import time; time.sleep(60)'],
        tmp_path / 'timeout.log', timeout=0.25, heartbeat=0.1,
    )
    assert status == 'TIMEOUT'
    assert elapsed < 5
