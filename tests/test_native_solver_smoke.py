"""Standalone Linux/Windows VSPAERO diagnostic; not the production API workflow.

The input is an isolated, single-Wake, zero-control-group G103A diagnostic
prepared from a prior successful native Linux run. It is NOT an accepted VV
case or a substitute for the reviewed 20-Wake GEO-03-G3 run.
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'tests/fixtures/native_solver_smoke'


@pytest.mark.native_solver
@pytest.mark.slow
def test_native_solver_smoke(vspaero_binary, native_workspace):
    for source in FIXTURE.iterdir():
        shutil.copy2(source, native_workspace / source.name)
    log = native_workspace / 'solver_stdout.log'
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '1'
    with log.open('wb') as output:
        # Same syntax as the validated native standalone smoke. This does not
        # alter the Python API execution path of any OpenVSP/VV production code.
        try:
            run = subprocess.run(
                [str(vspaero_binary), '-omp', '1', 'smoke'],
                cwd=native_workspace, env=env, stdin=subprocess.DEVNULL,
                stdout=output, stderr=subprocess.STDOUT,
                timeout=240, check=False,
            )
        except subprocess.TimeoutExpired:
            pytest.fail(f'Native solver exceeded 240s; log: {log}', pytrace=False)
    assert run.returncode == 0, f'VSPAERO returncode={run.returncode}; log: {log}'
    for ext in ('.polar', '.history', '.lod', '.adb'):
        result = native_workspace / f'smoke{ext}'
        assert result.is_file() and result.stat().st_size > 0, f'Missing native output: {result}'
    from openvsppy.vv.VVPostprocess import read_history
    # Keep the detailed .history parser exercised on actual native output.
    parsed = read_history(native_workspace / 'smoke.history')
    assert parsed
