"""Shared fixtures for OpenVSP/VSPAERO tests.

The checkout's src/ directory is imported directly; editable installation is optional.
These fixtures do not automatically execute native solvers.
"""
import importlib
import os
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src'
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@pytest.fixture
def native_openvsp():
    """Fail explicitly if a selected native test cannot load the OpenVSP API."""
    try:
        vsp = importlib.import_module('openvsp')
    except (ImportError, OSError) as exc:
        pytest.fail(
            'OpenVSP Python API unavailable (native test was explicitly selected): '
            f'{exc}. Install the matching OpenVSP Python extension and its native libraries.',
            pytrace=False,
        )
    if not hasattr(vsp, 'GetVSPVersion'):
        pytest.fail('Imported openvsp module is not the OpenVSP native API', pytrace=False)
    version = str(vsp.GetVSPVersion())
    if '3.50.5' not in version:
        pytest.fail(f'OpenVSP 3.50.5 required, found {version!r}', pytrace=False)
    return vsp


@pytest.fixture
def vspaero_binary():
    """Locate a runnable solver for deliberately selected native tests."""
    binary = os.environ.get('VSPAERO_BINARY')
    if not binary:
        binary = shutil.which('vspaero') or shutil.which('vspaero.exe')
    if not binary:
        pytest.fail('VSPAERO_BINARY must point to the vspaero(.exe) executable '
                    '(or the executable must be in PATH)', pytrace=False)
    path = Path(binary).expanduser().resolve()
    if not path.is_file():
        pytest.fail(f'VSPAERO executable not found: {path}', pytrace=False)
    return path


@pytest.fixture
def native_vspaero_api(native_openvsp, vspaero_binary):
    """Use the official OpenVSP API's solver location; never run VSPAERO directly."""
    vsp = native_openvsp
    if not vsp.SetVSPAEROPath(str(vspaero_binary.parent)):
        pytest.fail(f'OpenVSP SetVSPAEROPath failed for {vspaero_binary.parent}', pytrace=False)
    if Path(vsp.GetVSPAEROPath()).resolve() != vspaero_binary.parent:
        pytest.fail('OpenVSP solver directory does not match requested VSPAERO_BINARY', pytrace=False)
    return vsp


@pytest.fixture
def native_workspace(tmp_path, monkeypatch):
    """Ensure native output defaults cannot write into the checkout."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture
def isolated_model(native_workspace):
    """Copy model and its sibling inputs to per-test scratch space."""
    def copy(source):
        source = Path(source).resolve()
        target = native_workspace / source.parent.name
        if target.exists():
            raise ValueError(f'Model directory already prepared: {target}')
        shutil.copytree(source.parent, target)
        return target / source.name
    return copy
