"""Shared fixtures for isolated, offline V&V tests.

The checkout's src/ directory is imported directly; Editable Install is optional.
Never call the native VSPAERO executable from these fixtures.
"""
import json
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

from openvsppy.vv.VVDatabase import VVDatabase


@pytest.fixture
def study_factory(request):
    """Return a factory for independent Study folders with a fresh SQLite DB."""
    def create(folder):
        folder = Path(folder)
        (folder / '00_plan').mkdir(parents=True)
        (folder / '01_database').mkdir(exist_ok=True)
        (folder / '03_run_raw').mkdir(exist_ok=True)
        shutil.copy2(ROOT / 'docs/vv/VSPAERO_3.50.5_VV_Plan_v1.1.md', folder / '00_plan')
        for name in ('schema.sql', 'test_catalog.csv'):
            shutil.copy2(ROOT / 'vv_study/01_database' / name, folder / '01_database')
        settings = json.loads((ROOT / 'vv_study/study.json').read_text(encoding='utf-8'))
        settings['plan']['file'] = '00_plan/VSPAERO_3.50.5_VV_Plan_v1.1.md'
        (folder / 'study.json').write_text(json.dumps(settings, ensure_ascii=False), encoding='utf-8')
        (folder / 'manifest.json').write_text('{}', encoding='utf-8')
        db = VVDatabase(folder / '01_database/vv.sqlite')
        db.initialize(folder / 'study.json', folder / '01_database/schema.sql')
        request.addfinalizer(db.close)
        return folder, db
    return create


@pytest.fixture
def study(study_factory, tmp_path):
    return study_factory(tmp_path / 'vv_study')


@pytest.fixture
def configured_case(study):
    root,db=study
    src=ROOT/'examples/models/CRM-HS/CRM-HS.vsp3'
    db.add_geometry('VSP3505-VV-001',{'geometry_id':'G3','stage':'G3','vsp3_source':str(src),
           'physical_model':{'wing':'Thin','fuselage':'Thick'},'representation':{'set_indices':[3,4]}})
    cond={'condition_id':'FC-001','mach':0.6,'re_cref':1e7,'alpha_deg':3.,'beta_deg':0.,
          'controls_deg':{'elevator':0.,'aileron':0.,'rudder':0.},'trim_mode':'FIXED_CONTROL',
          'sref':10.,'bref':12.,'cref':2.,'moment_x':0.,'moment_y':0.,'moment_z':0.,
          'reference_length_unit':'m','axes_convention_id':'VSPAERO_NATIVE'}
    db.add_condition('VSP3505-VV-001',cond)
    db.add_setting('VSP3505-VV-001',{'setting_id':'SET1','analysis_mode':'STEADY_SWEEP','mesh_level':'MEDIUM',
                                   'requested':{'thin_set_index':3,'thick_set_index':4,
                                                'vspaero':{'WakeNumIter':5,'RootWakeNodes':8}}})
    db.resolve_case('GRID-01-MEDIUM','G3','FC-001','SET1')
    return root,db,cond


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
