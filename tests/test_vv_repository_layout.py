"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys

from openvsppy.vv.VVDatabase import REPO_ROOT
from openvsppy.vv.VVDatabase import digest
from openvsppy.vv.VVDatabase import execution_fingerprint
from openvsppy.vv.VVDatabase import load_study

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs/vv'


def test_vv_document_set_is_canonical():
    expected = {
        "VSPAERO_3.50.5_VV_Plan_v1.1.md",
        "VSPAERO_VV_DB_Program_Spec_v1.0.1.md",
        "VSPAERO_3.50.5_VV_Evaluation_Report_Template_v0.1.md",
        "Test_Record_Template.md",
        "DB_Report_Field_Mapping.md",
        "Methods_and_References.md",
        "VV_CLI_Guide.md",
        "VV_Implementation_Status.md",
        "CHANGELOG.md",
    }
    present = {p.name for p in DOCS.rglob("*.md")}
    assert expected <= present
    assert all(p.parent == DOCS for p in DOCS.rglob("*.md"))
    assert not (DOCS / "history").exists()


def test_plan_still_matches_reviewed_study_sha256():
    study = json.loads((ROOT / "vv_study" / "study.json").read_text(encoding="utf-8"))
    path = (ROOT / "vv_study" / study["plan"]["file"]).resolve()
    assert path == DOCS / "VSPAERO_3.50.5_VV_Plan_v1.1.md"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == study["plan"]["sha256"]


def test_report_template_keeps_thirteen_chapter_headers():
    template = (DOCS / "VSPAERO_3.50.5_VV_Evaluation_Report_Template_v0.1.md").read_text(encoding="utf-8")
    assert [int(n) for n in re.findall(r"(?m)^## (\d+)\.", template)] == list(range(1, 14))


def test_readme_vv_links_resolve():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    links = re.findall(r"\[[^\]]+\]\((docs/vv/[^)#]+)(?:#[^)]*)?\)", text)
    assert links
    for link in links:
        assert (ROOT / link).is_file(), link


def test_study_docs_and_package_layout():
    study_file, study, plan, catalog = load_study(REPO_ROOT/'vv_study/study.json')
    assert study['runtime']['openvsp_version_required'] == '3.50.5'
    assert plan.parent == REPO_ROOT/'docs/vv'
    assert plan.is_file() and digest(plan) == study['plan']['sha256']
    assert catalog == REPO_ROOT/'vv_study/01_database/test_catalog.csv'
    assert (REPO_ROOT/'src/openvsppy/vv/worker.py').is_file()
    assert (REPO_ROOT/'src/openvsppy/vspaero/AnalysisVSPAERO.py').is_file()
    assert not (REPO_ROOT/'src/VVRunner.py').exists()
    assert not (REPO_ROOT/'vv/__main__.py').exists()
    assert len(execution_fingerprint(REPO_ROOT)) == 64


def test_checkout_cli_and_worker_from_independent_working_directory(tmp_path):
    env = os.environ.copy()
    env['PYTHONPATH'] = str(REPO_ROOT/'src') + os.pathsep + env.get('PYTHONPATH', '')
    for argv in (('openvsppy.vv', '--help'), ('openvsppy.vv.worker', '--help')):
        proc = subprocess.run([sys.executable, '-m', *argv], cwd=tmp_path,
                              env=env, text=True, capture_output=True, timeout=30)
        assert proc.returncode == 0, proc.stderr
        assert 'usage:' in proc.stdout


def test_notebooks_bootstrap_src_without_editable_install():
    notebooks = sorted(REPO_ROOT.rglob('*.ipynb'))
    assert notebooks, 'Expected at least one example Notebook'
    for p in notebooks:
        nb = json.loads(p.read_text(encoding='utf-8'))
        code = nb['cells'][0]
        assert code['cell_type'] == 'code'
        assert 'repository-bootstrap' in code['metadata']['tags']
        source = ''.join(code['source'])
        assert 'pip install -e' in source  # The comment says installation is unnecessary.
        assert "sys.path.insert(0, source_dir)" in source
        assert all(c.get('outputs', []) == [] for c in nb['cells'] if c['cell_type']=='code')
        assert all(c.get('execution_count') is None for c in nb['cells'] if c['cell_type']=='code')
        assert not any('from src.' in ''.join(c['source']) for c in nb['cells'] if c['cell_type']=='code')
