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
        "plan.md",
        "database_spec.md",
        "evaluation_report_template.md",
        "test_record_template.md",
        "db_report_field_mapping.md",
        "methods_and_references.md",
        "cli_guide.md",
        "implementation_status.md",
        "CHANGELOG.md",
    }
    present = {p.name for p in DOCS.rglob("*.md")}
    assert expected == present
    assert all(p.parent == DOCS for p in DOCS.rglob("*.md"))
    assert not (DOCS / "history").exists()
    assert not (ROOT / "tools/compare_g103a_direct.py").exists()
    assert (ROOT / "tools/run_native_pytest.py").is_file()


def test_plan_still_matches_reviewed_study_sha256():
    study = json.loads((ROOT / "vv_study" / "study.json").read_text(encoding="utf-8"))
    path = (ROOT / "vv_study" / study["plan"]["file"]).resolve()
    assert path == DOCS / "plan.md"
    assert study["plan"]["version"] == "1.1"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == study["plan"]["sha256"]


def test_report_template_keeps_thirteen_chapter_headers():
    template = (DOCS / "evaluation_report_template.md").read_text(encoding="utf-8")
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


def test_reference_evidence_matches_manifest_bytes():
    """Canonical native evidence must survive ZIP/Git round-trips byte-for-byte."""
    evidence = ROOT / 'vv_study/02_reference/20261031_api_vs_direct'
    manifest = json.loads((evidence / 'EVIDENCE_MANIFEST.json').read_text(encoding='utf-8'))
    for relative, expected in manifest['files'].items():
        path = evidence / relative
        assert path.is_file(), relative
        assert path.stat().st_size == expected['size_bytes'], relative
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected['sha256'], relative


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
