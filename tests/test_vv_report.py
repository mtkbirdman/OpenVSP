"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path
import json
import re

from openvsppy.vv.VVReport import write_evaluation_report

ROOT = Path(__file__).resolve().parents[1]


def test_formal_13_chapter_report_and_evidence(study,tmp_path):
    _,db=study
    out=tmp_path/'report.md'
    write_evaluation_report(db,out)
    content=out.read_text(encoding='utf-8')
    import re
    assert len(re.findall(r'^## \d+\. ',content,re.M))==13
    assert 'ACC-01' in content and 'INSUFFICIENT_EVIDENCE' in content
    assert content.count('《未設定')>0


def test_new_cli_assess_report_export(study,tmp_path,capsys):
    from openvsppy.vv.__main__ import main
    root,db=study
    assert main(['--study',str(root/'study.json'),'--db',str(db.path),'assess-all',
                 '--output',str(tmp_path/'inventory.json')])==0
    result=json.loads(capsys.readouterr().out)
    assert result['test_ids']==25 and len(result['tests'])==25
    assert {t['engineering_verdict'] for t in result['tests']}=={'INSUFFICIENT_EVIDENCE'}
    report=tmp_path/'report.md'
    assert main(['--study',str(root/'study.json'),'--db',str(db.path),'report','--output',str(report)])==0
    assert 'ACC-01' in report.read_text(encoding='utf-8')
    assert 'INSUFFICIENT_EVIDENCE' in report.read_text(encoding='utf-8')
    assert main(['--study',str(root/'study.json'),'--db',str(db.path),'export',
                 '--output',str(tmp_path/'independent')])==0
    assert (tmp_path/'independent/tables/planned_case.csv').exists()
