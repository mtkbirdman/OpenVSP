"""Feature-focused V&V regression tests; offline checks do not establish native solver success."""
from __future__ import annotations

from pathlib import Path
import json
import os
import subprocess
import sys

import pytest

from openvsppy.vv.VVRunner import wait_for_worker

ROOT = Path(__file__).resolve().parents[1]


def test_cli_prepare_register_and_run_test_guard(study,tmp_path,capsys, configured_case):
    from openvsppy.vv.__main__ import main
    root,db,cond=configured_case
    template=json.loads((ROOT/'docs/vv/templates/case.json').read_text())
    template['geometry']['vsp3_source']=str(ROOT/'examples/models/CRM-HS/CRM-HS.vsp3')
    template['flight_condition']=cond
    template['numerical_setting']['requested']={
        'thin_set_index':3,'thick_set_index':4,
        'vspaero':{'WakeNumIter':5,'RootWakeNodes':8,'FixedWakeFlag':False}}
    baseline=tmp_path/'baseline.json';baseline.write_text(json.dumps(template))
    assert main(['--study',str(root/'study.json'),'--db',str(db.path),
                 'prepare','--test','WAKE-01','--base',str(baseline),'--dry-run'])==0
    assert not (root/'generated/WAKE-01').exists()
    capsys.readouterr()
    assert main(['--study',str(root/'study.json'),'--db',str(db.path),
                 'prepare','--test','WAKE-01','--base',str(baseline),'--register'])==0
    content=json.loads(capsys.readouterr().out)
    assert len(content['cases'])==3
    assert {c['registration'] for c in content['cases']}=={'REGISTERED'}
    with pytest.raises(ValueError,match='Only --jobs 1'):
        main(['--study',str(root/'study.json'),'--db',str(db.path),'run-test',
              '--test','WAKE-01','--binary',str(ROOT/'README.md'),'--jobs','4'])


def test_monitor_reads_wake_and_heartbeat_during_worker_execution(tmp_path):
    native = tmp_path / 'vspaero_api_stdout.log'
    script = (
        "import pathlib,sys,time\n"
        "p=pathlib.Path(sys.argv[1])\n"
        "for n in (1,5,10,20):\n"
        "    with p.open('ab') as log:\n"
        "        log.write(f'Wake Iter: {n} / 20 ... GMRES Iter: 0\\n'.encode())\n"
        "    time.sleep(0.08)\n"
    )
    cmd = [sys.executable, '-c', script, str(native)]
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    events = []
    try:
        assert wait_for_worker(proc, cmd, 5, native, progress_callback=events.append,
                               verbose=True, heartbeat_seconds=0.06, poll_seconds=0.02) == 0
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
    assert any(x['event'] == 'WAKE' and x['iteration'] == 20 for x in events)
    assert any(x['event'] == 'HEARTBEAT' for x in events)
    assert [x['iteration'] for x in events if x['event'] == 'WAKE'] == sorted(
        {x['iteration'] for x in events if x['event'] == 'WAKE'})


def test_monitor_missing_log_and_broken_display_do_not_fail(tmp_path):
    cmd = [sys.executable, '-c', 'import time;time.sleep(0.15)']
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    def bad_display(event):
        raise BrokenPipeError('simulated terminal disconnect')
    try:
        assert wait_for_worker(proc, cmd, 3, tmp_path / 'missing.log',
            progress_callback=bad_display, verbose=True,
            heartbeat_seconds=0.04, poll_seconds=0.02) == 0
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()


def test_monitor_timeout_is_not_extended_by_heartbeat(tmp_path):
    cmd = [sys.executable, '-c', 'import time;time.sleep(5)']
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        with pytest.raises(subprocess.TimeoutExpired):
            wait_for_worker(proc, cmd, 0.12, tmp_path / 'missing.log',
                progress_callback=lambda _: None, verbose=True,
                heartbeat_seconds=0.03, poll_seconds=0.02)
    finally:
        proc.kill()
        proc.wait(timeout=3)


def test_execute_all_cli_stdout_json_stderr_progress_and_jsonl(tmp_path, study_factory):
    study, db = study_factory(tmp_path / 'vv_study')
    db_path = db.path
    manifest = study / 'manifest.json'

    env = os.environ.copy()
    env['PYTHONPATH'] = str(ROOT / 'src')
    cmd = [sys.executable, '-m', 'openvsppy.vv', '--study', str(study/'study.json'),
           '--db', str(db_path), 'execute-all', '--config', str(manifest),
           '--no-run', '--verbose']
    result = subprocess.run(cmd, cwd=tmp_path, env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload['status'] == 'REVIEW_REPORT_GENERATED'
    assert payload['case_count'] == 68
    assert '[PREPARE]' in result.stderr
    assert '[REPORT]' in result.stderr
    assert '[SUMMARY]' in result.stderr
    assert '[VSPAERO]' not in result.stderr  # --no-run is not native success
    events = [json.loads(x) for x in (Path(payload['output_dir'])/'progress.jsonl')
              .read_text(encoding='utf-8').splitlines()]
    assert any(x.get('phase') == 'SUMMARY' for x in events)
    assert any(x.get('phase') == 'REPORT' for x in events)


def test_execute_all_no_binary_records_each_blocked_case_and_prints_summary(tmp_path, study_factory):
    study, db = study_factory(tmp_path / 'study')

    env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'src')
    cmd=[sys.executable,'-m','openvsppy.vv','--study',str(study/'study.json'),
        '--db',str(study/'01_database/vv.sqlite'),'execute-all',
        '--config',str(study/'manifest.json')]
    done=subprocess.run(cmd,cwd=tmp_path,env=env,capture_output=True,text=True,timeout=60)
    assert done.returncode==0,done.stderr
    report=json.loads(done.stdout)
    assert report['blocked_or_failed']==68
    assert '[CASE] 1/68' in done.stderr
    assert 'BLOCKED' in done.stderr
    assert '[SUMMARY] 68 Cases; 0 succeeded/reused; 68 blocked/failed' in done.stderr
    events=[json.loads(line) for line in (Path(report['output_dir'])/'progress.jsonl').read_text().splitlines()]
    finished=[e for e in events if e.get('phase')=='RUN' and e.get('event')=='FINISHED' and e.get('case_id')]
    assert len(finished)==68
    assert all(e['status']=='BLOCKED' for e in finished)
    assert report['status']=='REVIEW_REPORT_GENERATED'  # Not a claim that Cases passed
