"""Offline tests for report factual autofill; no native OpenVSP dependencies."""
from __future__ import annotations
import copy
import hashlib
import json
import re
from pathlib import Path

import pytest
from openvsppy.vv.VVPostprocess import collect_study_results
from openvsppy.vv.VVDatabase import digest
from openvsppy.vv.VVReport import write_evaluation_report
from openvsppy.vv.VVReportFacts import (
    _equivalent, _get, _json, _fmt, _record_settings, _selected,
    collect_report_facts, write_report_facts,
)


def _snapshot(root, db):
    s=collect_study_results(db,root)
    s['charts']=[]
    return s


def _sample_run():
    return {
        'run_id':'RUN-TEST','requested_snapshot_json':json.dumps({
            'setting':{'requested_json':json.dumps({'vspaero':{
                'WakeNumIter':20,'RootWakeNodes':8,'FixedWakeFlag':False,
                'ImplicitWake':False,'Symmetry':False},'thin_set_index':4,'thick_set_index':3})},
            'condition':{'sref':17.,'bref':12.,'cref':2.,'moment_x':0.,'moment_y':0.,'moment_z':0.}}),
        'effective_snapshot_json':json.dumps({'Symmetry':0,'ImplicitWake':0,
               'ThinGeomSet':4,'GeomSet':3,'Sref':17.,'bref':12.,'cref':2.}),
        'generated_snapshot_json':json.dumps({'solver_settings':{
              'WakeIters':20,'NumWakeNodes':8,'Symmetry':0,'ImplicitWake':0,
              'Sref':17.,'Bref':12.,'Cref':2.}}),
    }


def test_report_sections_1_to_4_are_filled_and_no_markers(study,tmp_path):
    root,db=study
    snapshot=_snapshot(root,db)
    before=digest(db.path)
    dest=tmp_path/'report.md'
    write_evaluation_report(db,dest,results=snapshot)
    text=dest.read_text(encoding='utf-8')
    assert len(re.findall(r'^## \d+\. ',text,re.M))==13
    for section in ('1.1','1.2','1.3','2.1','2.2','2.3','3.1','3.2','3.3',
                    '3.4','3.5','4.1','4.2','4.3','4.4'):
        assert re.search(r'^### '+re.escape(section)+r' ',text,re.M), section
    assert '<!-- VV:' not in text
    assert (tmp_path/'report_field_provenance.json').is_file()
    assert len(list((tmp_path/'test_reports').glob('*.md')))==25
    assert digest(db.path)==before
    assert 'NOT_EVALUATED' in text and 'INSUFFICIENT_EVIDENCE' in text


def test_audit_json_and_csv_have_traceable_field_ids(study,tmp_path):
    root,db=study
    f=collect_report_facts(db,_snapshot(root,db))
    write_report_facts(f,tmp_path)
    provenance=json.loads((tmp_path/'report_field_provenance.json').read_text(encoding='utf-8'))
    ids={x['field_id'] for x in provenance['fields']}
    assert '1.case_count' in ids
    assert any(x.startswith('section.4.') for x in ids)
    for name in ('report_field_audit.csv','report_settings_audit.csv',
                 'report_conditions_audit.csv','report_environment_audit.csv'):
        assert (tmp_path/'tables'/name).is_file()


def test_selected_snapshot_validates_case_and_denies_fabricated_run(study):
    root,db=study
    data=_snapshot(root,db)
    assert len(_selected(db,data))==len(data['cases'])
    wrong=copy.deepcopy(data)
    wrong['cases'][0]['run_id']='RUN-FAKE'
    wrong['cases'][0]['status']='SUCCEEDED'
    with pytest.raises(ValueError,match='Unmatched selected Run'):
        _selected(db,wrong)


def test_unspecified_profile_keeps_manual_baseline(study):
    root,db=study
    facts=collect_report_facts(db,_snapshot(root,db))
    assert '基準Case未確定' in facts['sections']['4.3']
    assert '未認定' in facts['sections']['2.3']
    assert 'PASS' not in facts['sections']['3.4'].replace('SUCCEEDEDをEngineering PASSへ変換しません','')


def test_mismatched_profile_rejected(study):
    root,db=study
    with pytest.raises(ValueError,match='Study ID mismatch'):
        collect_report_facts(db,_snapshot(root,db),{
            'schema_version':'REPORT_PROFILE_V1','study_id':'NOT-THIS-STUDY'})


def test_unknown_explicit_baseline_rejected(study):
    root,db=study
    with pytest.raises(ValueError,match='selected successful Case'):
        collect_report_facts(db,_snapshot(root,db),{
            'schema_version':'REPORT_PROFILE_V1','study_id':'VSP3505-VV-001',
            'baseline_case_id':'UNEXISTING-CASE'})


def test_requested_effective_generated_agree_and_fixed_wake_is_request_only():
    fields=[]
    entry={'case':{'case_id':'GEO-03-G3'},'run':_sample_run()}
    rows=_record_settings(entry,fields)
    values={r[0]:r for r in rows}
    assert values['Wake反復数'][1]=='20' and values['Wake反復数'][3]=='20'
    assert values['Wake反復数'][4]=='CROSSCHECKED'
    assert values['Wake節点数'][4]=='CROSSCHECKED'
    assert values['Fixed Wake'][1]=='false' and values['Fixed Wake'][3]=='—'
    assert values['Fixed Wake'][4]=='RECORDED'
    assert values['Implicit Wake'][4]=='CROSSCHECKED'
    assert values['Symmetry'][4]=='CROSSCHECKED'
    assert values['Thick Set'][4]=='CROSSCHECKED'
    assert values['Thin Set'][4]=='CROSSCHECKED'
    assert values['Sref [m²]'][4]=='CROSSCHECKED'
    assert any(x['source_values'] for x in fields)


def test_conflicting_settings_preserve_both_values():
    entry={'case':{'case_id':'GEO-03-G3'},'run':_sample_run()}
    generated=json.loads(entry['run']['generated_snapshot_json'])
    generated['solver_settings']['WakeIters']=10
    entry['run']['generated_snapshot_json']=json.dumps(generated)
    fields=[]; rows=_record_settings(entry,fields)
    result=next(x for x in rows if x[0]=='Wake反復数')
    assert result[1]=='20' and result[3]=='10' and result[4]=='CONFLICT'
    assert any(x['state']=='CONFLICT' and len(x['source_values'])==2 for x in fields)


@pytest.mark.parametrize('a,b,expected',[
    (False,0,True),(True,1,True),(False,1,False),
    (0.0,0,True),(0,1,False), (1.0,1.0+1e-10,True),
    (1.0,1.1,False),('A','A',True),('A','B',False),
])
def test_record_value_comparison(a,b,expected):
    assert _equivalent(a,b,boolean=isinstance(a,bool))==expected


def test_null_zero_false_not_conflated():
    assert _fmt(0)=='0'
    assert _fmt(False)=='false'
    assert _fmt(None)=='DB未記録'
    assert _get({'x':0},('x',))==0
    assert _get({'x':False},('x',)) is False


def test_legacy_missing_generated_tess_is_not_fabricated():
    entry={'case':{'case_id':'GRID-01-MEDIUM'},'run':_sample_run()}
    rows=_record_settings(entry,[])
    assert next(x for x in rows if x[0]=='Wing Tess_W')[4]=='NOT_RECORDED'
    assert next(x for x in rows if x[0]=='Wing SectTess_U')[4]=='NOT_RECORDED'


def test_report_profile_is_not_engineering_approval(study,tmp_path):
    root,db=study
    profile={'schema_version':'REPORT_PROFILE_V1','study_id':'VSP3505-VV-001',
        'document_control':{'classification':'internal','approval_status':'DRAFT',
                            'author':'Reviewer','reviewer':None,'approver':None}}
    write_evaluation_report(db,tmp_path/'draft.md',results=_snapshot(root,db),report_profile=profile)
    text=(tmp_path/'draft.md').read_text(encoding='utf-8')
    assert '| 情報区分 | internal |' in text
    assert 'Reviewer / 未設定 / 未設定' in text
    assert 'INSUFFICIENT_EVIDENCE' in text


def test_no_native_dependency_required():
    import sys
    from openvsppy.vv import VVReportFacts
    assert hasattr(VVReportFacts,'collect_report_facts')
    assert 'openvsp' not in VVReportFacts.__dict__


def test_all_provenance_field_states_valid(study):
    root,db=study
    facts=collect_report_facts(db,_snapshot(root,db))
    allowed={'RECORDED','CROSSCHECKED','CONFLICT','NOT_RECORDED','NOT_EXECUTED',
             'REVIEW_REQUIRED','MANUAL_DECISION','NOT_APPLICABLE'}
    assert {x['state'] for x in facts['fields']} <= allowed
    assert all(x['field_id'] and 'source_values' in x for x in facts['fields'])


def test_report_all_sections_gfm_width(study,tmp_path):
    from markdown_it import MarkdownIt
    root,db=study
    write_evaluation_report(db,tmp_path/'report.md',results=_snapshot(root,db))
    tokens=MarkdownIt('default').enable('table').parse((tmp_path/'report.md').read_text(encoding='utf-8'))
    assert sum(t.type=='table_open' for t in tokens)>=10
    assert sum(t.type=='td_open' for t in tokens)>25


def test_approved_and_provisional_rules_remain_distinct(study):
    root,db=study
    with db.transaction() as c:
        for rule_id,approval in (('RULE-TEMP','PROVISIONAL'),('RULE-APPROVED','APPROVED')):
            c.execute('''INSERT INTO acceptance_rule(rule_id,study_id,qoi,level,policy_json,
                approval,rationale,reviewer,approved_utc) VALUES (?,?,?,?,?,?,?,?,?)''',
                (rule_id,'VSP3505-VV-001','CLi','ENGINEERING','{}',approval,
                 'Synthetic test only','TESTER' if approval=='APPROVED' else None,
                 '2026-11-04T00:00:00Z' if approval=='APPROVED' else None))
    facts=collect_report_facts(db,_snapshot(root,db))
    section=facts['sections']['3.4']
    assert 'RULE-TEMP' in section and 'PROVISIONAL' in section
    assert 'RULE-APPROVED' in section and 'APPROVED' in section
    assert 'SUCCEEDEDをEngineering PASSへ変換しません' in section


def test_full_report_provenance_only_generated_from_selected_data(study,tmp_path):
    root,db=study
    snapshot=_snapshot(root,db)
    db_hash=digest(db.path)
    output=tmp_path/'evaluation_report.md'
    write_evaluation_report(db,output,results=snapshot,database_snapshot=db.path)
    md=output.read_text(encoding='utf-8')
    assert db_hash in md
    assert '未記録' in md
    assert 'GEO-01-G1-G' not in md or 'CONFLICT' in md or '予定構成のみ' in md
    assert digest(db.path)==db_hash
