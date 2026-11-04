"""Read-only, provenance-aware Ch.1–4 report facts from a frozen V&V database.

No solver import, database updates, implicit approvals, or historical Git inference.
"""
from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

STATES = {'RECORDED', 'CROSSCHECKED', 'CONFLICT', 'NOT_RECORDED',
          'NOT_EXECUTED', 'REVIEW_REQUIRED', 'MANUAL_DECISION', 'NOT_APPLICABLE'}
EXTERNAL_KINDS = {'WIND_TUNNEL', 'PANEL', 'NS', 'THEORY'}
SETTINGS = (
    ('Wake反復数', ('setting','requested_json','vspaero','WakeNumIter'), None, ('solver_settings','WakeIters')),
    ('Wake節点数', ('setting','requested_json','vspaero','RootWakeNodes'), None, ('solver_settings','NumWakeNodes')),
    ('Implicit Wake', ('setting','requested_json','vspaero','ImplicitWake'), ('ImplicitWake',), ('solver_settings','ImplicitWake')),
    ('Fixed Wake', ('setting','requested_json','vspaero','FixedWakeFlag'), None, None),
    ('Symmetry', ('setting','requested_json','vspaero','Symmetry'), ('Symmetry',), ('solver_settings','Symmetry')),
    ('Reference Flag', None, ('RefFlag',), None),
    ('Thick Set', ('setting','requested_json','thick_set_index'), ('GeomSet',), None),
    ('Thin Set', ('setting','requested_json','thin_set_index'), ('ThinGeomSet',), None),
    ('Sref [m²]', ('condition','sref'), ('Sref',), ('solver_settings','Sref')),
    ('Bref [m]', ('condition','bref'), ('bref',), ('solver_settings','Bref')),
    ('Cref [m]', ('condition','cref'), ('cref',), ('solver_settings','Cref')),
    ('Xcg [m]', ('condition','moment_x'), ('Xcg',), ('solver_settings','X_cg')),
    ('Ycg [m]', ('condition','moment_y'), ('Ycg',), ('solver_settings','Y_cg')),
    ('Zcg [m]', ('condition','moment_z'), ('Zcg',), ('solver_settings','Z_cg')),
    ('Wing Tess_W', ('setting','requested_json','geom_parms','WingGeom','Tess_W'), ('geom_parms','WingGeom','Tess_W'), None),
    ('Wing SectTess_U', ('setting','requested_json','geom_parms','WingGeom','SectTess_U'), ('geom_parms','WingGeom','SectTess_U'), None),
)
BOOL_FIELDS = {'Implicit Wake', 'Fixed Wake', 'Symmetry'}


def _json(raw):
    if isinstance(raw, (dict, list)):
        return raw
    if not raw:
        return {}
    try:
        return json.loads(raw)
    except (TypeError, json.JSONDecodeError):
        return {}


def _get(obj, path):
    if not path:
        return None
    for part in path:
        if not isinstance(obj, dict):
            return None
        obj = obj.get(part)
    return obj


def _fmt(value):
    if value is None or value == '':
        return 'DB未記録'
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return format(value, '.7g') if math.isfinite(value) else 'DB未記録'
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _equivalent(a, b, boolean=False):
    if boolean:
        def as_bool(x):
            if type(x) is bool: return x
            if x in (0, 0.0, '0', 'false', 'False'): return False
            if x in (1, 1.0, '1', 'true', 'True'): return True
            return x
        return as_bool(a) == as_bool(b)
    if isinstance(a, (int,float)) and isinstance(b,(int,float)):
        return math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-9)
    return a == b


def _table(headers, rows):
    def safe(v):
        return _fmt(v).replace('\\', '\\\\').replace('|', '\\|').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(map(safe,headers))+' |',
                      '| '+' | '.join('---' for _ in headers)+' |']+
                     ['| '+' | '.join(map(safe,row))+' |' for row in rows])+'\n'


def _field(fields, identifier, value, state, source=None, **extra):
    if state not in STATES:
        raise ValueError(state)
    item={'field_id':identifier, 'display_value':value,'state':state,
          'source_values':source or []}
    item.update(extra)
    fields.append(item)
    return item


def _row(db, table, key, value):
    if value is None: return None
    return db.conn.execute(f'SELECT * FROM {table} WHERE {key}=?', (value,)).fetchone()


def _selected(db, snapshot):
    """Reject an unverified or fabricated selected Case/Run relationship."""
    if snapshot is None:
        return []
    selected=[]
    seen=set()
    for case in snapshot['cases']:
        cid=case['case_id']; rid=case.get('run_id')
        if cid in seen:
            raise ValueError('Duplicate snapshot Case: '+cid)
        seen.add(cid)
        planned=_row(db,'planned_case','case_id',cid)
        if planned is None or planned['test_id']!=case['test_id']:
            raise ValueError('Snapshot planned Case mismatch: '+cid)
        if bool(rid)!=(case.get('status')=='SUCCEEDED'):
            raise ValueError('Inconsistent selected status: '+cid)
        if not rid:
            selected.append({'case':case, 'planned':dict(planned), 'run':None})
            continue
        matching=db.conn.execute('''SELECT r.* FROM run r JOIN case_point cp ON cp.run_id=r.run_id
            JOIN case_run_link l ON l.point_id=cp.point_id
            WHERE l.case_id=? AND r.run_id=? AND cp.point_status='MATCHED' AND r.status='SUCCEEDED'
            LIMIT 1''',(cid,rid)).fetchone()
        if matching is None:
            raise ValueError('Unmatched selected Run: '+cid+' / '+rid)
        for q in case.get('quantities',[]):
            linked=db.conn.execute('''SELECT 1 FROM quantity q
                JOIN dataset d ON d.dataset_id=q.dataset_id
                JOIN case_point cp ON cp.point_id=d.point_id
                JOIN case_run_link l ON l.point_id=cp.point_id
                WHERE q.quantity_id=? AND d.dataset_id=? AND cp.run_id=?
                AND cp.point_status='MATCHED' AND l.case_id=? LIMIT 1''',
                (q['quantity_id'],q['dataset_id'],rid,cid)).fetchone()
            if linked is None:
                raise ValueError('Snapshot Quantity not linked to selected Case/Run/Point: '+cid+'/'+q['quantity_id'])
        selected.append({'case':case,'planned':dict(planned),'run':dict(matching)})
    n=db.conn.execute('SELECT COUNT(*) FROM planned_case').fetchone()[0]
    if n!=len(selected):
        raise ValueError('Snapshot Case population differs from DB')
    return selected


def _classify(entry):
    return 'BENCHMARK' if entry['case']['test_id'].startswith('VER-') else 'G103A'


def _record_settings(entry, fields):
    run=entry['run']
    if not run: return []
    rid=run['run_id']; cid=entry['case']['case_id']
    req=_json(run['requested_snapshot_json']);eff=_json(run['effective_snapshot_json']);gen=_json(run['generated_snapshot_json'])
    req['setting']=dict(req.get('setting') or {})
    req['setting']['requested_json']=_json(req['setting'].get('requested_json'))
    rows=[]
    for label,rp,ep,gp in SETTINGS:
        values=[]
        for source,path,root in (('Requested',rp,req),('Effective',ep,eff),('Generated',gp,gen)):
            found=_get(root,path)
            if found is not None and found != '':
                values.append((source,found, '.'.join(path)))
        if not values:
            state='NOT_RECORDED'
        elif len(values)==1:
            state='RECORDED'
        elif all(_equivalent(values[0][1],v[1],label in BOOL_FIELDS) for v in values[1:]):
            state='CROSSCHECKED'
        else:
            state='CONFLICT'
        out={name:next((_fmt(v) for s,v,_ in values if s==name),'—') for name in ('Requested','Effective','Generated')}
        item=_field(fields,'4.3.'+cid+'.'+label, '; '.join(f'{s}={_fmt(v)}' for s,v,_ in values)
                    if values else 'DB未記録', state,
                    [{'source':'run.'+{'Requested':'requested_snapshot_json','Effective':'effective_snapshot_json','Generated':'generated_snapshot_json'}[s]+'.'+path,
                      'value':v,'run_id':rid,'case_id':cid} for s,v,path in values],case_ids=[cid])
        rows.append((label,out['Requested'],out['Effective'],out['Generated'],state,rid))
    return rows


def collect_report_facts(db, snapshot, report_profile=None):
    """Collect data and provenance; no INSERT/UPDATE, no file or solver access."""
    study=db.conn.execute('SELECT * FROM study LIMIT 1').fetchone()
    if study is None: raise ValueError('No Study')
    profile=report_profile or {}
    if profile:
        if profile.get('schema_version') != 'REPORT_PROFILE_V1' or profile.get('study_id')!=study['study_id']:
            raise ValueError('Report profile schema or Study ID mismatch')
        if profile.get('document_control',{}).get('approval_status') not in (None,'DRAFT','PROPOSED','APPROVED'):
            raise ValueError('Unknown document status')
    fields=[]; warnings=[]; chapters={}; selected=_selected(db,snapshot)
    successes=[s for s in selected if s['run']]
    run_ids=sorted({s['run']['run_id'] for s in successes})
    db_runs=db.conn.execute('SELECT status,COUNT(*) n FROM run GROUP BY status').fetchall()
    run_stats={r['status']:r['n'] for r in db_runs}
    qty_selected=sum(len(s['case'].get('quantities',[])) for s in successes)
    qty_all=db.conn.execute('SELECT COUNT(*) FROM quantity').fetchone()[0]
    counts={table:db.conn.execute('SELECT COUNT(*) FROM '+table).fetchone()[0]
            for table in ('acceptance_rule','comparison','assessment')}
    ext=db.conn.execute('SELECT kind,usage_class,evidence_class,COUNT(*) n FROM dataset GROUP BY kind,usage_class,evidence_class ORDER BY kind').fetchall()
    ext=[dict(x) for x in ext]
    _field(fields,'1.case_count',len(selected),'RECORDED',[{'source':'planned_case.case_id'}])
    _field(fields,'1.selected_success',len(successes),'RECORDED',[{'source':'results_snapshot.cases[].run_id'}])
    _field(fields,'1.quantities',qty_selected,'RECORDED',[{'source':'results_snapshot.cases[].quantities'}])
    geometries=[dict(row) for row in db.conn.execute('SELECT * FROM geometry ORDER BY geometry_id')]
    geometry_rows=[]
    geom_map={g['geometry_id']:g for g in geometries}
    for g in geometries:
        model=_json(g['physical_model_json']).get('model',{})
        repr_data=_json(g['representation_json']);role=repr_data.get('case_role')
        if not g['geometry_id'].startswith('GEO-'):
            continue
        mismatches=[]
        if role in ('G1','G2','G3') and g['stage']!=role:
            mismatches.append(f'stage={g["stage"]} / case_role={role}')
        if g['geometry_id']=='GEO-04-THICK-G':
            if model.get('wing')=='Thin':
                mismatches.append('physical_model.wing=Thin / 実行RequestedではThick側の可能性：要確認')
        state='CONFLICT' if mismatches else 'RECORDED'
        if mismatches:
            warnings.append({'id':'OPEN-GEOM-'+g['geometry_id'],'state':'CONFLICT',
                             'message':'; '.join(mismatches),'source':g['geometry_id']})
        _field(fields,'2.2.geometry.'+g['geometry_id'],f'{g["stage"]} / {role}',state,
               [{'source':'geometry.stage','value':g['stage'],'geometry_id':g['geometry_id']},
                {'source':'geometry.representation_json.case_role','value':role,'geometry_id':g['geometry_id']}])
        if g['geometry_id'] in {s['run']['geometry_id'] for s in successes}:
            execution='採用成功Runあり'
        else:
            execution='予定構成のみ（採用成功Runなし）'
        geometry_rows.append((g['geometry_id'], g['stage'],role or 'DB未記録',
                 ', '.join(f'{k}={v}' for k,v in model.items() if k in ('wing','fuselage','horizontal_tail','vertical_tail')) or 'DB未記録',
                 g['geometry_equivalence_status'],execution,state))
    for s in successes:
        g=geom_map.get(s['run']['geometry_id'])
        if not g: continue
        if g['geometry_id']=='GEO-04-THICK-G':
            setting=_json(s['run']['requested_snapshot_json']).get('setting',{})
            member=_json(setting.get('requested_json',{})).get('geometry_membership',{})
            if 'WingGeom' in member.get('thick',[]):
                _field(fields,'2.2.GEO-04-THICK.membership','WingGeom=Thick (Requested)', 'CONFLICT',
                      [{'source':'run.requested_snapshot_json.setting.requested_json.geometry_membership.thick',
                        'value':member['thick'],'run_id':s['run']['run_id']}])
    # Separate planned (no success guarantee) from selected successful conditions.
    cond_ids=sorted({s['run']['condition_id'] for s in successes})
    condition_rows=[]; condition_audit=[]
    for condition_id in cond_ids:
        cond=_row(db,'flight_condition','condition_id',condition_id)
        if cond is None: raise ValueError('Selected Run condition missing: '+condition_id)
        matching=[s for s in successes if s['run']['condition_id']==condition_id]
        controls=_json(cond['controls_deg_json'])
        group='BENCHMARK' if all(_classify(x)=='BENCHMARK' for x in matching) else (
            'G103A' if all(_classify(x)=='G103A' for x in matching) else 'MIXED')
        if group=='MIXED':
            warnings.append({'id':'OPEN-CONDITION-'+condition_id,'state':'CONFLICT',
                             'message':'BENCHMARK and G103A share a Condition ID','source':condition_id})
        detail=(group,condition_id,_fmt(cond['mach']),_fmt(cond['re_cref']),
            _fmt(cond['alpha_deg']),_fmt(cond['beta_deg']),
            ', '.join(f'{k}={_fmt(v)}°' for k,v in sorted(controls.items())),
            cond['trim_mode'], ', '.join(x['case']['case_id'] for x in matching))
        condition_rows.append(detail)
        condition_audit.append({**dict(cond),'group':group,
                                'case_ids':[s['case']['case_id'] for s in matching],
                                'run_ids':[s['run']['run_id'] for s in matching]})
        _field(fields,'4.4.condition.'+condition_id,' / '.join(_fmt(x) for x in detail[2:8]),'RECORDED',
               [{'source':'flight_condition.'+k,'value':cond[k],'condition_id':condition_id}
                for k in ('mach','re_cref','alpha_deg','beta_deg','sref','bref','cref','axes_convention_id')])
    refs=[x for x in ext if x['kind'] in EXTERNAL_KINDS]
    statuses=[]
    for t in db.conn.execute('SELECT test_id,requirement FROM test_definition ORDER BY test_id'):
        cases=[x for x in selected if x['case']['test_id']==t['test_id']]
        statuses.append((t['test_id'], len(cases),sum(x['run'] is not None for x in cases)))
    manual=profile.get('document_control',{})
    baseline=profile.get('baseline_case_id')
    if baseline and not any(x['case']['case_id']==baseline and x['run'] for x in selected):
        raise ValueError('baseline_case_id is not a selected successful Case: '+baseline)
    baseline_entry=next((x for x in successes if x['case']['case_id']==baseline),None)
    candidate=next((x for x in successes if x['case']['case_id']=='GEO-03-G3'),None)
    if not baseline:
        _field(fields,'4.3.baseline','基準Case未確定。候補：'+
              ('GEO-03-G3' if candidate else 'なし'),'MANUAL_DECISION')
    settings_audit=[]
    per_case_settings={}
    for item in successes:
        run_settings=_record_settings(item,fields)
        per_case_settings[item['case']['case_id']]=run_settings
        for label,requested,effective,generated,state,rid in run_settings:
            settings_audit.append({'case_id':item['case']['case_id'],'run_id':rid,
                'setting':label,'requested':requested,'effective':effective,
                'generated':generated,'state':state})
    env_rows=[]
    for rid in run_ids:
        r=next(s['run'] for s in successes if s['run']['run_id']==rid)
        env=_json(r['environment_json'])
        env_rows.append({'run_id':rid,'binary_sha256':r['binary_sha256'],
            'runner_sha256':r['runner_sha256'],'ncpu':env.get('ncpu'),
            'git_commit':env.get('git_commit'),
            'execution_source_sha256':env.get('execution_source_sha256'),
            'postprocess_source_sha256':env.get('postprocess_source_sha256'),
            'binary_path':env.get('binary_path')})
    distinct_bins=sorted({r['binary_sha256'] for r in env_rows if r['binary_sha256']})
    distinct_exec=sorted({r['execution_source_sha256'] for r in env_rows if r['execution_source_sha256']})
    distinct_threads=sorted({str(r['ncpu']) for r in env_rows if r['ncpu'] is not None})
    distinct_git=sorted({r['git_commit'] for r in env_rows if r['git_commit']})
    # The report's actual source is its own frozen snapshot, not current Git HEAD.
    chapter1=['**採用Runに基づく実績概要（技術認定ではありません）**','',
         _table(('項目','実績／状態'),[
             ('機体','G103A Twin II Acro（登録形状の形態別照合が必要）'),
             ('Intended Use','未承認。Test Familyから用途候補を提示するのみ'),
             ('G103A実行済み離散条件',
                'Mach='+(','.join(sorted({_fmt(c['mach']) for c in condition_audit if c['group']=='G103A'})) or '記録なし')+
                '; α='+(','.join(sorted({_fmt(c['alpha_deg']) for c in condition_audit if c['group']=='G103A'})) or '記録なし')+'°'),
             ('代表取得QoI',', '.join(sorted({q['qoi'] for e in successes for q in e['case'].get('quantities',[])
                                             if q['qoi'] in ('CLi','CDi','CLiw','CDiw','CMiy')})) or 'なし'),
             ('計画Case / 採用成功Case',f'{len(selected)} / {len(successes)}'),
             ('全Run / 成功Run（履歴）',f'{sum(run_stats.values())} / {run_stats.get("SUCCEEDED",0)}'),
             ('採用Quantity / DB全Quantity',f'{qty_selected} / {qty_all}'),
             ('外部理論・風洞・他コード参照',f'{sum(x["n"] for x in refs)} Dataset（今回のDB）'),
             ('Comparison / Rule / Assessment',f'{counts["comparison"]} / {counts["acceptance_rule"]} / {counts["assessment"]}'),
             ('工学評価 / 使用認定','NOT_EVALUATED / INSUFFICIENT_EVIDENCE')]),
         '',f'選択成功Run数：{len(run_ids)}種類。以下の条件表は離散的な実施条件であり、連続包絡の認定ではありません。',
         '']
    missing=[x for x in selected if x['run'] is None]
    unresolved=[x for x in missing if x['case'].get('status')=='UNRESOLVED']
    not_succeeded=[x for x in missing if x not in unresolved]
    if not_succeeded:
        warnings.append({'id':'OPEN-NO-SUCCESSFUL-RUN','state':'NOT_EXECUTED',
            'message':f'{len(not_succeeded)} Caseに採用成功Runなし',
            'source':'results_snapshot.cases[].run_id',
            'case_ids':[x['case']['case_id'] for x in not_succeeded]})
    if unresolved:
        warnings.append({'id':'OPEN-UNRESOLVED','state':'NOT_EXECUTED',
            'message':f'{len(unresolved)} Caseが未解決',
            'source':'planned_case.execute_policy',
            'case_ids':[x['case']['case_id'] for x in unresolved]})
    if not refs: warnings.append({'id':'OPEN-REFERENCE','state':'REVIEW_REQUIRED','message':'独立外部参照Dataset未登録','source':'dataset.kind'})
    if not counts['acceptance_rule']: warnings.append({'id':'OPEN-RULE','state':'REVIEW_REQUIRED','message':'承認用Rule未登録','source':'acceptance_rule'})
    if not counts['assessment']: warnings.append({'id':'OPEN-ASSESSMENT','state':'MANUAL_DECISION','message':'Applicability Assessment未登録','source':'assessment'})
    if any(not x['git_commit'] for x in env_rows): warnings.append({'id':'OPEN-EXEC-GIT','state':'NOT_RECORDED','message':'一部Runで実行時Git Commit未記録','source':'run.environment_json.git_commit'})
    if not baseline: warnings.append({'id':'OPEN-BASELINE','state':'MANUAL_DECISION','message':'正式な基準Case未指定','source':'report_profile.baseline_case_id'})
    chapter1+=['### 未解決事項（自動診断ID。正式Action IDではありません）','',
               _table(('診断ID','状態','記録に基づく内容'),
                      [(x['id'],x['state'],x['message']) for x in warnings[:18]]),'',
               '全診断・出典：[report_field_provenance.json](report_field_provenance.json)。']
    chapters['1.1']='\n'.join(chapter1[:4])+'\n'
    chapters['1.2']=_table(('Decision ID','関連Test','採用成功Case','判定'),[
        (dec, ', '.join(families),
         sum(n for tid,_,n in statuses if any(tid.startswith(f+'-') for f in families)),
         'NOT_EVALUATED' if dec in ('DEC-01','DEC-02') else 'INSUFFICIENT_EVIDENCE')
        for dec,families in [('DEC-01',('GRID',)),('DEC-02',('WAKE',)),
            ('DEC-03',('GEO','VAL')),('DEC-04',('STA','CTRL')),('DEC-05',('DRAG',))]])
    chapters['1.3']='\n'.join(chapter1[4:])+'\n'
    chapters['1']='**採用Runベースの事実を下記1.1～1.3に自動入力します。解析成功と技術承認は別です。**'
    chapters['2.1']='**計画されたTest Familyに基づく用途候補（未承認）**：VER・GEO・GRID・WAKE・XCD・VAL・STA・CTRL・DRAG・ACC。必要精度と未達時影響は設計要求として別途定義が必要です。\n'
    chapters['2.2']=_table(('Geometry ID','DB Stage','Case Role','記録上の構成','同等性','実行状況','状態'),geometry_rows)+ '\n**注意**：Stage／Role等の不整合は自動解消しません。\n'
    envelope_rows=[]
    for group in ('G103A','BENCHMARK'):
        group_conditions=[x for x in condition_audit if x['group']==group]
        all_conditions=[]
        for record in db.conn.execute("SELECT fc.*,pc.test_id FROM flight_condition fc JOIN planned_case pc ON pc.condition_id=fc.condition_id"):
            cid=record['condition_id']
            if cid not in {r['condition_id'] for r in all_conditions}:
                if record['test_id'].startswith('VER-') == (group=='BENCHMARK'):
                    all_conditions.append(dict(record))
        for label,key in (('Mach','mach'),('Re (Cref)','re_cref'),('α [deg]','alpha_deg'),('β [deg]','beta_deg')):
            planned=', '.join(sorted({_fmt(x[key]) for x in all_conditions})) or 'DB未記録'
            executed=', '.join(sorted({_fmt(x[key]) for x in group_conditions})) or '実行実績なし'
            envelope_rows.append((group,label,planned,executed,'未認定'))
    chapters['2.3']=_table(('群','評価軸','計画上の離散値（参考）','採用成功の離散値','認定済み範囲'),envelope_rows)
    chapters['2.3']+='\n外部独立参照：'+('未登録' if not refs else str(sum(x['n'] for x in refs))+' Dataset（適用性未審査）')
    chapters['2.3']+='。**離散的な成功条件を連続した認定包絡と解釈しません。**\n'
    chapters['3.1']=_table(('Test Family','計画Test','計画Case','採用成功Case'),[
        (name,sum(tid.startswith(name+'-') for tid,_,_ in statuses),
         sum(n for tid,n,_ in statuses if tid.startswith(name+'-')),
         sum(ok for tid,_,ok in statuses if tid.startswith(name+'-')))
        for name in ('VER','GEO','GRID','WAKE','XCD','VAL','STA','CTRL','DRAG','ACC')])
    qgroups=defaultdict(set)
    for s in successes:
        for q in s['case'].get('quantities',[]):
            qgroups[q['qoi']].add((q.get('unit'),q.get('variant'),q.get('component')))
    chapters['3.2']=f'採用成功Runから取得したQoIは **{len(qgroups)}種類**、Quantityは **{qty_selected}件**です。'
    chapters['3.2']+=' 単位`1`は軸・符号・正規化の確認済みを意味しません。\n\n'
    chapters['3.2']+=_table(('代表QoI','Unit / Variant / Component（実績）','技術規約'),[
        (q,', '.join('/'.join(map(_fmt,v)) for v in sorted(vals,key=str)[:3]),'要技術確認')
        for q,vals in list(sorted(qgroups.items()))[:24]])
    qids=[q['quantity_id'] for e in successes for q in e['case'].get('quantities',[])]
    qmetadata=[]
    for qid in qids[:120]:
        r=_row(db,'quantity','quantity_id',qid)
        if r is not None and r['qoi'] in ('CLi','CDi','CLiw','CDiw','CMix','CMiy','CMiz'):
            qmetadata.append((r['qoi'],r['native_name'],r['native_unit'],
                r['canonical_unit'] or '未記録',r['convention_id'],
                r['normalization_method_id']))
    if qmetadata:
        chapters['3.2']+='\n**代表Quantityの原記録（係数定義の外部整合性は未確認）**\n\n'
        chapters['3.2']+=_table(('QoI','Native name','Native unit','Canonical unit','Convention ID','Norm method'),
                                list(dict.fromkeys(qmetadata))[:12])
    chapters['3.3']=f'Comparison登録：**{counts["comparison"]}件**。不確かさの推算・適用条件は未評価です。式の存在は計算実績を意味しません。\n'
    chapters['3.4']=f'Acceptance Rule：**{counts["acceptance_rule"]}件**、Assessment：**{counts["assessment"]}件**。SUCCEEDEDをEngineering PASSへ変換しません。\n'
    rules=[dict(x) for x in db.conn.execute('SELECT rule_id,qoi,level,approval,reviewer,approved_utc FROM acceptance_rule ORDER BY rule_id')]
    if rules:
        chapters['3.4']+='\n**登録Ruleの承認記録（適用条件は別途審査）**\n\n'
        chapters['3.4']+=_table(('Rule ID','QoI','評価階層','DB承認状態','Reviewer','承認日時'),[
            (r['rule_id'],r['qoi'],r['level'],r['approval'],r['reviewer'],r['approved_utc']) for r in rules])
    else:
        chapters['3.4']+='\nRule未登録（NOT_EVALUATED）。\n'
    assessments=[dict(x) for x in db.conn.execute('SELECT assessment_id,numerical_status,engineering_status,disposition,approved_utc FROM assessment ORDER BY assessment_id')]
    if assessments:
        chapters['3.4']+='\n**Assessmentの独立した状態**\n\n'
        chapters['3.4']+=_table(('Assessment ID','Numerical','Engineering','Applicability/Disposition','承認日時'),[
           (x['assessment_id'],x['numerical_status'],x['engineering_status'],x['disposition'],x['approved_utc'])
           for x in assessments])

    chapters['3.5']=_table(('Kind','Usage','Evidence Class','Dataset数'),[
        (r['kind'],r['usage_class'],r['evidence_class'] or '未記録',r['n']) for r in ext]) + \
        '\nWIND_TUNNEL／NS／PANEL／THEORYの独立参照を、EXTERNAL_CD0で代替しません。\n'
    report_fingerprint=snapshot.get('visualization_source_sha256') if snapshot else None
    chapters['4.1']=_table(('項目','実行時の記録（採用Run）','状態／注意'),[
        ('計画対象OpenVSP版',study['target_openvsp_version'],'実行時のAPI版を実測した証拠ではない'),
        ('実行VSPAERO binary SHA-256',', '.join(x[:16]+'…' for x in distinct_bins) or '未記録','完全値はJSONに収録'),
        ('実行コード指紋',', '.join(x[:16]+'…' for x in distinct_exec) or '未記録','現行Git HEADではない'),
        ('実行時Git Commit',', '.join(distinct_git) if distinct_git else 'DB未記録','nullを補完しない'),
        ('再生成側Visualization指紋',(report_fingerprint or 'DB未記録')[:16]+'…' if report_fingerprint else 'DB未記録','実行時コード指紋と別'),
        ('OpenVSP API実測版','DB未記録','要実機確認'),
        ('実行スレッド数',', '.join(distinct_threads) or 'DB未記録','ncpu記録'),
        ('Python / OS / CPU','DB未記録','過去環境を推定しない'),
        ('DBスキーマ版',study['schema_version'],'ソフトウェア版とは別'),
        ('計画書SHA-256',study['plan_sha256'][:16]+'…','DBの記録値／実ファイル照合は別'),
        ('Study設定SHA-256',study['config_sha256'][:16]+'…','DBの記録値')])+'\n実行コードと報告書生成コードの指紋は別途管理します。\n'
    # Include representative configuration only when baseline was explicitly chosen.
    chapters['4.2']=_table(('Geometry ID','Stage','Role','採用成功Run','同等性'),[
        (g['geometry_id'],g['stage'],_json(g['representation_json']).get('case_role','DB未記録'),
         ', '.join(s['run']['run_id'][:12]+'…' for s in successes if s['run']['geometry_id']==g['geometry_id']) or 'なし',
         g['geometry_equivalence_status']) for g in geometries if g['geometry_id'].startswith('GEO-')])
    set_rows=[]
    for entry in successes:
        if not entry['case']['test_id'].startswith('GEO-'):
            continue
        run=entry['run']; effective=_json(run['effective_snapshot_json'])
        req=_json(run['requested_snapshot_json']).get('setting',{})
        membership=_json(req.get('requested_json')).get('geometry_membership',{})
        sets=effective.get('geometry_sets',{})
        set_rows.append((entry['case']['case_id'],run['run_id'],
            'Thin '+_fmt(effective.get('ThinGeomSet'))+' / Thick '+_fmt(effective.get('GeomSet')),
            'Thin '+_fmt(membership.get('thin'))+' / Thick '+_fmt(membership.get('thick')),
            'Thin '+_fmt(sets.get('thin'))+' / Thick '+_fmt(sets.get('thick'))))
    chapters['4.2']+='\n**Geometry試験の採用成功RunにおけるSet構成（IDは未確認の部品名へ変換しない）**\n\n'
    chapters['4.2']+=_table(('Case','Run','Set Index (Effective)','Member (Requested)','Geometry ID (Effective)'),set_rows)
    meshrows=[]
    for entry in successes:
        if not entry['case']['test_id'].startswith(('GEO-','GRID-')): continue
        for key,value in sorted(entry['case'].get('mesh',{}).items()):
            if isinstance(value,(int,float)) and math.isfinite(value):
                meshrows.append((entry['case']['case_id'],key,_fmt(value),entry['run']['run_id']))
    if meshrows:
        chapters['4.2']+='\n**採用Runのメッシュ診断（設定Tessとは別）**\n\n'
        chapters['4.2']+=_table(('Case','診断指標','実績','Run ID'),meshrows[:15])
        if len(meshrows)>15:
            chapters['4.2']+=f'\n残り{len(meshrows)-15}件は試験別詳細レポートを参照。\n'
    chapters['4.2']+='\nメッシュ診断値は各試験レポートの「メッシュ・履歴」を参照。Tess入力値とは区別します。\n'

    displayed=baseline_entry or candidate
    if baseline_entry:
        chapters['4.3']=f'**指定された表示用代表Case：{baseline}**（Engineering Acceptanceではありません）\n\n'
    else:
        chapters['4.3']='**基準Case未確定（MANUAL_DECISION）**。候補：'+('`GEO-03-G3`' if candidate else '候補なし')+'。基準値を自動採用しません。\n'
    if displayed:
        chapters['4.3']+='\n**参考設定（未認定）**：'+displayed['case']['case_id']+'\n\n'
        chapters['4.3']+=_table(('設定項目','Requested','Effective','Generated','記録整合','Run ID'),
            per_case_settings[displayed['case']['case_id']])

    chapters['4.3']+='\n各成功Caseの全設定：[tables/report_settings_audit.csv](tables/report_settings_audit.csv)。Tess値が欠測の場合、計画書の開始値で補完しません。\n'
    chapters['4.4']=_table(('群','Condition ID','Mach','Re(Cref)','α [deg]','β [deg]','舵角 [deg]','Trim','採用Case'),condition_rows)+\
       '\n基準寸法・モーメント基準点・軸規約は[tables/report_conditions_audit.csv](tables/report_conditions_audit.csv)にCondition IDごとに記録。\n'
    for section, content in chapters.items():
        _field(fields,'section.'+section,'generated','RECORDED', [{'source':'frozen SQLite/results_snapshot','section':section}])
    return {'sections':chapters,'fields':fields,'warnings':warnings,
            'counts':{'selected_cases':len(selected),'selected_success':len(successes),
                      'selected_quantity':qty_selected,'all_quantity':qty_all,
                      'historical_runs':sum(run_stats.values()),'historical_success':run_stats.get('SUCCEEDED',0)},
            'condition_audit':condition_audit,'settings_audit':settings_audit,
            'environment_audit':env_rows,'profile':profile}


def write_report_facts(facts, directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    dest=directory/'report_field_provenance.json'
    dest.write_text(json.dumps({k:facts[k] for k in ('fields','warnings','counts','profile')},ensure_ascii=False,
                               indent=2,sort_keys=True)+'\n',encoding='utf-8')
    tables=directory/'tables';tables.mkdir(exist_ok=True)
    for name,rows in (('report_field_audit.csv',facts['fields']),
                      ('report_settings_audit.csv',facts['settings_audit']),
                      ('report_conditions_audit.csv',facts['condition_audit']),
                      ('report_environment_audit.csv',facts['environment_audit'])):
        path=tables/name
        keys=sorted({k for r in rows for k in r})
        with path.open('w',encoding='utf-8-sig',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader()
            for r in rows:
                writer.writerow({k:json.dumps(v,ensure_ascii=False,sort_keys=True) if isinstance(v,(dict,list)) else v
                                 for k,v in r.items()})
    return dest
