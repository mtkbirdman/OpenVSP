"""SQLite study, case, run and evidence store for the 3.50.5 V&V plan.

No OpenVSP import: metadata and offline post-processing work without a solver.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sqlite3
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

# Checkout root, independent of worker and notebook working directories.
REPO_ROOT = Path(__file__).resolve().parents[3]


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as file:
        for chunk in iter(lambda: file.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def canonical(obj):
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)


def signature(obj):
    return hashlib.sha256(canonical(obj).encode('utf-8')).hexdigest()


def resolve_ncpu(requested, override=None):
    """Resolve VSPAERO OpenMP thread count without modifying immutable Case inputs.

    Precedence: explicit CLI override, reviewed Case NCPU, legacy default 1.
    A valid value is a positive JSON integer (not bool or float).
    """
    def validate(value):
        if type(value) is not int or value < 1:
            raise ValueError('NCPU must be a positive integer (>= 1)')
        return value

    options=requested.get('vspaero', {})
    if not isinstance(options,dict):
        raise ValueError('VSPAERO requested settings must be a JSON object')
    configured = validate(options.get('NCPU', 1))
    return configured if override is None else validate(override)


def uid(prefix):
    return prefix + '-' + uuid.uuid4().hex[:16]


# Files capable of changing the physics execution, not reviewer/report formatting.
# Changing this manifest is itself tracked by its resulting digest.
EXECUTION_SOURCES = (
    'src/openvsppy/vv/VVRunner.py', 'src/openvsppy/vv/worker.py',
    'src/openvsppy/vv/VVDatabase.py',
    'src/openvsppy/vspaero/AnalysisVSPAERO.py',
    'src/openvsppy/vspaero/VSPAEROMesh.py',
    'src/openvsppy/vspaero/VSPAEROStab.py',
    'src/openvsppy/vv/VVNative.py', 'src/openvsppy/core/util.py',
    'src/openvsppy/vv/VVAutomation.py',
    'src/openvsppy/aerodynamics/AnalysisParasiteDrag.py',
)
POSTPROCESS_SOURCES = (
    'src/openvsppy/vv/VVPostprocess.py', 'src/openvsppy/vv/VVNative.py',
    'src/openvsppy/vspaero/VSPAEROStab.py',
)
VISUALIZATION_SOURCES = ('src/openvsppy/vv/VVPlots.py', 'src/openvsppy/vv/VVReport.py')


def source_fingerprint(repo_root, paths):
    repo_root = Path(repo_root)
    return signature({name: digest(repo_root / name) for name in paths})


def execution_fingerprint(repo_root):
    return source_fingerprint(repo_root, EXECUTION_SOURCES)


def postprocess_fingerprint(repo_root):
    return source_fingerprint(repo_root, POSTPROCESS_SOURCES)


def load_study(path):
    path = Path(path).resolve()
    d = json.loads(path.read_text(encoding='utf-8-sig'))
    if d.get('runtime', {}).get('openvsp_version_required') != '3.50.5':
        raise ValueError('Study must target OpenVSP 3.50.5.')
    plan = (path.parent / d['plan']['file']).resolve()
    if not plan.is_file() or digest(plan) != d['plan']['sha256']:
        raise ValueError('V&V plan file is absent or SHA-256 differs from study.json.')
    catalog = (path.parent / d['test_catalog_csv']).resolve()
    if not catalog.is_file():
        raise FileNotFoundError(catalog)
    return path, d, plan, catalog


# Test families are versioned in the plan/catalog, not inferred from model output.
CASE_ROLES = {
    'VER-01': ('COARSE', 'MEDIUM', 'FINE'), 'VER-02': ('COARSE', 'MEDIUM', 'FINE'),
    'GEO-01': ('G1',), 'GEO-02': ('G2',), 'GEO-03': ('G3',),
    'GEO-04': ('THIN', 'THICK'),
    'GRID-01': ('COARSE', 'MEDIUM', 'FINE'), 'GRID-02': ('COARSE', 'MEDIUM', 'FINE'),
    'GRID-03': ('BASE', 'LOCAL_A', 'LOCAL_B'),
    'WAKE-01': ('ITER_5', 'ITER_10', 'ITER_20'),
    'WAKE-02': ('NODES_8', 'NODES_12', 'NODES_16'),
    'WAKE-03': ('EXPLICIT', 'IMPLICIT'),
    'XCD-01': ('BASE',), 'XCD-02': ('BASE',),
    'VAL-01': ('ALPHA_NEG', 'ALPHA_0', 'ALPHA_POS'),
    'VAL-02': ('BETA_NEG', 'BETA_0', 'BETA_POS'),
    'STA-01': ('STEADY', 'ADJOINT'),
    'STA-02': ('MESH_BASE', 'MESH_FINE', 'WAKE_BASE', 'WAKE_FINE'),
    'STA-03': ('P', 'Q', 'R'),
    'CTRL-01': tuple(f'DELTA_{v:+d}' for v in (-10,-5,-2,0,2,5,10)),
    'CTRL-02': tuple(f'DELTA_{v:+d}' for v in (-10,-5,-2,0,2,5,10)),
    'CTRL-03': tuple(f'DELTA_{v:+d}' for v in (-10,-5,-2,0,2,5,10)),
    'DRAG-01': ('METHOD_A',), 'DRAG-02': ('METHOD_B',), 'ACC-01': (),
}

TRANSITIONS = {
    'PLANNED': {'PRECHECKED', 'PRECHECK_FAILED', 'CANCELLED'},
    'PRECHECKED': {'RUNNING', 'PRECHECK_FAILED', 'CANCELLED'},
    'RUNNING': {'CAPTURED', 'EXECUTION_FAILED', 'INTERRUPTED'},
    'CAPTURED': {'PARSED', 'OUTPUT_INVALID'},
    'PARSED': {'SUCCEEDED', 'OUTPUT_INVALID'},
}
FINAL = {'SUCCEEDED','PRECHECK_FAILED','EXECUTION_FAILED','OUTPUT_INVALID','INTERRUPTED','CANCELLED'}


class VVDatabase:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(self.path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute('PRAGMA foreign_keys=ON')
        self.conn.execute('PRAGMA busy_timeout=5000')

    def close(self):
        self.conn.close()

    @contextmanager
    def transaction(self):
        # Savepoints allow a complete Run import + final Case link to be atomic,
        # even though legacy insert helpers call transaction() themselves.
        if self.conn.in_transaction:
            name='vv_nested_' + uuid.uuid4().hex
            self.conn.execute('SAVEPOINT '+name)
            try:
                yield self.conn
                self.conn.execute('RELEASE SAVEPOINT '+name)
            except Exception:
                self.conn.execute('ROLLBACK TO SAVEPOINT '+name)
                self.conn.execute('RELEASE SAVEPOINT '+name)
                raise
        else:
            self.conn.execute('BEGIN IMMEDIATE')
            try:
                yield self.conn
                self.conn.commit()
            except Exception:
                self.conn.rollback()
                raise

    def row(self, table, key_column, value):
        if table not in ('study','planned_case','run','geometry','flight_condition','numerical_setting','dataset','quantity','test_definition'):
            raise ValueError('Disallowed table name')
        if key_column not in ('study_id','case_id','run_id','geometry_id','condition_id','setting_id','dataset_id','quantity_id','test_id'):
            raise ValueError('Disallowed key name')
        return self.conn.execute(f'SELECT * FROM {table} WHERE {key_column}=?', (value,)).fetchone()

    def initialize(self, study_file, schema_file):
        study_path, d, _, catalog = load_study(study_file)
        ddl = Path(schema_file).read_text(encoding='utf-8')
        if self.conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='study'").fetchone():
            raise ValueError('Already initialized; refusing to overwrite study.')
        self.conn.executescript(ddl)
        rows = list(csv.DictReader(catalog.open(encoding='utf-8-sig', newline='')))
        ids = {r['TestID'] for r in rows}
        if len(rows)!=25 or ids != set(CASE_ROLES):
            raise ValueError('Catalog must contain the 25 approved Test IDs exactly once.')
        with self.transaction() as c:
            c.execute('INSERT INTO study VALUES(?,?,?,?,?,?,?,?)',(
                d['study_id'], d['name'], d['plan']['version'], d['plan']['sha256'],
                d['runtime']['openvsp_version_required'], digest(study_path),
                d['schema_version'],utc()))
            for r in rows:
                c.execute('INSERT INTO test_definition VALUES(?,?,?,?,?,?,?)',(
                    r['TestID'],d['study_id'],r['Title'],r['Requirement'],
                    r['ExecutionMode'],r['TestID'].split('-')[0],r['EvaluationRule']))
                for role in CASE_ROLES[r['TestID']]:
                    c.execute('INSERT INTO planned_case VALUES(?,?,?,?,?,?,?,?)',(
                        r['TestID']+'-'+role,r['TestID'],role,None,None,None,
                        'UNRESOLVED',canonical({'requested_qoi':r['QoI'],'prerequisites':r['Prerequisites']})))
        return len(rows)

    def lock_evaluation_selection(self, test_id, selection, label):
        """Freeze an explicit case-to-Run selection; never silently update it."""
        if not label or not selection: raise ValueError('Selection label and at least one Run required')
        for case_id,run_id in selection.items():
            hit=self.conn.execute("""SELECT 1 FROM case_run_link l
                JOIN case_point p ON p.point_id=l.point_id JOIN run r ON r.run_id=p.run_id
                JOIN planned_case c ON c.case_id=l.case_id
                WHERE l.case_id=? AND c.test_id=? AND r.run_id=?
                  AND p.point_status='MATCHED' AND r.status='SUCCEEDED' """,
                (case_id,test_id,run_id)).fetchone()
            if hit is None:raise ValueError('Invalid selected Case/Run pair: '+case_id+' / '+run_id)
        snapshot=canonical(selection)
        existing=self.conn.execute('SELECT selection_json FROM evaluation_selection WHERE label=?',(label,)).fetchone()
        if existing:
            if existing[0]!=snapshot:raise ValueError('Selection already frozen with different Run IDs')
            return label
        with self.transaction() as c:
            c.execute('INSERT INTO evaluation_selection VALUES(?,?,?,?,?)',
                      (label,test_id,snapshot,signature(selection),utc()))
        return label

    def add_geometry(self, study_id, record):
        source=Path(record['vsp3_source']).resolve()
        if not source.is_file() or source.suffix.lower()!='.vsp3':
            raise ValueError('Geometry needs an existing .vsp3 file.')
        h=digest(source)
        prior=self.row('geometry','geometry_id',record['geometry_id'])
        if prior:
            if prior['source_vsp3_sha256']!=h or json.loads(prior['representation_json'])!=record.get('representation',{}) or json.loads(prior['physical_model_json']).get('model')!=record.get('physical_model',{}):
                raise ValueError('Geometry ID already has different immutable input contents')
            return
        with self.transaction() as c:
            c.execute('INSERT INTO geometry VALUES(?,?,?,?,?,?,?,?)',(
                record['geometry_id'],study_id,record['stage'],
                canonical({'source_file':str(source),'model':record.get('physical_model',{})}),
                canonical(record.get('representation',{})),h,
                'UNVERIFIED', record.get('remarks','')))

    def add_condition(self, study_id, record):
        for name in ('mach','re_cref','alpha_deg','beta_deg','sref','bref','cref','moment_x','moment_y','moment_z'):
            if name not in record or record[name] is None or not math.isfinite(float(record[name])):
                raise ValueError(f'Flight condition requires finite {name}.')
        if min(record[k] for k in ('re_cref','sref','bref','cref')) <= 0 or record['mach']<0:
            raise ValueError('Mach/Re/reference dimensions out of range.')
        fields=['mach','re_cref','alpha_deg','beta_deg','controls_deg','trim_mode','sref','bref','cref',
                'moment_x','moment_y','moment_z','reference_length_unit','axes_convention_id']
        snap={k:record[k] for k in fields}
        old=self.row('flight_condition','condition_id',record['condition_id'])
        if old:
            if old['condition_hash']!=signature(snap):
                raise ValueError('FlightCondition ID cannot be reused with changed values')
            return
        with self.transaction() as c:
            c.execute('INSERT INTO flight_condition VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(
                record['condition_id'],study_id,record['mach'],record['re_cref'],record['alpha_deg'],
                record['beta_deg'],canonical(record['controls_deg']),record['trim_mode'],
                record['sref'],record['bref'],record['cref'],record['moment_x'],record['moment_y'],
                record['moment_z'],record['reference_length_unit'],record['axes_convention_id'],signature(snap)))

    def add_setting(self, study_id, record):
        requested=record['requested']
        old=self.row('numerical_setting','setting_id',record['setting_id'])
        if old:
            if old['setting_hash']!=signature(record):
                raise ValueError('NumericalSetting ID cannot be reused with changed values')
            return
        with self.transaction() as c:
            c.execute('INSERT INTO numerical_setting VALUES(?,?,?,?,?,?,?)',(
                record['setting_id'],study_id,record['analysis_mode'], record['mesh_level'],
                canonical(requested),signature(record),'Imported from explicit case configuration'))

    def resolve_case(self, case_id, geometry_id, condition_id, setting_id, policy='EXECUTE'):
        if policy not in ('EXECUTE','REUSE','IMPORT','MANUAL'):
            raise ValueError('Invalid policy.')
        with self.transaction() as c:
            row=c.execute('SELECT * FROM planned_case WHERE case_id=?',(case_id,)).fetchone()
            if row is None: raise KeyError(case_id)
            if c.execute('SELECT 1 FROM case_run_link WHERE case_id=?',(case_id,)).fetchone():
                raise ValueError('Case already linked to an executed Run; cannot mutate Case definition')
            study=c.execute('SELECT study_id FROM test_definition WHERE test_id=?',(row['test_id'],)).fetchone()[0]
            for table,key,id_ in (('geometry','geometry_id',geometry_id),('flight_condition','condition_id',condition_id),('numerical_setting','setting_id',setting_id)):
                target=c.execute(f'SELECT study_id FROM {table} WHERE {key}=?',(id_,)).fetchone()
                if target is None or target[0]!=study: raise ValueError(f'{table} is missing or from different study.')
            c.execute('UPDATE planned_case SET geometry_id=?,condition_id=?,setting_id=?,execute_policy=? WHERE case_id=?',
                      (geometry_id,condition_id,setting_id,policy,case_id))

    def add_extra_planned_case(self, study_id, test_id, role):
        """Append a specifically configured sensitivity case, preserving legacy evidence.

        Existing 68 catalog cases are unchanged. SQLite backup precedes the first
        schema-data extension; a repeated prepare with the same ID is harmless.
        """
        case_id=f'{test_id}-{role}'
        if self.row('planned_case','case_id',case_id):return case_id
        if test_id!='GRID-03' or role!='NREF_1':
            raise ValueError('Additional planned role not explicitly supported')
        import datetime
        backup=self.path.with_name(self.path.name+'.pre_grid03_'+datetime.datetime.now(
            datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')+'.bak')
        with sqlite3.connect(backup) as destination:
            self.conn.backup(destination)
        with self.transaction() as connection:
            connection.execute('INSERT INTO planned_case VALUES(?,?,?,?,?,?,?,?)',
                (case_id,test_id,role,None,None,None,'UNRESOLVED',canonical({
                    'requested_qoi':'CLi,CMiy,CDiw,mesh_quality',
                    'prerequisites':'GRID-03 explicit NRef sensitivity level 1'})))
        return case_id

    def register_case_document(self, study_id, definition):
        """Register one immutable input triple then resolve a catalog Case."""
        self.add_geometry(study_id,definition['geometry'])
        self.add_condition(study_id,definition['flight_condition'])
        self.add_setting(study_id,definition['numerical_setting'])
        self.resolve_case(definition['case_id'],definition['geometry']['geometry_id'],
                          definition['flight_condition']['condition_id'],
                          definition['numerical_setting']['setting_id'],
                          definition.get('execute_policy','EXECUTE'))

    def register_run(self, case_id, study_root, binary_path, runner_path, requested, *, retry_of=None, ncpu=1):
        if type(ncpu) is not int or ncpu < 1:
            raise ValueError('NCPU must be a positive integer (>= 1)')
        row=self.row('planned_case','case_id',case_id)
        if not row or row['execute_policy']!='EXECUTE' or None in (row['geometry_id'],row['condition_id'],row['setting_id']):
            raise ValueError('Case not resolved for executable analysis.')
        binary=Path(binary_path).resolve()
        if not binary.is_file(): raise FileNotFoundError(binary)
        runner=Path(runner_path).resolve()
        if not runner.is_file(): raise FileNotFoundError(runner)
        runid=uid('RUN')
        root=Path(study_root).resolve()
        repo=REPO_ROOT
        code_hash=execution_fingerprint(repo)
        git_commit=None
        head=repo/'.git'/'HEAD'
        if head.is_file():
            ref=head.read_text(encoding='utf-8').strip()
            if ref.startswith('ref: '):
                target=repo/'.git'/ref[5:]
                if target.is_file():git_commit=target.read_text(encoding='utf-8').strip()
            else:git_commit=ref
        rundir=root/'03_run_raw'/runid
        rundir.mkdir(parents=True,exist_ok=False)
        st=self.conn.execute('SELECT study_id FROM test_definition WHERE test_id=?',(row['test_id'],)).fetchone()[0]
        f=signature({'geometry_id':row['geometry_id'],'condition_id':row['condition_id'],
                     'setting_id':row['setting_id'],'binary_sha256':digest(binary), 'execution_source_sha256':code_hash, 'requested':requested, 'ncpu':ncpu})
        with self.transaction() as c:
            c.execute('INSERT INTO run VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(
                runid,st,row['geometry_id'],row['condition_id'],row['setting_id'],'PLANNED',
                str(rundir.relative_to(root)),f,digest(binary),digest(runner),canonical({**requested,'runtime':{'ncpu':ncpu}}),
                None,None,canonical({'binary_path':str(binary),'ncpu':ncpu,'execution_source_sha256':code_hash,'postprocess_source_sha256':postprocess_fingerprint(REPO_ROOT),'git_commit':git_commit}),None,None,utc(),None,None,retry_of))
            self._event(c,runid,None,'PLANNED','CREATED','Run registered')
        return runid,rundir

    @staticmethod
    def _event(c,rid,prev,next_,code,msg,details=None):
        c.execute('INSERT INTO run_event(run_id,occurred_utc,previous_status,next_status,event_type,event_code,message,details_json) VALUES(?,?,?,?,?,?,?,?)',
                  (rid,utc(),prev,next_,'STATUS',code,msg,canonical(details or {})))

    def transition(self, runid, next_, *, code='OK', message='', effective=None, generated=None, exit_code=None):
        with self.transaction() as c:
            row=c.execute('SELECT status FROM run WHERE run_id=?',(runid,)).fetchone()
            if row is None: raise KeyError(runid)
            prev=row[0]
            if next_ not in TRANSITIONS.get(prev,set()):
                raise ValueError(f'Invalid Run transition: {prev} -> {next_}.')
            if next_=='SUCCEEDED':
                if not c.execute("SELECT 1 FROM case_point WHERE run_id=? AND point_status='MATCHED'",(runid,)).fetchone():
                    raise ValueError('Cannot succeed without a verified matching result point.')
            c.execute('UPDATE run SET status=?,effective_snapshot_json=COALESCE(?,effective_snapshot_json), generated_snapshot_json=COALESCE(?,generated_snapshot_json),exit_code=COALESCE(?,exit_code),error_code=?,started_utc=CASE WHEN ?=\'RUNNING\' THEN ? ELSE started_utc END,ended_utc=CASE WHEN ? IN (\'SUCCEEDED\',\'PRECHECK_FAILED\',\'EXECUTION_FAILED\',\'OUTPUT_INVALID\',\'INTERRUPTED\',\'CANCELLED\') THEN ? ELSE ended_utc END WHERE run_id=?',
                      (next_,None if effective is None else canonical(effective),None if generated is None else canonical(generated),exit_code,None if code=='OK' else code,next_,utc(),next_,utc(),runid))
            self._event(c,runid,prev,next_,code,message)

    def register_artifact(self, study_root, study_id, path, *, run_id=None,dataset_id=None,role='RAW'):
        root=Path(study_root).resolve(); path=Path(path).resolve()
        if not path.is_file() or not path.is_relative_to(root):
            raise ValueError('Artifact must be an existing file under study root.')
        rel=path.relative_to(root).as_posix()
        if run_id:
            r=self.row('run','run_id',run_id)
            if not r or not path.is_relative_to(root/r['run_dir']):
                raise ValueError('Run artifact must be inside its run directory.')
        with self.transaction() as c:
            found=c.execute('SELECT artifact_id,sha256 FROM artifact WHERE relative_path=?',(rel,)).fetchone()
            h=digest(path)
            if found:
                if found['sha256']!=h: raise ValueError('Artifact was modified after registration.')
                return found['artifact_id']
            id_=uid('ART')
            c.execute('INSERT INTO artifact VALUES(?,?,?,?,?,?,?,?,?,?,?)',(
                id_,study_id,run_id,dataset_id,role,rel,h,path.stat().st_size,None,utc(),'VERIFIED'))
            return id_

    def verify_artifacts(self, study_root):
        root=Path(study_root).resolve();result=[]
        for a in self.conn.execute('SELECT artifact_id,relative_path,sha256 FROM artifact'):
            p=(root/a['relative_path']).resolve()
            good=p.is_relative_to(root) and p.is_file() and digest(p)==a['sha256']
            result.append((a['relative_path'],'VERIFIED' if good else 'CHANGED_OR_MISSING'))
        return result
