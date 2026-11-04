"""Parse native VSPAERO files and evaluate evidence, without importing OpenVSP.

The native coefficient values are never overwritten. Missing uncertainty or
unapproved acceptance limits lead to INCONCLUSIVE, not PASS.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from .VVDatabase import canonical, digest, uid, utc, REPO_ROOT

POLAR_CONDITIONS = {'Beta':'beta_deg','Mach':'mach','AoA':'alpha_deg','Re/1e6':'re_million'}
AERO_QOI = {'CLi','CDi','CLiw','CDiw','CDo','CLo','CLtot','CDtot',
            'CLwtot','CDwtot','CMiy','CMix','CMiz','CMytot','CMxtot','CMztot',
            'CSi','CStot','L/D','E','StallFactor'}


def read_vspaero_settings(path):
    """Read generated .vspaero key/value lines without executing expressions."""
    values={}
    for line in Path(path).read_text(encoding='utf-8-sig', errors='replace').splitlines():
        if '=' not in line: continue
        k,v=(piece.strip() for piece in line.split('=',1))
        if not k or not v: continue
        try:
            num=float(v.replace('D','E'))
            values[k]=num if math.isfinite(num) else v
        except ValueError:
            values[k]=v
    return values


def read_polar(path):
    """Read whitespace-delimited VSPAERO .polar; preserve every native column.

    The file has a human-readable banner before the real column header; never
    assume that the first or second physical line is the header.
    """
    lines=Path(path).read_text(encoding='utf-8-sig',errors='replace').splitlines()
    header=None;rows=[]
    for line in lines:
        parts=line.split()
        if len(parts)>10 and parts[:3]==['Beta','Mach','AoA']:
            header=parts
            continue
        if header is None or not parts:
            continue
        if len(parts)!=len(header):
            if rows: break
            continue
        try:
            values=[float(x.replace('D','E')) for x in parts]
        except ValueError:
            if rows: break
            continue
        # A native derived ratio such as L/D may be NaN at zero lift.  Keep
        # valid coefficients; represent each nonfinite column as unavailable.
        # Flight-point identifiers must remain finite to establish provenance.
        row={name:(value if math.isfinite(value) else None)
             for name,value in zip(header,values)}
        invalid_conditions=[name for name in POLAR_CONDITIONS if name in row and row[name] is None]
        if invalid_conditions:
            raise ValueError('Nonfinite polar condition columns: '+str(sorted(invalid_conditions)))
        rows.append(row)
    if not header or not rows:
        raise ValueError('No complete VSPAERO polar header and data rows found.')
    if len(set(header))!=len(header):
        raise ValueError('Duplicate polar column names.')
    return rows


def condition_match(observed, target, *, atol_deg=1e-6, rel=1e-7):
    for obs,target_key in (('Mach','mach'),('AoA','alpha_deg'),('Beta','beta_deg')):
        if obs not in observed or abs(observed[obs]-float(target[target_key])) > (atol_deg if obs!='Mach' else rel):
            return False
    if 'Re/1e6' in observed:
        if not math.isclose(observed['Re/1e6']*1e6, float(target['re_cref']), rel_tol=rel, abs_tol=1.0):
            return False
    return True


def insert_dataset(db, *, study_id, kind, source, metadata=None, point_id=None, condition_id=None,
                   geometry_id=None, usage='REFERENCE', evidence=None, method='RAW', version='1'):
    ds=uid('DS')
    with db.transaction() as c:
        c.execute('INSERT INTO dataset VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(
            ds,study_id,kind,point_id,condition_id,geometry_id,usage,evidence,
            'NATIVE' if kind=='VSPAERO' else 'IMPORTED',str(source),canonical(metadata or {}),
            method,version,utc()))
    return ds


def insert_quantity(db, dataset, qoi, value, *, variant='DEFAULT',unit='1',native=None,
                    component='TOTAL',artifact=None,locator=None,method='NATIVE',status='AVAILABLE'):
    if status=='AVAILABLE' and (value is None or not math.isfinite(float(value))):
        raise ValueError(f'Invalid numerical value: {qoi}')
    if value is not None and not math.isfinite(float(value)):
        raise ValueError(f'Nonfinite native numerical value: {qoi}')
    ident=uid('Q')
    with db.transaction() as c:
        c.execute('INSERT INTO quantity VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(
            ident,dataset,qoi,component,variant,native or qoi,
            None if value is None else float(value),unit,
            None if status!='AVAILABLE' else float(value),unit if status=='AVAILABLE' else None,
            'VSPAERO_NATIVE' if method=='NATIVE' else 'SOURCE_DECLARED',method,artifact,locator,status))
    return ident


def import_polar(db, *, study_root, run_id, path, tol_deg=1e-6):
    """Store each native output point; only a single matching point may complete a run."""
    run=db.row('run','run_id',run_id)
    if run is None: raise KeyError(run_id)
    cond=db.row('flight_condition','condition_id',run['condition_id'])
    polar=read_polar(path)
    aid=db.register_artifact(study_root,run['study_id'],path,run_id=run_id,role='POLAR_RAW')
    if len(polar)!=1:
        raise ValueError('1 Run = 1 polar point; multiple output rows retained in artifact but rejected.')
    row=polar[0]
    matched=condition_match(row,cond,atol_deg=tol_deg)
    point=uid('PT')
    with db.transaction() as c:
        c.execute('INSERT INTO case_point VALUES(?,?,?,?,?,?)',(
            point,run_id,0,run['condition_id'],canonical(row),'MATCHED' if matched else 'MISMATCH'))
    ds=insert_dataset(db,study_id=run['study_id'],kind='VSPAERO',source='VSPAERO .polar',
                      point_id=point,condition_id=run['condition_id'],geometry_id=run['geometry_id'],
                      metadata={'native_columns':list(row),'artifact_id':aid},method='POLAR_NATIVE')
    unavailable=[]
    for col,val in row.items():
        if col in POLAR_CONDITIONS: continue
        if val is None:
            unavailable.append(col)
        insert_quantity(db,ds,col,val,native=col,artifact=aid,
                        locator=f'row=0,col={col}',
                        status='NOT_AVAILABLE' if val is None else 'AVAILABLE',
                        method='NONFINITE_POLAR' if val is None else 'NATIVE')
    return {'dataset_id':ds,'point_id':point,'matched':matched,
            'qoi_count':len(row)-len(POLAR_CONDITIONS),'nonfinite_columns':unavailable}


def import_stab(db, *, study_root, run_id, path, control_mapping=None,
                mode='STEADY', point_id=None, method='STAB_NATIVE', code_sha=None):
    """Import a stability-only Run without a spurious .polar prerequisite.

    Raw values are NATIVE_ONLY until conventions and normalization are approved.
    For finite-difference .stab, verify every native control-group label against
    an explicitly reviewed physical group name; do not infer group numbers.
    """
    from .VVNative import read_mode_derivatives, check_stability_condition
    header, rows, native_groups = read_mode_derivatives(path, mode)
    run=db.row('run','run_id',run_id)
    if run is None:raise KeyError(run_id)
    cond=db.row('flight_condition','condition_id',run['condition_id'])
    control_mapping=control_mapping or {}
    if mode=='STEADY':
        native_columns={variant for _,variant,_ in rows if variant.startswith('ConGrp_')}
        if native_columns != set(control_mapping):
            raise ValueError('Control group mappings not verified: '+str(sorted(native_columns-set(control_mapping))))
        if any(native_groups.get(col)!=name for col,name in control_mapping.items()):
            raise ValueError('Approved control group differs from native .stab group name')
        if any(not isinstance(name,str) or not name for name in control_mapping.values()):
            raise ValueError('Control mapping must have approved group names')
    check_stability_condition(header,cond)
    if point_id is None:
        points=list(db.conn.execute('SELECT point_id FROM case_point WHERE run_id=?',(run_id,)))
        if points:
            if len(points)!=1:raise ValueError('Ambiguous native Run point')
            point_id=points[0]['point_id']
        else:
            point_id=uid('PT')
            with db.transaction() as c:
                c.execute('INSERT INTO case_point VALUES(?,?,?,?,?,?)',
                    (point_id,run_id,0,run['condition_id'],canonical(header),'MATCHED'))
    aid=db.register_artifact(study_root,run['study_id'],path,run_id=run_id,role='STABILITY_RAW')
    meta={'native_header':header,'groups':native_groups,
          'approved_mapping':control_mapping,'source_mode':mode,'artifact_id':aid,
          'physical_convention':'NATIVE_NOT_REVIEWED'}
    if code_sha:meta['postprocess_source_sha256']=code_sha
    ds=insert_dataset(db,study_id=run['study_id'],kind='VSPAERO',
        source='VSPAERO '+Path(path).suffix,point_id=point_id,
        condition_id=run['condition_id'],geometry_id=run['geometry_id'],
        metadata=meta,method=method,version=code_sha[:12] if code_sha else '1')
    for coef,variant,value in rows:
        insert_quantity(db,ds,coef,value,variant=variant,unit='NATIVE_UNVERIFIED',
            native=f'{coef}/{variant}',artifact=aid,
            locator=f'mode={mode},coef={coef},column={variant}',
            method='UNVERIFIED_STABILITY',status='NATIVE_ONLY')
    return {'dataset_id':ds,'point_id':point_id,'matched':True,'quantity_count':len(rows)}


def import_reference_csv(db, *, study_id, csv_path, metadata):
    """Import long-form CSV; require explicit convention/conditions metadata.

    Required CSV columns: qoi,value,unit; optional component,variant.
    This never silently declares matching conditions or approval.
    """
    kind=metadata.get('kind')
    if kind not in ('PANEL','NS','WIND_TUNNEL','THEORY','EXTERNAL_CD0'):
        raise ValueError('Reference kind must be explicit.')
    for key in ('source_id','axes_convention_id','sref','bref','cref','mach','re_cref','alpha_deg','beta_deg','uncertainty_status','usage_class'):
        if key not in metadata or metadata[key] is None:
            raise ValueError('Missing reference metadata: '+key)
    if metadata['usage_class'] not in ('CALIBRATION','VALIDATION','REFERENCE','NOT_APPLICABLE'):
        raise ValueError('Unknown usage_class.')
    if kind!='WIND_TUNNEL' and metadata['usage_class']=='VALIDATION':
        raise ValueError('Nonexperimental reference cannot be classified as experimental validation.')
    with Path(csv_path).open(encoding='utf-8-sig',newline='') as f:
        data=list(csv.DictReader(f))
    if not data or not {'qoi','value','unit'} <= set(data[0]):
        raise ValueError('Expected long CSV columns qoi,value,unit.')
    for rec in data:
        if not rec['qoi'].strip() or not rec['unit'].strip():
            raise ValueError('Missing QoI/unit')
        if not math.isfinite(float(rec['value'])):
            raise ValueError('Reference contains nonfinite value')
    ds=insert_dataset(db,study_id=study_id,kind=kind,source=metadata['source_id'],metadata=metadata,
                      usage=metadata['usage_class'],evidence=metadata.get('evidence_class'),method='REFERENCE_CSV')
    for rec in data:
        if not rec['qoi'].strip() or not rec['unit'].strip(): raise ValueError('Missing QoI/unit')
        insert_quantity(db,ds,rec['qoi'].strip(),float(rec['value']),unit=rec['unit'].strip(),
                        native=rec['qoi'].strip(),variant=rec.get('variant','DEFAULT') or 'DEFAULT',
                        component=rec.get('component','TOTAL') or 'TOTAL',method='REFERENCE_CSV')
    return {'dataset_id':ds,'count':len(data)}


def drag_a(cdiw,cd0):
    return float(cdiw)+float(cd0)


def drag_b(cdiw,cdo,cdo_baseline,cd0):
    return float(cdiw)+float(cdo)-float(cdo_baseline)+float(cd0)


def create_drag_dataset(db, *, study_id, polar_dataset, cd0_dataset, cdo_baseline=None, baseline_metadata=None):
    """Store A; optionally B only when a documented baseline was supplied."""
    def quantity(dataset,qoi):
        r=db.conn.execute("SELECT * FROM quantity WHERE dataset_id=? AND qoi=? AND variant='DEFAULT' AND value_status='AVAILABLE'",(dataset,qoi)).fetchone()
        if r is None: raise ValueError(f'Missing {qoi} in {dataset}')
        return r
    base=db.row('dataset','dataset_id',polar_dataset)
    cd_ds=db.row('dataset','dataset_id',cd0_dataset)
    if not base or not cd_ds or base['kind']!='VSPAERO' or cd_ds['kind']!='EXTERNAL_CD0':
        raise ValueError('Expected VSPAERO and EXTERNAL_CD0 datasets.')
    cdiw=quantity(polar_dataset,'CDiw')
    cd0=quantity(cd0_dataset,'CD0_build')
    meta=json.loads(cd_ds['metadata_json'])
    if any(meta.get(k) is None for k in ('sref','mach','re_cref')):
        raise ValueError('CD0 reference area, Mach, Re are required.')
    condition=db.row('flight_condition','condition_id',base['condition_id'])
    if not math.isclose(float(meta['sref']),condition['sref'],rel_tol=1e-8) or not math.isclose(float(meta['mach']),condition['mach'],abs_tol=1e-6) or not math.isclose(float(meta['re_cref']),condition['re_cref'],rel_tol=1e-7):
        raise ValueError('CD0 condition/reference area does not match solver polar.')
    ds=insert_dataset(db,study_id=study_id,kind='DERIVED',source='External CD0 correction',
                      point_id=base['point_id'],condition_id=base['condition_id'],
                      geometry_id=base['geometry_id'],metadata={'polar_dataset':polar_dataset,'cd0_dataset':cd0_dataset},
                      method='CD0_CORRECTION',version='1')
    aid=insert_quantity(db,ds,'CD_A',drag_a(cdiw['canonical_value'],cd0['canonical_value']),method='CD0_A_1')
    with db.transaction() as c:
        for input_,role in ((cdiw,'CDiw'),(cd0,'CD0_build')):
            c.execute('INSERT INTO quantity_dependency VALUES(?,?,?)',(aid,input_['quantity_id'],role))
    bid=None
    if cdo_baseline is not None:
        if not baseline_metadata or not baseline_metadata.get('condition_id'):
            raise ValueError('Method B baseline condition provenance must be specified.')
        cdo=quantity(polar_dataset,'CDo')
        bid=insert_quantity(db,ds,'CD_B',drag_b(cdiw['canonical_value'],cdo['canonical_value'],cdo_baseline,cd0['canonical_value']),method='CD0_B_1')
        # A numeric documented baseline is stored in Dataset metadata by caller; no fabricated quantity IDs.
        with db.transaction() as c:
            c.execute('UPDATE dataset SET metadata_json=? WHERE dataset_id=?',(
                canonical({'polar_dataset':polar_dataset,'cd0_dataset':cd0_dataset,
                           'cdo_baseline':cdo_baseline,'baseline_metadata':baseline_metadata}),ds))
            for q,role in ((cdiw,'CDiw'),(cd0,'CD0_build'),(cdo,'CDo')):
                c.execute('INSERT INTO quantity_dependency VALUES(?,?,?)',(bid,q['quantity_id'],role))
    return {'dataset_id':ds,'CD_A_quantity_id':aid,'CD_B_quantity_id':bid}


def central_difference(plus, minus, delta_degrees):
    """Per-radian derivative from physical deflection +/-delta_degrees."""
    delta=float(delta_degrees)
    if not 0 < delta <= 30: raise ValueError('Physical deflection step outside 0..30 deg.')
    return (float(plus)-float(minus))/(2*math.radians(delta))


def compare_quantities(db, *, test_id, predicted_id, reference_id, rule_id=None,
                       comparison_kind='DIRECT', condition_verified=False, uncertainty=None):
    p=db.row('quantity','quantity_id',predicted_id)
    r=db.row('quantity','quantity_id',reference_id)
    if not p or not r or p['value_status']!='AVAILABLE' or r['value_status']!='AVAILABLE':
        raise ValueError('Both quantities must be available.')
    if p['qoi']!=r['qoi'] or p['canonical_unit']!=r['canonical_unit'] or p['component']!=r['component']:
        raise ValueError('Incompatible QoI, unit, or component.')
    pd=db.row('dataset','dataset_id',p['dataset_id'])
    rd=db.row('dataset','dataset_id',r['dataset_id'])
    if not pd or not rd:raise ValueError('Missing compared dataset')
    pm=json.loads(pd['metadata_json']);rm=json.loads(rd['metadata_json'])
    if pm.get('axes_convention_id') and rm.get('axes_convention_id') and pm['axes_convention_id']!=rm['axes_convention_id']:
        raise ValueError('Reference axes differ: documented conversion required')
    delta=p['canonical_value']-r['canonical_value']
    rel=delta/r['canonical_value'] if abs(r['canonical_value'])>1e-14 else None
    verdict='INCONCLUSIVE'
    notes='No approved numerical/engineering limit or physical condition equivalence.'
    if condition_verified and rule_id:
        rule=db.conn.execute('SELECT * FROM acceptance_rule WHERE rule_id=?',(rule_id,)).fetchone()
        if rule and rule['approval']=='APPROVED' and rule['qoi']==p['qoi']:
            policy=json.loads(rule['policy_json'])
            if policy.get('max_abs_difference') is not None:
                verdict='PASS' if abs(delta)<=float(policy['max_abs_difference']) else 'FAIL'
                notes='Approved absolute-difference threshold; check physical scope separately.'
    uncertainty=uncertainty or {}
    val=validation_uncertainty(numerical=uncertainty.get('numerical'),
        input_=uncertainty.get('input'),experimental=uncertainty.get('experimental'),
        independent=uncertainty.get('independent',False),standard=uncertainty.get('standard',False))
    # A repeated report/export should not multiply identical comparisons.
    for old in db.conn.execute('''SELECT * FROM comparison WHERE test_id=? AND
        predicted_quantity_id=? AND reference_quantity_id=? AND comparison_kind=?
        AND (rule_id IS ? OR rule_id=?)''',
        (test_id,predicted_id,reference_id,comparison_kind,rule_id,rule_id)):
        settings=json.loads(old['uncertainty_metadata_json'])
        if (settings.get('condition_verified')==condition_verified and
            settings.get('validation_uncertainty')==val):
            return {'comparison_id':old['comparison_id'],'delta':old['delta_value'],
                    'relative_delta':old['relative_delta'],'verdict':old['verdict'],
                    'uncertainty':val,'status':'REUSED_COMPARISON'}
    ident=uid('CMP')
    with db.transaction() as c:
        c.execute('INSERT INTO comparison VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(
            ident,test_id,predicted_id,reference_id,comparison_kind,delta,rel,
            val.get('standard_uncertainty'),None,
            canonical({'condition_verified':condition_verified,'validation_uncertainty':val}),rule_id,verdict,notes,utc()))
    return {'comparison_id':ident,'delta':delta,'relative_delta':rel,'verdict':verdict,'uncertainty':val}


def screen_convergence(values, tolerance=None):
    """Diagnostic only; three values ordered coarse/medium/fine."""
    if len(values)!=3 or any(not math.isfinite(float(v)) for v in values):
        raise ValueError('Exactly three finite grid values required.')
    c,m,f=map(float,values)
    dc=m-c;df=f-m
    trend='MONOTONIC' if dc*df>0 else 'OSCILLATORY' if dc*df<0 else 'FLAT_OR_MIXED'
    return {'coarse_medium':dc,'medium_fine':df,'trend':trend,
            'provisional_screen': 'NOT_EVALUATED' if tolerance is None else ('WITHIN_PROVISIONAL_LIMIT' if abs(df)<=tolerance else 'OUTSIDE_PROVISIONAL_LIMIT'),
            'engineering_verdict':'INCONCLUSIVE', 'gci':'NOT_ASSESSED'}


def read_history(path):
    """Read 3.50.5 .history columns by native named header, not offsets.

    Records from multiple Solver Case blocks are supported. The values are
    raw, not normalized or judged converged. Header mismatch is an error.
    """
    lines=Path(path).read_text(encoding='utf-8-sig',errors='replace').splitlines()
    rows=[];header=None;case_number=0
    for line in lines:
        if 'Solver Case:' in line:
            case_number+=1;header=None
            continue
        parts=line.split()
        if case_number==0 or not parts:continue
        if parts[0]=='Iter' and {'CLi','CLiw','CDi','CDiw','CMiy','L2','Residual'} <= set(parts):
            # The last three header tokens in native 3.50.5 are:
            # L2 Residual, Max Residual, Wall_Time (each residual has 2 tokens).
            # Normalize the compound tokens explicitly, no column indices.
            head=' '.join(parts).replace('L2 Residual','L2_Residual').replace('Max Residual','Max_Residual').split()
            if len(head)!=len(set(head)):raise ValueError('Duplicate native history headings')
            header=head
            continue
        if header is None or not parts[0].isdigit():continue
        if len(parts)!=len(header):raise ValueError('Native history row/header column count mismatch')
        try: values=[float(x.replace('D','E')) for x in parts]
        except ValueError as exc:raise ValueError('Nonnumeric native history row') from exc
        if any(not math.isfinite(x) for x in values):raise ValueError('Nonfinite history value')
        row=dict(zip(header,values))
        if int(row['Iter'])!=row['Iter']:raise ValueError('Noninteger iteration')
        row.update({'case':case_number,'iteration':int(row['Iter']),
                    'l2_residual':row['L2_Residual'],
                    'max_residual':row['Max_Residual'],
                    'wall_time':row['Wall_Time']})
        rows.append(row)
    if case_number and not rows:raise ValueError('History has Solver Case but no valid data rows')
    return rows


def richardson_gci(values, grid_widths, safety_factor=1.25):
    """Three-grid asymptotic estimator; refuse invalid or unproven assumptions.

    Widths are physical/representative positive h, ordered coarse -> fine.
    A similar ratio is a prerequisite; this estimate alone does not establish
    asymptotic convergence, solution verification or acceptance.
    """
    if len(values)!=3 or len(grid_widths)!=3:
        raise ValueError('Three QoI values and physical representative grid widths required')
    v=list(map(float,values));h=list(map(float,grid_widths))
    if not all(math.isfinite(x) for x in v+h) or not h[0]>h[1]>h[2]>0:
        return {'status':'NOT_ASSESSED','reason':'Widths must be positive, strictly refined and finite'}
    r21=h[0]/h[1];r32=h[1]/h[2]
    if min(r21,r32)<1.25 or abs(r21-r32)/max(r21,r32)>0.02:
        return {'status':'NOT_ASSESSED','reason':'Three widths do not have a sufficiently uniform refinement ratio'}
    d21=v[0]-v[1];d32=v[1]-v[2]
    if d21*d32<=0 or abs(d32)<1e-14:
        return {'status':'NOT_ASSESSED','reason':'Nonmonotone/flat differences or zero denominator'}
    p=math.log(abs(d21/d32))/math.log(r32)
    if not 0<p<=12 or not math.isfinite(p):
        return {'status':'NOT_ASSESSED','reason':'Nonpositive or implausible observed order'}
    r=r32
    extrapolated=v[2]+(v[2]-v[1])/(r**p-1)
    fine_error=abs(v[2]-v[1])/(r**p-1)
    coarse_error=abs(v[1]-v[0])/(r**p-1)
    # Absolute form, no division by near-zero QoI.
    return {'status':'CANDIDATE_ONLY','order':p,'extrapolated':extrapolated,
            'gci_fine_absolute':safety_factor*fine_error,
            'gci_medium_absolute':safety_factor*coarse_error,
            'refinement_ratio':r,'asymptotic_verified':False,
            'reason':'GCI estimates only; asymptotic range, grid equivalence and solver convergence require independent review'}


def validation_uncertainty(*, numerical, input_, experimental, independent=False, standard=False):
    """Combine standard uncertainties only with explicit independence review."""
    if not independent or not standard:
        return {'status':'NOT_ASSESSED','reason':'Component independence/standard uncertainty not verified'}
    values=[numerical,input_,experimental]
    if any(x is None or not math.isfinite(float(x)) or float(x)<0 for x in values):
        return {'status':'NOT_ASSESSED','reason':'Missing/nonfinite/negative uncertainty component'}
    return {'status':'COMPUTED','standard_uncertainty':math.sqrt(sum(float(x)**2 for x in values)),
            'components':{'numerical':numerical,'input':input_,'experimental':experimental},
            'interpretation':'Standard uncertainty of independent components; not an automatic validation PASS'}


def store_assessment(db, *, study_id, test_id, qoi, evidence, reason,
                     case_id=None, rule_ids=None, numerical='NOT_EVALUATED',
                     engineering='NOT_EVALUATED', disposition='INSUFFICIENT_EVIDENCE'):
    """Persist explicit, unapproved evidence snapshots. No hidden mutable selection."""
    ident=uid('ASM')
    if numerical=='PASS' or engineering=='PASS':
        raise ValueError('Automated evidence inventory cannot approve an engineering result')
    with db.transaction() as c:
        c.execute('INSERT INTO assessment VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',(
            ident,study_id,test_id,case_id,qoi,canonical({}),numerical,
            engineering,disposition,canonical(rule_ids or []),canonical(evidence),
            canonical({}),reason,None,None,utc()))
    return ident


def evaluate_control_sweep(values, qoi):
    """Extract finite-deflection increments; enforce physical deflection input.

    values: dict physical angle (degrees) -> coefficient at same conditions.
    No control convention or derivative signs are assumed to be approved.
    """
    samples={float(k):float(v) for k,v in values.items()}
    result={'qoi':qoi,'samples_deg':samples,'status':'INSUFFICIENT_EVIDENCE'}
    derivatives={}
    for step in (2.,5.,10.):
        if step in samples and -step in samples:
            derivatives[str(int(step))]=central_difference(samples[step],samples[-step],step)
    result['symmetric_per_rad']=derivatives
    if '2' in derivatives and '5' in derivatives:
        result['slope_difference_per_rad']=derivatives['5']-derivatives['2']
    if 0. in samples:
        result['finite_deflection_increments']={str(a):v-samples[0.] for a,v in sorted(samples.items())}
    return result


def read_mesh_stats(stdout_path):
    """Capture explicitly named solver mesh counts, with missing fields omitted."""
    import re
    data=Path(stdout_path).read_text(encoding='utf-8',errors='replace')
    patterns={
        'triangles_total':r'Pure tri mesh has\s+(\d+)\s+tris',
        'kutta_nodes':r'Number of Kutta Nodes:\s*(\d+)',
        'wake_trailing_nodes':r'NumberOfWakeTrailingNodes_:\s*(\d+)',
        'wake_nodes_total':r'NumberOfWakeNodes:\s*(\d+)',
        'wake_loops':r'NumberOfWakeLoops:\s*(\d+)',
    }
    result={}
    for name,pattern in patterns.items():
        found=re.findall(pattern,data,flags=re.IGNORECASE)
        if found:result[name]=int(found[-1])
    return result


def compare_reference_dataset(db, *, test_id, predicted_dataset, reference_dataset,
                              conditions_reviewed=False, rule_id=None, uncertainty=None):
    """Compare compatible native quantities with an explicit independent source.

    No coordinate conversion/interpolation is invented. Every pair receives
    its own comparison record; missing requirements are blockers, not PASS.
    """
    expected={'VER-01':'THEORY','VER-02':'THEORY','XCD-01':'PANEL',
              'XCD-02':'NS','VAL-01':'WIND_TUNNEL','VAL-02':'WIND_TUNNEL'}
    if test_id not in expected:raise ValueError('Reference comparison undefined for '+test_id)
    pds=db.row('dataset','dataset_id',predicted_dataset)
    rds=db.row('dataset','dataset_id',reference_dataset)
    if not pds or not rds or pds['kind']!='VSPAERO' or rds['kind']!=expected[test_id]:
        raise ValueError('Wrong source kind for comparison')
    if not pds['condition_id']:raise ValueError('Missing VSPAERO condition ID')
    cond=db.row('flight_condition','condition_id',pds['condition_id'])
    meta=json.loads(rds['metadata_json'])
    for key in ('mach','re_cref','alpha_deg','beta_deg','sref','bref','cref'):
        if key not in meta or not math.isclose(float(meta[key]),float(cond[key]),rel_tol=1e-7,abs_tol=1e-7):
            raise ValueError('Reference flight/reference condition mismatch: '+key)
    if meta.get('axes_convention_id')!=cond['axes_convention_id']:
        raise ValueError('Reference axes convention mismatch; reviewed transform needed')
    if not conditions_reviewed or not meta.get('geometry_equivalence_reviewed',False):
        raise ValueError('Explicit review of geometry and axis/normalization conventions required')
    if test_id.startswith('VAL') and rds['usage_class']!='VALIDATION':
        raise ValueError('Wind tunnel data reserved for calibration cannot be reused for validation')
    observed={}
    for q in db.conn.execute("SELECT * FROM quantity WHERE dataset_id=? AND value_status='AVAILABLE'",(predicted_dataset,)):
        observed[(q['qoi'],q['component'],q['canonical_unit'])]=q
    results=[]
    for ref in db.conn.execute("SELECT * FROM quantity WHERE dataset_id=? AND value_status='AVAILABLE'",(reference_dataset,)):
        native=observed.get((ref['qoi'],ref['component'],ref['canonical_unit']))
        if native:
            results.append(compare_quantities(db,test_id=test_id,predicted_id=native['quantity_id'],
                reference_id=ref['quantity_id'],rule_id=rule_id,
                condition_verified=True,uncertainty=uncertainty))
    if not results:raise ValueError('No matched QoI, component and canonical units')
    return results


def reprocess_run(db, *, study_root, run_id):
    """Append an immutable native re-read dataset without rerunning VSPAERO.

    Native .polar, .stab, .adjoint.stab and P/Q/R formats are all supported.
    Repeat reprocessing of the same postprocessor fingerprint is idempotent.
    """
    from .VVDatabase import postprocess_fingerprint
    from .VVNative import required_native_files, read_mode_derivatives, check_stability_condition
    root=Path(study_root).resolve();repo=REPO_ROOT
    code=postprocess_fingerprint(repo)
    run=db.row('run','run_id',run_id)
    if not run or run['status']!='SUCCEEDED':raise ValueError('Only succeeded Runs can be reprocessed')
    for a in db.conn.execute('SELECT * FROM artifact WHERE run_id=?',(run_id,)):
        f=root/a['relative_path']
        if not f.is_file() or digest(f)!=a['sha256']:
            raise ValueError('Immutable raw artifact changed: '+a['relative_path'])
    setting=db.row('numerical_setting','setting_id',run['setting_id'])
    mode=setting['analysis_mode']
    method='POLAR_REPROCESS' if mode=='STEADY_SWEEP' else 'STABILITY_REPROCESS'
    prior=db.conn.execute("""SELECT dataset_id,metadata_json FROM dataset WHERE method_id=?
        AND point_id IN (SELECT point_id FROM case_point WHERE run_id=?)""",(method,run_id)).fetchall()
    for item in prior:
        if json.loads(item['metadata_json']).get('postprocess_source_sha256')==code:
            return {'dataset_id':item['dataset_id'],'status':'REUSED_POSTPROCESS'}
    point=db.conn.execute("SELECT * FROM case_point WHERE run_id=? AND point_status='MATCHED'",(run_id,)).fetchone()
    if not point:raise ValueError('Matched source point missing')
    raw=root/run['run_dir']/'raw'
    suffix=required_native_files(mode)[0]
    path=raw/('model_effective'+suffix)
    if not path.is_file():
        # Old individually imported Runs can use an arbitrary basename. Only
        # permit a unique, checksum-registered native source within this Run.
        candidates=[f for f in raw.iterdir() if f.is_file() and f.name.endswith(suffix)]
        if len(candidates)!=1:
            raise ValueError('Expected one immutable native source: '+suffix)
        path=candidates[0]
    cond=db.row('flight_condition','condition_id',run['condition_id'])
    aid=db.conn.execute('SELECT artifact_id FROM artifact WHERE relative_path=?',
        (path.relative_to(root).as_posix(),)).fetchone()
    if not aid:raise ValueError('Unregistered native source artifact')
    if mode=='STEADY_SWEEP':
        rows=read_polar(path)
        if len(rows)!=1 or not condition_match(rows[0],cond):
            raise ValueError('Polar condition inconsistent with source Run')
        metadata={'raw_artifact_id':aid[0],'postprocess_source_sha256':code,'source_mode':mode}
        values=[(key,'DEFAULT',value) for key,value in rows[0].items() if key not in POLAR_CONDITIONS]
    else:
        header,values,groups=read_mode_derivatives(path,mode)
        check_stability_condition(header,cond)
        metadata={'raw_artifact_id':aid[0],'postprocess_source_sha256':code,
                  'source_mode':mode,'native_header':header,'groups':groups,
                  'physical_convention':'NATIVE_NOT_REVIEWED'}
        if mode=='STEADY':
            original=db.conn.execute("""SELECT metadata_json FROM dataset WHERE point_id=?
                AND method_id='STAB_NATIVE' ORDER BY rowid LIMIT 1""",(point['point_id'],)).fetchone()
            if not original:raise ValueError('Missing reviewed original control group mapping')
            mapping=json.loads(original['metadata_json']).get('approved_mapping',{})
            native_columns={variant for _,variant,_ in values if variant.startswith('ConGrp_')}
            if native_columns!=set(mapping) or any(groups.get(c)!=name for c,name in mapping.items()):
                raise ValueError('Original reviewed control-group mapping no longer matches')
            metadata['approved_mapping']=mapping
    with db.transaction():
        ds=insert_dataset(db,study_id=run['study_id'],kind='VSPAERO',
            source='Reprocessed native '+suffix,point_id=point['point_id'],
            condition_id=run['condition_id'],geometry_id=run['geometry_id'],
            metadata=metadata,method=method,version=code[:12])
        for key,variant,value in values:
            status='NOT_AVAILABLE' if value is None else ('AVAILABLE' if mode=='STEADY_SWEEP' else 'NATIVE_ONLY')
            insert_quantity(db,ds,key,value,variant=variant,
                unit='1' if mode=='STEADY_SWEEP' else 'NATIVE_UNVERIFIED',
                native=f'{key}/{variant}',artifact=aid[0],
                locator=f'mode={mode},coef={key},column={variant}',
                method='NATIVE' if mode=='STEADY_SWEEP' else 'UNVERIFIED_STABILITY',status=status)
    return {'dataset_id':ds,'status':'REPROCESSED'}


def collect_study_results(db, study_root, selected_runs=None):
    """Freeze an explicit Case/Run/Dataset snapshot for review, without approval.

    Preserve the native names of all stability derivatives. Native unverified
    values may be plotted under their original names, never silently converted
    to flight-dynamics conventions.
    """
    root=Path(study_root).resolve()
    selected_runs=selected_runs or {}
    unknown=set(selected_runs)-{row['case_id'] for row in db.conn.execute('SELECT case_id FROM planned_case')}
    if unknown:raise ValueError('Unknown selected Case IDs: '+str(sorted(unknown)))
    rows=[]
    for case in db.conn.execute('SELECT * FROM planned_case ORDER BY test_id,role'):
        cid=case['case_id']
        result={'test_id':case['test_id'],'case_id':cid,'role':case['role'],
                'policy':case['execute_policy'],'status':'UNRESOLVED','run_id':None,
                'datasets':[],'quantities':[],'history':[],'mesh':{},'issues':[]}
        if cid in selected_runs:
            chosen=db.conn.execute('''SELECT r.* FROM run r
                JOIN case_point p ON p.run_id=r.run_id
                JOIN case_run_link l ON l.point_id=p.point_id
                WHERE l.case_id=? AND r.run_id=? AND p.point_status='MATCHED'
                    AND r.status='SUCCEEDED' ''',(cid,selected_runs[cid])).fetchone()
            if not chosen:raise ValueError(f'Invalid selected successful Run for {cid}: {selected_runs[cid]}')
        else:
            chosen=db.conn.execute('''SELECT r.* FROM run r
                JOIN case_point p ON p.run_id=r.run_id
                JOIN case_run_link l ON l.point_id=p.point_id
                WHERE l.case_id=? AND p.point_status='MATCHED' AND r.status='SUCCEEDED'
                ORDER BY r.rowid DESC LIMIT 1''',(cid,)).fetchone()
        if chosen:
            result['run_id']=chosen['run_id']
            result['status']='SUCCEEDED'
            result['condition_id']=chosen['condition_id']
            condition=db.row('flight_condition','condition_id',chosen['condition_id'])
            result['flight_condition']={name:condition[name] for name in
                ('mach','re_cref','alpha_deg','beta_deg','sref','bref','cref')}
            result['flight_condition']['controls_deg']=json.loads(condition['controls_deg_json'])
            # Only points explicitly linked to THIS Case are admissible.
            # A multi-point solver Run can serve multiple planned Cases.
            points=db.conn.execute('''SELECT pt.point_id FROM case_point pt
                JOIN case_run_link link ON link.point_id=pt.point_id
                WHERE pt.run_id=? AND link.case_id=? AND pt.point_status=?''',
                (chosen['run_id'],cid,'MATCHED')).fetchall()
            for point in points:
                for ds in db.conn.execute('''SELECT * FROM dataset WHERE point_id=?
                    ORDER BY rowid''',(point['point_id'],)):
                    result['datasets'].append({'dataset_id':ds['dataset_id'],
                        'method_id':ds['method_id'],'source':ds['source_description']})
                    for q in db.conn.execute('''SELECT * FROM quantity WHERE dataset_id=?
                        ORDER BY qoi,variant''',(ds['dataset_id'],)):
                        # Native unverified derivatives retain their own names and unit.
                        value=q['canonical_value'] if q['value_status']=='AVAILABLE' else q['native_value']
                        if value is None:continue
                        result['quantities'].append({'dataset_id':ds['dataset_id'],
                            'quantity_id':q['quantity_id'],'qoi':q['qoi'],
                            'variant':q['variant'],'component':q['component'],
                            'value':value,'unit':q['canonical_unit'] or q['native_unit'],
                            'status':q['value_status']})
            raw=root/chosen['run_dir']/'raw'
            histories=list(raw.glob('*.history')) if raw.is_dir() else []
            if len(histories)==1:
                try:result['history']=read_history(histories[0])
                except (OSError,ValueError) as exc:result['issues'].append(f'history: {exc}')
            log=raw/'vspaero_api_stdout.log'
            if not log.is_file():log=root/chosen['run_dir']/'logs'/'stdout.log'
            if log.is_file():result['mesh']=read_mesh_stats(log)
        else:
            result['status']='NO_SUCCESSFUL_RUN' if case['geometry_id'] else 'UNRESOLVED'
        rows.append(result)
    from .VVDatabase import CASE_ROLES
    # Plot and difference order is physical Test order, not SQL alphabetical order.
    def case_order(item):
        if item['test_id'].startswith('CTRL-'):
            return (item['test_id'],float(item['role'].split('_')[-1]))
        return (item['test_id'],CASE_ROLES[item['test_id']].index(item['role'])
                if item['role'] in CASE_ROLES[item['test_id']]
                else len(CASE_ROLES[item['test_id']]))
    rows.sort(key=case_order)
    # A report may only display comparisons whose predicted Quantity belongs
    # to the explicitly selected successful Run, not historical alternatives.
    selected_quantities={q['quantity_id'] for case in rows for q in case['quantities']}
    selected_datasets={d['dataset_id'] for case in rows for d in case['datasets']}
    comparisons=[]
    excluded_comparison_ids=[]
    for c in db.conn.execute('''SELECT c.*, p.native_value predicted_value,
            r.native_value reference_value,p.native_unit predicted_unit,r.native_unit reference_unit,
            p.qoi qoi,p.dataset_id predicted_dataset,r.dataset_id reference_dataset
            FROM comparison c JOIN quantity p ON p.quantity_id=c.predicted_quantity_id
            JOIN quantity r ON r.quantity_id=c.reference_quantity_id
            ORDER BY c.test_id,c.comparison_id'''):
        if c['predicted_quantity_id'] not in selected_quantities:
            excluded_comparison_ids.append(c['comparison_id'])
            continue
        comparisons.append({'test_id':c['test_id'],'comparison_id':c['comparison_id'],
            'qoi':c['qoi'],'predicted':c['predicted_value'],'reference':c['reference_value'],
            'delta':c['delta_value'],'unit':c['predicted_unit'],
            'predicted_dataset':c['predicted_dataset'],'reference_dataset':c['reference_dataset']})
    # External CD0 corrections may be attached to a source Case in another
    # Test family. Attribute them to DRAG-01/02 by the *stored method*, not by
    # inventing additional solver Runs for a non-solver correction.
    drag_polars=[]
    for ds in db.conn.execute("SELECT * FROM dataset WHERE kind='DERIVED' AND method_id='CD0_CORRECTION' ORDER BY rowid"):
        meta=json.loads(ds['metadata_json'])
        polar_id=meta.get('polar_dataset')
        if not polar_id or polar_id not in selected_datasets:continue
        source=db.row('dataset','dataset_id',polar_id)
        if not source:continue
        lift=db.conn.execute("""SELECT qoi,canonical_value FROM quantity
            WHERE dataset_id=? AND variant='DEFAULT' AND component='TOTAL'
              AND value_status='AVAILABLE' AND qoi IN ('CLi','CLiw','CLtot')
            ORDER BY CASE qoi WHEN 'CLi' THEN 0 WHEN 'CLiw' THEN 1 ELSE 2 END LIMIT 1""",(polar_id,)).fetchone()
        base=db.conn.execute("SELECT canonical_value FROM quantity WHERE dataset_id=? AND qoi='CDiw' AND value_status='AVAILABLE'",(polar_id,)).fetchone()
        for q in db.conn.execute("SELECT * FROM quantity WHERE dataset_id=? AND qoi IN ('CD_A','CD_B')",(ds['dataset_id'],)):
            # Method-B datasets contain a duplicate auxiliary CD_A; show only B
            # when an independent method-A dataset exists for the same point.
            if q['qoi']=='CD_A' and 'cdo_baseline' in meta:
                other=db.conn.execute("SELECT COUNT(*) FROM dataset WHERE method_id='CD0_CORRECTION' AND point_id=?",
                                      (ds['point_id'],)).fetchone()[0]
                if other>1:continue
            drag_polars.append({'test_id':'DRAG-01' if q['qoi']=='CD_A' else 'DRAG-02',
                'dataset_id':ds['dataset_id'],'polar_dataset_id':polar_id,
                'cd0_dataset_id':meta.get('cd0_dataset'),'quantity_id':q['quantity_id'],
                'qoi':q['qoi'],'cd':q['canonical_value'],
                'cd_native_wake':base['canonical_value'] if base else None,
                'cl':lift['canonical_value'] if lift else None,
                'cl_qoi':lift['qoi'] if lift else None,
                'condition_id':source['condition_id']})
    return {'source':'SQLite Run/Point/Dataset/Quantity snapshot',
            'selection_policy':'EXPLICIT_WHERE_PROVIDED_ELSE_LATEST_SUCCESSFUL_FROZEN_IN_SNAPSHOT',
            'cases':rows,'comparisons':comparisons,'drag_polars':drag_polars,
            'excluded_historical_comparison_ids':excluded_comparison_ids}


def evaluate_test(db, study_id, test_id, *, selection=None, label=None):
    """Legacy per-Test numerical screening, separate from CLI and acceptance.

    The integrated execute-all workflow uses collect_study_results instead.
    This command preserves the earlier explicit Run-selection semantics.
    """
    from .VVDatabase import CASE_ROLES
    from .VVAutomation import case_difference
    roles={'GRID-01':('COARSE','MEDIUM','FINE'),
           'GRID-02':('COARSE','MEDIUM','FINE'),
           'WAKE-01':('ITER_5','ITER_10','ITER_20'),
           'WAKE-02':('NODES_8','NODES_12','NODES_16')}
    if label:
        if not selection:raise ValueError('Frozen evaluation requires an explicit Run selection')
        db.lock_evaluation_selection(test_id,selection,label)
    if test_id.startswith('CTRL-'):
        if selection is None:raise ValueError('CTRL evaluation requires an explicit Run selection')
        samples={}
        for role in CASE_ROLES[test_id]:
            case_id=test_id+'-'+role
            if case_id not in selection:raise ValueError('Missing selected Run: '+case_id)
            ds=db.conn.execute('''SELECT ds.dataset_id FROM case_run_link l
                JOIN case_point pt ON pt.point_id=l.point_id
                JOIN run r ON r.run_id=pt.run_id JOIN dataset ds ON ds.point_id=pt.point_id
                WHERE l.case_id=? AND r.run_id=? AND r.status='SUCCEEDED'
                  AND ds.method_id IN ('POLAR_NATIVE','POLAR_REPROCESS')
                ORDER BY ds.rowid DESC LIMIT 1''',(case_id,selection[case_id])).fetchone()
            if ds is None:raise ValueError('No matching polar dataset: '+case_id)
            values=db.conn.execute('''SELECT qoi,canonical_value FROM quantity WHERE dataset_id=?
                AND value_status='AVAILABLE' AND component='TOTAL' AND variant='DEFAULT' ''',
                (ds['dataset_id'],)).fetchall()
            samples[int(role.split('_')[1])]={q['qoi']:q['canonical_value'] for q in values}
        qois={'CTRL-01':('CMiy','CLi'),'CTRL-02':('CMix','CMiz'),
              'CTRL-03':('CMiz','CSi','CMix')}[test_id]
        result={q:evaluate_control_sweep({a:v[q] for a,v in samples.items() if q in v},q)
                for q in qois}
        evidence=[{'case_id':cid,'run_id':rid} for cid,rid in selection.items()]
        aid=store_assessment(db,study_id=study_id,test_id=test_id,qoi='CONTROL_SENSITIVITY',
            evidence=evidence,reason='Numerical screening only; no approved engineering thresholds')
        return {'test_id':test_id,'assessment_id':aid,'results':result}
    if test_id not in roles:
        return {'test_id':test_id,'verdict':'NOT_EVALUATED',
                'reason':'Dedicated numerical screening not available for this Test ID'}
    if selection is not None and not isinstance(selection,dict):
        raise ValueError('Selection must map Case ID to Run ID')
    equivalence=case_difference(db,test_id)
    samples=[];sources={}
    for role in roles[test_id]:
        cid=test_id+'-'+role
        if selection is not None and cid not in selection:
            raise ValueError('Missing explicit selected Run: '+cid)
        selected=db.conn.execute('''SELECT r.run_id,ds.dataset_id FROM case_run_link l
            JOIN case_point pt ON pt.point_id=l.point_id JOIN run r ON r.run_id=pt.run_id
            JOIN dataset ds ON ds.point_id=pt.point_id WHERE l.case_id=? AND r.status='SUCCEEDED'
              AND (? IS NULL OR r.run_id=?) AND pt.point_status='MATCHED'
              AND ds.kind='VSPAERO' AND ds.method_id IN ('POLAR_NATIVE','POLAR_REPROCESS')
            ORDER BY r.rowid DESC,ds.rowid DESC LIMIT 1''',
            (cid,selection.get(cid) if selection is not None else None,
             selection.get(cid) if selection is not None else None)).fetchone()
        if selected is None:
            sources[cid]=None;samples.append({});continue
        values=db.conn.execute('''SELECT qoi,canonical_value FROM quantity
            WHERE dataset_id=? AND variant='DEFAULT' AND component='TOTAL'
              AND value_status='AVAILABLE' AND qoi IN ('CLi','CMiy','CDiw','CLiw','CDi')''',
            (selected['dataset_id'],)).fetchall()
        sources[cid]={'run_id':selected['run_id'],'dataset_id':selected['dataset_id']}
        samples.append({v['qoi']:v['canonical_value'] for v in values})
    limits={'CLi':0.005,'CMiy':0.002,'CDiw':0.0005,'CLiw':None,'CDi':None}
    outcome={}
    for qoi,limit in limits.items():
        if equivalence:
            outcome[qoi]={'verdict':'INCONCLUSIVE','reason':'Case equivalence not verified: '+'; '.join(equivalence)}
        elif all(qoi in s for s in samples):
            outcome[qoi]=screen_convergence([s[qoi] for s in samples],limit)
        else:
            outcome[qoi]={'verdict':'INCONCLUSIVE','reason':'Missing successful Run/QoI'}
    return {'test_id':test_id,'status':'SCREENING_ONLY',
            'equivalence_errors':equivalence,'sources':sources,'results':outcome}
