"""OpenVSP 3.50.5 native result contracts, as written by vspaero.C.

No solver execution, interpretation of flight-dynamics axes, or silent fallback.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

MODE_FILES = {
    'STEADY_SWEEP': ('.polar', '.history', '.lod', '.adb'),
    'STEADY': ('.stab',),
    'ADJOINT': ('.adjoint.stab',),
    'P': ('.pstab',), 'Q': ('.qstab',), 'R': ('.rstab',),
}
COEFFICIENTS = ('CFx','CFy','CFz','CMx','CMy','CMz','CL','CD','CS','CMl','CMm','CMn')


def required_native_files(mode):
    try:
        return MODE_FILES[mode]
    except KeyError:
        raise ValueError('Unsupported 3.50.5 analysis mode: ' + str(mode)) from None


def native_condition_header(path):
    """Read the scalar part of WriteCaseHeader; require actual case identifiers."""
    result = {}
    for line in Path(path).read_text(encoding='utf-8-sig', errors='replace').splitlines():
        if line.lstrip().startswith('Case'):
            break
        match = re.match(r'^\s*(Sref|Bref|Cref|Mach|AoA|Beta|Vinf|Rho)_?\s+([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\b', line)
        if match:
            result[match[1]] = float(match[2].replace('D','E'))
    if not all(k in result and math.isfinite(result[k]) for k in ('Mach','AoA','Beta','Sref','Bref','Cref')):
        raise ValueError('Incomplete/nonfinite native stability case header: ' + str(path))
    return result


def check_stability_condition(header, condition):
    fields = {'Mach':'mach','AoA':'alpha_deg','Beta':'beta_deg',
              'Sref':'sref','Bref':'bref','Cref':'cref'}
    for observed, requested in fields.items():
        expected = float(condition[requested])
        if not math.isclose(header[observed], expected, rel_tol=2e-6, abs_tol=2e-6):
            raise ValueError(f'Native stability condition mismatch: {observed}={header[observed]} != {expected}')


def read_mode_derivatives(path, mode):
    """Return (header, [(coefficient, variant, value)], native control groups).

    Derivatives are native-only. The field names and table variant are retained;
    no axis/phase/sign convention is inferred from a numeric value.
    """
    if mode == 'STEADY_SWEEP':
        raise ValueError('A steady polar is not a stability derivative file')
    header = native_condition_header(path)
    if mode == 'STEADY':
        from ..vspaero.VSPAEROStab import read_vspaero_stab
        stab = read_vspaero_stab(path)
        rows = [(str(coef), str(col), float(value))
                for coef, data in stab.derivatives.iterrows()
                for col, value in data.items()]
        groups = dict(stab.control_groups)
        names = {coef for coef,_,_ in rows}
        if names != set(COEFFICIENTS) or len(rows)!=len(COEFFICIENTS)*len(stab.derivatives.columns):
            raise ValueError('Incomplete native steady stability derivative matrix')
    elif mode == 'ADJOINT':
        lines = Path(path).read_text(encoding='utf-8-sig', errors='replace').splitlines()
        rows = []
        groups = {}  # The adjoint numerical columns do not prove physical CS group names.
        for i, line in enumerate(lines):
            parts = line.split()
            if len(parts) < 9 or parts[0] != 'Coef' or parts[1] not in ('Total','Inviscid','Viscous'):
                continue
            variant = parts[1]
            columns = parts[2:]
            if columns[:7] != ['Alpha','Beta','p','q','r','Mach','U']:
                raise ValueError('Unknown adjoint derivative columns: ' + str(columns))
            table = []
            for following in lines[i+1:]:
                items = following.split()
                if not items or items[0].startswith('#'):
                    continue
                if items[0] not in COEFFICIENTS:
                    if table:
                        break
                    continue
                if len(items) != len(columns)+2:
                    raise ValueError('Malformed adjoint derivative row: ' + following.strip())
                table.extend((items[0], variant+'_'+column, float(value.replace('D','E')))
                             for column, value in zip(['Base']+columns, items[1:]))
            if len(table) != len(COEFFICIENTS)*(len(columns)+1):
                raise ValueError('Incomplete adjoint derivative table: ' + variant)
            rows.extend(table)
        if not rows:
            raise ValueError('No native adjoint derivative table')
    else:
        expected = {'P':'_p', 'Q':'_(q + alpha_dot)', 'R':'_(r - beta_dot)'}[mode]
        lines = Path(path).read_text(encoding='utf-8-sig', errors='replace').splitlines()
        rows = []
        groups = {}
        for line in lines:
            match = re.match(r'^\s*([A-Za-z][A-Za-z0-9_]*(?:_\([^)]*\))?)\s+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)\s*$', line)
            if not match or not match[1].endswith(expected):
                continue
            coeff = match[1][:-len(expected)]
            if coeff in COEFFICIENTS:
                rows.append((match[1], 'NATIVE', float(match[2].replace('D','E'))))
        if len(rows) != len(COEFFICIENTS) or {name[:-len(expected)] for name,_,_ in rows} != set(COEFFICIENTS):
            raise ValueError(f'Expected 12 distinct {mode} derivatives in native file, got {len(rows)}')
    keys = [(name, variant) for name, variant, _ in rows]
    if not rows or len(set(keys)) != len(keys):
        raise ValueError('Missing or duplicate native stability derivatives')
    if any(not math.isfinite(value) for _, _, value in rows):
        raise ValueError('Nonfinite native stability derivative')
    return header, rows, groups


def classify_control_surfaces(surfaces, selected_geometry_names, all_geometry_names):
    """Match native CS names to parent geometry, retaining approved group gains.

    OpenVSP names a CS as '<GeometryName>_SurfN_<label>'. The longest
    matching geometry name avoids a prefix collision (e.g. Wing/WingTip).
    Unknown/misnamed surfaces must fail instead of being silently disabled.
    """
    present, absent = [], []
    for surface in surfaces:
        parents=[name for name in all_geometry_names
                 if str(surface).startswith(name+'_Surf')]
        if len(parents)!=1:
            raise ValueError('Cannot resolve Control Surface parent geometry: '+str(surface))
        (present if parents[0] in selected_geometry_names else absent).append(str(surface))
    return present, absent


def inspect_native_mesh(geometry_path, key_path):
    """Screen native VSPGEOM v3 topology before spawning the solver.

    A mesh with thousands of vertices and only seven loops, as observed in
    GRID-03-NREF_1, cannot be accepted as a valid aerodynamic surface mesh.
    This is a structural preflight, not a GCI/numerical-convergence verdict.
    """
    from pathlib import Path
    geom=Path(geometry_path)
    key=Path(key_path)
    report={'status':'UNVERIFIED','geometry_file':geom.name,'key_file':key.name,
            'nodes':None,'loops':None,'vkey_parts':None,'issues':[]}
    try:
        with geom.open('r',encoding='utf-8',errors='replace') as stream:
            marker=stream.readline().strip()
            version_line=stream.readline().strip()
            dimensions=stream.readline().split()
        if not marker.lower().startswith('# vspgeom v3') or len(dimensions)<2:
            report['issues'].append('Unrecognized VSPGEOM v3 header')
        else:
            nodes,loops=map(int,dimensions[:2])
            report['nodes'],report['loops']=nodes,loops
            if nodes<=0 or loops<=0:
                report['issues'].append('Nonpositive mesh node/loop count')
            elif nodes>=100 and loops<max(20,nodes//50):
                report['issues'].append(
                    f'Implausibly few surface loops: {loops} loops for {nodes} nodes')
    except (OSError,ValueError,IndexError) as exc:
        report['issues'].append(f'Cannot parse VSPGEOM header: {exc}')
    try:
        if key.is_file():
            lines=key.read_text(encoding='utf-8',errors='replace').splitlines()
            if len(lines)>2:
                report['vkey_parts']=int(lines[2].strip())
                if report['vkey_parts']<=0:
                    report['issues'].append('VKEY contains no surface parts')
        elif report['nodes'] is not None:
            report['issues'].append('VKEY file missing')
    except (OSError,ValueError) as exc:
        report['issues'].append(f'Cannot read VKEY part count: {exc}')
    report['status']='INVALID' if report['issues'] else ('PASS' if report['nodes'] is not None else 'UNVERIFIED')
    return report


def diagnose_solver_log(path):
    """Classify observed native log stages; do not fabricate an exit code."""
    from pathlib import Path
    log=Path(path)
    text=log.read_text(encoding='utf-8',errors='replace') if log.is_file() else ''
    stages=(('ADJOINT_LINEAR_SYSTEM','Solving adjoint linear systems'),
            ('WAKE_SOLVE','Wake Iter:'),
            ('WAKE_GRID_INITIALIZED','Updating wake grids'),
            ('MESH_COMPONENT_MAPPING','Renumbering surfaces'),
            ('MESH_READ','NumNodes:'))
    stage=next((name for name,marker in stages if marker.lower() in text.lower()),'NOT_OBSERVED')
    if 'not associated with any component' in text:
        classification='INVALID_MESH_COMPONENT'
    elif 'Could not find control surface' in text:
        classification='INVALID_CONTROL_SURFACE'
    elif stage=='WAKE_GRID_INITIALIZED':
        classification='SOLVER_STOPPED_AFTER_WAKE_GRID_UNEXPLAINED'
    elif stage=='ADJOINT_LINEAR_SYSTEM':
        classification='ADJOINT_NOT_COMPLETED'
    else:
        classification='NATIVE_OUTPUT_MISSING_UNEXPLAINED'
    return {'classification':classification,'last_stage':stage,
            'native_solver_exit_code':'NOT_EXPOSED_BY_PYTHON_API',
            'log_file':str(log),'log_bytes':log.stat().st_size if log.is_file() else 0}


def classify_worker_error(stderr_text, native_diagnosis):
    """Prioritize observed Python preflight exceptions over absent solver output."""
    if 'VSPAERO mesh preflight failed:' in stderr_text:
        return 'MESH_PREFLIGHT_FAILED'
    if 'Partial Control Group geometry selection:' in stderr_text:
        return 'PARTIAL_CONTROL_GROUP'
    if 'Nonzero deflection for omitted Control Group:' in stderr_text:
        return 'OMITTED_CONTROL_GROUP_WITH_DEFLECTION'
    if 'Cannot resolve Control Surface parent geometry:' in stderr_text:
        return 'UNKNOWN_CONTROL_SURFACE_PARENT'
    return native_diagnosis['classification']
