# https://openvsp.org/pyapi_docs/latest/openvsp.html
# https://github.com/OpenVSP/OpenVSP/blob/fc4a97b75c92399fecf74b21f2a57087f89e12b7/examples/scripts/TestAnalysisVSPAERO.vspscript

import os
import sys
from contextlib import contextmanager

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import openvsp as vsp

from .ISAspecification import *

@contextmanager
def suppress_stdout():

    # 標準出力を一時的に/dev/nullにリダイレクト
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

def vsp_sweep(vsp, alpha, mach, reynolds=[1e6], verbose=1):

    if verbose:
        print('\n-> Calculate alpha & mach sweep analysis\n')

    # //==== Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File ====//
    
    # Set defaults
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)

    # List inputs, type, and current values
    if verbose:
        print(compgeom_name)
        vsp.PrintAnalysisInputs(compgeom_name)
        print('')

    # Execute
    if verbose:
        print('\tExecuting...')
    compgeom_resid = vsp.ExecAnalysis(compgeom_name)
    if verbose:
        print('\tCOMPLETE')

    # Get & Display Results
    if verbose:
        vsp.PrintResults(compgeom_resid)
        print('')

    # //==== Analysis: VSPAero Sweep ====//

    # Set defaults
    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))
    if 'UnsteadyType' in analysis_inputs:
        vsp.SetIntAnalysisInput(analysis_name, 'UnsteadyType', [vsp.STABILITY_OFF], 0)

    # Reference geometry set
    geom_set = [0]
    vsp.SetIntAnalysisInput(analysis_name, 'GeomSet', geom_set, 0)
    ref_flag = [1]
    vsp.SetIntAnalysisInput(analysis_name, 'RefFlag', ref_flag, 0)
    wid = vsp.FindGeomsWithName('WingGeom')
    vsp.SetStringAnalysisInput(analysis_name, 'WingID', wid, 0)

    # Freestream Parameters
    alpha_npts = [len(alpha)]
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [alpha[0]], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaEnd", [alpha[-1]], 0)
    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", alpha_npts, 0)
    mach_npts = [len(mach)]
    vsp.SetDoubleAnalysisInput(analysis_name, 'MachStart', [mach[0]], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, 'MachEnd', [mach[-1]], 0)
    vsp.SetIntAnalysisInput(analysis_name, 'MachNpts', mach_npts, 0)
    reynolds_npts = [len(reynolds)]
    vsp.SetDoubleAnalysisInput(analysis_name, 'ReCref', [reynolds[0]], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, 'ReCrefEnd', [reynolds[-1]], 0)
    vsp.SetIntAnalysisInput(analysis_name, 'ReCrefNpts', reynolds_npts, 0)
    vsp.Update()

    # List inputs, type, and current values
    if verbose:
        print(analysis_name)
        vsp.PrintAnalysisInputs(analysis_name)
        print('')

    # Execute
    if verbose:
        print('\tExecuting...')
        rid = vsp.ExecAnalysis(analysis_name)
        print('\tCOMPLETE')
    else:
        with suppress_stdout():
            rid = vsp.ExecAnalysis('VSPAEROSweep')

    # Get & Display Results
    # vsp.PrintResults(rid)

    return rid

def MyFindParm(cs_group_container_id, parm_name):

    # Find param_id from cs_group_container_id and parm_name
    parm_ids = vsp.FindContainerParmIDs( cs_group_container_id )
    for parm_id in parm_ids:
        if parm_name == vsp.GetParmName(parm_id):
            return parm_id
    return ''



# G103A stability-derivative preflight specification.
# This validator is intentionally not generic: it checks the naming and set
# conventions used by the G103A OpenVSP model before running stability analyses.
G103A_EXPECTED_GEOMS = {
    'FuselageGeom': 'FUSELAGE',
    'WingGeom': 'WING',
    'HTailGeom': 'WING',
    'VTailGeom': 'WING',
}

G103A_EXPECTED_SUBSURFACES = {
    'WingGeom': 'AILERON',
    'HTailGeom': 'ELEVATOR',
    'VTailGeom': 'RUDDER',
}

G103A_EXPECTED_SETS = {
    'ThickGeom': ['FuselageGeom'],
    'ThinGeom': ['WingGeom', 'HTailGeom', 'VTailGeom'],
}

G103A_EXPECTED_CONTROL_GROUPS = {
    'AILERON_GROUP': {
        'geom_name': 'WingGeom',
        'subsurface_name': 'AILERON',
        'expected_gains': [1.0, 1.0],
    },
    'ELEVATOR_GROUP': {
        'geom_name': 'HTailGeom',
        'subsurface_name': 'ELEVATOR',
        'expected_gains': [1.0, -1.0],
    },
    'RUDDER_GROUP': {
        'geom_name': 'VTailGeom',
        'subsurface_name': 'RUDDER',
        'expected_gains': [1.0],
    },
}

G103A_REF_GEOM_NAME = 'WingGeom'


def validate_vsp3_for_stability_derivatives(vsp3_path, *, verbose=1):
    """
    Validate whether a G103A-style .vsp3 model is structurally ready for
    VSPAERO stability-derivative calculation.

    This is a static preflight check. It reads the .vsp3 file and validates
    geometry names, subsurfaces, thick/thin sets, VSPAERO symmetry,
    control-surface groups, gains, reference values, and moment-reference
    coordinates.

    It does not run VSPAEROComputeGeometry, VSPAEROSweep, STABILITY_DEFAULT,
    or parse .stab output.

    Parameters
    ----------
    vsp3_path : str or os.PathLike
        Path to the .vsp3 file to read.
    verbose : int or bool, optional
        0: no stdout, 1: major progress, 2: detailed validation summary.

    Returns
    -------
    dict
        Validation report with passed/errors/warnings/infos and summaries.
    """

    report = {
        'passed': False,
        'errors': [],
        'warnings': [],
        'infos': [],
        'geom_summary': {},
        'subsurface_summary': {},
        'symmetry_summary': {},
        'set_summary': {},
        'control_group_summary': {},
        'vspaero_settings_summary': {},
    }

    def add(kind, code, message, context=None):
        report[kind].append({
            'code': code,
            'message': message,
            'context': context or {},
        })

    def vprint(level, message):
        if verbose and int(verbose) >= level:
            print(message)

    def normalize_type_name(type_name):
        return str(type_name).strip().upper()

    def is_finite_number(value):
        try:
            return np.isfinite(float(value))
        except (TypeError, ValueError):
            return False

    def get_container_parm_value(container_id, parm_name):
        parm_id = MyFindParm(container_id, parm_name)
        if not parm_id:
            return None, ''
        return vsp.GetParmVal(parm_id), parm_id

    vsp3_path = os.fspath(vsp3_path)
    if not os.path.isfile(vsp3_path):
        add('errors', 'FILE_NOT_FOUND', 'The specified .vsp3 file was not found.', {'vsp3_path': vsp3_path})
        return report

    vprint(1, f'\n-> Validate G103A stability-derivative preflight: {vsp3_path}')

    # Read the model into a clean OpenVSP session. This function intentionally
    # validates the file as a standalone model, not as an insert into the
    # currently loaded vehicle.
    try:
        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(vsp3_path)
        vsp.Update()
    except Exception as err:
        add('errors', 'READ_VSP3_FAILED', 'OpenVSP failed to read the .vsp3 file.', {'error': repr(err)})
        return report

    add('infos', 'READ_VSP3_OK', 'The .vsp3 file was read by OpenVSP.', {'vsp3_path': vsp3_path})

    # 1. Required Geoms.
    geom_ids = {}
    vprint(1, '  Checking required Geoms...')
    for geom_name, expected_type in G103A_EXPECTED_GEOMS.items():
        ids = list(vsp.FindGeomsWithName(geom_name))
        summary = {
            'found': bool(ids),
            'geom_ids': ids,
            'geom_id': ids[0] if len(ids) == 1 else None,
            'type_name': None,
            'expected_type': expected_type,
            'passed': False,
        }

        if len(ids) == 0:
            add('errors', 'MISSING_GEOM', f"Required Geom '{geom_name}' was not found.", {'geom_name': geom_name})
        elif len(ids) > 1:
            add('errors', 'DUPLICATE_GEOM', f"Required Geom '{geom_name}' is ambiguous because multiple Geoms have the same name.", {'geom_name': geom_name, 'geom_ids': ids})
        else:
            geom_id = ids[0]
            geom_ids[geom_name] = geom_id
            type_name = normalize_type_name(vsp.GetGeomTypeName(geom_id))
            summary['type_name'] = type_name
            summary['passed'] = type_name == expected_type
            if type_name != expected_type:
                add('errors', 'GEOM_TYPE_MISMATCH', f"Geom '{geom_name}' has type '{type_name}', expected '{expected_type}'.", {'geom_name': geom_name, 'geom_id': geom_id})

        report['geom_summary'][geom_name] = summary

    # 2. Required control-surface Subsurfaces.
    subsurface_ids = {}
    vprint(1, '  Checking required control-surface Subsurfaces...')
    expected_control_type = int(getattr(vsp, 'SS_CONTROL', 3))
    for geom_name, subsurface_name in G103A_EXPECTED_SUBSURFACES.items():
        summary = {
            'found': False,
            'geom_name': geom_name,
            'subsurface_name': subsurface_name,
            'subsurface_id': None,
            'subsurface_type': None,
            'expected_type': 'SS_CONTROL',
            'passed': False,
        }

        geom_id = geom_ids.get(geom_name)
        if not geom_id:
            add('errors', 'SUBSURFACE_PARENT_MISSING', f"Cannot check Subsurface '{subsurface_name}' because parent Geom '{geom_name}' is missing.", {'geom_name': geom_name})
            report['subsurface_summary'][geom_name] = summary
            continue

        sub_ids = list(vsp.GetSubSurfIDVec(geom_id))
        matches = [sid for sid in sub_ids if vsp.GetSubSurfName(sid) == subsurface_name]
        if len(matches) == 0:
            add('errors', 'MISSING_SUBSURFACE', f"Required Subsurface '{subsurface_name}' was not found on '{geom_name}'.", {'geom_name': geom_name, 'available_subsurfaces': [vsp.GetSubSurfName(sid) for sid in sub_ids]})
        elif len(matches) > 1:
            add('errors', 'DUPLICATE_SUBSURFACE', f"Subsurface '{subsurface_name}' appears multiple times on '{geom_name}'.", {'geom_name': geom_name, 'subsurface_ids': matches})
        else:
            subsurface_id = matches[0]
            subsurface_type = int(vsp.GetSubSurfType(subsurface_id))
            subsurface_ids[(geom_name, subsurface_name)] = subsurface_id
            summary.update({
                'found': True,
                'subsurface_id': subsurface_id,
                'subsurface_type': subsurface_type,
                'passed': subsurface_type == expected_control_type,
            })
            if subsurface_type != expected_control_type:
                add('errors', 'SUBSURFACE_TYPE_MISMATCH', f"Subsurface '{subsurface_name}' on '{geom_name}' is not SS_CONTROL.", {'geom_name': geom_name, 'subsurface_id': subsurface_id, 'subsurface_type': subsurface_type})

        report['subsurface_summary'][geom_name] = summary

    # 4. ThickGeom / ThinGeom sets.
    vprint(1, '  Checking ThickGeom / ThinGeom sets...')
    for set_name, expected_geom_names in G103A_EXPECTED_SETS.items():
        expected_names = set(expected_geom_names)
        summary = {
            'found': False,
            'set_name': set_name,
            'set_index': None,
            'expected_geom_names': expected_geom_names,
            'actual_geom_names': [],
            'missing_geom_names': [],
            'unexpected_geom_names': [],
            'passed': False,
        }

        try:
            set_index = int(vsp.GetSetIndex(set_name))
        except Exception:
            set_index = -1

        if set_index < 0:
            add('errors', 'MISSING_SET', f"Required set '{set_name}' was not found.", {'set_name': set_name})
            report['set_summary'][set_name] = summary
            continue

        try:
            set_geom_ids = list(vsp.GetGeomSet(set_name))
        except Exception:
            set_geom_ids = list(vsp.GetGeomSetAtIndex(set_index))

        actual_names = [vsp.GetGeomName(gid) for gid in set_geom_ids]
        actual_name_set = set(actual_names)
        missing = sorted(expected_names - actual_name_set)
        unexpected = sorted(actual_name_set - expected_names)

        summary.update({
            'found': True,
            'set_index': set_index,
            'actual_geom_names': actual_names,
            'missing_geom_names': missing,
            'unexpected_geom_names': unexpected,
            'passed': not missing and not unexpected,
        })

        if missing:
            add('errors', 'SET_MISSING_GEOM', f"Set '{set_name}' is missing required Geoms.", {'set_name': set_name, 'missing_geom_names': missing})
        if unexpected:
            add('errors', 'SET_HAS_UNEXPECTED_GEOM', f"Set '{set_name}' contains unexpected Geoms.", {'set_name': set_name, 'unexpected_geom_names': unexpected})

        report['set_summary'][set_name] = summary

    # 5. VSPAERO settings container, reference values, and set indices.
    vprint(1, '  Checking VSPAERO settings...')
    vspaero_settings_id = vsp.FindContainer('VSPAEROSettings', 0)
    settings_summary = {
        'container_id': vspaero_settings_id,
        'Sref': None,
        'bref': None,
        'cref': None,
        'Xcg': None,
        'Ycg': None,
        'Zcg': None,
        'GeomSet': None,
        'ThinGeomSet': None,
        'Symmetry': None,
        'expected_GeomSet': report['set_summary'].get('ThickGeom', {}).get('set_index'),
        'expected_ThinGeomSet': report['set_summary'].get('ThinGeom', {}).get('set_index'),
        'passed': False,
    }

    symmetry_summary = {
        'source': 'VSPAEROSettings',
        'container_id': vspaero_settings_id,
        'parm_name': 'Symmetry',
        'parm_id': '',
        'value': None,
        'expected_value': 0.0,
        'passed': False,
    }

    if not vspaero_settings_id:
        add('errors', 'MISSING_VSPAERO_SETTINGS', "The 'VSPAEROSettings' container was not found.")
    else:
        for parm_name in ['Sref', 'bref', 'cref', 'Xcg', 'Ycg', 'Zcg', 'GeomSet', 'ThinGeomSet']:
            value, parm_id = get_container_parm_value(vspaero_settings_id, parm_name)
            settings_summary[parm_name] = value
            if value is None:
                add('errors', 'MISSING_VSPAERO_PARM', f"VSPAERO setting '{parm_name}' was not found.", {'parm_name': parm_name})

        symmetry_value, symmetry_parm_id = get_container_parm_value(vspaero_settings_id, 'Symmetry')
        symmetry_summary['parm_id'] = symmetry_parm_id
        symmetry_summary['value'] = symmetry_value
        settings_summary['Symmetry'] = symmetry_value

        if symmetry_parm_id == '':
            add('errors', 'MISSING_VSPAERO_SYMMETRY', "VSPAERO setting 'Symmetry' was not found.", {'parm_name': 'Symmetry'})
        elif int(round(symmetry_value)) != 0:
            add('errors', 'VSPAERO_XZ_SYMMETRY_ENABLED', 'VSPAERO Settings Symmetry must be 0 for stability-derivative calculation.', symmetry_summary)
        else:
            symmetry_summary['passed'] = True

        for parm_name in ['Sref', 'bref', 'cref']:
            value = settings_summary[parm_name]
            if value is not None and (not is_finite_number(value) or float(value) <= 0):
                add('errors', 'INVALID_REFERENCE_VALUE', f"VSPAERO reference value '{parm_name}' must be positive and finite.", {'parm_name': parm_name, 'value': value})

        for parm_name in ['Xcg', 'Ycg', 'Zcg']:
            value = settings_summary[parm_name]
            if value is not None and not is_finite_number(value):
                add('errors', 'INVALID_MOMENT_REFERENCE', f"VSPAERO moment reference '{parm_name}' must be finite.", {'parm_name': parm_name, 'value': value})

        if settings_summary['GeomSet'] is not None and settings_summary['expected_GeomSet'] is not None:
            if int(round(settings_summary['GeomSet'])) != int(settings_summary['expected_GeomSet']):
                add('errors', 'VSPAERO_GEOMSET_MISMATCH', "VSPAERO GeomSet is not the 'ThickGeom' set.", settings_summary)

        if settings_summary['ThinGeomSet'] is not None and settings_summary['expected_ThinGeomSet'] is not None:
            if int(round(settings_summary['ThinGeomSet'])) != int(settings_summary['expected_ThinGeomSet']):
                add('errors', 'VSPAERO_THINGEOMSET_MISMATCH', "VSPAERO ThinGeomSet is not the 'ThinGeom' set.", settings_summary)

    report['symmetry_summary'] = symmetry_summary
    report['vspaero_settings_summary'] = settings_summary

    # 6. VSPAERO control-surface groups and gains.
    vprint(1, '  Checking VSPAERO control-surface groups and gains...')
    try:
        control_group_names = [vsp.GetVSPAEROControlGroupName(i) for i in range(vsp.GetNumControlSurfaceGroups())]
    except Exception as err:
        control_group_names = []
        add('errors', 'CONTROL_GROUP_LIST_FAILED', 'Failed to read VSPAERO control-surface groups.', {'error': repr(err)})

    for group_name, expected in G103A_EXPECTED_CONTROL_GROUPS.items():
        geom_name = expected['geom_name']
        subsurface_name = expected['subsurface_name']
        expected_gains = expected['expected_gains']
        summary = {
            'found': group_name in control_group_names,
            'group_name': group_name,
            'group_index': None,
            'geom_name': geom_name,
            'subsurface_name': subsurface_name,
            'active_control_surfaces': [],
            'expected_gains': expected_gains,
            'actual_gains': [],
            'passed': False,
        }

        if group_name not in control_group_names:
            add('errors', 'MISSING_CONTROL_GROUP', f"Required control-surface group '{group_name}' was not found.", {'available_groups': control_group_names})
            report['control_group_summary'][group_name] = summary
            continue

        group_index = control_group_names.index(group_name)
        summary['group_index'] = group_index
        active_names = list(vsp.GetActiveCSNameVec(group_index))
        summary['active_control_surfaces'] = active_names
        if not active_names:
            add('errors', 'CONTROL_GROUP_INACTIVE', f"Control-surface group '{group_name}' has no active control surfaces.", {'group_name': group_name})
        elif not any(geom_name in name and subsurface_name in name for name in active_names):
            add('warnings', 'CONTROL_GROUP_ACTIVE_NAME_UNCLEAR', f"Control-surface group '{group_name}' is active, but its active names do not clearly include both the expected Geom and Subsurface names.", {'group_name': group_name, 'active_control_surfaces': active_names})

        subsurface_id = subsurface_ids.get((geom_name, subsurface_name))
        if vspaero_settings_id and subsurface_id:
            for gain_index, expected_gain in enumerate(expected_gains):
                parm_name = f'Surf_{subsurface_id}_{gain_index}_Gain'
                value, parm_id = get_container_parm_value(vspaero_settings_id, parm_name)
                summary['actual_gains'].append(value)
                if value is None:
                    add('errors', 'MISSING_CONTROL_GAIN', f"Gain parm '{parm_name}' was not found for group '{group_name}'.", {'group_name': group_name, 'parm_name': parm_name})
                elif abs(float(value) - float(expected_gain)) > 1e-6:
                    add('errors', 'CONTROL_GAIN_MISMATCH', f"Control-surface group '{group_name}' has unexpected gain.", {'group_name': group_name, 'parm_name': parm_name, 'expected_gain': expected_gain, 'actual_gain': value})

        summary['passed'] = (
            summary['found']
            and bool(summary['active_control_surfaces'])
            and len(summary['actual_gains']) == len(expected_gains)
            and all(value is not None and abs(float(value) - float(expected)) <= 1e-6 for value, expected in zip(summary['actual_gains'], expected_gains))
        )
        report['control_group_summary'][group_name] = summary

    settings_summary['passed'] = not any(item['code'].startswith('VSPAERO_') or item['code'].startswith('INVALID_') or item['code'].startswith('MISSING_VSPAERO_') for item in report['errors'])
    report['passed'] = len(report['errors']) == 0

    if verbose:
        status = 'PASSED' if report['passed'] else 'FAILED'
        print(f'  Validation {status}: {len(report["errors"])} error(s), {len(report["warnings"])} warning(s)')
        if int(verbose) >= 2:
            for error in report['errors']:
                print(f"    ERROR   {error['code']}: {error['message']}")
            for warning in report['warnings']:
                print(f"    WARNING {warning['code']}: {warning['message']}")

    return report

def set_control_surface(geom_name, deflection, cs_group_name, sub_id=0, gains=(1,1), verbose=1):

    # Set control surface group
    current_group_names = [vsp.GetVSPAEROControlGroupName(group_index) for group_index in range(vsp.GetNumControlSurfaceGroups())]
    if cs_group_name in current_group_names:
        group_index = current_group_names.index(cs_group_name)
    else:
        group_index = vsp.CreateVSPAEROControlSurfaceGroup()
        vsp.SetVSPAEROControlGroupName(cs_group_name, group_index)

    # Add subsurface to control surface group
    cs_name_vec = vsp.GetAvailableCSNameVec(group_index)
    selcted = [1 + i for i, cs_name in enumerate(cs_name_vec) if geom_name in cs_name]
    vsp.AddSelectedToCSGroup(selcted, group_index)
    
    # Set gain (deferential)
    if verbose:
        print('\n', cs_group_name, ':', vsp.GetActiveCSNameVec(group_index))
        print('\t', f"{'container_id':14s}", f"{'parm_name':24s}", f"{'group_name':22s}", f"{'parm_id':14s}", f"{'value':8s}")
    geom_id = vsp.FindGeomsWithName(geom_name)[0]
    group_name = 'ControlSurfaceGroup_' + str(group_index)
    cs_group_container_id = vsp.FindContainer('VSPAEROSettings', 0)    
    for i, gain in enumerate(gains):
        parm_name = 'Surf_' + vsp.GetSubSurf(geom_id, sub_id) + '_' + str(i) + '_Gain'
        parm_id = MyFindParm(cs_group_container_id, parm_name)
        # parm_id = vsp.FindParm(cs_group_container_id, parm_name, group_name) # なぜか反応しない...
        if parm_id != '':
            vsp.SetParmVal(parm_id, gain)
            if verbose:
                print('\t', f'{cs_group_container_id:14s}', f'{parm_name:24s}', f'{group_name:22s}', f'{parm_id:14s}', f'{gain:8.3f}')
        else:
            if verbose:
                print('\t', 'Failed to find ' + parm_name)

    # Set defrection angle
    parm_name = 'DeflectionAngle'
    parm_id = vsp.FindParm(cs_group_container_id, parm_name, group_name)
    vsp.SetParmVal(parm_id, deflection)
    vsp.Update()

    # Check if deflection angle successfully applied
    if deflection == vsp.GetParmVal(parm_id):
        if verbose:
            print('\t', f'{cs_group_container_id:14s}', f'{parm_name:24s}', f'{group_name:22s}', f'{parm_id:14s}', f'{deflection:8.3f}')
    else:
        if verbose:
            print('Failed to set deflection angle')
        exit()

    return parm_id

def get_polar_result(result_ids):

    # Loop through the results associated with the given result_ids
    for result_id in vsp.GetStringResults(result_ids, 'ResultsVec'):
        # Check if the result name is 'VSPAERO_Polar'
        if vsp.GetResultsName(result_id) == 'VSPAERO_Polar':
            data, columns = [], []
            # Loop through all available data names for this result_id
            for data_name in vsp.GetAllDataNames(result_id):
                # Get double results for the data_name and append if available
                double_results = vsp.GetDoubleResults(result_id, data_name, 0)
                if double_results:
                    data.append(double_results)
                    columns.append(data_name)
            # If data was found, return it as a pandas DataFrame
            if data:
                return pd.DataFrame(np.array(data).T, columns=columns)
    
    # Return an empty DataFrame if no data is found
    return pd.DataFrame()

def vsp_stability_derivatives(vsp3_path, *, alpha=2.0, mach=0.1, reynolds=4.4e6, verbose=1):
    """
    Run VSPAERO steady 6DOF stability-derivative analysis for a validated
    G103A-style .vsp3 model.

    This function assumes validate_vsp3_for_stability_derivatives() has already
    passed. It does not repair geometry, control-surface groups, gains, or sets.
    It reads the .vsp3 file, runs VSPAEROComputeGeometry with the fixed
    ThickGeom/ThinGeom sets, runs VSPAEROSweep with STABILITY_DEFAULT at a
    single alpha/Mach/Reynolds point, and returns the VSPAERO_Stab result as a
    pandas DataFrame together with useful result metadata.

    The moment-reference position is always taken from the .vsp3 file. This
    function intentionally has no cg argument.
    """

    result = {
        'passed': False,
        'errors': [],
        'warnings': [],
        'infos': [],
        'wrapper_result_id': '',
        'stab_result_id': '',
        'result_names': [],
        'data_names': [],
        'derivatives': pd.DataFrame(),
    }

    def add(kind, code, message, context=None):
        result[kind].append({
            'code': code,
            'message': message,
            'context': context or {},
        })

    def vprint(level, message):
        if verbose and int(verbose) >= level:
            print(message)

    def set_input_if_available(analysis_name, input_names, setter, input_name, values):
        if input_name in input_names:
            setter(analysis_name, input_name, values, 0)
            return True
        return False

    vsp3_path = os.fspath(vsp3_path)
    if not os.path.isfile(vsp3_path):
        add('errors', 'FILE_NOT_FOUND', 'The specified .vsp3 file was not found.', {'vsp3_path': vsp3_path})
        return result

    vprint(1, f'\n-> Calculate G103A VSPAERO stability derivatives: {vsp3_path}')

    try:
        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(vsp3_path)
        vsp.Update()
    except Exception as err:
        add('errors', 'READ_VSP3_FAILED', 'OpenVSP failed to read the .vsp3 file.', {'error': repr(err)})
        return result

    thick_set = int(vsp.GetSetIndex('ThickGeom'))
    thin_set = int(vsp.GetSetIndex('ThinGeom'))
    if thick_set < 0:
        add('errors', 'MISSING_SET', "Required set 'ThickGeom' was not found.")
    if thin_set < 0:
        add('errors', 'MISSING_SET', "Required set 'ThinGeom' was not found.")
    if result['errors']:
        return result

    wing_ids = list(vsp.FindGeomsWithName(G103A_REF_GEOM_NAME))
    if len(wing_ids) != 1:
        add('errors', 'REF_GEOM_NOT_FOUND', f"Reference Geom '{G103A_REF_GEOM_NAME}' must exist exactly once.", {'geom_ids': wing_ids})
        return result
    wing_id = wing_ids[0]

    # VSPAERO geometry preparation.  The set names are fixed by the G103A
    # convention, but the indexes are read from the model so the code remains
    # understandable and does not depend on magic numbers such as 3 and 4.
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)
    compgeom_inputs = set(vsp.GetAnalysisInputNames(compgeom_name))
    set_input_if_available(compgeom_name, compgeom_inputs, vsp.SetIntAnalysisInput, 'GeomSet', [thick_set])
    set_input_if_available(compgeom_name, compgeom_inputs, vsp.SetIntAnalysisInput, 'ThinGeomSet', [thin_set])

    if verbose and int(verbose) >= 3:
        print(compgeom_name)
        vsp.PrintAnalysisInputs(compgeom_name)
        print('')

    vprint(1, '  Executing VSPAEROComputeGeometry...')
    if verbose:
        compgeom_result_id = vsp.ExecAnalysis(compgeom_name)
    else:
        with suppress_stdout():
            compgeom_result_id = vsp.ExecAnalysis(compgeom_name)
    result['infos'].append({
        'code': 'COMPUTE_GEOMETRY_COMPLETE',
        'message': 'VSPAEROComputeGeometry completed.',
        'context': {'result_id': compgeom_result_id, 'thick_set': thick_set, 'thin_set': thin_set},
    })

    # Steady 6DOF stability analysis at one flight condition.
    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))

    set_input_if_available(analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'GeomSet', [thick_set])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'ThinGeomSet', [thin_set])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'RefFlag', [1])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetStringAnalysisInput, 'WingID', [wing_id])

    set_input_if_available(analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'AlphaStart', [alpha])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'AlphaEnd', [alpha])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'AlphaNpts', [1])

    set_input_if_available(analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'MachStart', [mach])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'MachEnd', [mach])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'MachNpts', [1])

    set_input_if_available(analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'ReCref', [reynolds])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'ReCrefEnd', [reynolds])
    set_input_if_available(analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'ReCrefNpts', [1])

    if 'UnsteadyType' in analysis_inputs:
        vsp.SetIntAnalysisInput(analysis_name, 'UnsteadyType', [vsp.STABILITY_DEFAULT], 0)
    else:
        add('errors', 'MISSING_UNSTEADYTYPE_INPUT', "VSPAEROSweep does not expose the 'UnsteadyType' input in this OpenVSP version.")
        return result

    vsp.Update()

    if verbose and int(verbose) >= 3:
        print(analysis_name)
        vsp.PrintAnalysisInputs(analysis_name)
        print('')

    vprint(1, '  Executing VSPAEROSweep with STABILITY_DEFAULT...')
    if verbose:
        wrapper_result_id = vsp.ExecAnalysis(analysis_name)
    else:
        with suppress_stdout():
            wrapper_result_id = vsp.ExecAnalysis(analysis_name)
    result['wrapper_result_id'] = wrapper_result_id

    # VSPAEROSweep normally returns a wrapper result with ResultsVec children.
    # Search those children first, then fall back to FindLatestResultsID so this
    # remains usable across minor OpenVSP result-wrapper differences.
    child_result_ids = []
    try:
        child_result_ids = list(vsp.GetStringResults(wrapper_result_id, 'ResultsVec'))
    except Exception:
        child_result_ids = []

    for result_id in child_result_ids:
        name = vsp.GetResultsName(result_id)
        result['result_names'].append(name)
        if name == 'VSPAERO_Stab':
            result['stab_result_id'] = result_id

    if not result['stab_result_id']:
        latest_stab_id = vsp.FindLatestResultsID('VSPAERO_Stab')
        if latest_stab_id:
            result['stab_result_id'] = latest_stab_id
            if 'VSPAERO_Stab' not in result['result_names']:
                result['result_names'].append('VSPAERO_Stab')

    if not result['stab_result_id']:
        add('errors', 'MISSING_VSPAERO_STAB', "VSPAERO_Stab was not found after STABILITY_DEFAULT analysis.", {'result_names': result['result_names']})
        return result

    data, columns = [], []
    data_names = list(vsp.GetAllDataNames(result['stab_result_id']))
    result['data_names'] = data_names
    for data_name in data_names:
        values = vsp.GetDoubleResults(result['stab_result_id'], data_name, 0)
        if values:
            data.append(values)
            columns.append(data_name)

    if data:
        result['derivatives'] = pd.DataFrame(np.array(data).T, columns=columns)
    else:
        add('errors', 'EMPTY_VSPAERO_STAB', 'VSPAERO_Stab was found, but no numeric derivative data was readable.')
        return result

    result['passed'] = len(result['errors']) == 0
    if verbose:
        status = 'PASSED' if result['passed'] else 'FAILED'
        print(f'  Stability derivative calculation {status}: {len(result["errors"])} error(s), {len(result["warnings"])} warning(s)')
        if int(verbose) >= 2:
            print('  VSPAERO_Stab data names:')
            print('   ', ', '.join(columns))

    return result

def _get_trim(x, vsp, alpha, mach, reynolds, CMy_tol=5e-4):

    # Set the elevator deflection and control surface parameters
    _ = set_control_surface(geom_name='HTailGeom', deflection=x, cs_group_name='ELEVATOR_GROUP', gains=(1,-1), verbose=0)
    
    # Perform a sweep analysis with the given alpha and mach
    result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[mach], reynolds=[reynolds], verbose=0)
    
    # Retrieve polar results from the sweep analysis
    df = get_polar_result(result_ids)
    
    # Get the CMy (moment coefficient about the y-axis) value
    if list(df.columns):

        CMy = df['CMy'].values[0]
        f = np.abs(CMy)
    
        # Print the current deflection (x) and CMy value
        print(f'\rx = {x:8.4f} ', f'CMy = {CMy:10.6f}', end='')
        
        # If the absolute value of CMy is smaller than the tolerance, raise StopIteration to exit
        if f < CMy_tol:  
            raise StopIteration(x)
        
        # Return the absolute value of CMy to be minimized
        return f
    
    else:
        raise StopIteration(999)

def vsp_trimed_sweep(vsp, alpha_list, Weight, altitude=0, dT=0):

    g = get_gravity() # [m/s^2]
    density = get_density(altitude=altitude, dT=dT) # [kg/m^3]
    Sref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'Sref')[0] # [m^2]
    cref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'cref')[0] # [m]

    # Create an empty DataFrame to store trimmed polar results
    trimed_polar = pd.DataFrame()
    
    # Iterate over the alpha angles
    for alpha in alpha_list:

        result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[0], verbose=0)
        df = get_polar_result(result_ids)
        velocity = np.sqrt((2*Weight*g)/(density*Sref*df['CL'].values[0]))
        mach = velocity_to_mach(velocity=velocity, dT=0, altitude=0)
        reynolds = velocity_to_reynolds(velocity=mach_to_velocity(mach, altitude=altitude, dT=dT), length=cref, altitude=altitude, dT=dT)

        try:
            # Minimize the absolute CMy to find the optimal elevator deflection (x) for the given alpha
            res = minimize_scalar(fun=_get_trim, args=(vsp, alpha, mach, reynolds), method='bounded', bounds=(-20, 20))
            de = res.x
        except StopIteration as err:
            # If StopIteration is raised, use the optimal deflection (err.value) and print 'Done'
            if err.value == 999:
                print('\tError\t', 'Alpha 'f'{alpha:8.3f}, velocity {velocity:8.3f}, mach {mach:8.3f}, reynolds {reynolds:8.3e}')
                de = None
            else:
                print('\tDone\t', 'Alpha 'f'{alpha:8.3f}, velocity {velocity:8.3f}, mach {mach:8.3f}, reynolds {reynolds:8.3e}')
                de = err.value
        
        if de:
            # Set the elevator deflection with the optimized value
            _ = set_control_surface(geom_name='HTailGeom', deflection=de, cs_group_name='ELEVATOR_GROUP', gains=(1,-1), verbose=0)
            
            # Perform another sweep with the optimized elevator deflection
            result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[mach], reynolds=[reynolds], verbose=0)
            
            # Retrieve polar results and add the deflection value to the results DataFrame
            df = get_polar_result(result_ids)
            df['de'] = de
            
            # Append the results to the trimmed polar DataFrame
            trimed_polar = pd.concat([trimed_polar, df])
    
    trimed_polar['gamma'] = np.arctan(1/trimed_polar['L_D'])
    trimed_polar['Velocity'] = np.sqrt((2*Weight*g)/(density*Sref*trimed_polar['CL']*np.cos(trimed_polar['gamma'])))
    trimed_polar['Vx'] = trimed_polar['Velocity']*np.cos(trimed_polar['gamma'])*(3.6)
    trimed_polar['Vz'] = trimed_polar['Velocity']*np.sin(trimed_polar['gamma'])

    return trimed_polar

def make_CDo_correction(vsp, trimed_polar, Weight, CDpCL=0.0065, thickness=0.12, interference_factor=1, xTr=(0,0), altitude=0, dT=0):

    # Function to calculate the frictional drag coefficient in turbulent regions
    def Cf_turbulance(reynolds):
        return 0.455/(np.log10(reynolds)**2.58)

    # Function to calculate the frictional drag coefficient in the laminar flow regime
    def Cf_laminar(reynolds):
        return 1.32824 / np.sqrt(reynolds)
    
    g = get_gravity()
    Sref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'Sref')[0] # Get wing area [m^2]
    density = get_density(altitude=altitude, dT=dT) # Get density based on altitude [kg/m^3]
    reynolds = trimed_polar['Re_1e6'].values * 1e6 # Get Reynolds number

    # Correct drag based on percentage of laminar flow area
    CDo = 0
    for laminar_persent in xTr:
        reynolds_laminar = np.maximum(reynolds * laminar_persent, 1e3)  # Reynolds number in laminar flow region
        CDo += Cf_turbulance(reynolds) - Cf_turbulance(reynolds_laminar) * laminar_persent + Cf_laminar(reynolds_laminar) * laminar_persent

    # form factor
    form_factor = 1 + 2 * thickness + 60 * (thickness ** 4)

    # Corrected CD0, CDtotal, and lift-drag ratio are calculated and added to the data frame
    trimed_polar['CDo_corr'] = CDo * form_factor * interference_factor + CDpCL * (trimed_polar['CL'].values) ** 2
    trimed_polar['CDtot_corr'] = trimed_polar['CDo_corr'].values + trimed_polar['CDi'].values
    trimed_polar['L_D_corr'] = trimed_polar['CL'].values / trimed_polar['CDtot_corr'].values

    # Calculate angle of attack from modified lift-drag ratio
    trimed_polar['gamma'] = np.arctan(1/trimed_polar['L_D_corr'].values)

    # Calculate velocity
    trimed_polar['Velocity'] = np.sqrt((2*Weight*g)/(density*Sref*trimed_polar['CL'].values*np.cos(trimed_polar['gamma'].values)))

    # Calculates horizontal velocity (Vx) and vertical velocity (Vz) and adds them to the data frame
    trimed_polar['Vx'] = trimed_polar['Velocity'].values*np.cos(trimed_polar['gamma'])*(3.6)  # horizontal velocity [km/h]
    trimed_polar['Vz'] = trimed_polar['Velocity'].values*np.sin(trimed_polar['gamma'].values)  # vertical velocity [m/s]

    return trimed_polar  # Returns corrected data

def vsp_sweep_wig(vsp, alpha, mach, reynolds, height, AnalysisMethod=0, verbose=1):

    if verbose:
        print('\n-> Calculate alpha & mach sweep analysis\n')

    # //==== Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File ====//
    
    # Set defaults
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)
    if AnalysisMethod:
        vsp.SetIntAnalysisInput(compgeom_name, 'AnalysisMethod', [AnalysisMethod], 0)

    # List inputs, type, and current values
    if verbose:
        print(compgeom_name)
        vsp.PrintAnalysisInputs(compgeom_name)
        print('')

    # Execute
    if verbose:
        print('\tExecuting...')
    compgeom_resid = vsp.ExecAnalysis(compgeom_name)
    if verbose:
        print('\tCOMPLETE')

    # Get & Display Results
    if verbose:
        vsp.PrintResults(compgeom_resid)
        print('')

    # //==== Analysis: VSPAero Sweep ====//

    # Set defaults
    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))
    if 'UnsteadyType' in analysis_inputs:
        vsp.SetIntAnalysisInput(analysis_name, 'UnsteadyType', [vsp.STABILITY_OFF], 0)

    # Reference geometry set
    geom_set = [0]
    vsp.SetIntAnalysisInput(analysis_name, 'GeomSet', geom_set, 0)
    ref_flag = [1]
    vsp.SetIntAnalysisInput(analysis_name, 'RefFlag', ref_flag, 0)
    wid = vsp.FindGeomsWithName('WingGeom')
    vsp.SetStringAnalysisInput(analysis_name, 'WingID', wid, 0)

    df = pd.DataFrame()
    for re in reynolds:
        for ma in mach:
            for he in height:
                for al in alpha:
                    # Freestream Parameters
                    vsp.SetDoubleAnalysisInput(analysis_name, "AlphaStart", [al], 0)
                    vsp.SetIntAnalysisInput(analysis_name, "AlphaNpts", [1], 0)
                    vsp.SetDoubleAnalysisInput(analysis_name, 'MachStart', [ma], 0)
                    vsp.SetIntAnalysisInput(analysis_name, 'MachNpts', [1], 0)
                    vsp.SetDoubleAnalysisInput(analysis_name, 'ReCref', [re], 0)
                    vsp.SetIntAnalysisInput(analysis_name, 'ReCrefNpts', [1], 0)
                    vsp.Update()

                    # Ground Effect
                    vsp.SetIntAnalysisInput(analysis_name, 'GroundEffectToggle', [1], 0)
                    vsp.SetDoubleAnalysisInput(analysis_name, 'GroundEffect', [he], 0)

                    # List inputs, type, and current values
                    if verbose:
                        print(analysis_name)
                        vsp.PrintAnalysisInputs(analysis_name)
                        print('')

                    # Execute
                    if verbose:
                        print('\tExecuting...')
                        rid = vsp.ExecAnalysis(analysis_name)
                        print('\tCOMPLETE')
                    else:
                        with suppress_stdout():
                            rid = vsp.ExecAnalysis('VSPAEROSweep')
                    
                    tmp = get_polar_result(rid)
                    bref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'bref')[0] # [m]
                    tmp['Alpha'], tmp['Height'], tmp['H_bref'] = al, he, he/bref
                    df = pd.concat([df, tmp])

    return df