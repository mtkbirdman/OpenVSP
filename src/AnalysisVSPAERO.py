# AnalysisVSPAERO.py
import os
import time
import numpy as np
import pandas as pd
from scipy.optimize import root_scalar
from scipy.optimize import minimize_scalar
import openvsp as vsp

from .util import (
    analysis_duration_seconds,
    get_container_parm_value,
    find_one_geom,
    set_control_surface,
    results_dataframe,
    set_analysis_input_if_available,
    suppress_stdout,
)

from .ISAspecification import *

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
        'expected_gains': [-1.0],
    },
}

G103A_REF_GEOM_NAME = 'WingGeom'

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
    vsp.SetDoubleAnalysisInput(analysis_name, 'AlphaStart', [alpha[0]], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, 'AlphaEnd', [alpha[-1]], 0)
    vsp.SetIntAnalysisInput(analysis_name, 'AlphaNpts', alpha_npts, 0)

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

def vsp_trimed_sweep(vsp, alpha_list, Weight, altitude=0, dT=0):
        
    def _get_trim(x, vsp, alpha, mach, reynolds, CMy_tol=5e-4):
        # Set the elevator deflection and control surface parameters
        _ = set_control_surface(vsp, geom_name='HTailGeom', deflection=x, cs_group_name='ELEVATOR_GROUP', gains=(1, -1), verbose=0)

        # Perform a sweep analysis with the given alpha and mach
        result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[mach], reynolds=[reynolds], verbose=0)

        # Retrieve polar results from the sweep analysis
        df = results_dataframe(vsp, result_ids, 'VSPAERO Polar')

        # Get the CMytot (moment coefficient about the y-axis) value
        if list(df.columns):
            CMytot = df['CMytot'].values[0]
            f = np.abs(CMytot)
            # Print the current deflection (x) and CMytot value
            print(f'\rx = {x:8.4f} ', f'CMytot = {CMytot:10.6f}', end='')
            # If the absolute value of CMytot is smaller than the tolerance, raise StopIteration to exit
            if f < CMy_tol:
                raise StopIteration(x)
            # Return the absolute value of CMytot to be minimized
            return f
        else:
            raise StopIteration(999)

    g = get_gravity()  # [m/s*2]
    density = get_density(altitude=altitude, dT=dT)  # [kg/m*3]
    Sref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'Sref')[0]  # [m*2]
    cref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'cref')[0]  # [m]

    # Create an empty DataFrame to store trimmed polar results
    trimed_polar = pd.DataFrame()

    # Iterate over the alpha angles
    for alpha in alpha_list:
        result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[0], verbose=0)
        df = results_dataframe(vsp, result_ids, 'VSPAERO Polar')
        velocity = np.sqrt((2 * Weight * g) / (density * Sref * df['CL'].values[0]))
        mach = velocity_to_mach(velocity=velocity, dT=0, altitude=0)
        reynolds = velocity_to_reynolds(velocity=mach_to_velocity(mach, altitude=altitude, dT=dT), length=cref, altitude=altitude, dT=dT)

        try:
            # Minimize the absolute CMytot to find the optimal elevator deflection (x) for the given alpha
            res = minimize_scalar(fun=_get_trim, args=(vsp, alpha, mach, reynolds), method='bounded', bounds=(-20, 20))
            de = res.x
        except StopIteration as err:
            # If StopIteration is raised, use the optimal deflection (err.value) and print 'Done'
            if err.value == 999:
                print('\tError\t', 'Alpha ' f'{alpha:8.3f}, velocity {velocity:8.3f}, mach {mach:8.3f}, reynolds {reynolds:8.3e}')
                de = None
            else:
                print('\tDone\t', 'Alpha ' f'{alpha:8.3f}, velocity {velocity:8.3f}, mach {mach:8.3f}, reynolds {reynolds:8.3e}')
                de = err.value

        if de:
            # Set the elevator deflection with the optimized value
            _ = set_control_surface(vsp, geom_name='HTailGeom', deflection=de, cs_group_name='ELEVATOR_GROUP', gains=(1, -1), verbose=0)
            # Perform another sweep with the optimized elevator deflection
            result_ids = vsp_sweep(vsp=vsp, alpha=[alpha], mach=[mach], reynolds=[reynolds], verbose=0)
            # Retrieve polar results and add the deflection value to the results DataFrame
            df = results_dataframe(vsp, result_ids, 'VSPAERO Polar')
            df['de'] = de
            # Append the results to the trimmed polar DataFrame
            trimed_polar = pd.concat([trimed_polar, df])

    trimed_polar['gamma'] = np.arctan(1 / trimed_polar['L_D'])
    trimed_polar['Velocity'] = np.sqrt((2 * Weight * g) / (density * Sref * trimed_polar['CL'] * np.cos(trimed_polar['gamma'])))
    trimed_polar['Vx'] = trimed_polar['Velocity'] * np.cos(trimed_polar['gamma']) * 3.6
    trimed_polar['vz'] = trimed_polar['Velocity'] * np.sin(trimed_polar['gamma'])
    return trimed_polar

def vsp_stability_derivatives(
    vsp3_path,
    *,
    alpha=2.0,
    mach=0.1,
    reynolds=4.4e6,
    verbose=1,
    vspaero_verbose=0,
    ncpu=None,
    wake_num_iter=None,
    fixed_wake_flag=None,
    redirect_file='',
    stop_before_run=False,
):
    """
    Run VSPAERO steady 6DOF stability-derivative analysis for a validated
    G103A-style .vsp3 model.

    This function keeps the execution path intentionally direct:
    read .vsp3 -> VSPAEROComputeGeometry -> VSPAEROSweep/STABILITY_DEFAULT
    -> VSPAERO_Stab extraction.

    Parameters added for long Vv-Gamma sweeps
    -----------------------------------------
    ncpu : int or None
        If provided, set VSPAEROSweep/NCPU. OpenVSP passes this to VSPAERO as
        the OpenMP thread count when the solver build supports OpenMP.
    wake_num_iter : int or None
        If provided, set VSPAEROSweep/WakeNumIter, which is written to the
        .vspaero setup file as WakeIters.
    fixed_wake_flag : bool or None
        If provided, set VSPAEROSweep/FixedWakeFlag. This can change the wake
        model and must be treated as an accuracy/speed trade-off.
    redirect_file : str or None
        If not None, set VSPAEROSweep/RedirectFile. Use '' to suppress solver
        stdout, 'stdout' to show it, or a file path to capture it.
    stop_before_run : bool
        If True, write VSPAERO input files and stop before the solver run. This
        is a diagnostic mode; no VSPAERO_Stab result is expected.
    vspaero_verbose : int or bool
        Controls OpenVSP/VSPAERO-internal detail only. The public verbose
        argument remains for human-readable progress.
    """
    import time

    result = {
        'passed': False,
        'errors': [],
        'warnings': [],
        'infos': [],
        'wrapper_result_id': '',
        'compute_geometry_result_id': '',
        'stab_result_id': '',
        'result_names': [],
        'data_names': [],
        'derivatives': pd.DataFrame(),
        'timing': {},
        'stopped_before_run': bool(stop_before_run),
        'vspaero_settings': {
            'ncpu': ncpu,
            'wake_num_iter': wake_num_iter,
            'fixed_wake_flag': fixed_wake_flag,
            'redirect_file': redirect_file,
            'stop_before_run': bool(stop_before_run),
        },
    }

    total_start = time.perf_counter()

    def add(kind, code, message, context=None):
        result[kind].append({
            'code': code,
            'message': message,
            'context': context or {},
        })

    def vprint(level, message):
        if verbose and int(verbose) >= level:
            print(message, flush=True)

    vsp3_path = os.fspath(vsp3_path)
    if not os.path.isfile(vsp3_path):
        add('errors', 'FILE_NOT_FOUND', 'The specified .vsp3 file was not found.', {'vsp3_path': vsp3_path})
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    vprint(1, f'\n-> Calculate G103A VSPAERO stability derivatives: {vsp3_path}')

    read_start = time.perf_counter()
    try:
        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(vsp3_path)
        vsp.Update()
    except Exception as err:
        add('errors', 'READ_VSP3_FAILED', 'OpenVSP failed to read the .vsp3 file.', {'error': repr(err)})
        result['timing']['read_vsp3_elapsed_s'] = time.perf_counter() - read_start
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result
    result['timing']['read_vsp3_elapsed_s'] = time.perf_counter() - read_start

    thick_set = int(vsp.GetSetIndex('ThickGeom'))
    thin_set = int(vsp.GetSetIndex('ThinGeom'))
    if thick_set < 0:
        add('errors', 'MISSING_SET', "Required set 'ThickGeom' was not found.")
    if thin_set < 0:
        add('errors', 'MISSING_SET', "Required set 'ThinGeom' was not found.")
    if result['errors']:
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    wing_ids = list(vsp.FindGeomsWithName(G103A_REF_GEOM_NAME))
    if len(wing_ids) != 1:
        add('errors', 'REF_GEOM_NOT_FOUND', f"Reference Geom '{G103A_REF_GEOM_NAME}' must exist exactly once.", {'geom_ids': wing_ids})
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result
    wing_id = wing_ids[0]

    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)
    compgeom_inputs = set(vsp.GetAnalysisInputNames(compgeom_name))
    set_analysis_input_if_available(vsp, compgeom_name, compgeom_inputs, vsp.SetIntAnalysisInput, 'GeomSet', [thick_set])
    set_analysis_input_if_available(vsp, compgeom_name, compgeom_inputs, vsp.SetIntAnalysisInput, 'ThinGeomSet', [thin_set])

    vprint(1, ' Executing VSPAEROComputeGeometry...')
    compgeom_start = time.perf_counter()
    try:
        if vspaero_verbose:
            compgeom_result_id = vsp.ExecAnalysis(compgeom_name)
        else:
            with suppress_stdout():
                compgeom_result_id = vsp.ExecAnalysis(compgeom_name)
    except Exception as err:
        add('errors', 'VSPAERO_COMPUTE_GEOMETRY_FAILED', 'VSPAEROComputeGeometry failed.', {'error': repr(err)})
        result['timing']['compute_geometry_elapsed_s'] = time.perf_counter() - compgeom_start
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    result['compute_geometry_result_id'] = compgeom_result_id
    result['timing']['compute_geometry_elapsed_s'] = time.perf_counter() - compgeom_start
    result['timing']['compute_geometry_analysis_duration_s'] = analysis_duration_seconds(vsp, compgeom_result_id)

    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'GeomSet', [thick_set])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'ThinGeomSet', [thin_set])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'RefFlag', [1])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetStringAnalysisInput, 'RefGeomID', [wing_id])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'AlphaStart', [alpha])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'AlphaEnd', [alpha])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'AlphaNpts', [1])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'MachStart', [mach])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'MachEnd', [mach])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'MachNpts', [1])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'ReCref', [reynolds])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetDoubleAnalysisInput, 'ReCrefEnd', [reynolds])
    set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'ReCrefNpts', [1])

    if 'UnsteadyType' in analysis_inputs:
        vsp.SetIntAnalysisInput(analysis_name, 'UnsteadyType', [vsp.STABILITY_DEFAULT], 0)
    else:
        add('errors', 'MISSING_UNSTEADYTYPE_INPUT', "VSPAEROSweep does not expose the 'UnsteadyType' input in this OpenVSP version.")
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    # Explicit speed/diagnostic settings. Each setting is only applied when the
    # current OpenVSP build exposes the corresponding Analysis input.
    if ncpu is not None:
        if int(ncpu) < 1:
            add('errors', 'INVALID_NCPU', 'ncpu must be a positive integer.', {'ncpu': ncpu})
        elif not set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'NCPU', [int(ncpu)]):
            add('warnings', 'MISSING_NCPU_INPUT', "VSPAEROSweep does not expose the 'NCPU' input.")

    if wake_num_iter is not None:
        if int(wake_num_iter) < 0:
            add('errors', 'INVALID_WAKE_NUM_ITER', 'wake_num_iter must be zero or a positive integer.', {'wake_num_iter': wake_num_iter})
        elif not set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'WakeNumIter', [int(wake_num_iter)]):
            add('warnings', 'MISSING_WAKE_NUM_ITER_INPUT', "VSPAEROSweep does not expose the 'WakeNumIter' input.")

    if fixed_wake_flag is not None:
        if not set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'FixedWakeFlag', [1 if fixed_wake_flag else 0]):
            add('warnings', 'MISSING_FIXED_WAKE_FLAG_INPUT', "VSPAEROSweep does not expose the 'FixedWakeFlag' input.")

    if redirect_file is not None:
        if not set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetStringAnalysisInput, 'RedirectFile', [str(redirect_file)]):
            add('warnings', 'MISSING_REDIRECT_FILE_INPUT', "VSPAEROSweep does not expose the 'RedirectFile' input.")

    if stop_before_run:
        if not set_analysis_input_if_available(vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput, 'StopBeforeRun', [1]):
            add('warnings', 'MISSING_STOP_BEFORE_RUN_INPUT', "VSPAEROSweep does not expose the 'StopBeforeRun' input.")

    if result['errors']:
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    vprint(
        1,
        'Executing VSPAEROSweep with STABILITY_DEFAULT'
        # f' (NCPU={ncpu}, WakeNumIter={wake_num_iter}, RedirectFile={redirect_file!r})...',
    )

    sweep_start = time.perf_counter()
    try:
        if vspaero_verbose:
            wrapper_result_id = vsp.ExecAnalysis(analysis_name)
        else:
            with suppress_stdout():
                wrapper_result_id = vsp.ExecAnalysis(analysis_name)
    except Exception as err:
        add('errors', 'VSPAERO_SWEEP_FAILED', 'VSPAEROSweep failed.', {'error': repr(err)})
        result['timing']['vspaero_sweep_elapsed_s'] = time.perf_counter() - sweep_start
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    result['wrapper_result_id'] = wrapper_result_id
    result['timing']['vspaero_sweep_elapsed_s'] = time.perf_counter() - sweep_start
    result['timing']['vspaero_sweep_analysis_duration_s'] = analysis_duration_seconds(vsp, wrapper_result_id)

    if stop_before_run:
        add('infos', 'STOPPED_BEFORE_RUN', 'VSPAEROSweep was stopped before solver execution by request.')
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    extraction_start = time.perf_counter()
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
        add('errors', 'MISSING_VSPAERO_STAB', 'VSPAERO_Stab was not found after STABILITY_DEFAULT analysis.', {'result_names': result['result_names']})
        result['timing']['result_extraction_elapsed_s'] = time.perf_counter() - extraction_start
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
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
        result['timing']['result_extraction_elapsed_s'] = time.perf_counter() - extraction_start
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    result['timing']['result_extraction_elapsed_s'] = time.perf_counter() - extraction_start
    result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
    result['passed'] = len(result['errors']) == 0

    if verbose:
        status = 'PASSED' if result['passed'] else 'FAILED'
        print(f' Stability derivative calculation {status}: {len(result["errors"])} error(s), {len(result["warnings"])} warning(s)')
        if int(verbose) >= 2:
            timing = result['timing']
            print(
                ' Timing: '
                f"compute_geometry={timing.get('compute_geometry_elapsed_s', float('nan')):.1f} s, "
                f"vspaero_sweep={timing.get('vspaero_sweep_elapsed_s', float('nan')):.1f} s, "
                f"result_extraction={timing.get('result_extraction_elapsed_s', float('nan')):.1f} s, "
                f"total={timing.get('total_elapsed_s', float('nan')):.1f} s"
            )
    if vspaero_verbose and int(vspaero_verbose) >= 2:
        print(' VSPAERO_Stab data names:')
        print(' ', ', '.join(columns))

    return result

def vsp_sweep_wig(
    vsp,
    alpha_deg,
    mach,
    reynolds,
    height,
    *,
    analysis_method=None,
    trimmed=False,
    elevator_bounds=(-20.0, 20.0),
    cmy_tolerance=5e-4,
    max_trim_evaluations=30,
    verbose=1,
):
    """Run steady VSPAERO ground-effect polar or pitch-trimmed polar cases.

    Parameters
    ----------
    vsp
        Imported ``openvsp`` module with a model already loaded.
    alpha_deg, mach, reynolds, height
        Sequences of angle of attack [deg], Mach number, Reynolds number based
        on ``cref``, and ground-effect height.

        ``height`` is the height above the ground of the CG / moment-reference
        point defined by ``Xcg``, ``Ycg``, and ``Zcg`` in VSPAERO Settings.
        The same CG height is retained for every requested angle of attack.
    analysis_method : int or None, optional
        Legacy fallback used only when the installed OpenVSP build does not
        expose both ``GeomSet`` and ``ThinGeomSet``.
    trimmed : bool, optional
        False returns the polar at the elevator deflection already stored in
        the loaded model. True keeps each requested angle of attack fixed and
        finds the elevator deflection that satisfies ``CMytot = 0``.
    elevator_bounds : tuple[float, float], optional
        Lower and upper elevator-deflection limits [deg] used to bracket the
        pitch-trim solution.
    cmy_tolerance : float, optional
        Maximum accepted absolute value of the final ``CMytot``.
    max_trim_evaluations : int, optional
        Maximum number of VSPAERO evaluations allowed for one trimmed case,
        including the final verification run.
    verbose : int or bool, optional
        0: no progress output.
        1: model summary and one final result summary per requested case.
        2: level 1 plus every VSPAERO evaluation, analysis input, and result.

    Returns
    -------
    pandas.DataFrame
        One row for each requested angle/Mach/Reynolds/ground condition. The
        ground-effect-off reference case is included once for every
        angle/Mach/Reynolds combination. Its ``CGHeight`` columns are NaN.

        When ``trimmed=True``, each row is the final verified ``CMytot = 0``
        result at the specified angle of attack. Required-CL interpolation is
        intentionally outside this function.
    """

    error_mgr = vsp.ErrorMgrSingleton.getInstance()

    def pop_openvsp_errors():
        errors = []
        while error_mgr.GetNumTotalErrors() > 0:
            error = error_mgr.PopLastError()
            errors.append(error.GetErrorString())
        return errors

    alpha_deg = list(alpha_deg)
    mach = list(mach)
    reynolds = list(reynolds)
    height = list(height)
    if not alpha_deg or not mach or not reynolds or not height:
        raise ValueError('alpha_deg, mach, reynolds, and height must not be empty.')

    if trimmed:
        elevator_lower = float(elevator_bounds[0])
        elevator_upper = float(elevator_bounds[1])
        if elevator_lower >= elevator_upper:
            raise ValueError('elevator_bounds must satisfy lower < upper.')
        if cmy_tolerance <= 0.0:
            raise ValueError('cmy_tolerance must be positive.')
        if int(max_trim_evaluations) < 4:
            raise ValueError('max_trim_evaluations must be at least 4.')

    vsp.Update()
    wing_id = find_one_geom(vsp, G103A_REF_GEOM_NAME)
    vsp3_path = vsp.GetVSPFileName()
    openvsp_version = vsp.GetVSPVersion()

    settings_id = vsp.FindContainer('VSPAEROSettings', 0)
    settings = {}
    if settings_id:
        for name in ('Sref', 'bref', 'cref', 'Xcg', 'Ycg', 'Zcg', 'Symmetry', 'RefFlag'):
            value, parm_id = get_container_parm_value(vsp, settings_id, name)
            if parm_id:
                settings[name] = value

    # Locate the existing elevator group once. Polar-only calculations do not
    # require it, but its current deflection is recorded when available.
    elevator_parm_id = ''
    initial_elevator_deg = np.nan
    elevator_group_names = [
        vsp.GetVSPAEROControlGroupName(group_index)
        for group_index in range(vsp.GetNumControlSurfaceGroups())
    ]
    elevator_group_matches = [
        group_index
        for group_index, group_name in enumerate(elevator_group_names)
        if group_name == 'ELEVATOR_GROUP'
    ]

    if len(elevator_group_matches) == 1 and settings_id:
        elevator_group_index = elevator_group_matches[0]
        elevator_group_name = f'ControlSurfaceGroup_{elevator_group_index}'
        elevator_parm_id = vsp.FindParm(
            settings_id,
            'DeflectionAngle',
            elevator_group_name,
        )
        if elevator_parm_id and str(elevator_parm_id).upper() != 'NONE':
            initial_elevator_deg = float(vsp.GetParmVal(elevator_parm_id))
        else:
            elevator_parm_id = ''

    if trimmed:
        if len(elevator_group_matches) != 1:
            raise RuntimeError(
                "trimmed=True requires exactly one existing 'ELEVATOR_GROUP'. "
                f'found={len(elevator_group_matches)}'
            )
        active_elevator_surfaces = list(
            vsp.GetActiveCSNameVec(elevator_group_matches[0])
        )
        if not active_elevator_surfaces:
            raise RuntimeError("'ELEVATOR_GROUP' has no active control surfaces.")
        if not elevator_parm_id:
            raise RuntimeError(
                "Could not find 'DeflectionAngle' for the existing "
                "'ELEVATOR_GROUP'."
            )
        if not elevator_lower <= initial_elevator_deg <= elevator_upper:
            raise ValueError(
                'The elevator deflection stored in the loaded model must lie '
                f'within elevator_bounds. initial={initial_elevator_deg}, '
                f'bounds={elevator_bounds}'
            )

    # ------------------------------------------------------------------
    # 1. Compute the mixed thick/thin VSPAERO geometry once.
    # ------------------------------------------------------------------
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)
    compgeom_inputs = set(vsp.GetAnalysisInputNames(compgeom_name))

    use_saved_thick_thin = {'GeomSet', 'ThinGeomSet'} <= compgeom_inputs
    if use_saved_thick_thin:
        thick_set = int(vsp.GetIntAnalysisInput(compgeom_name, 'GeomSet')[0])
        thin_set = int(vsp.GetIntAnalysisInput(compgeom_name, 'ThinGeomSet')[0])
        set_source = 'loaded .vsp3 VSPAERO settings'
    else:
        thick_set = 0
        thin_set = None
        set_source = 'legacy fallback: GeomSet=0'
        if 'GeomSet' in compgeom_inputs:
            vsp.SetIntAnalysisInput(compgeom_name, 'GeomSet', [thick_set], 0)
        if analysis_method is not None and 'AnalysisMethod' in compgeom_inputs:
            vsp.SetIntAnalysisInput(
                compgeom_name,
                'AnalysisMethod',
                [int(analysis_method)],
                0,
            )

    set_none = int(getattr(vsp, 'SET_NONE', -1))
    thick_geom_names = []
    if thick_set != set_none:
        thick_geom_names = [
            vsp.GetGeomName(geom_id)
            for geom_id in vsp.GetGeomSetAtIndex(thick_set)
        ]
    thin_geom_names = []
    if thin_set is not None and thin_set != set_none:
        thin_geom_names = [
            vsp.GetGeomName(geom_id)
            for geom_id in vsp.GetGeomSetAtIndex(thin_set)
        ]

    ground_condition_count = 1 + len(height)
    requested_case_count = (
        len(alpha_deg)
        * len(mach)
        * len(reynolds)
        * ground_condition_count
    )

    if verbose:
        mode = 'trimmed polar (fixed alpha, CMytot=0)' if trimmed else 'polar'
        print('\n-> VSPAERO ground-effect sweep', flush=True)
        print(f' Mode          : {mode}', flush=True)
        print(f' OpenVSP       : {openvsp_version}', flush=True)
        print(f' VSP3          : {vsp3_path or "<unsaved model>"}', flush=True)
        print(f' Reference wing: {G103A_REF_GEOM_NAME} ({wing_id})', flush=True)
        print(f' Geometry sets : {set_source}', flush=True)
        print(f'   Thick [{thick_set}]: {thick_geom_names or ["<none>"]}', flush=True)
        print(f'   Thin  [{thin_set}]: {thin_geom_names or ["<none>"]}', flush=True)
        if settings:
            print(
                ' VSPAERO refs  : '
                + ', '.join(f'{name}={value:.6g}' for name, value in settings.items()),
                flush=True,
            )
        if trimmed:
            print(
                f' Elevator trim : initial={initial_elevator_deg:.6g} deg, '
                f'bounds=({elevator_lower:.6g}, {elevator_upper:.6g}) deg, '
                f'|CMytot|<={cmy_tolerance:.3g}',
                flush=True,
            )
        print(
            f' Cases         : {requested_case_count} '
            f'({len(alpha_deg)} alpha x {len(mach)} Mach x '
            f'{len(reynolds)} Re x {ground_condition_count} ground conditions)',
            flush=True,
        )

    if verbose and int(verbose) >= 2:
        print(f'\n{compgeom_name} inputs')
        vsp.PrintAnalysisInputs(compgeom_name)

    pop_openvsp_errors()
    compgeom_start = time.perf_counter()
    if verbose and int(verbose) >= 2:
        compgeom_result_id = vsp.ExecAnalysis(compgeom_name)
    else:
        with suppress_stdout():
            compgeom_result_id = vsp.ExecAnalysis(compgeom_name)
    compgeom_errors = pop_openvsp_errors()
    compgeom_elapsed = time.perf_counter() - compgeom_start

    compgeom_result_name = (
        vsp.GetResultsName(compgeom_result_id)
        if compgeom_result_id
        else ''
    )
    vspgeom_files = (
        list(vsp.GetStringResults(compgeom_result_id, 'VSPGeomFileName'))
        if compgeom_result_id
        else []
    )
    if (
        not compgeom_result_id
        or compgeom_result_name != 'VSPAERO_Geom'
        or not vspgeom_files
        or not os.path.isfile(vspgeom_files[0])
    ):
        raise RuntimeError(
            'VSPAEROComputeGeometry did not produce a usable VSPAERO_Geom result.\n'
            f'  result_id={compgeom_result_id!r}\n'
            f'  result_name={compgeom_result_name!r}\n'
            f'  vspgeom_files={vspgeom_files!r}\n'
            f'  OpenVSP errors={compgeom_errors or ["<none reported>"]}\n'
            'Likely causes: the selected Thick/Thin sets contain no usable geometry; '
            'the geometry mesh could not be generated or written; the working directory '
            'is not writable; or the installed OpenVSP version does not support the '
            'requested analysis inputs.'
        )

    if verbose:
        print(
            f' Geometry ready : {vspgeom_files[0]} '
            f'({compgeom_elapsed:.2f} s)',
            flush=True,
        )
    if verbose and int(verbose) >= 2:
        vsp.PrintResults(compgeom_result_id)

    # ------------------------------------------------------------------
    # 2. Configure a one-point VSPAERO analysis reused by polar and trim runs.
    # ------------------------------------------------------------------
    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))

    required_inputs = {
        'GeomSet', 'RefFlag', 'WingID',
        'AlphaStart', 'AlphaEnd', 'AlphaNpts',
        'MachStart', 'MachEnd', 'MachNpts',
        'ReCref', 'ReCrefEnd', 'ReCrefNpts',
        'GroundEffectToggle', 'GroundEffect',
    }
    missing_inputs = sorted(required_inputs - analysis_inputs)
    if missing_inputs:
        raise RuntimeError(
            f'{analysis_name} is missing required inputs: {missing_inputs}. '
            'The installed OpenVSP API is not compatible with this implementation.'
        )

    if use_saved_thick_thin and 'ThinGeomSet' in analysis_inputs:
        vsp.SetIntAnalysisInput(analysis_name, 'GeomSet', [thick_set], 0)
        vsp.SetIntAnalysisInput(analysis_name, 'ThinGeomSet', [thin_set], 0)
    else:
        vsp.SetIntAnalysisInput(analysis_name, 'GeomSet', [0], 0)

    if 'UnsteadyType' in analysis_inputs:
        vsp.SetIntAnalysisInput(
            analysis_name,
            'UnsteadyType',
            [vsp.STABILITY_OFF],
            0,
        )
    if 'RedirectFile' in analysis_inputs:
        redirect_file = 'stdout' if verbose and int(verbose) >= 2 else ''
        vsp.SetStringAnalysisInput(
            analysis_name,
            'RedirectFile',
            [redirect_file],
            0,
        )

    component_ref = int(getattr(vsp, 'COMPONENT_REF', 1))
    vsp.SetIntAnalysisInput(analysis_name, 'RefFlag', [component_ref], 0)
    vsp.SetStringAnalysisInput(analysis_name, 'WingID', [wing_id], 0)

    bref = float(vsp.GetDoubleAnalysisInput(analysis_name, 'bref')[0])
    if bref <= 0.0:
        raise ValueError(f'VSPAERO reference span bref must be positive. bref={bref}')

    def run_single_case(
        alpha,
        mach_value,
        reynolds_value,
        ground_effect_enabled,
        cg_height,
        elevator_deg=None,
        evaluation_label='',
    ):
        if elevator_parm_id and elevator_deg is not None:
            vsp.SetParmVal(elevator_parm_id, float(elevator_deg))

        vsp.SetDoubleAnalysisInput(analysis_name, 'AlphaStart', [float(alpha)], 0)
        vsp.SetDoubleAnalysisInput(analysis_name, 'AlphaEnd', [float(alpha)], 0)
        vsp.SetIntAnalysisInput(analysis_name, 'AlphaNpts', [1], 0)
        vsp.SetDoubleAnalysisInput(analysis_name, 'MachStart', [float(mach_value)], 0)
        vsp.SetDoubleAnalysisInput(analysis_name, 'MachEnd', [float(mach_value)], 0)
        vsp.SetIntAnalysisInput(analysis_name, 'MachNpts', [1], 0)
        vsp.SetDoubleAnalysisInput(analysis_name, 'ReCref', [float(reynolds_value)], 0)
        vsp.SetDoubleAnalysisInput(analysis_name, 'ReCrefEnd', [float(reynolds_value)], 0)
        vsp.SetIntAnalysisInput(analysis_name, 'ReCrefNpts', [1], 0)
        vsp.SetIntAnalysisInput(
            analysis_name,
            'GroundEffectToggle',
            [1 if ground_effect_enabled else 0],
            0,
        )
        vsp.SetDoubleAnalysisInput(
            analysis_name,
            'GroundEffect',
            [float(cg_height) if ground_effect_enabled else 0.0],
            0,
        )
        vsp.Update()

        if verbose and int(verbose) >= 2:
            elevator_text = (
                '<unchanged>' if elevator_deg is None else f'{float(elevator_deg):.8g}'
            )
            print(
                f'   evaluation {evaluation_label}: '
                f'alpha={float(alpha):.8g}, Mach={float(mach_value):.8g}, '
                f'Re={float(reynolds_value):.8g}, '
                f'ground={"ON" if ground_effect_enabled else "OFF"}, '
                f'CGHeight={cg_height if ground_effect_enabled else np.nan}, '
                f'elevator={elevator_text}',
                flush=True,
            )
            vsp.PrintAnalysisInputs(analysis_name)

        pop_openvsp_errors()
        evaluation_start = time.perf_counter()
        if verbose and int(verbose) >= 2:
            wrapper_result_id = vsp.ExecAnalysis(analysis_name)
        else:
            with suppress_stdout():
                wrapper_result_id = vsp.ExecAnalysis(analysis_name)
        api_errors = pop_openvsp_errors()
        elapsed = time.perf_counter() - evaluation_start

        child_result_ids = (
            list(vsp.GetStringResults(wrapper_result_id, 'ResultsVec'))
            if wrapper_result_id
            else []
        )
        child_result_names = [
            vsp.GetResultsName(result_id)
            for result_id in child_result_ids
        ]
        case_result = (
            results_dataframe(vsp, wrapper_result_id, 'VSPAERO_Polar')
            if wrapper_result_id
            else pd.DataFrame()
        )

        if case_result.empty:
            raise RuntimeError(
                'VSPAERO returned no VSPAERO_Polar result.\n'
                f'  alpha={alpha}, Mach={mach_value}, Re={reynolds_value}\n'
                f'  ground_effect_enabled={ground_effect_enabled}, '
                f'CGHeight={cg_height if ground_effect_enabled else None}\n'
                f'  elevator_deg={elevator_deg}\n'
                f'  wrapper_result_id={wrapper_result_id!r}\n'
                f'  child_results={child_result_names!r}\n'
                f'  OpenVSP errors={api_errors or ["<none reported>"]}\n'
                'Likely causes: the VSPAERO executable/path is unavailable; '
                'the .vspgeom or setup file is invalid; the selected Thick/Thin '
                'sets are inconsistent; the ground plane intersects or is too '
                'close to the model; the solver terminated before writing output; '
                'or the output files could not be read.'
            )
        if len(case_result.index) != 1:
            raise RuntimeError(
                'A one-point VSPAERO analysis returned an unexpected number of rows. '
                f'rows={len(case_result.index)}'
            )

        duration = analysis_duration_seconds(vsp, wrapper_result_id)
        duration_seconds = elapsed if duration is None else duration

        if verbose and int(verbose) >= 2:
            vsp.PrintResults(wrapper_result_id)
            print(case_result.to_string(index=False), flush=True)

        return case_result.copy(), float(duration_seconds)

    # ------------------------------------------------------------------
    # 3. Run the ground-effect-off reference and each requested CG height.
    # ------------------------------------------------------------------
    ground_conditions = [(False, None)] + [(True, float(value)) for value in height]
    case_frames = []
    case_index = 0

    try:
        for reynolds_value in reynolds:
            for mach_value in mach:
                for ground_effect_enabled, cg_height in ground_conditions:
                    for alpha in alpha_deg:
                        case_index += 1
                        ground_text = (
                            f'ON, CGHeight={cg_height:.8g}'
                            if ground_effect_enabled
                            else 'OFF'
                        )
                        if verbose:
                            print(
                                f' [{case_index:>3}/{requested_case_count}] '
                                f'alpha={float(alpha):.8g} deg, '
                                f'Mach={float(mach_value):.8g}, '
                                f'Re={float(reynolds_value):.8g}, '
                                f'ground={ground_text}, '
                                f'trimmed={bool(trimmed)}',
                                flush=True,
                            )

                        if not trimmed:
                            case_result, case_duration = run_single_case(
                                alpha,
                                mach_value,
                                reynolds_value,
                                ground_effect_enabled,
                                cg_height,
                                elevator_deg=None,
                                evaluation_label='polar',
                            )
                            elevator_deg = initial_elevator_deg
                            cmy_residual = (
                                float(case_result['CMytot'].iloc[0])
                                if 'CMytot' in case_result.columns
                                else np.nan
                            )
                            trim_evaluations = 0
                        else:
                            evaluation_count = 0
                            evaluation_cache = {}
                            case_duration = 0.0
                            solver_evaluation_limit = int(max_trim_evaluations) - 1

                            def evaluate_cmy(elevator_deg):
                                nonlocal evaluation_count, case_duration
                                key = round(float(elevator_deg), 12)
                                if key in evaluation_cache:
                                    return evaluation_cache[key][0]
                                if evaluation_count >= solver_evaluation_limit:
                                    raise RuntimeError(
                                        'Pitch-trim evaluation limit reached before '
                                        'the final verification run. '
                                        f'max_trim_evaluations={max_trim_evaluations}'
                                    )

                                result, duration = run_single_case(
                                    alpha,
                                    mach_value,
                                    reynolds_value,
                                    ground_effect_enabled,
                                    cg_height,
                                    elevator_deg=float(elevator_deg),
                                    evaluation_label=f'trim {evaluation_count + 1}',
                                )
                                if 'CMytot' not in result.columns:
                                    raise RuntimeError(
                                        "VSPAERO_Polar does not contain the required 'CMytot' column."
                                    )
                                cmy = float(result['CMytot'].iloc[0])
                                if not np.isfinite(cmy):
                                    raise RuntimeError(
                                        f'VSPAERO returned a non-finite CMytot value: {cmy}'
                                    )

                                evaluation_count += 1
                                case_duration += duration
                                evaluation_cache[key] = (cmy, result)
                                return cmy

                            initial_cmy = evaluate_cmy(initial_elevator_deg)
                            trim_elevator_deg = float(initial_elevator_deg)

                            if abs(initial_cmy) > cmy_tolerance:
                                lower_cmy = evaluate_cmy(elevator_lower)
                                upper_cmy = evaluate_cmy(elevator_upper)

                                exact_root = None
                                for candidate_deg, candidate_cmy in (
                                    (elevator_lower, lower_cmy),
                                    (elevator_upper, upper_cmy),
                                ):
                                    if abs(candidate_cmy) <= cmy_tolerance:
                                        exact_root = float(candidate_deg)
                                        break

                                if exact_root is not None:
                                    trim_elevator_deg = exact_root
                                else:
                                    brackets = []
                                    if lower_cmy * initial_cmy < 0.0:
                                        brackets.append((elevator_lower, initial_elevator_deg))
                                    if initial_cmy * upper_cmy < 0.0:
                                        brackets.append((initial_elevator_deg, elevator_upper))

                                    if not brackets:
                                        raise RuntimeError(
                                            'No CMytot=0 elevator-trim root was bracketed within '
                                            f'elevator_bounds={elevator_bounds}.\n'
                                            f'  alpha={alpha}, Mach={mach_value}, '
                                            f'Re={reynolds_value}\n'
                                            f'  ground_effect_enabled={ground_effect_enabled}, '
                                            f'CGHeight={cg_height}\n'
                                            f'  CMytot(lower={elevator_lower})={lower_cmy}\n'
                                            f'  CMytot(initial={initial_elevator_deg})={initial_cmy}\n'
                                            f'  CMytot(upper={elevator_upper})={upper_cmy}'
                                        )

                                    bracket = min(
                                        brackets,
                                        key=lambda values: (
                                            abs(values[1] - values[0]),
                                            abs(0.5 * (values[0] + values[1]) - initial_elevator_deg),
                                        ),
                                    )
                                    remaining_evaluations = (
                                        solver_evaluation_limit - evaluation_count
                                    )
                                    if remaining_evaluations <= 0:
                                        raise RuntimeError(
                                            'No VSPAERO evaluations remain for the '
                                            'pitch-trim root solve.'
                                        )

                                    root_result = root_scalar(
                                        evaluate_cmy,
                                        bracket=bracket,
                                        method='brentq',
                                        xtol=1e-4,
                                        rtol=1e-8,
                                        maxiter=max(1, remaining_evaluations),
                                    )
                                    if not root_result.converged:
                                        raise RuntimeError(
                                            'The pitch-trim root solver did not converge. '
                                            f'flag={root_result.flag!r}'
                                        )
                                    trim_elevator_deg = float(root_result.root)

                            if evaluation_count >= int(max_trim_evaluations):
                                raise RuntimeError(
                                    'No VSPAERO evaluation remains for final trim verification.'
                                )

                            case_result, final_duration = run_single_case(
                                alpha,
                                mach_value,
                                reynolds_value,
                                ground_effect_enabled,
                                cg_height,
                                elevator_deg=trim_elevator_deg,
                                evaluation_label='trim verification',
                            )
                            evaluation_count += 1
                            case_duration += final_duration
                            cmy_residual = float(case_result['CMytot'].iloc[0])
                            if abs(cmy_residual) > cmy_tolerance:
                                raise RuntimeError(
                                    'Pitch trim failed final verification.\n'
                                    f'  alpha={alpha}, Mach={mach_value}, '
                                    f'Re={reynolds_value}\n'
                                    f'  ground_effect_enabled={ground_effect_enabled}, '
                                    f'CGHeight={cg_height}\n'
                                    f'  elevator_deg={trim_elevator_deg}\n'
                                    f'  CMytot={cmy_residual}\n'
                                    f'  tolerance={cmy_tolerance}'
                                )

                            elevator_deg = trim_elevator_deg
                            trim_evaluations = evaluation_count

                        case_result.insert(0, 'Trimmed', bool(trimmed))
                        case_result.insert(1, 'GroundEffectEnabled', bool(ground_effect_enabled))
                        case_result.insert(2, 'InputAlpha_deg', float(alpha))
                        case_result.insert(3, 'InputMach', float(mach_value))
                        case_result.insert(4, 'InputReCref', float(reynolds_value))
                        case_result.insert(
                            5,
                            'CGHeight',
                            float(cg_height) if ground_effect_enabled else np.nan,
                        )
                        case_result.insert(
                            6,
                            'CGHeight_bref',
                            float(cg_height) / bref if ground_effect_enabled else np.nan,
                        )
                        case_result.insert(7, 'Elevator_deg', float(elevator_deg))
                        case_result.insert(8, 'CMyResidual', cmy_residual)
                        case_result.insert(
                            9,
                            'TrimFunctionEvaluations',
                            int(trim_evaluations),
                        )
                        case_frames.append(case_result)

                        if verbose:
                            summary_columns = [
                                name
                                for name in (
                                    'Elevator_deg', 'CL', 'CDtot', 'CDi',
                                    'CD0', 'L_D', 'CMytot', 'CMyResidual',
                                )
                                if name in case_result.columns
                            ]
                            print(f'   complete: {case_duration:.2f} s', flush=True)
                            print(
                                case_result[summary_columns].to_string(index=False),
                                flush=True,
                            )
    finally:
        if trimmed and elevator_parm_id:
            vsp.SetParmVal(elevator_parm_id, float(initial_elevator_deg))
            vsp.Update()

    return pd.concat(case_frames, ignore_index=True)

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

    vprint(1, ' Checking required Geoms...')
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
            add(
                'errors',
                'DUPLICATE_GEOM',
                f"Required Geom '{geom_name}' is ambiguous because multiple Geoms have the same name.",
                {'geom_name': geom_name, 'geom_ids': ids},
            )
        else:
            geom_id = ids[0]
            geom_ids[geom_name] = geom_id
            type_name = normalize_type_name(vsp.GetGeomTypeName(geom_id))
            summary['type_name'] = type_name
            summary['passed'] = type_name == expected_type
            if type_name != expected_type:
                add(
                    'errors',
                    'GEOM_TYPE_MISMATCH',
                    f"Geom '{geom_name}' has type '{type_name}', expected '{expected_type}'.",
                    {'geom_name': geom_name, 'geom_id': geom_id},
                )
        report['geom_summary'][geom_name] = summary

    # 2. Required control-surface Subsurfaces.
    subsurface_ids = {}

    vprint(1, ' Checking required control-surface Subsurfaces...')
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
            add(
                'errors',
                'SUBSURFACE_PARENT_MISSING',
                f"Cannot check Subsurface '{subsurface_name}' because parent Geom '{geom_name}' is missing.",
                {'geom_name': geom_name},
            )
            report['subsurface_summary'][geom_name] = summary
            continue

        sub_ids = list(vsp.GetSubSurfIDVec(geom_id))
        matches = [sid for sid in sub_ids if vsp.GetSubSurfName(sid) == subsurface_name]
        if len(matches) == 0:
            add(
                'errors',
                'MISSING_SUBSURFACE',
                f"Required Subsurface '{subsurface_name}' was not found on '{geom_name}'.",
                {'geom_name': geom_name, 'available_subsurfaces': [vsp.GetSubSurfName(sid) for sid in sub_ids]},
            )
        elif len(matches) > 1:
            add(
                'errors',
                'DUPLICATE_SUBSURFACE',
                f"Subsurface '{subsurface_name}' appears multiple times on '{geom_name}'.",
                {'geom_name': geom_name, 'subsurface_ids': matches},
            )
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
                add(
                    'errors',
                    'SUBSURFACE_TYPE_MISMATCH',
                    f"Subsurface '{subsurface_name}' on '{geom_name}' is not SS_CONTROL.",
                    {'geom_name': geom_name, 'subsurface_id': subsurface_id, 'subsurface_type': subsurface_type},
                )
        report['subsurface_summary'][geom_name] = summary

    # 4. ThickGeom / ThinGeom sets.
    vprint(1, ' Checking ThickGeom / ThinGeom sets...')
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
    vprint(1, ' Checking VSPAERO settings...')
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
            value, parm_id = get_container_parm_value(vsp, vspaero_settings_id, parm_name)
            settings_summary[parm_name] = value
            if value is None:
                add('errors', 'MISSING_VSPAERO_PARM', f"VSPAERO setting '{parm_name}' was not found.", {'parm_name': parm_name})

        symmetry_value, symmetry_parm_id = get_container_parm_value(vsp, vspaero_settings_id, 'Symmetry')
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
    vprint(1, ' Checking VSPAERO control-surface groups and gains...')
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
            add(
                'warnings',
                'CONTROL_GROUP_ACTIVE_NAME_UNCLEAR',
                f"Control-surface group '{group_name}' is active, but its active names do not clearly include both the expected Geom and Subsurface names.",
                {'group_name': group_name, 'active_control_surfaces': active_names},
            )

        subsurface_id = subsurface_ids.get((geom_name, subsurface_name))
        if vspaero_settings_id and subsurface_id:
            for gain_index, expected_gain in enumerate(expected_gains):
                parm_name = f'Surf_{subsurface_id}_{gain_index}_Gain'
                value, parm_id = get_container_parm_value(vsp, vspaero_settings_id, parm_name)
                summary['actual_gains'].append(value)
                if value is None:
                    add('errors', 'MISSING_CONTROL_GAIN', f"Gain parm '{parm_name}' was not found for group '{group_name}'.", {'group_name': group_name, 'parm_name': parm_name})
                elif abs(float(value) - float(expected_gain)) > 1e-6:
                    add(
                        'errors',
                        'CONTROL_GAIN_MISMATCH',
                        f"Control-surface group '{group_name}' has unexpected gain.",
                        {'group_name': group_name, 'parm_name': parm_name, 'expected_gain': expected_gain, 'actual_gain': value},
                    )
        summary['passed'] = (
            summary['found']
            and bool(summary['active_control_surfaces'])
            and len(summary['actual_gains']) == len(expected_gains)
            and all(value is not None and abs(float(value) - float(expected)) <= 1e-6 for value, expected in zip(summary['actual_gains'], expected_gains))
        )
        report['control_group_summary'][group_name] = summary

    settings_summary['passed'] = not any(
        item['code'].startswith('VSPAERO_') or item['code'].startswith('INVALID_') or item['code'].startswith('MISSING_VSPAERO_')
        for item in report['errors']
    )
    report['passed'] = len(report['errors']) == 0

    if verbose:
        status = 'PASSED' if report['passed'] else 'FAILED'
        print(f' Validation {status}: {len(report["errors"])} error(s), {len(report["warnings"])} warning(s)')
        if int(verbose) >= 2:
            for error in report['errors']:
                print(f"    ERROR {error['code']}: {error['message']}")
            for warning in report['warnings']:
                print(f"    WARNING {warning['code']}: {warning['message']}")

    return report

def make_CDo_correction(
    vsp, 
    trimed_polar, 
    Weight, 
    CDpCL=0.0065, 
    thickness=0.12, 
    interference_factor=1, 
    xTr=(0, 0), 
    altitude=0, 
    dT=0
):

    # Function to calculate the frictional drag coefficient in turbulent regions
    def Cf_turbulance(reynolds):
        return 0.455 / (np.log10(reynolds) ** 2.58)

    # Function to calculate the frictional drag coefficient in the laminar flow regime
    def Cf_laminar(reynolds):
        return 1.32824 / np.sqrt(reynolds)

    g = get_gravity()
    Sref = vsp.GetDoubleAnalysisInput('VSPAEROSweep', 'Sref')[0]  # Get wing area [m2]
    density = get_density(altitude=altitude, dT=dT)  # Get density based on altitude [kg/m*3]
    reynolds = trimed_polar['Re_1e6'].values * 1e6  # Get Reynolds number

    # Correct drag based on percentage of laminar flow area
    CDo = 0
    for laminar_persent in xTr:
        reynolds_laminar = np.maximum(reynolds * laminar_persent, 1e3)  # Reynolds number in laminar flow region
        CDo += (
            Cf_turbulance(reynolds)
            - Cf_turbulance(reynolds_laminar) * laminar_persent
            + Cf_laminar(reynolds_laminar) * laminar_persent
        )

    # form factor
    form_factor = 1 + 2 * thickness + 60 * (thickness ** 4)

    # Corrected CD0, CDtotal, and lift-drag ratio are calculated and added to the data frame
    trimed_polar['CDo_corr'] = CDo * form_factor * interference_factor + CDpCL * (trimed_polar['CL'].values) ** 2
    trimed_polar['CDtot_corr'] = trimed_polar['CDo_corr'].values + trimed_polar['CDi'].values
    trimed_polar['L_D_corr'] = trimed_polar['CL'].values / trimed_polar['CDtot_corr'].values

    # Calculate angle of attack from modified lift-drag ratio
    trimed_polar['gamma'] = np.arctan(1 / trimed_polar['L_D_corr'].values)

    # Calculate velocity
    trimed_polar['Velocity'] = np.sqrt((2 * Weight * g) / (density * Sref * trimed_polar['CL'].values * np.cos(trimed_polar['gamma'].values)))

    # Calculates horizontal velocity (Vx) and vertical velocity (Vz) and adds them to the data frame
    trimed_polar['Vx'] = trimed_polar['Velocity'].values * np.cos(trimed_polar['gamma'].values) * 3.6  # horizontal velocity [km/h]
    trimed_polar['vz'] = trimed_polar['Velocity'].values * np.sin(trimed_polar['gamma'].values)  # vertical velocity [m/s]
    return trimed_polar  # Returns corrected data
