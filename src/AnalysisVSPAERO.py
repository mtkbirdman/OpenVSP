# AnalysisVSPAERO.py
import os
import time
import warnings
import numpy as np
import pandas as pd
import openvsp as vsp

from .util import (
    analysis_duration_seconds,
    get_container_parm_value,
    find_one_geom,
    find_vspaero_settings_container,
    results_dataframe,
    set_analysis_input_if_available,
    suppress_stdout,
)

from .ISAspecification import *

# Default reference-geometry name retained by legacy trim workflows.
DEFAULT_REF_GEOM_NAME = 'WingGeom'

def vsp_sweep(
    vsp,
    alpha,
    mach,
    reynolds=[1e6],
    verbose=1,
    *,
    ncpu=None,
    wake_num_iter=None,
    wake_num_nodes=None,
    fixed_wake_flag=None,
    thick_geom_set=None,
    thin_geom_set=None,
):
    """Run an OpenVSP VSPAERO angle-of-attack and Mach-number sweep.

    Parameters
    ----------
    vsp : module or object
        OpenVSP Python API object used to configure and execute analyses.
    alpha : sequence of float
        Angle-of-attack sweep values in degrees. The first and last values
        define ``AlphaStart`` and ``AlphaEnd``.
    mach : sequence of float
        Mach-number sweep values. The first and last values define
        ``MachStart`` and ``MachEnd``.
    reynolds : sequence of float, optional
        Reynolds-number sweep values based on VSPAERO reference chord.
    verbose : int or bool, optional
        Print OpenVSP analysis inputs, progress, and results when truthy.
    ncpu : int or None, optional
        Value assigned to ``VSPAEROSweep/NCPU``. ``None`` leaves the OpenVSP
        default unchanged.
    wake_num_iter : int or None, optional
        Value assigned to ``VSPAEROSweep/WakeNumIter``. Valid explicit values
        are 3 through 255. ``None`` leaves the OpenVSP default unchanged.
    wake_num_nodes : int or None, optional
        Value assigned to ``VSPAEROSweep/NumWakeNodes``. ``None`` leaves the
        OpenVSP default unchanged.
    fixed_wake_flag : bool or None, optional
        Value assigned to ``VSPAEROSweep/FixedWakeFlag``. When enabled, the
        fixed-wake setting takes precedence over ``wake_num_iter``. ``None``
        leaves the OpenVSP default unchanged.
    thick_geom_set : int or None, optional
        Thick-geometry set index used by both ``VSPAEROComputeGeometry`` and
        ``VSPAEROSweep``. ``None`` preserves the selection saved in the loaded
        model. Pass ``vsp.SET_NONE`` to request no thick geometry explicitly.
    thin_geom_set : int or None, optional
        Thin-geometry set index used by both ``VSPAEROComputeGeometry`` and
        ``VSPAEROSweep``. ``None`` preserves the selection saved in the loaded
        model. Pass ``vsp.SET_NONE`` to request no thin geometry explicitly.

    Returns
    -------
    str
        OpenVSP result identifier returned by ``VSPAEROSweep``.

    Raises
    ------
    ValueError
        If ``ncpu`` is less than one, ``wake_num_iter`` is outside the
        supported range 3 through 255, or ``wake_num_nodes`` is negative.
    RuntimeError
        If the current OpenVSP build does not expose an explicitly requested
        execution-setting input.

    Notes
    -----
    This function performs one ``VSPAEROComputeGeometry`` call followed by one
    native VSPAERO sweep. It does not perform elevator trimming.
    """
    if verbose:
        print('\n-> Calculate alpha & mach sweep analysis\n')

    # //==== Analysis: VSPAero Compute Geometry to Create Vortex Lattice DegenGeom File ====//
    # Use the thick/thin selection saved in the loaded model unless the caller
    # explicitly overrides it.  ComputeGeometry and Sweep always receive the
    # same effective selection.
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)
    compgeom_inputs = set(vsp.GetAnalysisInputNames(compgeom_name))
    if 'GeomSet' not in compgeom_inputs:
        raise RuntimeError(
            "VSPAEROComputeGeometry does not expose the 'GeomSet' input."
        )

    saved_thick_set = int(vsp.GetIntAnalysisInput(compgeom_name, 'GeomSet')[0])
    saved_thin_set = (
        int(vsp.GetIntAnalysisInput(compgeom_name, 'ThinGeomSet')[0])
        if 'ThinGeomSet' in compgeom_inputs
        else None
    )
    if thin_geom_set is not None and 'ThinGeomSet' not in compgeom_inputs:
        raise RuntimeError(
            "VSPAEROComputeGeometry does not expose the 'ThinGeomSet' input."
        )

    thick_set = saved_thick_set if thick_geom_set is None else int(thick_geom_set)
    thin_set = saved_thin_set if thin_geom_set is None else int(thin_geom_set)
    set_none = int(getattr(vsp, 'SET_NONE', -1))
    if thick_set == set_none and (thin_set is None or thin_set == set_none):
        raise RuntimeError('The requested VSPAERO settings select no geometry.')

    vsp.SetIntAnalysisInput(compgeom_name, 'GeomSet', [thick_set], 0)
    if thin_set is not None:
        vsp.SetIntAnalysisInput(compgeom_name, 'ThinGeomSet', [thin_set], 0)

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

    if ncpu is not None:
        ncpu = int(ncpu)
        if ncpu < 1:
            raise ValueError('ncpu must be a positive integer.')
        if 'NCPU' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'NCPU' input."
            )
        vsp.SetIntAnalysisInput(analysis_name, 'NCPU', [ncpu], 0)

    if wake_num_iter is not None:
        wake_num_iter = int(wake_num_iter)
        if not 3 <= wake_num_iter <= 255:
            raise ValueError('wake_num_iter must be between 3 and 255.')
        if 'WakeNumIter' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'WakeNumIter' input."
            )
        vsp.SetIntAnalysisInput(
            analysis_name,
            'WakeNumIter',
            [wake_num_iter],
            0,
        )

    if wake_num_nodes is not None:
        wake_num_nodes = int(wake_num_nodes)
        if wake_num_nodes < 0:
            raise ValueError('wake_num_nodes must be zero or greater.')
        if 'NumWakeNodes' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'NumWakeNodes' input."
            )
        vsp.SetIntAnalysisInput(
            analysis_name,
            'NumWakeNodes',
            [wake_num_nodes],
            0,
        )

    if fixed_wake_flag is not None:
        if 'FixedWakeFlag' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'FixedWakeFlag' input."
            )
        vsp.SetIntAnalysisInput(
            analysis_name,
            'FixedWakeFlag',
            [1 if fixed_wake_flag else 0],
            0,
        )

    # Keep the sweep on the exact same effective thick/thin selection used to
    # generate the VSPAERO geometry.  Reference dimensions and reference-wing
    # choice remain the values saved in the model.
    set_analysis_input_if_available(
        vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput,
        'GeomSet', [thick_set],
    )
    if thin_set is not None:
        set_analysis_input_if_available(
            vsp, analysis_name, analysis_inputs, vsp.SetIntAnalysisInput,
            'ThinGeomSet', [thin_set],
        )

    effective_ncpu = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'NCPU')[0])
        if 'NCPU' in analysis_inputs
        else None
    )
    effective_wake_num_iter = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'WakeNumIter')[0])
        if 'WakeNumIter' in analysis_inputs
        else None
    )
    effective_wake_num_nodes = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'NumWakeNodes')[0])
        if 'NumWakeNodes' in analysis_inputs
        else None
    )
    effective_fixed_wake_flag = (
        bool(vsp.GetIntAnalysisInput(analysis_name, 'FixedWakeFlag')[0])
        if 'FixedWakeFlag' in analysis_inputs
        else None
    )
    if wake_num_iter is not None and effective_fixed_wake_flag:
        warnings.warn(
            'wake_num_iter was specified while FixedWakeFlag is enabled; '
            'the fixed-wake setting takes precedence.',
            RuntimeWarning,
            stacklevel=2,
        )

    if 'UnsteadyType' in analysis_inputs:
        vsp.SetIntAnalysisInput(analysis_name, 'UnsteadyType', [vsp.STABILITY_OFF], 0)

    if verbose:
        print(
            ' VSPAERO run : '
            f'NCPU={effective_ncpu}, '
            f'WakeNumIter={effective_wake_num_iter}, '
            f'NumWakeNodes={effective_wake_num_nodes}, '
            f'FixedWakeFlag={effective_fixed_wake_flag}'
        )

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

def _trim_elevator_at_fixed_alpha(
    evaluate_cmy,
    initial_elevator_deg,
    *,
    previous_cmy_slope=None,
    elevator_bounds=(-20.0, 20.0),
    cmy_tolerance=5e-4,
    max_trim_evaluations=30,
    trim_step_deg=1.0,
):
    r"""Find the elevator deflection that satisfies the pitch-trim residual.

    The normal path uses the trimmed elevator deflection and elevator
    effectiveness inherited from the preceding case.  Every VSPAERO result is
    checked immediately, and the calculation stops at the first evaluated
    point satisfying ``abs(CMytot) <= cmy_tolerance``.  If continuation and
    open secant corrections do not converge, the function establishes a
    sign-changing bracket and continues with a safeguarded secant method.

    Parameters
    ----------
    evaluate_cmy : callable
        Function called as ``evaluate_cmy(elevator_deg, evaluation_number)``.
        It must return ``(cmy, result, duration_seconds)``, where ``cmy`` is
        ``CMytot``, ``result`` is the one-point polar result, and
        ``duration_seconds`` is the VSPAERO execution duration.
    initial_elevator_deg : float
        Initial elevator-deflection estimate in degrees.  The calling workflow
        should normally pass the trimmed deflection from the preceding case.
    previous_cmy_slope : float or None, optional
        Elevator effectiveness inherited from the preceding case,
        :math:`\partial C_{My}/\partial\delta_e`, in ``1/deg``.  A missing,
        non-finite, or nearly zero value disables the modified-Newton step.
    elevator_bounds : tuple of float, optional
        Inclusive lower and upper elevator-deflection limits in degrees.
    cmy_tolerance : float, optional
        Accepted absolute ``CMytot`` residual.
    max_trim_evaluations : int, optional
        Maximum number of unique calls to ``evaluate_cmy`` for one trim case.
    trim_step_deg : float, optional
        Elevator increment in degrees used for the first local-slope probe and
        for bracket expansion.  The inherited-slope Newton correction is
        limited to twice this value so a stale slope cannot cause a large jump.

    Returns
    -------
    result : pandas.DataFrame
        One-point VSPAERO polar at the accepted elevator deflection.
    elevator_deg : float
        Accepted elevator deflection in degrees.
    cmy_residual : float
        ``CMytot`` at ``elevator_deg``.
    cmy_slope : float or None
        Estimated local :math:`\partial C_{My}/\partial\delta_e` in ``1/deg``.
    evaluation_count : int
        Number of unique calls made to ``evaluate_cmy``.
    total_duration : float
        Sum of the reported VSPAERO durations in seconds.

    Raises
    ------
    ValueError
        If the bounds, tolerance, evaluation limit, or probe step are invalid.
    RuntimeError
        If VSPAERO returns a non-finite pitching moment, the evaluation limit
        is reached, no sign-changing bracket can be established inside the
        elevator limits, or a narrow bracket does not contain an evaluated
        point satisfying ``cmy_tolerance``.

    Notes
    -----
    Elevator deflection is one scalar trim variable.  Bracket expansion tests
    lower and higher values of that same variable; it does not represent
    separate left and right elevator controls.

    For nearly linear elevator effectiveness, the first case normally needs
    three evaluations: the initial point, one slope probe, and one secant root
    prediction.  A continued case normally needs one evaluation when the
    inherited deflection already passes, or two evaluations when one
    modified-Newton correction is required.
    """

    lower, upper = (float(value) for value in elevator_bounds)
    if lower >= upper:
        raise ValueError('elevator_bounds must satisfy lower < upper.')
    if cmy_tolerance <= 0.0:
        raise ValueError('cmy_tolerance must be positive.')
    if int(max_trim_evaluations) < 1:
        raise ValueError('max_trim_evaluations must be positive.')
    if trim_step_deg <= 0.0:
        raise ValueError('trim_step_deg must be positive.')

    initial_elevator_deg = float(initial_elevator_deg)
    if not np.isfinite(initial_elevator_deg):
        raise RuntimeError('The initial elevator deflection is not finite.')
    initial_elevator_deg = min(max(initial_elevator_deg, lower), upper)

    slope_epsilon = 1.0e-10
    elevator_epsilon = 1.0e-9
    max_correction_deg = 2.0 * float(trim_step_deg)
    evaluation_limit = int(max_trim_evaluations)

    inherited_slope = None
    if previous_cmy_slope is not None:
        candidate = float(previous_cmy_slope)
        if np.isfinite(candidate) and abs(candidate) > slope_epsilon:
            inherited_slope = candidate

    cache = {}
    evaluation_count = 0
    total_duration = 0.0

    print('trim: ', end='', flush=True)

    def evaluate(elevator_deg):
        """Run and cache one VSPAERO point."""
        nonlocal evaluation_count, total_duration

        elevator_deg = min(max(float(elevator_deg), lower), upper)
        key = round(elevator_deg, 10)
        if key in cache:
            return cache[key]
        if evaluation_count >= evaluation_limit:
            raise RuntimeError(
                'Pitch-trim evaluation limit reached. '
                f'max_trim_evaluations={evaluation_limit}'
            )

        cmy, result, duration = evaluate_cmy(
            elevator_deg,
            evaluation_count + 1,
        )
        cmy = float(cmy)
        if not np.isfinite(cmy):
            raise RuntimeError(
                f'VSPAERO returned a non-finite CMytot value: {cmy}'
            )

        evaluation_count += 1
        total_duration += float(duration)
        cache[key] = {
            'elevator_deg': elevator_deg,
            'cmy': cmy,
            'result': result,
        }
        print('*', end='', flush=True)
        return cache[key]

    def slope_between(first, second):
        delta_elevator = second['elevator_deg'] - first['elevator_deg']
        if abs(delta_elevator) <= elevator_epsilon:
            return None
        slope = (second['cmy'] - first['cmy']) / delta_elevator
        if not np.isfinite(slope) or abs(slope) <= slope_epsilon:
            return None
        return float(slope)

    def find_bracket():
        entries = sorted(cache.values(), key=lambda entry: entry['elevator_deg'])
        for lower_entry, upper_entry in zip(entries, entries[1:]):
            if lower_entry['cmy'] * upper_entry['cmy'] < 0.0:
                return lower_entry, upper_entry
        return None

    def finish(final, fallback_slope=None):
        slope_candidates = []
        for other in cache.values():
            if other is final:
                continue
            slope = slope_between(other, final)
            if slope is None:
                continue
            distance = abs(other['elevator_deg'] - final['elevator_deg'])
            slope_candidates.append(
                (abs(distance - float(trim_step_deg)), slope)
            )

        final_slope = None
        if slope_candidates:
            final_slope = min(slope_candidates, key=lambda item: item[0])[1]
        elif fallback_slope is not None:
            candidate = float(fallback_slope)
            if np.isfinite(candidate) and abs(candidate) > slope_epsilon:
                final_slope = candidate

        print('', flush=True)
        return (
            final['result'].copy(),
            final['elevator_deg'],
            final['cmy'],
            final_slope,
            evaluation_count,
            total_duration,
        )

    first = evaluate(initial_elevator_deg)
    if abs(first['cmy']) <= cmy_tolerance:
        return finish(first, inherited_slope)

    # Establish a second current-case point.  A continued case first tries one
    # limited modified-Newton correction.  The first case instead measures the
    # local elevator effectiveness with one nearby deflection.
    second_elevator_deg = None
    if inherited_slope is not None:
        correction = -first['cmy'] / inherited_slope
        correction = min(max(correction, -max_correction_deg), max_correction_deg)
        candidate = min(max(first['elevator_deg'] + correction, lower), upper)
        if abs(candidate - first['elevator_deg']) > elevator_epsilon:
            second_elevator_deg = candidate

    if second_elevator_deg is None:
        higher = first['elevator_deg'] + float(trim_step_deg)
        lower_probe = first['elevator_deg'] - float(trim_step_deg)
        if higher <= upper:
            second_elevator_deg = higher
        elif lower_probe >= lower:
            second_elevator_deg = lower_probe
        else:
            raise RuntimeError(
                'No distinct elevator point is available inside '
                f'elevator_bounds={elevator_bounds}.'
            )

    second = evaluate(second_elevator_deg)
    latest_slope = slope_between(first, second) or inherited_slope
    if abs(second['cmy']) <= cmy_tolerance:
        return finish(second, latest_slope)

    # Nearly linear cases should finish after one or two current-case secant
    # corrections.  Each correction is limited so a poor local slope cannot
    # jump across most of the allowed elevator range.
    previous = first
    current = second
    for _ in range(2):
        current_slope = slope_between(previous, current)
        if current_slope is None:
            break
        latest_slope = current_slope

        correction = -current['cmy'] / current_slope
        candidate = min(max(current['elevator_deg'] + correction, lower), upper)

        bracket = find_bracket()
        if bracket is not None:
            bracket_lower, bracket_upper = bracket
            if not (
                bracket_lower['elevator_deg']
                < candidate
                < bracket_upper['elevator_deg']
            ):
                candidate = 0.5 * (
                    bracket_lower['elevator_deg']
                    + bracket_upper['elevator_deg']
                )

        if (
            not np.isfinite(candidate)
            or abs(candidate - current['elevator_deg']) <= elevator_epsilon
            or round(candidate, 10) in cache
        ):
            break

        predicted = evaluate(candidate)
        if abs(predicted['cmy']) <= cmy_tolerance:
            return finish(predicted, latest_slope)
        previous, current = current, predicted

    # If the open corrections did not converge, establish a sign-changing
    # bracket around the best evaluated point.  When a usable slope is known,
    # test the deflection direction predicted to reduce CMytot first.
    bracket = find_bracket()
    expansion_step = float(trim_step_deg)
    while bracket is None and evaluation_count < evaluation_limit:
        best = min(cache.values(), key=lambda entry: abs(entry['cmy']))
        direction = 1.0
        if latest_slope is not None and abs(latest_slope) > slope_epsilon:
            direction = 1.0 if -best['cmy'] / latest_slope >= 0.0 else -1.0

        candidate_offsets = (
            direction * expansion_step,
            -direction * expansion_step,
        )
        evaluated_new_point = False
        for offset in candidate_offsets:
            candidate = min(max(best['elevator_deg'] + offset, lower), upper)
            if (
                abs(candidate - best['elevator_deg']) <= elevator_epsilon
                or round(candidate, 10) in cache
            ):
                continue

            entry = evaluate(candidate)
            evaluated_new_point = True
            if abs(entry['cmy']) <= cmy_tolerance:
                return finish(entry, latest_slope)

            bracket = find_bracket()
            if bracket is not None or evaluation_count >= evaluation_limit:
                break

        if bracket is not None:
            break
        if not evaluated_new_point:
            break
        expansion_step *= 2.0

    if bracket is None:
        best = min(cache.values(), key=lambda entry: abs(entry['cmy']))
        sampled = ', '.join(
            f"{entry['elevator_deg']:.8g}:{entry['cmy']:.8g}"
            for entry in sorted(
                cache.values(),
                key=lambda entry: entry['elevator_deg'],
            )
        )
        print('', flush=True)
        raise RuntimeError(
            'No evaluated elevator deflections bracket CMytot=0 within '
            f'elevator_bounds={elevator_bounds}.\n'
            f"  best_elevator_deg={best['elevator_deg']}\n"
            f"  best_CMytot={best['cmy']}\n"
            f'  tolerance={cmy_tolerance}\n'
            f'  evaluations={evaluation_count}\n'
            f'  sampled elevator_deg:CMytot = {sampled}'
        )

    # Safe path.  Use a secant interpolation while it remains well inside the
    # bracket; otherwise use the midpoint.  Only an actually evaluated point
    # satisfying the CMytot residual is accepted as the trim result.
    bracket_lower, bracket_upper = bracket
    while evaluation_count < evaluation_limit:
        width = bracket_upper['elevator_deg'] - bracket_lower['elevator_deg']
        if width <= elevator_epsilon:
            break

        denominator = bracket_upper['cmy'] - bracket_lower['cmy']
        if abs(denominator) > slope_epsilon:
            candidate = (
                bracket_lower['elevator_deg'] * bracket_upper['cmy']
                - bracket_upper['elevator_deg'] * bracket_lower['cmy']
            ) / denominator
        else:
            candidate = np.nan

        fraction = (
            (candidate - bracket_lower['elevator_deg']) / width
            if np.isfinite(candidate)
            else np.nan
        )
        if not np.isfinite(fraction) or fraction <= 0.1 or fraction >= 0.9:
            candidate = 0.5 * (
                bracket_lower['elevator_deg']
                + bracket_upper['elevator_deg']
            )

        if round(float(candidate), 10) in cache:
            candidate = 0.5 * (
                bracket_lower['elevator_deg']
                + bracket_upper['elevator_deg']
            )
        if round(float(candidate), 10) in cache:
            break

        entry = evaluate(candidate)
        if abs(entry['cmy']) <= cmy_tolerance:
            bracket_slope = slope_between(bracket_lower, bracket_upper)
            return finish(entry, bracket_slope or latest_slope)

        if bracket_lower['cmy'] * entry['cmy'] < 0.0:
            bracket_upper = entry
        else:
            bracket_lower = entry

    best = min(cache.values(), key=lambda entry: abs(entry['cmy']))
    sampled = ', '.join(
        f"{entry['elevator_deg']:.8g}:{entry['cmy']:.8g}"
        for entry in sorted(
            cache.values(),
            key=lambda entry: entry['elevator_deg'],
        )
    )
    print('', flush=True)
    raise RuntimeError(
        'Pitch trim did not produce an evaluated point satisfying the '
        'CMytot residual tolerance.\n'
        f"  best_elevator_deg={best['elevator_deg']}\n"
        f"  best_CMytot={best['cmy']}\n"
        f'  tolerance={cmy_tolerance}\n'
        f'  evaluations={evaluation_count}\n'
        f"  final_bracket=({bracket_lower['elevator_deg']}, "
        f"{bracket_upper['elevator_deg']})\n"
        f'  sampled elevator_deg:CMytot = {sampled}'
    )

def vsp_trimmed_sweep(
    vsp,
    alpha_list,
    Weight,
    altitude=0,
    dT=0,
    *,
    elevator_bounds=(-20.0, 20.0),
    cmy_tolerance=5e-4,
    max_trim_evaluations=30,
    trim_step_deg=1.0,
    analysis_method=None,
    ncpu=None,
    wake_num_iter=None,
    fixed_wake_flag=None,
    verbose=1,
):
    r"""Calculate a fixed-alpha, pitch-trimmed glider polar.

    For each requested angle of attack, the function first obtains the lift
    coefficient used to determine flight speed, Mach number, and Reynolds
    number. It then varies the existing ``ELEVATOR_GROUP`` deflection until
    ``CMytot`` satisfies the requested tolerance. VSPAERO geometry is built
    once for the complete sweep.

    Parameters
    ----------
    vsp : module or object
        OpenVSP Python API object used to configure and execute analyses.
    alpha_list : sequence of float
        Angles of attack in degrees.
    Weight : float
        Aircraft mass in kilograms. The legacy parameter name is retained for
        API compatibility.
    altitude : float, optional
        Geopotential altitude in metres used by the ISA functions.
    dT : float, optional
        ISA temperature offset in kelvin.
    elevator_bounds : tuple of float, optional
        Inclusive lower and upper elevator-deflection limits in degrees.
    cmy_tolerance : float, optional
        Accepted absolute ``CMytot`` residual.
    max_trim_evaluations : int, optional
        Maximum VSPAERO evaluations permitted for each trim case.
    trim_step_deg : float, optional
        Elevator increment used to measure local effectiveness when no reliable
        slope is inherited, and as the initial fallback bracket spacing.
    analysis_method : int or None, optional
        Legacy ``VSPAEROComputeGeometry/AnalysisMethod`` override. It is used
        only when the loaded OpenVSP build does not expose separate thick and
        thin geometry-set inputs.
    ncpu : int or None, optional
        Value assigned to ``VSPAEROSweep/NCPU``. ``None`` leaves the OpenVSP
        default unchanged.
    wake_num_iter : int or None, optional
        Value assigned to ``VSPAEROSweep/WakeNumIter``. Valid explicit values
        are 3 through 255.
    fixed_wake_flag : bool or None, optional
        Value assigned to ``VSPAEROSweep/FixedWakeFlag``. When enabled, the
        fixed-wake setting takes precedence over ``wake_num_iter``.
    verbose : int or bool, optional
        ``0`` suppresses progress, ``1`` prints case summaries, and values of
        ``2`` or greater also print detailed OpenVSP analysis information.

    Returns
    -------
    pandas.DataFrame
        Concatenated trimmed-polar rows. Diagnostic columns include
        ``Elevator_deg``, ``CMyResidual``, ``CMyElevatorSlope_per_deg``,
        ``TrimFunctionEvaluations``, and execution durations.

    Raises
    ------
    ValueError
        If input sequences, VSPAERO reference values, or numerical trim
        settings are invalid.
    RuntimeError
        If the required reference wing or elevator group is unavailable,
        ``ELEVATOR_GROUP`` is disabled in ``Control Group Angles``, an
        OpenVSP analysis fails, or a trim root cannot be obtained.

    Notes
    -----
    Elevator deflection and :math:`\partial C_{My}/\partial\delta_e` are
    continued from one angle-of-attack case to the next. The inherited slope
    normally reduces later trim cases to one or two VSPAERO evaluations.
    """

    alpha_list = list(alpha_list)
    if not alpha_list:
        raise ValueError('alpha_list must not be empty.')

    error_mgr = vsp.ErrorMgrSingleton.getInstance()

    def pop_openvsp_errors():
        errors = []
        while error_mgr.GetNumTotalErrors() > 0:
            error = error_mgr.PopLastError()
            errors.append(error.GetErrorString())
        return errors

    # Resolve the existing elevator once. During trim only DeflectionAngle is
    # changed; the control-surface group and gains are not rebuilt.
    vsp.Update()
    wing_id = find_one_geom(vsp, DEFAULT_REF_GEOM_NAME)
    settings_id = find_vspaero_settings_container(vsp)
    elevator_group_names = [
        vsp.GetVSPAEROControlGroupName(group_index)
        for group_index in range(vsp.GetNumControlSurfaceGroups())
    ]
    elevator_group_matches = [
        group_index
        for group_index, group_name in enumerate(elevator_group_names)
        if group_name == 'ELEVATOR_GROUP'
    ]
    if len(elevator_group_matches) != 1:
        raise RuntimeError(
            "Pitch trim requires exactly one existing 'ELEVATOR_GROUP'. "
            f'found={len(elevator_group_matches)}'
        )
    elevator_group_index = elevator_group_matches[0]
    if not list(vsp.GetActiveCSNameVec(elevator_group_index)):
        raise RuntimeError("'ELEVATOR_GROUP' has no active control surfaces.")

    elevator_group_parm_group = (
        f'ControlSurfaceGroup_{elevator_group_index}'
    )
    elevator_active_parm_id = vsp.FindParm(
        settings_id,
        'ActiveFlag',
        elevator_group_parm_group,
    )
    if (
        not elevator_active_parm_id
        or str(elevator_active_parm_id).upper() == 'NONE'
    ):
        raise RuntimeError(
            "Could not find 'ActiveFlag' for the existing "
            "'ELEVATOR_GROUP'."
        )
    if not bool(vsp.GetParmVal(elevator_active_parm_id)):
        raise RuntimeError(
            "'ELEVATOR_GROUP' is disabled in VSPAERO "
            "'Control Group Angles'. Enable it before running pitch trim."
        )

    elevator_parm_id = vsp.FindParm(
        settings_id,
        'DeflectionAngle',
        elevator_group_parm_group,
    )
    if not elevator_parm_id or str(elevator_parm_id).upper() == 'NONE':
        raise RuntimeError(
            "Could not find 'DeflectionAngle' for the existing "
            "'ELEVATOR_GROUP'."
        )
    original_elevator_deg = float(vsp.GetParmVal(elevator_parm_id))

    # Build the mixed Thick/Thin geometry once for the complete trimmed sweep.
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)
    compgeom_inputs = set(vsp.GetAnalysisInputNames(compgeom_name))
    use_saved_thick_thin = {'GeomSet', 'ThinGeomSet'} <= compgeom_inputs
    if use_saved_thick_thin:
        thick_set = int(vsp.GetIntAnalysisInput(compgeom_name, 'GeomSet')[0])
        thin_set = int(vsp.GetIntAnalysisInput(compgeom_name, 'ThinGeomSet')[0])
    else:
        thick_set = 0
        thin_set = None
        if 'GeomSet' in compgeom_inputs:
            vsp.SetIntAnalysisInput(compgeom_name, 'GeomSet', [thick_set], 0)
        if analysis_method is not None and 'AnalysisMethod' in compgeom_inputs:
            vsp.SetIntAnalysisInput(
                compgeom_name,
                'AnalysisMethod',
                [int(analysis_method)],
                0,
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
            'VSPAEROComputeGeometry did not produce a usable '
            'VSPAERO_Geom result.\n'
            f'  result_id={compgeom_result_id!r}\n'
            f'  result_name={compgeom_result_name!r}\n'
            f'  vspgeom_files={vspgeom_files!r}\n'
            f'  OpenVSP errors={compgeom_errors or ["<none reported>"]}'
        )

    # Configure one reusable one-point polar analysis.
    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))

    if ncpu is not None:
        ncpu = int(ncpu)
        if ncpu < 1:
            raise ValueError('ncpu must be a positive integer.')
        if 'NCPU' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'NCPU' input."
            )
        vsp.SetIntAnalysisInput(analysis_name, 'NCPU', [ncpu], 0)

    if wake_num_iter is not None:
        wake_num_iter = int(wake_num_iter)
        if not 3 <= wake_num_iter <= 255:
            raise ValueError('wake_num_iter must be between 3 and 255.')
        if 'WakeNumIter' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'WakeNumIter' input."
            )
        vsp.SetIntAnalysisInput(
            analysis_name,
            'WakeNumIter',
            [wake_num_iter],
            0,
        )

    if fixed_wake_flag is not None:
        if 'FixedWakeFlag' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'FixedWakeFlag' input."
            )
        vsp.SetIntAnalysisInput(
            analysis_name,
            'FixedWakeFlag',
            [1 if fixed_wake_flag else 0],
            0,
        )

    effective_ncpu = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'NCPU')[0])
        if 'NCPU' in analysis_inputs
        else None
    )
    effective_wake_num_iter = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'WakeNumIter')[0])
        if 'WakeNumIter' in analysis_inputs
        else None
    )
    effective_fixed_wake_flag = (
        bool(vsp.GetIntAnalysisInput(analysis_name, 'FixedWakeFlag')[0])
        if 'FixedWakeFlag' in analysis_inputs
        else None
    )
    if wake_num_iter is not None and effective_fixed_wake_flag:
        warnings.warn(
            'wake_num_iter was specified while FixedWakeFlag is enabled; '
            'the fixed-wake setting takes precedence.',
            RuntimeWarning,
            stacklevel=2,
        )

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
            f'{analysis_name} is missing required inputs: {missing_inputs}.'
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

    Sref = float(vsp.GetDoubleAnalysisInput(analysis_name, 'Sref')[0])
    cref = float(vsp.GetDoubleAnalysisInput(analysis_name, 'cref')[0])
    if Sref <= 0.0 or cref <= 0.0:
        raise ValueError(
            'VSPAERO reference values must be positive. '
            f'Sref={Sref}, cref={cref}'
        )

    def run_polar_point(
        alpha,
        mach_value,
        reynolds_value,
        elevator_deg,
        evaluation_label,
    ):
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
        vsp.SetIntAnalysisInput(analysis_name, 'GroundEffectToggle', [0], 0)
        vsp.SetDoubleAnalysisInput(analysis_name, 'GroundEffect', [0.0], 0)
        vsp.Update()

        if verbose and int(verbose) >= 2:
            print(
                f'   evaluation {evaluation_label}: '
                f'alpha={float(alpha):.8g}, Mach={float(mach_value):.8g}, '
                f'Re={float(reynolds_value):.8g}, '
                f'elevator={float(elevator_deg):.8g}',
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

        polar = pd.DataFrame()
        for result_name in ('VSPAERO_Polar', 'VSPAERO Polar'):
            polar = (
                results_dataframe(vsp, wrapper_result_id, result_name)
                if wrapper_result_id
                else pd.DataFrame()
            )
            if not polar.empty:
                break

        if polar.empty:
            raise RuntimeError(
                'VSPAERO returned no polar result.\n'
                f'  alpha={alpha}, Mach={mach_value}, Re={reynolds_value}\n'
                f'  elevator_deg={elevator_deg}\n'
                f'  child_results={child_result_names!r}\n'
                f'  OpenVSP errors={api_errors or ["<none reported>"]}'
            )
        if len(polar.index) != 1:
            raise RuntimeError(
                'A one-point VSPAERO analysis returned an unexpected number '
                f'of rows: {len(polar.index)}'
            )

        if 'CLtot' in polar.columns and 'CL' not in polar.columns:
            polar['CL'] = polar['CLtot']
        if 'CL' in polar.columns and 'CLtot' not in polar.columns:
            polar['CLtot'] = polar['CL']
        if 'CMytot' in polar.columns and 'CMy' not in polar.columns:
            polar['CMy'] = polar['CMytot']
        if 'CMy' in polar.columns and 'CMytot' not in polar.columns:
            polar['CMytot'] = polar['CMy']

        duration = analysis_duration_seconds(vsp, wrapper_result_id)
        duration_seconds = elapsed if duration is None else duration

        if verbose and int(verbose) >= 2:
            vsp.PrintResults(wrapper_result_id)
            print(polar.to_string(index=False), flush=True)

        return polar.copy(), float(duration_seconds)

    g = get_gravity()
    density = get_density(altitude=altitude, dT=dT)
    previous_trim_deg = original_elevator_deg
    previous_cmy_slope = None
    frames = []

    if verbose:
        print('\n-> VSPAERO trimmed sweep', flush=True)
        print(
            f' Geometry ready : {vspgeom_files[0]} '
            f'({compgeom_elapsed:.2f} s)',
            flush=True,
        )
        print(f' Cases          : {len(alpha_list)}', flush=True)
        print(
            ' VSPAERO run    : '
            f'NCPU={effective_ncpu}, '
            f'WakeNumIter={effective_wake_num_iter}, '
            f'FixedWakeFlag={effective_fixed_wake_flag}',
            flush=True,
        )

    try:
        for case_index, alpha in enumerate(alpha_list, start=1):
            if verbose:
                print(
                    f' [{case_index:>3}/{len(alpha_list)}] '
                    f'alpha={float(alpha):.8g} deg',
                    flush=True,
                )

            # The speed-defining point always uses the elevator state that was
            # present when this function was called.
            initial_polar, initial_duration = run_polar_point(
                alpha,
                0.0,
                1.0e6,
                original_elevator_deg,
                'speed condition',
            )
            if 'CL' not in initial_polar.columns:
                raise RuntimeError(
                    "VSPAERO polar does not contain 'CL' or 'CLtot'."
                )
            initial_cl = float(initial_polar['CL'].iloc[0])
            if not np.isfinite(initial_cl) or initial_cl <= 0.0:
                raise RuntimeError(
                    'A positive finite CL is required to calculate velocity. '
                    f'alpha={alpha}, CL={initial_cl}'
                )

            velocity = np.sqrt(
                (2.0 * float(Weight) * g)
                / (density * Sref * initial_cl)
            )
            mach_value = velocity_to_mach(
                velocity=velocity,
                dT=dT,
                altitude=altitude,
            )
            reynolds_value = velocity_to_reynolds(
                velocity=velocity,
                length=cref,
                altitude=altitude,
                dT=dT,
            )

            def evaluate_cmy(elevator_deg, evaluation_number):
                result, duration = run_polar_point(
                    alpha,
                    mach_value,
                    reynolds_value,
                    elevator_deg,
                    f'trim {evaluation_number}',
                )
                if 'CMytot' not in result.columns:
                    raise RuntimeError(
                        "VSPAERO polar does not contain 'CMytot' or 'CMy'."
                    )
                return float(result['CMytot'].iloc[0]), result, duration

            (
                result,
                elevator_deg,
                cmy_residual,
                cmy_slope,
                trim_evaluations,
                trim_duration,
            ) = _trim_elevator_at_fixed_alpha(
                evaluate_cmy,
                previous_trim_deg,
                previous_cmy_slope=previous_cmy_slope,
                elevator_bounds=elevator_bounds,
                cmy_tolerance=cmy_tolerance,
                max_trim_evaluations=max_trim_evaluations,
                trim_step_deg=trim_step_deg,
            )
            previous_trim_deg = elevator_deg
            previous_cmy_slope = cmy_slope

            result.insert(0, 'InputAlpha_deg', float(alpha))
            result.insert(1, 'InputMach', float(mach_value))
            result.insert(2, 'InputReCref', float(reynolds_value))
            result.insert(3, 'Elevator_deg', float(elevator_deg))
            result.insert(4, 'CMyResidual', float(cmy_residual))
            result.insert(
                5,
                'CMyElevatorSlope_per_deg',
                float(cmy_slope) if cmy_slope is not None else np.nan,
            )
            result.insert(6, 'TrimFunctionEvaluations', int(trim_evaluations))
            result.insert(7, 'NCPU', effective_ncpu)
            result.insert(8, 'WakeNumIter', effective_wake_num_iter)
            result.insert(9, 'FixedWakeFlag', effective_fixed_wake_flag)
            result.insert(10, 'SpeedConditionDuration_s', float(initial_duration))
            result.insert(11, 'TrimDuration_s', float(trim_duration))
            result['de'] = float(elevator_deg)
            frames.append(result)

            if verbose:
                print(
                    f'   complete: elevator={elevator_deg:.6g} deg, '
                    f'CMytot={cmy_residual:.3e}, '
                    f'dCMy/dde={cmy_slope if cmy_slope is not None else np.nan:.6g} 1/deg, '
                    f'evaluations={trim_evaluations}, '
                    f'V={velocity:.6g}, Mach={mach_value:.6g}, '
                    f'Re={reynolds_value:.6g}',
                    flush=True,
                )
    finally:
        vsp.SetParmVal(elevator_parm_id, original_elevator_deg)
        vsp.Update()

    trimmed_polar = pd.concat(frames, ignore_index=True)
    trimmed_polar['gamma'] = np.arctan(1.0 / trimmed_polar['L_D'])
    trimmed_polar['Velocity'] = np.sqrt(
        (2.0 * float(Weight) * g)
        / (
            density
            * Sref
            * trimmed_polar['CL']
            * np.cos(trimmed_polar['gamma'])
        )
    )
    trimmed_polar['Vx'] = (
        trimmed_polar['Velocity']
        * np.cos(trimmed_polar['gamma'])
        * 3.6
    )
    trimmed_polar['vz'] = (
        trimmed_polar['Velocity']
        * np.sin(trimmed_polar['gamma'])
    )
    return trimmed_polar

def vsp_stability_derivatives(
    vsp3_path,
    *,
    alpha=2.0,
    mach=0.1,
    reynolds=4.4e6,
    stability_type=vsp.STABILITY_DEFAULT,
    verbose=1,
    vspaero_verbose=0,
    ncpu=None,
    wake_num_iter=None,
    wake_num_nodes=None,
    fixed_wake_flag=None,
    redirect_file='',
    stop_before_run=False,
):
    r"""Run a steady VSPAERO stability-derivative analysis.

    The function loads the ``.vsp3`` model, uses the model's saved VSPAERO
    thick/thin geometry selection and reference settings, executes
    ``VSPAEROComputeGeometry``, then runs ``VSPAEROSweep`` in either the normal
    steady finite-difference or steady adjoint stability mode.

    No particular Geom names, Set names, or control-surface group names are
    required.  A model may therefore use only thin geometry, only thick
    geometry, or a mixed thick/thin representation as configured in the
    ``.vsp3`` file.

    Parameters
    ----------
    vsp3_path : str or os.PathLike
        Path to the OpenVSP ``.vsp3`` model.
    alpha : float, optional
        Analysis angle of attack in degrees.
    mach : float, optional
        Analysis Mach number.
    reynolds : float, optional
        Reynolds number based on the saved VSPAERO reference chord.
    stability_type : int, optional
        Steady stability-derivative method. Use ``vsp.STABILITY_DEFAULT`` for
        the normal finite-difference calculation or ``vsp.STABILITY_ADJOINT``
        for the steady adjoint calculation.
    verbose : int or bool, optional
        ``0`` suppresses progress, ``1`` prints major progress and the final
        status, and values of ``2`` or greater also print timing details.
    vspaero_verbose : int or bool, optional
        Controls OpenVSP and VSPAERO internal output. A false value suppresses
        solver output where supported.
    ncpu : int or None, optional
        Value assigned to ``VSPAEROSweep/NCPU``. ``None`` leaves the OpenVSP
        default unchanged.
    wake_num_iter : int or None, optional
        Value assigned to ``VSPAEROSweep/WakeNumIter``. Valid explicit values
        are 3 through 255.
    wake_num_nodes : int or None, optional
        Value assigned to ``VSPAEROSweep/NumWakeNodes``. ``None`` leaves the
        loaded/default value unchanged.
    fixed_wake_flag : bool or None, optional
        Value assigned to ``VSPAEROSweep/FixedWakeFlag``. When enabled, the
        fixed-wake setting takes precedence over ``wake_num_iter``.
    redirect_file : str or None, optional
        Value assigned to ``VSPAEROSweep/RedirectFile``. Use an empty string to
        suppress redirected output, ``"stdout"`` to display it, or a path to
        capture it. ``None`` leaves the setting unchanged.
    stop_before_run : bool, optional
        If ``True``, request generation of VSPAERO input files without running
        the solver. No ``VSPAERO_Stab`` result is then expected.

    Returns
    -------
    dict
        Diagnostic report containing ``passed``, ``errors``, ``warnings``,
        ``infos``, OpenVSP result identifiers, the stability-derivative
        ``pandas.DataFrame``, timing data, and the effective VSPAERO settings.
    """
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
            'stability_type': stability_type,
            'ncpu': ncpu,
            'wake_num_iter': wake_num_iter,
            'wake_num_nodes': wake_num_nodes,
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

    stability_types = {
        int(vsp.STABILITY_DEFAULT): ('STABILITY_DEFAULT', '.stab'),
        int(vsp.STABILITY_ADJOINT): ('STABILITY_ADJOINT', '.adjoint.stab'),
    }
    stability_type = int(stability_type)
    if stability_type not in stability_types:
        add(
            'errors',
            'INVALID_STABILITY_TYPE',
            'stability_type must be vsp.STABILITY_DEFAULT or '
            'vsp.STABILITY_ADJOINT.',
            {'stability_type': stability_type},
        )
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    stability_name, stab_file_suffix = stability_types[stability_type]
    result['vspaero_settings'].update({
        'stability_type': stability_type,
        'stability_name': stability_name,
        'stab_file_suffix': stab_file_suffix,
    })

    vsp3_path = os.fspath(vsp3_path)
    if not os.path.isfile(vsp3_path):
        add(
            'errors',
            'FILE_NOT_FOUND',
            'The specified .vsp3 file was not found.',
            {'vsp3_path': vsp3_path},
        )
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    vprint(1, f'\n-> Calculate VSPAERO stability derivatives: {vsp3_path}')

    read_start = time.perf_counter()
    try:
        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(vsp3_path)
        vsp.Update()
    except Exception as err:
        add(
            'errors',
            'READ_VSP3_FAILED',
            'OpenVSP failed to read the .vsp3 file.',
            {'error': repr(err)},
        )
        result['timing']['read_vsp3_elapsed_s'] = time.perf_counter() - read_start
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result
    result['timing']['read_vsp3_elapsed_s'] = time.perf_counter() - read_start

    # Use the geometry selection saved in the model.  SET_NONE is a valid
    # selection for either side; only selecting no geometry on both sides is
    # rejected before VSPAEROComputeGeometry.
    compgeom_name = 'VSPAEROComputeGeometry'
    vsp.SetAnalysisInputDefaults(compgeom_name)
    compgeom_inputs = set(vsp.GetAnalysisInputNames(compgeom_name))
    if 'GeomSet' not in compgeom_inputs:
        add(
            'errors',
            'MISSING_GEOMSET_INPUT',
            "VSPAEROComputeGeometry does not expose the 'GeomSet' input.",
        )
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    thick_set = int(vsp.GetIntAnalysisInput(compgeom_name, 'GeomSet')[0])
    thin_set = (
        int(vsp.GetIntAnalysisInput(compgeom_name, 'ThinGeomSet')[0])
        if 'ThinGeomSet' in compgeom_inputs
        else None
    )
    set_none = int(getattr(vsp, 'SET_NONE', -1))
    if thick_set == set_none and (thin_set is None or thin_set == set_none):
        add(
            'errors',
            'NO_AERO_GEOMETRY',
            'The saved VSPAERO settings select no thick or thin geometry.',
            {'GeomSet': thick_set, 'ThinGeomSet': thin_set},
        )
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    result['vspaero_settings'].update({
        'effective_geom_set': thick_set,
        'effective_thin_geom_set': thin_set,
    })

    vprint(
        1,
        ' Executing VSPAEROComputeGeometry '
        f'(GeomSet={thick_set}, ThinGeomSet={thin_set})...',
    )
    compgeom_start = time.perf_counter()
    try:
        if vspaero_verbose:
            compgeom_result_id = vsp.ExecAnalysis(compgeom_name)
        else:
            with suppress_stdout():
                compgeom_result_id = vsp.ExecAnalysis(compgeom_name)
    except Exception as err:
        add(
            'errors',
            'VSPAERO_COMPUTE_GEOMETRY_FAILED',
            'VSPAEROComputeGeometry failed.',
            {'error': repr(err)},
        )
        result['timing']['compute_geometry_elapsed_s'] = time.perf_counter() - compgeom_start
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    result['compute_geometry_result_id'] = compgeom_result_id
    result['timing']['compute_geometry_elapsed_s'] = time.perf_counter() - compgeom_start
    result['timing']['compute_geometry_analysis_duration_s'] = analysis_duration_seconds(vsp, compgeom_result_id)

    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))

    # Keep VSPAEROSweep on the exact same saved thick/thin selection used to
    # generate the VSPAERO geometry.  Reference dimensions and reference-wing
    # choice are intentionally left at the model's saved values.
    set_analysis_input_if_available(
        vsp,
        analysis_name,
        analysis_inputs,
        vsp.SetIntAnalysisInput,
        'GeomSet',
        [thick_set],
    )
    if thin_set is not None:
        set_analysis_input_if_available(
            vsp,
            analysis_name,
            analysis_inputs,
            vsp.SetIntAnalysisInput,
            'ThinGeomSet',
            [thin_set],
        )

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
        vsp.SetIntAnalysisInput(analysis_name, 'UnsteadyType', [stability_type], 0)
    else:
        add(
            'errors',
            'MISSING_UNSTEADYTYPE_INPUT',
            "VSPAEROSweep does not expose the 'UnsteadyType' input in this OpenVSP version.",
        )
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    if ncpu is not None:
        if int(ncpu) < 1:
            add('errors', 'INVALID_NCPU', 'ncpu must be a positive integer.', {'ncpu': ncpu})
        elif not set_analysis_input_if_available(
            vsp,
            analysis_name,
            analysis_inputs,
            vsp.SetIntAnalysisInput,
            'NCPU',
            [int(ncpu)],
        ):
            add('warnings', 'MISSING_NCPU_INPUT', "VSPAEROSweep does not expose the 'NCPU' input.")

    if wake_num_iter is not None:
        wake_num_iter = int(wake_num_iter)
        if not 3 <= wake_num_iter <= 255:
            add(
                'errors',
                'INVALID_WAKE_NUM_ITER',
                'wake_num_iter must be between 3 and 255.',
                {'wake_num_iter': wake_num_iter},
            )
        elif not set_analysis_input_if_available(
            vsp,
            analysis_name,
            analysis_inputs,
            vsp.SetIntAnalysisInput,
            'WakeNumIter',
            [wake_num_iter],
        ):
            add('warnings', 'MISSING_WAKE_NUM_ITER_INPUT', "VSPAEROSweep does not expose the 'WakeNumIter' input.")

    if wake_num_nodes is not None:
        wake_num_nodes = int(wake_num_nodes)
        if wake_num_nodes < 0:
            add(
                'errors',
                'INVALID_WAKE_NUM_NODES',
                'wake_num_nodes must be zero or greater.',
                {'wake_num_nodes': wake_num_nodes},
            )
        elif not set_analysis_input_if_available(
            vsp,
            analysis_name,
            analysis_inputs,
            vsp.SetIntAnalysisInput,
            'NumWakeNodes',
            [wake_num_nodes],
        ):
            add('warnings', 'MISSING_WAKE_NUM_NODES_INPUT', "VSPAEROSweep does not expose the 'NumWakeNodes' input.")

    if fixed_wake_flag is not None:
        if not set_analysis_input_if_available(
            vsp,
            analysis_name,
            analysis_inputs,
            vsp.SetIntAnalysisInput,
            'FixedWakeFlag',
            [1 if fixed_wake_flag else 0],
        ):
            add('warnings', 'MISSING_FIXED_WAKE_FLAG_INPUT', "VSPAEROSweep does not expose the 'FixedWakeFlag' input.")

    if redirect_file is not None:
        if not set_analysis_input_if_available(
            vsp,
            analysis_name,
            analysis_inputs,
            vsp.SetStringAnalysisInput,
            'RedirectFile',
            [str(redirect_file)],
        ):
            add('warnings', 'MISSING_REDIRECT_FILE_INPUT', "VSPAEROSweep does not expose the 'RedirectFile' input.")

    if stop_before_run:
        if not set_analysis_input_if_available(
            vsp,
            analysis_name,
            analysis_inputs,
            vsp.SetIntAnalysisInput,
            'StopBeforeRun',
            [1],
        ):
            add('warnings', 'MISSING_STOP_BEFORE_RUN_INPUT', "VSPAEROSweep does not expose the 'StopBeforeRun' input.")

    if result['errors']:
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    effective_stability_type = int(vsp.GetIntAnalysisInput(analysis_name, 'UnsteadyType')[0])
    effective_stability_name, effective_stab_file_suffix = stability_types.get(
        effective_stability_type,
        (f'UNKNOWN_{effective_stability_type}', ''),
    )
    effective_ncpu = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'NCPU')[0])
        if 'NCPU' in analysis_inputs
        else None
    )
    effective_wake_num_iter = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'WakeNumIter')[0])
        if 'WakeNumIter' in analysis_inputs
        else None
    )
    effective_wake_num_nodes = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'NumWakeNodes')[0])
        if 'NumWakeNodes' in analysis_inputs
        else None
    )
    effective_fixed_wake_flag = (
        bool(vsp.GetIntAnalysisInput(analysis_name, 'FixedWakeFlag')[0])
        if 'FixedWakeFlag' in analysis_inputs
        else None
    )
    effective_ref_flag = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'RefFlag')[0])
        if 'RefFlag' in analysis_inputs
        else None
    )
    try:
        effective_ref_geom_id = vsp.GetVSPAERORefWingID()
    except Exception:
        effective_ref_geom_id = ''

    result['vspaero_settings'].update({
        'effective_stability_type': effective_stability_type,
        'effective_stability_name': effective_stability_name,
        'effective_stab_file_suffix': effective_stab_file_suffix,
        'effective_ncpu': effective_ncpu,
        'effective_wake_num_iter': effective_wake_num_iter,
        'effective_wake_num_nodes': effective_wake_num_nodes,
        'effective_fixed_wake_flag': effective_fixed_wake_flag,
        'effective_ref_flag': effective_ref_flag,
        'effective_ref_geom_id': effective_ref_geom_id,
    })

    if wake_num_iter is not None and effective_fixed_wake_flag:
        add(
            'warnings',
            'WAKE_ITERATIONS_DISABLED_BY_FIXED_WAKE',
            'wake_num_iter was specified while FixedWakeFlag is enabled; '
            'the fixed-wake setting takes precedence.',
        )

    vprint(
        1,
        f' Executing VSPAEROSweep with {stability_name} '
        f'(NCPU={effective_ncpu}, '
        f'WakeNumIter={effective_wake_num_iter}, '
        f'NumWakeNodes={effective_wake_num_nodes}, '
        f'FixedWakeFlag={effective_fixed_wake_flag})',
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
        add(
            'errors',
            'MISSING_VSPAERO_STAB',
            f'VSPAERO_Stab was not found after {stability_name} analysis.',
            {'result_names': result['result_names']},
        )
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
        add(
            'errors',
            'EMPTY_VSPAERO_STAB',
            'VSPAERO_Stab was found, but no numeric derivative data was readable.',
        )
        result['timing']['result_extraction_elapsed_s'] = time.perf_counter() - extraction_start
        result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
        return result

    result['timing']['result_extraction_elapsed_s'] = time.perf_counter() - extraction_start
    result['timing']['total_elapsed_s'] = time.perf_counter() - total_start
    result['passed'] = len(result['errors']) == 0

    if verbose:
        status = 'PASSED' if result['passed'] else 'FAILED'
        print(
            f' Stability derivative calculation {status}: '
            f'{len(result["errors"])} error(s), '
            f'{len(result["warnings"])} warning(s)'
        )
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
    trim_step_deg=1.0,
    ncpu=None,
    wake_num_iter=None,
    fixed_wake_flag=None,
    verbose=1,
):
    r"""Run VSPAERO ground-effect polar or pitch-trimmed cases.

    A ground-effect-off reference is added automatically for every requested
    angle, Mach number, and Reynolds number. Each value in ``height`` is passed
    to ``VSPAEROSweep/GroundEffect`` as the height of the VSPAERO CG / moment
    reference point above the ground. With ``trimmed=True``, the flow condition
    and height remain fixed while the elevator is varied until ``CMytot`` is
    within tolerance.

    Parameters
    ----------
    vsp : module or object
        OpenVSP Python API object used to configure and execute analyses.
    alpha_deg : sequence of float
        Angles of attack in degrees.
    mach : sequence of float
        Mach numbers.
    reynolds : sequence of float
        Reynolds numbers based on VSPAERO reference chord.
    height : sequence of float
        Positive CG / moment-reference heights above the ground in model length
        units. Do not add a large sentinel value for the out-of-ground-effect
        case; that reference is generated automatically.
    analysis_method : int or None, optional
        Legacy ``VSPAEROComputeGeometry/AnalysisMethod`` override. It is used
        only when separate thick and thin geometry-set inputs are unavailable.
    trimmed : bool, optional
        If ``True``, solve ``CMytot = 0`` by varying the existing elevator
        control group at each fixed aerodynamic condition.
    elevator_bounds : tuple of float, optional
        Inclusive lower and upper elevator-deflection limits in degrees.
    cmy_tolerance : float, optional
        Accepted absolute ``CMytot`` residual.
    max_trim_evaluations : int, optional
        Maximum VSPAERO evaluations permitted for each trim case.
    trim_step_deg : float, optional
        Elevator increment used to measure local effectiveness when no reliable
        slope is inherited, and as the initial fallback bracket spacing.
    ncpu : int or None, optional
        Value assigned to ``VSPAEROSweep/NCPU``. ``None`` leaves the OpenVSP
        default unchanged.
    wake_num_iter : int or None, optional
        Value assigned to ``VSPAEROSweep/WakeNumIter``. Valid explicit values
        are 3 through 255.
    fixed_wake_flag : bool or None, optional
        Value assigned to ``VSPAEROSweep/FixedWakeFlag``. When enabled, the
        fixed-wake setting takes precedence over ``wake_num_iter``.
    verbose : int or bool, optional
        ``0`` suppresses progress, ``1`` prints case summaries, and values of
        ``2`` or greater also print detailed OpenVSP analysis information.

    Returns
    -------
    pandas.DataFrame
        Concatenated VSPAERO polar rows. Trimmed results include
        ``Elevator_deg``, ``CMyResidual``, ``CMyElevatorSlope_per_deg``, and
        ``TrimFunctionEvaluations``. Ground-effect metadata include
        ``GroundEffectEnabled``, ``CGHeight``, and ``CGHeight_bref``.

    Raises
    ------
    ValueError
        If any required input sequence is empty, a reference value is invalid,
        or a numerical setting is invalid.
    RuntimeError
        If the required geometry, analysis input, elevator group, VSPAERO
        result, or trim root is unavailable.

    Notes
    -----
    For trimmed cases, the internal calculation order is out of ground effect,
    followed by finite CG heights from largest to smallest.  This makes the
    continuation path approach the ground gradually.  Returned rows are sorted
    back to the caller's original ``height`` order.  Both elevator deflection
    and :math:`\partial C_{My}/\partial\delta_e` are continued through the
    internal sequence.
    """

    alpha_deg = list(alpha_deg)
    mach = list(mach)
    reynolds = list(reynolds)
    height = list(height)
    if not alpha_deg or not mach or not reynolds or not height:
        raise ValueError(
            'alpha_deg, mach, reynolds, and height must not be empty.'
        )

    error_mgr = vsp.ErrorMgrSingleton.getInstance()

    def pop_openvsp_errors():
        errors = []
        while error_mgr.GetNumTotalErrors() > 0:
            error = error_mgr.PopLastError()
            errors.append(error.GetErrorString())
        return errors

    # Read model information and resolve the existing elevator once.
    vsp.Update()
    wing_id = find_one_geom(vsp, DEFAULT_REF_GEOM_NAME)
    vsp3_path = vsp.GetVSPFileName()
    openvsp_version = vsp.GetVSPVersion()

    settings_id = find_vspaero_settings_container(vsp)
    settings = {}
    if settings_id:
        for name in (
            'Sref', 'bref', 'cref', 'Xcg', 'Ycg', 'Zcg',
            'Symmetry', 'RefFlag',
        ):
            value, parm_id = get_container_parm_value(vsp, settings_id, name)
            if parm_id:
                settings[name] = value

    elevator_parm_id = ''
    original_elevator_deg = np.nan
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
        parm_id = vsp.FindParm(
            settings_id,
            'DeflectionAngle',
            f'ControlSurfaceGroup_{elevator_group_index}',
        )
        if parm_id and str(parm_id).upper() != 'NONE':
            elevator_parm_id = parm_id
            original_elevator_deg = float(vsp.GetParmVal(parm_id))

    if trimmed:
        if len(elevator_group_matches) != 1:
            raise RuntimeError(
                "Pitch trim requires exactly one existing 'ELEVATOR_GROUP'. "
                f'found={len(elevator_group_matches)}'
            )
        if not list(vsp.GetActiveCSNameVec(elevator_group_matches[0])):
            raise RuntimeError("'ELEVATOR_GROUP' has no active control surfaces.")
        if not elevator_parm_id:
            raise RuntimeError(
                "Could not find 'DeflectionAngle' for the existing "
                "'ELEVATOR_GROUP'."
            )

    # Build the mixed Thick/Thin geometry once for all polar and trim cases.
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
            'VSPAEROComputeGeometry did not produce a usable '
            'VSPAERO_Geom result.\n'
            f'  result_id={compgeom_result_id!r}\n'
            f'  result_name={compgeom_result_name!r}\n'
            f'  vspgeom_files={vspgeom_files!r}\n'
            f'  OpenVSP errors={compgeom_errors or ["<none reported>"]}'
        )

    # Configure one reusable one-point VSPAERO sweep.
    analysis_name = 'VSPAEROSweep'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))

    if ncpu is not None:
        ncpu = int(ncpu)
        if ncpu < 1:
            raise ValueError('ncpu must be a positive integer.')
        if 'NCPU' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'NCPU' input."
            )
        vsp.SetIntAnalysisInput(analysis_name, 'NCPU', [ncpu], 0)

    if wake_num_iter is not None:
        wake_num_iter = int(wake_num_iter)
        if not 3 <= wake_num_iter <= 255:
            raise ValueError('wake_num_iter must be between 3 and 255.')
        if 'WakeNumIter' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'WakeNumIter' input."
            )
        vsp.SetIntAnalysisInput(
            analysis_name,
            'WakeNumIter',
            [wake_num_iter],
            0,
        )

    if fixed_wake_flag is not None:
        if 'FixedWakeFlag' not in analysis_inputs:
            raise RuntimeError(
                "VSPAEROSweep does not expose the 'FixedWakeFlag' input."
            )
        vsp.SetIntAnalysisInput(
            analysis_name,
            'FixedWakeFlag',
            [1 if fixed_wake_flag else 0],
            0,
        )

    effective_ncpu = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'NCPU')[0])
        if 'NCPU' in analysis_inputs
        else None
    )
    effective_wake_num_iter = (
        int(vsp.GetIntAnalysisInput(analysis_name, 'WakeNumIter')[0])
        if 'WakeNumIter' in analysis_inputs
        else None
    )
    effective_fixed_wake_flag = (
        bool(vsp.GetIntAnalysisInput(analysis_name, 'FixedWakeFlag')[0])
        if 'FixedWakeFlag' in analysis_inputs
        else None
    )
    if wake_num_iter is not None and effective_fixed_wake_flag:
        warnings.warn(
            'wake_num_iter was specified while FixedWakeFlag is enabled; '
            'the fixed-wake setting takes precedence.',
            RuntimeWarning,
            stacklevel=2,
        )

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
            f'{analysis_name} is missing required inputs: {missing_inputs}.'
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
        raise ValueError(
            f'VSPAERO reference span bref must be positive. bref={bref}'
        )

    def run_polar_point(
        alpha,
        mach_value,
        reynolds_value,
        ground_effect_enabled,
        cg_height,
        elevator_deg,
        evaluation_label,
    ):
        if elevator_deg is not None:
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
            print(
                f'   evaluation {evaluation_label}: '
                f'alpha={float(alpha):.8g}, Mach={float(mach_value):.8g}, '
                f'Re={float(reynolds_value):.8g}, '
                f'ground={"ON" if ground_effect_enabled else "OFF"}, '
                f'CGHeight={cg_height if ground_effect_enabled else np.nan}, '
                f'elevator={elevator_deg}',
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

        polar = pd.DataFrame()
        for result_name in ('VSPAERO_Polar', 'VSPAERO Polar'):
            polar = (
                results_dataframe(vsp, wrapper_result_id, result_name)
                if wrapper_result_id
                else pd.DataFrame()
            )
            if not polar.empty:
                break

        if polar.empty:
            raise RuntimeError(
                'VSPAERO returned no polar result.\n'
                f'  alpha={alpha}, Mach={mach_value}, Re={reynolds_value}\n'
                f'  ground_effect_enabled={ground_effect_enabled}, '
                f'CGHeight={cg_height if ground_effect_enabled else None}\n'
                f'  elevator_deg={elevator_deg}\n'
                f'  child_results={child_result_names!r}\n'
                f'  OpenVSP errors={api_errors or ["<none reported>"]}'
            )
        if len(polar.index) != 1:
            raise RuntimeError(
                'A one-point VSPAERO analysis returned an unexpected number '
                f'of rows: {len(polar.index)}'
            )

        if 'CLtot' in polar.columns and 'CL' not in polar.columns:
            polar['CL'] = polar['CLtot']
        if 'CL' in polar.columns and 'CLtot' not in polar.columns:
            polar['CLtot'] = polar['CL']
        if 'CMytot' in polar.columns and 'CMy' not in polar.columns:
            polar['CMy'] = polar['CMytot']
        if 'CMy' in polar.columns and 'CMytot' not in polar.columns:
            polar['CMytot'] = polar['CMy']

        duration = analysis_duration_seconds(vsp, wrapper_result_id)
        duration_seconds = elapsed if duration is None else duration

        if verbose and int(verbose) >= 2:
            vsp.PrintResults(wrapper_result_id)
            print(polar.to_string(index=False), flush=True)

        return polar.copy(), float(duration_seconds)

    # Keep the public output in the caller's height order.  Trimmed cases are
    # calculated internally from ground-effect off to the largest finite height
    # and then toward the ground, which makes continuation between cases more
    # gradual without changing the returned row order.
    ground_conditions = [(0, False, None)] + [
        (output_index, True, float(value))
        for output_index, value in enumerate(height, start=1)
    ]
    calculation_ground_conditions = ground_conditions
    if trimmed:
        calculation_ground_conditions = [ground_conditions[0]] + sorted(
            ground_conditions[1:],
            key=lambda condition: condition[2],
            reverse=True,
        )

    requested_case_count = (
        len(alpha_deg)
        * len(mach)
        * len(reynolds)
        * len(ground_conditions)
    )

    if verbose:
        mode = 'trimmed polar (fixed alpha, CMytot=0)' if trimmed else 'polar'
        print('\n-> VSPAERO ground-effect sweep', flush=True)
        print(f' Mode          : {mode}', flush=True)
        print(f' OpenVSP       : {openvsp_version}', flush=True)
        print(f' VSP3          : {vsp3_path or "<unsaved model>"}', flush=True)
        print(
            f' Reference wing: {DEFAULT_REF_GEOM_NAME} ({wing_id})',
            flush=True,
        )
        print(f' Geometry sets : {set_source}', flush=True)
        print(
            f'   Thick [{thick_set}]: {thick_geom_names or ["<none>"]}',
            flush=True,
        )
        print(
            f'   Thin  [{thin_set}]: {thin_geom_names or ["<none>"]}',
            flush=True,
        )
        if settings:
            print(
                ' VSPAERO refs  : '
                + ', '.join(
                    f'{name}={value:.6g}'
                    for name, value in settings.items()
                ),
                flush=True,
            )
        if trimmed:
            print(
                f' Elevator trim : initial={original_elevator_deg:.6g} '
                f'deg, bounds={elevator_bounds}, '
                f'|CMytot|<={cmy_tolerance:.3g}',
                flush=True,
            )
        print(
            f' Geometry ready : {vspgeom_files[0]} '
            f'({compgeom_elapsed:.2f} s)',
            flush=True,
        )
        print(f' Cases         : {requested_case_count}', flush=True)
        print(
            ' VSPAERO run   : '
            f'NCPU={effective_ncpu}, '
            f'WakeNumIter={effective_wake_num_iter}, '
            f'FixedWakeFlag={effective_fixed_wake_flag}',
            flush=True,
        )

    case_frames = []
    case_index = 0
    case_group_index = 0

    try:
        # Each angle starts out of ground effect.  Trimmed finite-height cases
        # then proceed from the largest CG height toward the ground so the
        # inherited elevator deflection and slope change gradually.
        for reynolds_value in reynolds:
            for mach_value in mach:
                for alpha in alpha_deg:
                    output_group_index = case_group_index
                    case_group_index += 1
                    next_initial_elevator_deg = original_elevator_deg
                    next_cmy_slope = None
                    for (
                        ground_output_index,
                        ground_effect_enabled,
                        cg_height,
                    ) in calculation_ground_conditions:
                        case_index += 1
                        output_order = (
                            output_group_index * len(ground_conditions)
                            + ground_output_index
                        )
                        if verbose:
                            ground_text = (
                                f'ON, CGHeight={cg_height:.8g}'
                                if ground_effect_enabled
                                else 'OFF'
                            )
                            print(
                                f' [{case_index:>3}/{requested_case_count}] '
                                f'alpha={float(alpha):.8g} deg, '
                                f'Mach={float(mach_value):.8g}, '
                                f'Re={float(reynolds_value):.8g}, '
                                f'ground={ground_text}, '
                                f'trimmed={bool(trimmed)}',
                                flush=True,
                            )

                        if trimmed:
                            def evaluate_cmy(elevator_deg, evaluation_number):
                                result, duration = run_polar_point(
                                    alpha,
                                    mach_value,
                                    reynolds_value,
                                    ground_effect_enabled,
                                    cg_height,
                                    elevator_deg,
                                    f'trim {evaluation_number}',
                                )
                                if 'CMytot' not in result.columns:
                                    raise RuntimeError(
                                        "VSPAERO polar does not contain "
                                        "'CMytot' or 'CMy'."
                                    )
                                return (
                                    float(result['CMytot'].iloc[0]),
                                    result,
                                    duration,
                                )

                            (
                                case_result,
                                elevator_deg,
                                cmy_residual,
                                cmy_slope,
                                trim_evaluations,
                                case_duration,
                            ) = _trim_elevator_at_fixed_alpha(
                                evaluate_cmy,
                                next_initial_elevator_deg,
                                previous_cmy_slope=next_cmy_slope,
                                elevator_bounds=elevator_bounds,
                                cmy_tolerance=cmy_tolerance,
                                max_trim_evaluations=max_trim_evaluations,
                                trim_step_deg=trim_step_deg,
                            )
                            next_initial_elevator_deg = elevator_deg
                            next_cmy_slope = cmy_slope
                        else:
                            case_result, case_duration = run_polar_point(
                                alpha,
                                mach_value,
                                reynolds_value,
                                ground_effect_enabled,
                                cg_height,
                                None,
                                'polar',
                            )
                            elevator_deg = original_elevator_deg
                            cmy_residual = (
                                float(case_result['CMytot'].iloc[0])
                                if 'CMytot' in case_result.columns
                                else np.nan
                            )
                            cmy_slope = np.nan
                            trim_evaluations = 0

                        case_result.insert(0, 'Trimmed', bool(trimmed))
                        case_result.insert(
                            1,
                            'GroundEffectEnabled',
                            bool(ground_effect_enabled),
                        )
                        case_result.insert(2, 'InputAlpha_deg', float(alpha))
                        case_result.insert(3, 'InputMach', float(mach_value))
                        case_result.insert(
                            4,
                            'InputReCref',
                            float(reynolds_value),
                        )
                        case_result.insert(5, 'NCPU', effective_ncpu)
                        case_result.insert(
                            6,
                            'WakeNumIter',
                            effective_wake_num_iter,
                        )
                        case_result.insert(
                            7,
                            'FixedWakeFlag',
                            effective_fixed_wake_flag,
                        )
                        case_result.insert(
                            8,
                            'CGHeight',
                            float(cg_height)
                            if ground_effect_enabled
                            else np.nan,
                        )
                        case_result.insert(
                            9,
                            'CGHeight_bref',
                            float(cg_height) / bref
                            if ground_effect_enabled
                            else np.nan,
                        )
                        case_result.insert(
                            10,
                            'Elevator_deg',
                            float(elevator_deg),
                        )
                        case_result.insert(
                            11,
                            'CMyResidual',
                            float(cmy_residual),
                        )
                        case_result.insert(
                            12,
                            'CMyElevatorSlope_per_deg',
                            float(cmy_slope)
                            if cmy_slope is not None
                            else np.nan,
                        )
                        case_result.insert(
                            13,
                            'TrimFunctionEvaluations',
                            int(trim_evaluations),
                        )
                        case_result.insert(
                            14,
                            'CaseDuration_s',
                            float(case_duration),
                        )
                        case_frames.append((output_order, case_result))

                        if verbose:
                            summary_columns = [
                                name
                                for name in (
                                    'Elevator_deg', 'CL', 'CDtot', 'CDi',
                                    'CD0', 'L_D', 'CMytot', 'CMyResidual',
                                    'CMyElevatorSlope_per_deg',
                                    'TrimFunctionEvaluations',
                                )
                                if name in case_result.columns
                            ]
                            print(
                                f'   complete: {case_duration:.2f} s',
                                flush=True,
                            )
                            print(
                                case_result[summary_columns].to_string(index=False),
                                flush=True,
                            )
    finally:
        if elevator_parm_id and np.isfinite(original_elevator_deg):
            vsp.SetParmVal(elevator_parm_id, original_elevator_deg)
            vsp.Update()

    case_frames.sort(key=lambda item: item[0])
    return pd.concat(
        [case_result for _, case_result in case_frames],
        ignore_index=True,
    )

def validate_vsp3_for_stability_derivatives(vsp3_path, *, verbose=1):
    """Validate saved VSPAERO settings for stability-derivative analysis.

    The check is intentionally model-agnostic.  It does not require specific
    Geom names, Set names, Subsurface names, or control-surface groups.  The
    model's saved ``GeomSet`` and ``ThinGeomSet`` selections are treated as the
    aerodynamic model definition, including valid all-thin or all-thick cases.

    The function reads the ``.vsp3`` file but does not execute VSPAERO.
    """
    report = {
        'passed': False,
        'errors': [],
        'warnings': [],
        'infos': [],
        'set_summary': {},
        'symmetry_summary': {},
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

    def is_finite_number(value):
        try:
            return np.isfinite(float(value))
        except (TypeError, ValueError):
            return False

    vsp3_path = os.fspath(vsp3_path)
    if not os.path.isfile(vsp3_path):
        add(
            'errors',
            'FILE_NOT_FOUND',
            'The specified .vsp3 file was not found.',
            {'vsp3_path': vsp3_path},
        )
        return report

    vprint(1, f'\n-> Validate VSPAERO stability-derivative preflight: {vsp3_path}')

    try:
        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(vsp3_path)
        vsp.Update()
    except Exception as err:
        add(
            'errors',
            'READ_VSP3_FAILED',
            'OpenVSP failed to read the .vsp3 file.',
            {'error': repr(err)},
        )
        return report

    add('infos', 'READ_VSP3_OK', 'The .vsp3 file was read by OpenVSP.', {'vsp3_path': vsp3_path})

    # Read the same saved geometry selection that VSPAEROComputeGeometry will
    # use.  Set names are deliberately irrelevant here.
    vprint(1, ' Checking saved VSPAERO geometry selection...')
    compgeom_name = 'VSPAEROComputeGeometry'
    try:
        vsp.SetAnalysisInputDefaults(compgeom_name)
        compgeom_inputs = set(vsp.GetAnalysisInputNames(compgeom_name))
    except Exception as err:
        add(
            'errors',
            'VSPAERO_COMPUTE_GEOMETRY_INPUTS_FAILED',
            'Could not read VSPAEROComputeGeometry inputs.',
            {'error': repr(err)},
        )
        return report

    if 'GeomSet' not in compgeom_inputs:
        add(
            'errors',
            'MISSING_GEOMSET_INPUT',
            "VSPAEROComputeGeometry does not expose the 'GeomSet' input.",
        )
        return report

    thick_set = int(vsp.GetIntAnalysisInput(compgeom_name, 'GeomSet')[0])
    thin_set = (
        int(vsp.GetIntAnalysisInput(compgeom_name, 'ThinGeomSet')[0])
        if 'ThinGeomSet' in compgeom_inputs
        else None
    )
    set_none = int(getattr(vsp, 'SET_NONE', -1))

    selected_geom_ids = set()
    for label, set_index in (('GeomSet', thick_set), ('ThinGeomSet', thin_set)):
        summary = {
            'set_index': set_index,
            'selected': set_index is not None and set_index != set_none,
            'geom_ids': [],
            'geom_names': [],
        }
        if summary['selected']:
            try:
                geom_ids = list(vsp.GetGeomSetAtIndex(int(set_index)))
            except Exception as err:
                add(
                    'errors',
                    'GEOM_SET_READ_FAILED',
                    f'Could not read the selected {label}.',
                    {'set_index': set_index, 'error': repr(err)},
                )
                geom_ids = []
            summary['geom_ids'] = geom_ids
            summary['geom_names'] = [vsp.GetGeomName(geom_id) for geom_id in geom_ids]
            selected_geom_ids.update(geom_ids)
        report['set_summary'][label] = summary

    if not selected_geom_ids:
        add(
            'errors',
            'NO_AERO_GEOMETRY',
            'The saved VSPAERO geometry selection contains no Geoms.',
            {'GeomSet': thick_set, 'ThinGeomSet': thin_set},
        )
    elif not report['set_summary']['GeomSet']['geom_ids']:
        add(
            'infos',
            'NO_THICK_GEOMETRY',
            'No thick geometry is selected; the stability analysis will use the selected thin geometry only.',
        )

    # Validate only the reference quantities and symmetry state required by the
    # stability calculation.  The actual geometry/control naming is not part of
    # analysis readiness.
    vprint(1, ' Checking VSPAERO reference settings...')
    settings_id = find_vspaero_settings_container(vsp)
    settings_summary = {
        'container_id': settings_id,
        'Sref': None,
        'bref': None,
        'cref': None,
        'Xcg': None,
        'Ycg': None,
        'Zcg': None,
        'RefFlag': None,
        'RefGeomID': '',
        'GeomSet': thick_set,
        'ThinGeomSet': thin_set,
        'Symmetry': None,
    }

    if not settings_id:
        add('errors', 'MISSING_VSPAERO_SETTINGS', "The 'VSPAEROSettings' container was not found.")
    else:
        for parm_name in ('Sref', 'bref', 'cref', 'Xcg', 'Ycg', 'Zcg', 'RefFlag', 'Symmetry'):
            value, parm_id = get_container_parm_value(vsp, settings_id, parm_name)
            settings_summary[parm_name] = value
            if parm_id == '':
                add(
                    'errors',
                    'MISSING_VSPAERO_PARM',
                    f"VSPAERO setting '{parm_name}' was not found.",
                    {'parm_name': parm_name},
                )

        for parm_name in ('Sref', 'bref', 'cref'):
            value = settings_summary[parm_name]
            if value is not None and (not is_finite_number(value) or float(value) <= 0.0):
                add(
                    'errors',
                    'INVALID_REFERENCE_VALUE',
                    f"VSPAERO reference value '{parm_name}' must be positive and finite.",
                    {'parm_name': parm_name, 'value': value},
                )

        for parm_name in ('Xcg', 'Ycg', 'Zcg'):
            value = settings_summary[parm_name]
            if value is not None and not is_finite_number(value):
                add(
                    'errors',
                    'INVALID_MOMENT_REFERENCE',
                    f"VSPAERO moment reference '{parm_name}' must be finite.",
                    {'parm_name': parm_name, 'value': value},
                )

        symmetry_value = settings_summary['Symmetry']
        if symmetry_value is not None and int(round(float(symmetry_value))) != 0:
            add(
                'errors',
                'VSPAERO_XZ_SYMMETRY_ENABLED',
                'VSPAERO Settings Symmetry must be 0 for full stability-derivative calculation.',
                {'Symmetry': symmetry_value},
            )

        try:
            ref_geom_id = vsp.GetVSPAERORefWingID()
        except Exception:
            ref_geom_id = ''
        settings_summary['RefGeomID'] = ref_geom_id

        component_ref = int(getattr(vsp, 'COMPONENT_REF', 1))
        ref_flag = settings_summary['RefFlag']
        if ref_flag is not None and int(round(float(ref_flag))) == component_ref:
            if not ref_geom_id or str(ref_geom_id).upper() == 'NONE':
                add(
                    'errors',
                    'MISSING_REFERENCE_GEOM',
                    'VSPAERO uses component reference dimensions, but no reference wing is selected.',
                )

    report['symmetry_summary'] = {
        'value': settings_summary['Symmetry'],
        'expected_value': 0.0,
        'passed': (
            settings_summary['Symmetry'] is not None
            and int(round(float(settings_summary['Symmetry']))) == 0
        ),
    }
    report['vspaero_settings_summary'] = settings_summary

    # Control groups are optional for stability derivatives.  Record what is
    # present without turning missing aileron/elevator/rudder conventions into
    # analysis errors.
    vprint(1, ' Checking optional VSPAERO control-surface groups...')
    try:
        for group_index in range(vsp.GetNumControlSurfaceGroups()):
            group_name = vsp.GetVSPAEROControlGroupName(group_index)
            report['control_group_summary'][group_name] = {
                'group_index': group_index,
                'active_control_surfaces': list(vsp.GetActiveCSNameVec(group_index)),
            }
    except Exception as err:
        add(
            'warnings',
            'CONTROL_GROUP_LIST_FAILED',
            'Could not read VSPAERO control-surface groups.',
            {'error': repr(err)},
        )

    report['passed'] = len(report['errors']) == 0

    if verbose:
        status = 'PASSED' if report['passed'] else 'FAILED'
        print(
            f' Validation {status}: {len(report["errors"])} error(s), '
            f'{len(report["warnings"])} warning(s)'
        )
        if int(verbose) >= 2:
            print(
                '    Geometry selection: '
                f'GeomSet={thick_set} '
                f'({len(report["set_summary"]["GeomSet"]["geom_ids"])} Geom(s)), '
                f'ThinGeomSet={thin_set} '
                f'({len(report["set_summary"]["ThinGeomSet"]["geom_ids"])} Geom(s))'
            )
            for error in report['errors']:
                print(f"    ERROR {error['code']}: {error['message']}")
            for warning in report['warnings']:
                print(f"    WARNING {warning['code']}: {warning['message']}")
            for info in report['infos']:
                print(f"    INFO {info['code']}: {info['message']}")

    return report

def make_CDo_correction(
    vsp, 
    trimmed_polar, 
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
    reynolds = trimmed_polar['Re_1e6'].values * 1e6  # Get Reynolds number

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
    trimmed_polar['CDo_corr'] = CDo * form_factor * interference_factor + CDpCL * (trimmed_polar['CL'].values) ** 2
    trimmed_polar['CDtot_corr'] = trimmed_polar['CDo_corr'].values + trimmed_polar['CDi'].values
    trimmed_polar['L_D_corr'] = trimmed_polar['CL'].values / trimmed_polar['CDtot_corr'].values

    # Calculate angle of attack from modified lift-drag ratio
    trimmed_polar['gamma'] = np.arctan(1 / trimmed_polar['L_D_corr'].values)

    # Calculate velocity
    trimmed_polar['Velocity'] = np.sqrt((2 * Weight * g) / (density * Sref * trimmed_polar['CL'].values * np.cos(trimmed_polar['gamma'].values)))

    # Calculates horizontal velocity (Vx) and vertical velocity (Vz) and adds them to the data frame
    trimmed_polar['Vx'] = trimmed_polar['Velocity'].values * np.cos(trimmed_polar['gamma'].values) * 3.6  # horizontal velocity [km/h]
    trimmed_polar['vz'] = trimmed_polar['Velocity'].values * np.sin(trimmed_polar['gamma'].values)  # vertical velocity [m/s]
    return trimmed_polar  # Returns corrected data

