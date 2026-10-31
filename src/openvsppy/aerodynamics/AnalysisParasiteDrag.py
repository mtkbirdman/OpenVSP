"""OpenVSP ParasiteDrag analysis and VSPAERO polar correction.

The functions in this module keep the two analyses separate:

* ``vsp_parasite_drag`` runs OpenVSP's Parasite Drag Tool at specified flight
  speeds using the ParasiteDrag settings saved in the loaded ``.vsp3`` model.
* ``apply_parasite_drag_to_polar`` combines those parasite-drag results with a
  VSPAERO polar using ``CL = CLi`` and ``CD = CD0 + CDi``.

Component skin-friction equations, form factors, laminar percentages,
interference factors, roughness, and excrescences remain model settings in
OpenVSP rather than becoming a second set of Python arguments.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pandas as pd

from ..core.ISAspecification import mach_to_velocity
from ..core.util import (
    analysis_duration_seconds,
    find_vspaero_settings_container,
    get_container_parm_value,
    suppress_stdout,
)


def vsp_parasite_drag(
    vsp,
    velocity_mps,
    altitude_m=0.0,
    dT=0.0,
    *,
    sref=None,
    geom_set=None,
    verbose=1,
):
    """Run OpenVSP ParasiteDrag for one or more true airspeeds.

    Parameters
    ----------
    vsp : module or object
        OpenVSP Python API object with the desired ``.vsp3`` already loaded.
    velocity_mps : float or sequence of float
        True airspeed in m/s.  Each value produces one ParasiteDrag result row.
    altitude_m : float, optional
        Geometric altitude in metres for the US Standard Atmosphere 1976 model.
    dT : float, optional
        Temperature offset from the standard atmosphere in K.  Temperature
        differences in K and deg C have the same numeric value.
    sref : float or None, optional
        Manual reference area.  The numeric area must use the square of the
        ParasiteDrag ``LengthUnit`` saved in the loaded model.  ``None`` keeps
        the reference-area setting saved in the model.
    geom_set : int or None, optional
        ParasiteDrag geometry-set index.  ``None`` keeps the set saved in the
        loaded model.
    verbose : int or bool, optional
        ``0`` suppresses progress, ``1`` prints one line per condition, and
        values of ``2`` or greater also print the OpenVSP analysis inputs.

    Returns
    -------
    pandas.DataFrame
        One row per requested speed with ``Velocity_mps``, ``Altitude_m``,
        ``Mach``, ``Re_L``, ``Sref``, ``CD0``, and ``AnalysisDuration_s``.

    Notes
    -----
    The function intentionally leaves component Parasite Drag properties and
    friction/form-factor selections in the loaded model.  Geometry is
    recomputed for the first condition and reused for subsequent speeds.
    """

    velocities = np.asarray(velocity_mps, dtype=float)
    if velocities.ndim == 0:
        velocities = velocities.reshape(1)
    if velocities.ndim != 1 or velocities.size == 0:
        raise ValueError('velocity_mps must be a scalar or a non-empty 1D sequence.')
    if not np.isfinite(velocities).all() or np.any(velocities <= 0.0):
        raise ValueError('velocity_mps must contain only positive finite values.')

    altitude_m = float(altitude_m)
    dT = float(dT)
    if sref is not None:
        sref = float(sref)
        if not np.isfinite(sref) or sref <= 0.0:
            raise ValueError('sref must be a positive finite value.')

    analysis_name = 'ParasiteDrag'
    vsp.SetAnalysisInputDefaults(analysis_name)
    analysis_inputs = set(vsp.GetAnalysisInputNames(analysis_name))

    required_inputs = {
        'VelocityUnit',
        'AltLengthUnit',
        'TempUnit',
        'FreestreamPropChoice',
        'Vinf',
        'Altitude',
        'DeltaTemp',
        'RecomputeGeom',
    }
    if sref is not None:
        required_inputs.update({'RefFlag', 'Sref'})
    if geom_set is not None:
        required_inputs.add('GeomSet')

    missing_inputs = sorted(required_inputs - analysis_inputs)
    if missing_inputs:
        raise RuntimeError(
            f'{analysis_name} is missing required inputs: {missing_inputs}.'
        )

    vsp.SetIntAnalysisInput(
        analysis_name,
        'VelocityUnit',
        [int(vsp.V_UNIT_M_S)],
        0,
    )
    vsp.SetIntAnalysisInput(
        analysis_name,
        'AltLengthUnit',
        [int(vsp.PD_UNITS_METRIC)],
        0,
    )
    vsp.SetIntAnalysisInput(
        analysis_name,
        'TempUnit',
        [int(vsp.TEMP_UNIT_K)],
        0,
    )
    vsp.SetIntAnalysisInput(
        analysis_name,
        'FreestreamPropChoice',
        [int(vsp.ATMOS_TYPE_US_STANDARD_1976)],
        0,
    )

    if geom_set is not None:
        vsp.SetIntAnalysisInput(analysis_name, 'GeomSet', [int(geom_set)], 0)

    if sref is not None:
        vsp.SetIntAnalysisInput(
            analysis_name,
            'RefFlag',
            [int(vsp.MANUAL_REF)],
            0,
        )
        vsp.SetDoubleAnalysisInput(analysis_name, 'Sref', [sref], 0)

    if 'FileName' in analysis_inputs:
        vsp.SetStringAnalysisInput(analysis_name, 'FileName', [os.devnull], 0)

    vsp.SetDoubleAnalysisInput(analysis_name, 'Altitude', [altitude_m], 0)
    vsp.SetDoubleAnalysisInput(analysis_name, 'DeltaTemp', [dT], 0)

    if verbose and int(verbose) >= 2:
        print(f'\n{analysis_name} inputs')
        vsp.PrintAnalysisInputs(analysis_name)

    rows = []
    result_fields = {
        'Velocity_mps': 'FC_Vinf',
        'Altitude_m': 'FC_Alt',
        'Mach': 'FC_Mach',
        'Re_L': 'Re_L',
        'Sref': 'FC_Sref',
        'CD0': 'Total_CD_Total',
    }

    for index, velocity in enumerate(velocities):
        vsp.SetDoubleAnalysisInput(analysis_name, 'Vinf', [float(velocity)], 0)
        vsp.SetIntAnalysisInput(
            analysis_name,
            'RecomputeGeom',
            [1 if index == 0 else 0],
            0,
        )

        start = time.perf_counter()
        if verbose and int(verbose) >= 2:
            result_id = vsp.ExecAnalysis(analysis_name)
        else:
            with suppress_stdout():
                result_id = vsp.ExecAnalysis(analysis_name)
        elapsed = time.perf_counter() - start

        if not result_id:
            raise RuntimeError(
                'ParasiteDrag returned no result ID for '
                f'V={float(velocity):.8g} m/s.'
            )

        row = {}
        missing_results = []
        for column_name, result_name in result_fields.items():
            values = vsp.GetDoubleResults(result_id, result_name, 0)
            if not values:
                missing_results.append(result_name)
                continue
            row[column_name] = float(values[0])

        if missing_results:
            raise RuntimeError(
                'ParasiteDrag result is missing required fields: '
                f'{missing_results}.'
            )

        duration = analysis_duration_seconds(vsp, result_id)
        row['AnalysisDuration_s'] = (
            float(elapsed) if duration is None else float(duration)
        )
        rows.append(row)

        if verbose:
            print(
                ' ParasiteDrag: '
                f'V={row["Velocity_mps"]:.6g} m/s, '
                f'Mach={row["Mach"]:.6g}, '
                f'CD0={row["CD0"]:.8g}',
                flush=True,
            )

    return pd.DataFrame(rows)


def apply_parasite_drag_to_polar(
    vsp,
    polar,
    altitude_m=0.0,
    dT=0.0,
    *,
    sref=None,
    geom_set=None,
    verbose=1,
):
    """Replace a VSPAERO polar's lift/drag model with CLi and ParasiteDrag CD0.

    The corrected coefficients are defined as

    ``CL = CLi``

    ``CD = CD0 + CDi``

    where ``CD0`` is calculated by OpenVSP ParasiteDrag at the Mach number of
    each input row.  The input DataFrame is not modified.

    Parameters
    ----------
    vsp : module or object
        OpenVSP Python API object with the same model used for the VSPAERO polar.
    polar : pandas.DataFrame
        VSPAERO polar containing ``CLi``, ``CDi``, and either ``InputMach`` or
        ``Mach``.
    altitude_m : float, optional
        Altitude in metres used for the ParasiteDrag atmosphere calculation.
    dT : float, optional
        Temperature offset from standard atmosphere in K.
    sref : float or None, optional
        Reference area used by ParasiteDrag.  ``None`` reads ``Sref`` from the
        loaded model's VSPAERO settings so ``CD0`` and ``CDi`` use the same
        reference-area value.
    geom_set : int or None, optional
        ParasiteDrag geometry-set index.  ``None`` keeps the model setting.
    verbose : int or bool, optional
        Passed through to :func:`vsp_parasite_drag`.

    Returns
    -------
    pandas.DataFrame
        Copy of ``polar`` with ``CL``, ``CD0``, ``CD``, and ``L_D`` replaced or
        added according to the correction above.  If the input already has a
        ``CD0`` column, its original values are retained as ``CD0_VSPAERO``.

    Notes
    -----
    This is an aerodynamic-coefficient correction at the existing VSPAERO
    operating points.  It does not iterate Mach/Reynolds number after replacing
    the lift and drag coefficients, and therefore is not a coupled fixed-weight
    glide-performance solver.
    """

    if not isinstance(polar, pd.DataFrame):
        raise TypeError('polar must be a pandas.DataFrame.')
    if polar.empty:
        return polar.copy()

    required_columns = {'CLi', 'CDi'}
    missing_columns = sorted(required_columns - set(polar.columns))
    if missing_columns:
        raise ValueError(f'polar is missing required columns: {missing_columns}.')

    if 'InputMach' in polar.columns:
        mach_column = 'InputMach'
    elif 'Mach' in polar.columns:
        mach_column = 'Mach'
    else:
        raise ValueError("polar must contain either 'InputMach' or 'Mach'.")

    if sref is None:
        settings_id = find_vspaero_settings_container(vsp)
        if not settings_id:
            raise RuntimeError('Could not find the loaded model VSPAERO settings.')
        sref, sref_parm_id = get_container_parm_value(vsp, settings_id, 'Sref')
        if not sref_parm_id or sref is None or float(sref) <= 0.0:
            raise RuntimeError('Could not read a positive VSPAERO Sref.')
    sref = float(sref)

    mach = polar[mach_column].to_numpy(dtype=float)
    if not np.isfinite(mach).all() or np.any(mach < 0.0):
        raise ValueError(f'{mach_column} must contain only non-negative finite values.')

    velocity_mps = np.array(
        [mach_to_velocity(value, altitude=altitude_m, dT=dT) for value in mach],
        dtype=float,
    )
    parasite_drag = vsp_parasite_drag(
        vsp,
        velocity_mps,
        altitude_m=altitude_m,
        dT=dT,
        sref=sref,
        geom_set=geom_set,
        verbose=verbose,
    )

    if len(parasite_drag) != len(polar):
        raise RuntimeError(
            'ParasiteDrag returned a different number of rows than the VSPAERO polar.'
        )

    corrected = polar.copy()
    if 'CD0' in corrected.columns:
        corrected['CD0_VSPAERO'] = corrected['CD0'].to_numpy(copy=True)

    corrected['CL'] = corrected['CLi'].to_numpy(dtype=float)
    corrected['CD0'] = parasite_drag['CD0'].to_numpy(dtype=float)
    corrected['CD'] = corrected['CD0'].to_numpy() + corrected['CDi'].to_numpy(dtype=float)
    corrected['L_D'] = corrected['CL'].to_numpy() / corrected['CD'].to_numpy()

    return corrected
