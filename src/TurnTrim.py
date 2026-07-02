"""
Steady turn trim solvers based on a VSPAERO .stab linear aero model.

For a powered steady level turn, the solver uses the ten turn parameters

    V, alpha, beta, phi, theta, Omega, delta_e, delta_a, delta_r, T

and solves the seven steady-level-turn equations after any three parameters
are fixed by the caller.

For an unpowered gliding turn, T is fixed to zero.  The solver uses the nine
turn parameters

    V, alpha, beta, phi, theta, Omega, delta_e, delta_a, delta_r

and solves the six force/moment balance equations after any three parameters
are fixed by the caller.  The height constraint is not imposed; the resulting
z_dot, h_dot, and sink_rate are returned as derived values.

Internal units are SI-style and radians unless the caller intentionally uses a
consistent non-SI unit system:

    V        : m/s, or the same length/time unit used with rho, mass, Sref
    alpha    : rad
    beta     : rad
    phi      : rad
    theta    : rad
    Omega    : rad/s, or rad/Tunit consistent with V
    delta_*  : rad
    T        : N, or the force unit consistent with rho, V, Sref
    mass     : kg, or the mass unit consistent with rho and g
    rho      : kg/m^3, or density consistent with V, Sref, mass, g

If rho=None is passed to solve_steady_level_turn() or solve_steady_gliding_turn(),
the Rho value in the .stab file is used.  This is convenient when the .stab file and the caller's mass,
length, and time units are consistent.  It is not a unit conversion.

Coordinate conventions used in this module:

    Solver body axes : +x forward, +y right, +z down
    OpenVSP body axes: +X aft, +Y right, +Z up

The .stab CL/CD/CS coefficients are converted back to the solver body-axis
CX/CY/CZ using alpha and beta.  The .stab CMl/CMm/CMn columns are used as
Cl/Cm/Cn.  They are body-axis moments with aviation sign convention, not true
wind-axis moments.  The p/q/r derivatives are reduced-rate derivatives:

    p_hat = p * Bref / (2 V)
    q_hat = q * Cref / (2 V)
    r_hat = r * Bref / (2 V)

Control derivatives are with respect to VSPAERO control-group deflection in
radians.  The control group names in the .stab case table are used to map
ConGrp_1, ConGrp_2, ... to delta_a, delta_e, and delta_r.
"""

from __future__ import annotations

import math
import os
import re
from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

TURN_PARAMETERS = [
    "V",
    "alpha",
    "beta",
    "phi",
    "theta",
    "Omega",
    "delta_e",
    "delta_a",
    "delta_r",
    "T",
]

GLIDING_TURN_PARAMETERS = [
    "V",
    "alpha",
    "beta",
    "phi",
    "theta",
    "Omega",
    "delta_e",
    "delta_a",
    "delta_r",
]

AERO_COEFFICIENTS = ["CL", "CD", "CS", "CMl", "CMm", "CMn"]
BASE_AERO_COLUMNS = ["CFx", "CFy", "CFz", "CMx", "CMy", "CMz", "CL", "CD", "CS", "CMl", "CMm", "CMn"]
STABILITY_CASE_NAMES = {"Base_Aero", "Alpha", "Beta", "Roll__Rate", "Pitch_Rate", "Yaw___Rate", "Mach"}
STAB_DERIVATIVE_COLUMNS = [
    "Alpha",
    "Beta",
    "p",
    "q",
    "r",
    "Mach",
    "U",
    "ConGrp_1",
    "ConGrp_2",
    "ConGrp_3",
]

CONTROL_NAME_HINTS = {
    "delta_a": "AILERON",
    "delta_e": "ELEVATOR",
    "delta_r": "RUDDER",
}
CONTROL_COLUMNS_FALLBACK = {
    "delta_a": "ConGrp_1",
    "delta_e": "ConGrp_2",
    "delta_r": "ConGrp_3",
}

@dataclass
class VSPAEROStab:
    """Parsed values from a VSPAERO .stab file."""

    path: str
    references: dict
    base_condition: dict
    base_aero: dict
    derivatives: pd.DataFrame
    control_groups: dict

    @property
    def Sref(self) -> float:
        return float(self.references["Sref"])

    @property
    def Cref(self) -> float:
        return float(self.references["Cref"])

    @property
    def Bref(self) -> float:
        return float(self.references["Bref"])

    @property
    def alpha0(self) -> float:
        return math.radians(float(self.base_condition.get("AoA", 0.0)))

    @property
    def beta0(self) -> float:
        return math.radians(float(self.base_condition.get("Beta", 0.0)))

    @property
    def V0(self) -> float:
        return float(self.base_condition.get("Vinf", 0.0) or 0.0)

    @property
    def rho0(self) -> float | None:
        value = self.base_condition.get("Rho")
        return None if value is None else float(value)

def _clean_stab_name(name: str) -> str:
    name = name.strip()
    return name[:-1] if name.endswith("_") else name

def read_vspaero_stab(stab_path: str | os.PathLike) -> VSPAEROStab:
    """
    Read the parts of a VSPAERO .stab file needed by solve_steady_level_turn().

    This parser is intentionally narrow: it reads the standard scalar header,
    the Base_Aero row, the control-group case rows, and the derivative table
    headed by "Coef Base Aero".
    """

    stab_path = os.fspath(stab_path)
    with open(stab_path, "r", encoding="utf-8", errors="ignore") as file:
        lines = file.readlines()

    scalar_values = {}
    for line in lines:
        if line.lstrip().startswith("Case"):
            break
        match = re.match(r"^\s*([A-Za-z][A-Za-z0-9_]*_?)\s+([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\b", line)
        if match:
            scalar_values[_clean_stab_name(match.group(1))] = float(match.group(2))

    references = {
        "Sref": scalar_values["Sref"],
        "Cref": scalar_values["Cref"],
        "Bref": scalar_values["Bref"],
        "Xcg": scalar_values.get("Xcg"),
        "Ycg": scalar_values.get("Ycg"),
        "Zcg": scalar_values.get("Zcg"),
    }
    base_condition = {
        "Mach": scalar_values.get("Mach"),
        "AoA": scalar_values.get("AoA", 0.0),
        "Beta": scalar_values.get("Beta", 0.0),
        "Rho": scalar_values.get("Rho"),
        "Vinf": scalar_values.get("Vinf", 0.0),
    }

    base_aero = None
    derivative_start = None
    control_group_names = []
    for index, line in enumerate(lines):
        if line.lstrip().startswith("Base_Aero"):
            parts = line.split()
            base_aero = {name: float(value) for name, value in zip(BASE_AERO_COLUMNS, parts[3:15])}

        parts = line.split()
        if len(parts) >= 15 and parts[0] not in STABILITY_CASE_NAMES and not parts[0].startswith("#"):
            try:
                float(parts[1])
            except ValueError:
                pass
            else:
                control_group_names.append(parts[0])

        if line.lstrip().startswith("Coef") and "Alpha" in line and "ConGrp_1" in line:
            derivative_start = index + 4
            break

    if base_aero is None:
        raise ValueError("Base_Aero row was not found in the .stab file.")
    if derivative_start is None:
        raise ValueError("Derivative table was not found in the .stab file.")

    derivative_rows = []
    derivative_column_names = None
    for line in lines:
        if line.lstrip().startswith("Coef") and "Alpha" in line and "ConGrp_1" in line:
            header_parts = line.split()
            derivative_column_names = header_parts[2:]
            break

    if derivative_column_names is None:
        derivative_column_names = STAB_DERIVATIVE_COLUMNS

    for line in lines[derivative_start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < len(derivative_column_names) + 2:
            break
        coef_name = parts[0]
        if coef_name not in BASE_AERO_COLUMNS:
            break
        values = [float(value) for value in parts[1:len(derivative_column_names) + 2]]
        derivative_rows.append([coef_name] + values)

    columns = ["Coef", "Base"] + derivative_column_names
    derivatives = pd.DataFrame(derivative_rows, columns=columns).set_index("Coef")
    control_groups = {f"ConGrp_{index + 1}": name for index, name in enumerate(control_group_names)}

    return VSPAEROStab(
        path=stab_path,
        references=references,
        base_condition=base_condition,
        base_aero=base_aero,
        derivatives=derivatives,
        control_groups=control_groups,
    )

def _resolve_rho(stab: VSPAEROStab, rho: float | None) -> float:
    if rho is not None:
        return float(rho)
    if stab.rho0 is None:
        raise ValueError("rho=None was passed, but the .stab file does not contain Rho.")
    return float(stab.rho0)

def _control_columns_from_stab(stab: VSPAEROStab, control_map: Mapping[str, str] | None = None) -> dict:
    """
    Map delta_a/e/r to ConGrp_* columns.

    control_map values may be either ConGrp_* names or control-group names from
    the .stab case table.  Without an explicit map, group names containing
    AILERON, ELEVATOR, and RUDDER are used, with the old ConGrp_1/2/3 convention
    as a fallback for G103A-style files.
    """

    available_columns = set(stab.derivatives.columns)
    group_to_column = {group_name: column for column, group_name in stab.control_groups.items()}
    result = {}

    if control_map:
        for delta_name, target in control_map.items():
            if target in available_columns:
                result[delta_name] = target
            elif target in group_to_column:
                result[delta_name] = group_to_column[target]
        return result

    for delta_name, hint in CONTROL_NAME_HINTS.items():
        for column, group_name in stab.control_groups.items():
            if hint in group_name.upper():
                result[delta_name] = column
                break

    for delta_name, fallback_column in CONTROL_COLUMNS_FALLBACK.items():
        if delta_name not in result and fallback_column in available_columns:
            result[delta_name] = fallback_column

    return result

def _default_initial_guess(stab: VSPAEROStab, fixed: Mapping[str, float], mass: float, rho: float, g: float) -> dict:
    V = float(fixed.get("V", stab.V0 if stab.V0 else 30.0))
    alpha = float(fixed.get("alpha", stab.alpha0))
    beta = float(fixed.get("beta", stab.beta0))
    phi = float(fixed.get("phi", 0.0))
    theta = float(fixed.get("theta", alpha))
    Omega = float(fixed.get("Omega", 0.0))
    qbar = 0.5 * rho * V * V
    T = qbar * stab.Sref * float(stab.base_aero.get("CD", 0.02))

    return {
        "V": V,
        "alpha": alpha,
        "beta": beta,
        "phi": phi,
        "theta": theta,
        "Omega": Omega,
        "delta_e": 0.0,
        "delta_a": 0.0,
        "delta_r": 0.0,
        "T": T if math.isfinite(T) else 0.0,
    }

def _complete_parameters(fixed: Mapping[str, float], unknown_names: list[str], unknown_values: np.ndarray) -> dict:
    params = {name: float(value) for name, value in fixed.items()}
    params.update({name: float(value) for name, value in zip(unknown_names, unknown_values)})
    return params

def _steady_turn_derived(params: Mapping[str, float], stab: VSPAEROStab, rho: float) -> dict:
    V = params["V"]
    alpha = params["alpha"]
    beta = params["beta"]
    phi = params["phi"]
    theta = params["theta"]
    Omega = params["Omega"]

    u = V * math.cos(alpha) * math.cos(beta)
    v = V * math.sin(beta)
    w = V * math.sin(alpha) * math.cos(beta)

    p = -Omega * math.sin(theta)
    q = Omega * math.sin(phi) * math.cos(theta)
    r = Omega * math.cos(phi) * math.cos(theta)

    p_hat = p * stab.Bref / (2.0 * V)
    q_hat = q * stab.Cref / (2.0 * V)
    r_hat = r * stab.Bref / (2.0 * V)

    qbar = 0.5 * rho * V * V
    turn_radius = math.inf if abs(Omega) < 1e-12 else V / Omega
    z_dot = (
        -u * math.sin(theta)
        + v * math.sin(phi) * math.cos(theta)
        + w * math.cos(phi) * math.cos(theta)
    )
    h_dot = -z_dot

    return {
        "u": u,
        "v": v,
        "w": w,
        "p": p,
        "q": q,
        "r": r,
        "p_hat": p_hat,
        "q_hat": q_hat,
        "r_hat": r_hat,
        "qbar": qbar,
        "turn_radius": turn_radius,
        "z_dot": z_dot,
        "h_dot": h_dot,
        "sink_rate": z_dot,
    }

def _linear_aero_coefficients(
    params: Mapping[str, float],
    derived: Mapping[str, float],
    stab: VSPAEROStab,
    control_columns: Mapping[str, str],
) -> dict:
    u_increment = 0.0
    if abs(stab.V0) > 1e-12:
        u_increment = (params["V"] - stab.V0) / stab.V0

    increments = {
        "Alpha": params["alpha"] - stab.alpha0,
        "Beta": params["beta"] - stab.beta0,
        "p": derived["p_hat"],
        "q": derived["q_hat"],
        "r": derived["r_hat"],
        "Mach": 0.0,
        "U": u_increment,
    }
    for delta_name, column_name in control_columns.items():
        increments[column_name] = params[delta_name]

    coefficients = {}
    for coef_name in AERO_COEFFICIENTS:
        row = stab.derivatives.loc[coef_name]
        value = float(row["Base"])
        for derivative_name, increment in increments.items():
            if derivative_name in row.index:
                value += float(row[derivative_name]) * float(increment)
        coefficients[coef_name] = value

    alpha = params["alpha"]
    beta = params["beta"]
    CL = coefficients["CL"]
    CD = coefficients["CD"]
    CS = coefficients["CS"]
    axial_opposite_forward = CD * math.cos(beta) + CS * math.sin(beta)

    coefficients["CX"] = -axial_opposite_forward * math.cos(alpha) + CL * math.sin(alpha)
    coefficients["CY"] = -CD * math.sin(beta) + CS * math.cos(beta)
    coefficients["CZ"] = -axial_opposite_forward * math.sin(alpha) - CL * math.cos(alpha)
    coefficients["Cl"] = coefficients["CMl"]
    coefficients["Cm"] = coefficients["CMm"]
    coefficients["Cn"] = coefficients["CMn"]
    return coefficients

def _steady_turn_residuals(
    params: Mapping[str, float],
    stab: VSPAEROStab,
    mass: float,
    rho: float,
    g: float,
    control_columns: Mapping[str, str],
) -> dict:
    V = params["V"]
    alpha = params["alpha"]
    beta = params["beta"]
    phi = params["phi"]
    theta = params["theta"]
    Omega = params["Omega"]
    T = params["T"]

    derived = _steady_turn_derived(params, stab, rho)
    coefficients = _linear_aero_coefficients(params, derived, stab, control_columns)

    qbar_s = derived["qbar"] * stab.Sref
    force_scale = qbar_s if abs(qbar_s) > 1e-12 else 1.0

    required_x = mass * (
        Omega * V * math.cos(theta) * (math.sin(phi) * math.sin(alpha) * math.cos(beta) - math.cos(phi) * math.sin(beta))
        + g * math.sin(theta)
    )
    required_y = mass * (
        Omega * V * math.cos(beta) * (math.cos(phi) * math.cos(theta) * math.cos(alpha) + math.sin(theta) * math.sin(alpha))
        - g * math.sin(phi) * math.cos(theta)
    )
    required_z = mass * (
        -Omega * V * math.sin(theta) * math.sin(beta)
        -Omega * V * math.sin(phi) * math.cos(theta) * math.cos(alpha) * math.cos(beta)
        -g * math.cos(phi) * math.cos(theta)
    )

    height_residual = (
        -math.cos(alpha) * math.cos(beta) * math.sin(theta)
        + math.sin(beta) * math.sin(phi) * math.cos(theta)
        + math.sin(alpha) * math.cos(beta) * math.cos(phi) * math.cos(theta)
    )

    return {
        "force_x": (T + qbar_s * coefficients["CX"] - required_x) / force_scale,
        "force_y": (qbar_s * coefficients["CY"] - required_y) / force_scale,
        "force_z": (qbar_s * coefficients["CZ"] - required_z) / force_scale,
        "moment_l": coefficients["Cl"],
        "moment_m": coefficients["Cm"],
        "moment_n": coefficients["Cn"],
        "height": height_residual,
    }

def solve_steady_level_turn(
    fixed: Mapping[str, float],
    stab_path: str | os.PathLike,
    mass: float,
    rho: float | None = None,
    initial_guess: Mapping[str, float] | None = None,
    *,
    g: float = 9.80665,
    bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    control_map: Mapping[str, str] | None = None,
    residual_tol: float = 1e-6,
    verbose: int | bool = 0,
) -> dict:
    """
    Solve a steady level turn by fixing any three of the ten turn parameters.

    Parameters
    ----------
    fixed
        Dict containing exactly three fixed parameters from TURN_PARAMETERS.
        Example: {"V": 30.0, "Omega": 0.08, "beta": 0.0}
    stab_path
        Path to a VSPAERO .stab file.
    mass
        Aircraft mass.  Use units consistent with rho, V, Sref, and g.
    rho
        Density used in the force-balance equations.  If None, the .stab Rho
        value is used as-is.  This does not convert units; the caller remains
        responsible for using a consistent unit system.
    initial_guess
        Optional guesses for unknown parameters.  Missing guesses are filled
        from the .stab base condition and simple neutral defaults.
    g
        Gravity acceleration.
    bounds
        Optional lower/upper bound dicts.  Only unknown parameter names matter.
    control_map
        Optional mapping from delta_a/e/r to a ConGrp_* column or a control
        group name in the .stab case table.  If omitted, AILERON/ELEVATOR/RUDDER
        group names are detected, with ConGrp_1/2/3 as fallback.
    residual_tol
        Maximum accepted absolute residual for passed=True.
    verbose
        If true, print the least_squares termination message.
    """

    fixed = dict(fixed)
    unknown_names = [name for name in TURN_PARAMETERS if name not in fixed]
    invalid_names = [name for name in fixed if name not in TURN_PARAMETERS]
    if invalid_names:
        raise ValueError(f"Unknown fixed parameter(s): {invalid_names}")
    if len(fixed) != 3:
        raise ValueError("fixed must contain exactly three turn parameters.")

    stab = read_vspaero_stab(stab_path)
    rho_used = _resolve_rho(stab, rho)
    control_columns = _control_columns_from_stab(stab, control_map)

    guess = _default_initial_guess(stab, fixed, mass, rho_used, g)
    if initial_guess:
        guess.update({name: float(value) for name, value in initial_guess.items() if name in TURN_PARAMETERS})
    guess.update({name: float(value) for name, value in fixed.items()})

    x0 = np.array([guess[name] for name in unknown_names], dtype=float)

    if bounds is None:
        lower = np.full(len(unknown_names), -np.inf)
        upper = np.full(len(unknown_names), np.inf)
    else:
        lower_dict, upper_dict = bounds
        lower = np.array([lower_dict.get(name, -np.inf) for name in unknown_names], dtype=float)
        upper = np.array([upper_dict.get(name, np.inf) for name in unknown_names], dtype=float)

    residual_names = ["force_x", "force_y", "force_z", "moment_l", "moment_m", "moment_n", "height"]

    def residual_vector(x: np.ndarray) -> np.ndarray:
        params = _complete_parameters(fixed, unknown_names, x)
        residuals = _steady_turn_residuals(params, stab, mass, rho_used, g, control_columns)
        return np.array([residuals[name] for name in residual_names], dtype=float)

    opt = least_squares(residual_vector, x0, bounds=(lower, upper))
    solution = _complete_parameters(fixed, unknown_names, opt.x)
    derived = _steady_turn_derived(solution, stab, rho_used)
    coefficients = _linear_aero_coefficients(solution, derived, stab, control_columns)
    residuals = _steady_turn_residuals(solution, stab, mass, rho_used, g, control_columns)

    max_abs_residual = max(abs(value) for value in residuals.values())
    passed = bool(opt.success and np.isfinite(max_abs_residual) and max_abs_residual <= residual_tol)

    if verbose:
        print(opt.message)
        print(f"max_abs_residual = {max_abs_residual:.6e}")

    return {
        "passed": passed,
        "message": opt.message,
        "fixed": fixed,
        "unknown_names": unknown_names,
        "solution": solution,
        "derived": derived,
        "coefficients": coefficients,
        "residuals": residuals,
        "max_abs_residual": max_abs_residual,
        "residual_tol": float(residual_tol),
        "cost": float(opt.cost),
        "optimality": float(opt.optimality),
        "nfev": int(opt.nfev),
        "rho": rho_used,
        "rho_source": ".stab" if rho is None else "argument",
        "control_columns": control_columns,
        "stab": {
            "path": stab.path,
            "references": stab.references,
            "base_condition": stab.base_condition,
            "base_aero": stab.base_aero,
            "control_groups": stab.control_groups,
        },
    }

def solve_steady_gliding_turn(
    fixed: Mapping[str, float],
    stab_path: str | os.PathLike,
    mass: float,
    rho: float | None = None,
    initial_guess: Mapping[str, float] | None = None,
    *,
    g: float = 9.80665,
    bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    control_map: Mapping[str, str] | None = None,
    residual_tol: float = 1e-6,
    verbose: int | bool = 0,
) -> dict:
    """
    Solve an unpowered steady gliding turn by fixing any three gliding parameters.

    This is the glider counterpart of solve_steady_level_turn().  The thrust is
    fixed internally as T = 0, and the height constraint is not included in the
    solved residual vector.  The resulting vertical motion is returned in
    derived["z_dot"], derived["h_dot"], and derived["sink_rate"].

    Parameters
    ----------
    fixed
        Dict containing exactly three fixed parameters from
        GLIDING_TURN_PARAMETERS.  T is not accepted because it is fixed to zero
        by the gliding-turn model.  Example: {"V": 30.0, "Omega": 0.08,
        "beta": 0.0}
    stab_path
        Path to a VSPAERO .stab file.
    mass
        Aircraft mass.  Use units consistent with rho, V, Sref, and g.
    rho
        Density used in the force-balance equations.  If None, the .stab Rho
        value is used as-is.  This does not convert units; the caller remains
        responsible for using a consistent unit system.
    initial_guess
        Optional guesses for unknown parameters.  Missing guesses are filled
        from the .stab base condition and simple neutral defaults.  For gliding
        turns, theta often needs a more deliberate initial guess than in a level
        turn because the height constraint is not imposed.
    g
        Gravity acceleration.
    bounds
        Optional lower/upper bound dicts.  Only unknown parameter names matter.
    control_map
        Optional mapping from delta_a/e/r to a ConGrp_* column or a control
        group name in the .stab case table.  If omitted, AILERON/ELEVATOR/RUDDER
        group names are detected, with ConGrp_1/2/3 as fallback.
    residual_tol
        Maximum accepted absolute residual for passed=True.
    verbose
        If true, print the least_squares termination message.
    """

    fixed = dict(fixed)
    unknown_names = [name for name in GLIDING_TURN_PARAMETERS if name not in fixed]
    invalid_names = [name for name in fixed if name not in GLIDING_TURN_PARAMETERS]
    if invalid_names:
        raise ValueError(f"Unknown fixed gliding-turn parameter(s): {invalid_names}")
    if len(fixed) != 3:
        raise ValueError("fixed must contain exactly three gliding-turn parameters.")

    stab = read_vspaero_stab(stab_path)
    rho_used = _resolve_rho(stab, rho)
    control_columns = _control_columns_from_stab(stab, control_map)

    guess = _default_initial_guess(stab, fixed, mass, rho_used, g)
    if initial_guess:
        guess.update({name: float(value) for name, value in initial_guess.items() if name in GLIDING_TURN_PARAMETERS})
    guess.update({name: float(value) for name, value in fixed.items()})
    guess["T"] = 0.0

    x0 = np.array([guess[name] for name in unknown_names], dtype=float)

    if bounds is None:
        lower = np.full(len(unknown_names), -np.inf)
        upper = np.full(len(unknown_names), np.inf)
    else:
        lower_dict, upper_dict = bounds
        lower = np.array([lower_dict.get(name, -np.inf) for name in unknown_names], dtype=float)
        upper = np.array([upper_dict.get(name, np.inf) for name in unknown_names], dtype=float)

    residual_names = ["force_x", "force_y", "force_z", "moment_l", "moment_m", "moment_n"]

    def residual_vector(x: np.ndarray) -> np.ndarray:
        params = _complete_parameters(fixed, unknown_names, x)
        params["T"] = 0.0
        residuals = _steady_turn_residuals(params, stab, mass, rho_used, g, control_columns)
        return np.array([residuals[name] for name in residual_names], dtype=float)

    opt = least_squares(residual_vector, x0, bounds=(lower, upper))
    solution = _complete_parameters(fixed, unknown_names, opt.x)
    solution["T"] = 0.0
    derived = _steady_turn_derived(solution, stab, rho_used)
    coefficients = _linear_aero_coefficients(solution, derived, stab, control_columns)
    all_residuals = _steady_turn_residuals(solution, stab, mass, rho_used, g, control_columns)
    residuals = {name: all_residuals[name] for name in residual_names}

    max_abs_residual = max(abs(value) for value in residuals.values())
    passed = bool(opt.success and np.isfinite(max_abs_residual) and max_abs_residual <= residual_tol)

    if verbose:
        print(opt.message)
        print(f"max_abs_residual = {max_abs_residual:.6e}")
        print(f"sink_rate = {derived['sink_rate']:.6e}")

    return {
        "passed": passed,
        "mode": "gliding_turn",
        "thrust_model": "T = 0",
        "message": opt.message,
        "fixed": fixed,
        "unknown_names": unknown_names,
        "solution": solution,
        "derived": derived,
        "coefficients": coefficients,
        "residuals": residuals,
        "height_residual": all_residuals["height"],
        "max_abs_residual": max_abs_residual,
        "residual_tol": float(residual_tol),
        "cost": float(opt.cost),
        "optimality": float(opt.optimality),
        "nfev": int(opt.nfev),
        "rho": rho_used,
        "rho_source": ".stab" if rho is None else "argument",
        "control_columns": control_columns,
        "stab": {
            "path": stab.path,
            "references": stab.references,
            "base_condition": stab.base_condition,
            "base_aero": stab.base_aero,
            "control_groups": stab.control_groups,
        },
    }

