"""
Steady turn trim solvers based on a VSPAERO .stab linear aero model.
For a powered steady level turn, the solver uses the ten turn parameters
    V, alpha, beta, phi, theta, Omega, delta_e, delta_a, delta_r, T
and solves the seven steady-level-turn equations after any three parameters
are fixed by the caller.
For an unpowered gliding turn, T is fixed to zero. The solver uses the nine
turn parameters
    V, alpha, beta, phi, theta, Omega, delta_e, delta_a, delta_r
and solves the six force/moment balance equations after any three parameters
are fixed by the caller. The height constraint is not imposed; the resulting
z_dot, h_dot, and sink_rate are returned as derived values.
Internal units are SI-style and radians unless the caller intentionally uses a
consistent non-SI unit system:
    V       : m/s, or the same length/time unit used with rho, mass, Sref
    alpha   : rad
    beta    : rad
    phi     : rad
    theta   : rad
    Omega   : rad/s, or rad/Tunit consistent with V
    delta_* : rad
    T       : N, or the force unit consistent with rho, V, Sref
    mass    : kg, or the mass unit consistent with rho and g
    rho     : kg/m^3, or density consistent with V, Sref, mass, g
If rho=None is passed to solve_steady_level_turn() or solve_steady_gliding_turn(),
the Rho value in the .stab file is used. This is convenient when the .stab file and the caller's mass,
length, and time units are consistent. It is not a unit conversion.
Coordinate conventions used in this module:
    Solver body axes : +x forward, +y right, +z down
    OpenVSP body axes: +X aft, +Y right, +Z up
The .stab CL/CD/CS coefficients are converted back to the solver body-axis
CX/CY/CZ using alpha and beta. The .stab CMl/CMm/CMn columns are used as
Cl/Cm/Cn. They are body-axis moments with aviation sign convention, not true
wind-axis moments. The p/q/r derivatives are reduced-rate derivatives:
    p_hat = p * Bref / (2 V)
    q_hat = q * Cref / (2 V)
    r_hat = r * Bref / (2 V)
Control derivatives are with respect to VSPAERO control-group deflection in
radians. The control group names in the .stab case table are used to map
ConGrp_1, ConGrp_2, ... to delta_a, delta_e, and delta_r.
"""
from __future__ import annotations
import math
import os
from typing import Mapping

import numpy as np
from scipy.optimize import least_squares

from .VSPAEROStab import (
    evaluate_stab_linear_aero,
    VSPAEROStab,
    read_vspaero_stab,
    resolve_control_columns_from_stab,
)

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


def _resolve_rho(stab: VSPAEROStab, rho: float | None) -> float:
    if rho is not None:
        return float(rho)
    if stab.rho0 is None:
        raise ValueError("rho=None was passed, but the .stab file does not contain Rho.")
    return float(stab.rho0)

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
    """Calculate velocities, angular rates, and flight-path geometry.

    The horizontal flight-path radius uses the horizontal speed, not the total
    airspeed.  In still air, the horizontal track angle differs from the heading
    by the constant ``track_heading_offset``.  Because that offset is constant
    in a steady turn, the track angle rotates at the same rate ``Omega``.
    """
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

    # These are the horizontal velocity components parallel and perpendicular
    # to the heading direction.  For arbitrary psi, the pair rotates by psi,
    # so its magnitude is the actual horizontal flight-path speed.
    track_forward_component = (
        u * math.cos(theta)
        + v * math.sin(phi) * math.sin(theta)
        + w * math.cos(phi) * math.sin(theta)
    )
    track_right_component = v * math.cos(phi) - w * math.sin(phi)
    horizontal_speed = math.hypot(track_forward_component, track_right_component)
    track_heading_offset = math.atan2(track_right_component, track_forward_component)
    horizontal_path_radius = (
        math.inf
        if abs(Omega) < 1.0e-12
        else horizontal_speed / abs(Omega)
    )

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
        "track_forward_component": track_forward_component,
        "track_right_component": track_right_component,
        "horizontal_speed": horizontal_speed,
        "track_heading_offset": track_heading_offset,
        "horizontal_path_radius": horizontal_path_radius,
        "z_dot": z_dot,
        "h_dot": h_dot,
        "sink_rate": z_dot,
    }

def _turn_aero_coefficients(
    params: Mapping[str, float],
    derived: Mapping[str, float],
    stab: VSPAEROStab,
    control_columns: Mapping[str, str],
) -> dict:
    """Build turn-state increments and evaluate the shared .stab aero model."""

    u_increment = 0.0
    if abs(stab.V0) > 1.0e-12:
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
    return evaluate_stab_linear_aero(
        stab,
        increments,
        alpha=params["alpha"],
        beta=params["beta"],
    )


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
    coefficients = _turn_aero_coefficients(params, derived, stab, control_columns)
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
        Density used in the force-balance equations. If None, the .stab Rho
        value is used as-is. This does not convert units; the caller remains
        responsible for using a consistent unit system.
    initial_guess
        Optional guesses for unknown parameters. Missing guesses are filled
        from the .stab base condition and simple neutral defaults.
    g
        Gravity acceleration.
    bounds
        Optional lower/upper bound dicts. Only unknown parameter names matter.
    control_map
        Optional mapping from delta_a/e/r to a ConGrp_* column or a control
        group name in the .stab case table. If omitted, AILERON/ELEVATOR/RUDDER
        group names are detected, with ConGrp_1/2/3 as fallback.
    residual_tol
        Maximum accepted absolute residual for passed=True.
    verbose
        If True, print the least_squares termination message.
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
    control_columns = resolve_control_columns_from_stab(stab, control_map)
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
    coefficients = _turn_aero_coefficients(solution, derived, stab, control_columns)
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
    solved residual vector. The resulting vertical motion is returned in
    derived["z_dot"], derived["h_dot"], and derived["sink_rate"].
    Parameters
    ----------
    fixed
        Dict containing exactly three fixed parameters from
        GLIDING_TURN_PARAMETERS.  T is not accepted because it is fixed to zero
        by the gliding-turn model. Example: {"V": 30.0, "Omega": 0.08,
        "beta": 0.0}
    stab_path
        Path to a VSPAERO .stab file.
    mass
        Aircraft mass. Use units consistent with rho, V, Sref, and g.
    rho
        Density used in the force-balance equations. If None, the .stab Rho
        value is used as-is. This does not convert units; the caller remains
        responsible for using a consistent unit system.
    initial_guess
        Optional guesses for unknown parameters. Missing guesses are filled
        from the .stab base condition and simple neutral defaults. For gliding
        turns, theta often needs a more deliberate initial guess than in a level
        turn because the height constraint is not imposed.
    g
        Gravity acceleration.
    bounds
        Optional lower/upper bound dicts. Only unknown parameter names matter.
    control_map
        Optional mapping from delta_a/e/r to a ConGrp_* column or a control
        group name in the .stab case table. If omitted, AILERON/ELEVATOR/RUDDER
        group names are detected, with ConGrp_1/2/3 as fallback.
    residual_tol
        Maximum accepted absolute residual for passed=True.
    verbose
        If True, print the least_squares termination message.
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
    control_columns = resolve_control_columns_from_stab(stab, control_map)
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
    coefficients = _turn_aero_coefficients(solution, derived, stab, control_columns)
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

def solve_rudder_limit_turn(
    stab_path: str | os.PathLike,
    mass: float,
    delta_r_max: float,
    *,
    mode: str = "gliding",
    V: float | None = None,
    delta_a: float = 0.0,
    rho: float | None = None,
    initial_guess: Mapping[str, float] | None = None,
    g: float = 9.80665,
    bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    control_map: Mapping[str, str] | None = None,
    residual_tol: float = 1e-6,
    verbose: int | bool = 0,
) -> dict:
    """Solve both rudder-limit turns and select the largest absolute bank angle.

    The solved endpoint cases are

        V = V0, delta_a = constant, delta_r = -delta_r_max
        V = V0, delta_a = constant, delta_r = +delta_r_max

    where V0 is read from the .stab file unless V is supplied explicitly.  Each
    endpoint is solved with the existing steady level- or gliding-turn solver.
    The accepted endpoint with the largest ``abs(phi)`` is returned through the
    ordinary top-level ``solution``, ``derived``, ``coefficients``, and
    ``residuals`` entries.  The two complete endpoint results remain available
    as ``negative_trim`` and ``positive_trim``.

    This is an exact numerical solution of the nonlinear steady-turn equations
    used by this module, with the .stab linear aerodynamic model.  It is not a
    closed-form analytic solution and does not prove monotonicity between zero
    rudder and the two rudder limits.
    """

    if mode not in {"level", "gliding"}:
        raise ValueError("mode must be 'level' or 'gliding'.")

    delta_r_max = float(delta_r_max)
    if not math.isfinite(delta_r_max) or delta_r_max <= 0.0:
        raise ValueError("delta_r_max must be positive and finite.")

    delta_a = float(delta_a)
    if not math.isfinite(delta_a):
        raise ValueError("delta_a must be finite.")

    stab = read_vspaero_stab(stab_path)
    V_used = float(stab.V0 if V is None else V)
    if not math.isfinite(V_used) or V_used <= 0.0:
        raise ValueError("V must be positive, or the .stab file must contain a positive Vinf.")

    solver = solve_steady_level_turn if mode == "level" else solve_steady_gliding_turn
    trim_verbose = 1 if verbose and int(verbose) >= 2 else 0
    endpoint_trims: dict[str, dict] = {}

    for side, delta_r in (("negative", -delta_r_max), ("positive", delta_r_max)):
        fixed = {"V": V_used, "delta_a": delta_a, "delta_r": delta_r}
        try:
            endpoint_trims[side] = solver(
                fixed,
                stab_path,
                mass,
                rho=rho,
                initial_guess=initial_guess,
                g=g,
                bounds=bounds,
                control_map=control_map,
                residual_tol=residual_tol,
                verbose=trim_verbose,
            )
        except Exception as exc:
            endpoint_trims[side] = {
                "passed": False,
                "message": repr(exc),
                "fixed": fixed,
                "solution": {},
                "derived": {},
                "coefficients": {},
                "residuals": {},
                "max_abs_residual": math.nan,
                "cost": math.nan,
                "nfev": 0,
            }

    accepted = [
        (side, trim)
        for side, trim in endpoint_trims.items()
        if bool(trim.get("passed", False))
        and math.isfinite(float((trim.get("solution", {}) or {}).get("phi", math.nan)))
    ]
    both_sides_passed = len(accepted) == 2

    if not accepted:
        return {
            "passed": False,
            "message": "Neither rudder-limit endpoint produced an accepted turn-trim solution.",
            "mode": f"{mode}_rudder_limit_turn",
            "rudder_limit_mode": mode,
            "rudder_limit_complete": False,
            "selected_side": "",
            "delta_r_max": delta_r_max,
            "limiting_delta_r": math.nan,
            "max_abs_phi": math.nan,
            "fixed_V": V_used,
            "fixed_delta_a": delta_a,
            "solution": {},
            "derived": {},
            "coefficients": {},
            "residuals": {},
            "max_abs_residual": math.nan,
            "cost": math.nan,
            "nfev": 0,
            "negative_trim": endpoint_trims["negative"],
            "positive_trim": endpoint_trims["positive"],
        }

    selected_side, selected_trim = max(
        accepted,
        key=lambda item: abs(float(item[1]["solution"]["phi"])),
    )
    result = dict(selected_trim)
    result["solver_message"] = selected_trim.get("message", "")
    result["message"] = (
        f"Selected the {selected_side} rudder-limit endpoint from "
        f"{len(accepted)} accepted endpoint solution(s)."
    )
    result.update(
        {
            "passed": True,
            "mode": f"{mode}_rudder_limit_turn",
            "rudder_limit_mode": mode,
            "rudder_limit_complete": both_sides_passed,
            "selected_side": selected_side,
            "delta_r_max": delta_r_max,
            "limiting_delta_r": float(selected_trim["solution"]["delta_r"]),
            "max_abs_phi": abs(float(selected_trim["solution"]["phi"])),
            "fixed_V": V_used,
            "fixed_delta_a": delta_a,
            "negative_trim": endpoint_trims["negative"],
            "positive_trim": endpoint_trims["positive"],
        }
    )
    return result

