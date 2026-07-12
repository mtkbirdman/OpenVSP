"""
Rudder roll-rate gain validation utilities.

This module compares quasi-steady rudder roll-rate gains from a VSPAERO .stab
file with finite-time gains from a 6DOF nonlinear rigid-body simulation using
the same .stab linear aerodynamic model.

Angles are radians internally.  The .stab control derivatives are assumed to be
per radian, and p/q/r derivatives are assumed to be reduced-rate derivatives:

    phat = p * Bref / (2 V)
    qhat = q * Cref / (2 V)
    rhat = r * Bref / (2 V)

The caller must pass mass and inertia values in a unit system consistent with
Sref, Cref, Bref, Vinf, Rho, and g used for the .stab case.

For the 6DOF equations, .stab CFx/CFy/CFz are not used directly as solver
body-axis forces.  The aerodynamic force coefficients are evaluated from
CL/CD/CS and converted to the solver body axes (+x forward, +y right, +z down)
using the instantaneous alpha and beta.
"""
from __future__ import annotations

import importlib
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from .TrimTurnSolver import (
    resolve_control_columns_from_stab,
    read_vspaero_stab,
)

def _stab_value(stab, coef: str, column: str) -> float:
    if coef not in stab.derivatives.index:
        raise KeyError(f"Coefficient row {coef!r} was not found in {stab.path}.")
    if column not in stab.derivatives.columns:
        raise KeyError(f"Derivative column {column!r} was not found in {stab.path}.")
    return float(stab.derivatives.loc[coef, column])

def _evaluate_stab_linear_coefficients(
    stab,
    increments: Mapping[str, float],
    *,
    alpha: float,
    beta: float,
) -> dict[str, float]:
    """
    Evaluate .stab linear aero and convert wind-axis forces to solver body axes.

    .stab CFx/CFy/CFz are retained only as raw diagnostic values. The 6DOF
    solver uses CX/CY/CZ built from CL/CD/CS for the solver body axes:
    +x forward, +y right, +z down.
    """
    coefficients: dict[str, float] = {}
    for coef_name in ["CL", "CD", "CS", "CMl", "CMm", "CMn"]:
        row = stab.derivatives.loc[coef_name]
        value = float(row["Base"])
        for derivative_name, increment in increments.items():
            if derivative_name in row.index:
                value += float(row[derivative_name]) * float(increment)
        coefficients[coef_name] = value

    for coef_name in ["CFx", "CFy", "CFz"]:
        if coef_name in stab.derivatives.index:
            row = stab.derivatives.loc[coef_name]
            value = float(row["Base"])
            for derivative_name, increment in increments.items():
                if derivative_name in row.index:
                    value += float(row[derivative_name]) * float(increment)
            coefficients[f"{coef_name}_raw"] = value

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

def _solver_axis_derivative_value(stab, coef: str, column: str) -> float:
    """Return a derivative for solver-axis coefficients, including derived CY."""
    alias = {"Cl": "CMl", "Cm": "CMm", "Cn": "CMn"}
    if coef in alias:
        return _stab_value(stab, alias[coef], column)
    if coef in {"CL", "CD", "CS", "CMl", "CMm", "CMn", "CFx", "CFy", "CFz"}:
        return _stab_value(stab, coef, column)
    if coef not in {"CX", "CY", "CZ"}:
        raise KeyError(f"Unsupported solver-axis coefficient {coef!r}.")

    eps = 1.0e-6
    increments: dict[str, float] = {column: eps}
    alpha = stab.alpha0 + (eps if column == "Alpha" else 0.0)
    beta = stab.beta0 + (eps if column == "Beta" else 0.0)
    base = _evaluate_stab_linear_coefficients(stab, {}, alpha=stab.alpha0, beta=stab.beta0)[coef]
    perturbed = _evaluate_stab_linear_coefficients(stab, increments, alpha=alpha, beta=beta)[coef]
    return float((perturbed - base) / eps)

def solver_axis_derivative_value_from_stab(stab, coef: str, column: str) -> float:
    """Return a .stab derivative in the solver-axis coefficient convention.

    This is the public wrapper used by higher-level post-processing code.
    For CY/CX/CZ it returns the derivative after converting CL/CD/CS to the
    solver body axes (+x forward, +y right, +z down). For Cl/Cm/Cn it maps to
    the VSPAERO CMl/CMm/CMn rows.
    """
    return _solver_axis_derivative_value(stab, coef, column)

def calculate_quasi_steady_rudder_roll_gain_from_stab(
    stab_path: str | Path,
    *,
    control_map: Mapping[str, str] | None = None,
    side_force_coef: str = "CY",
    roll_moment_coef: str = "CMl",
    yaw_moment_coef: str = "CMn",
) -> dict[str, Any]:
    """Calculate quasi-steady rudder roll-rate gains from a VSPAERO .stab file."""
    stab = read_vspaero_stab(stab_path)
    control_columns = resolve_control_columns_from_stab(stab, control_map)
    rudder_column = control_columns.get("delta_r")
    if rudder_column is None:
        raise KeyError("Could not resolve rudder control column. Pass control_map={'delta_r': 'ConGrp_N'} explicitly.")

    matrix = np.array(
        [
            [
                _solver_axis_derivative_value(stab, side_force_coef, "Beta"),
                _solver_axis_derivative_value(stab, side_force_coef, "p"),
                _solver_axis_derivative_value(stab, side_force_coef, "r"),
            ],
            [
                _solver_axis_derivative_value(stab, roll_moment_coef, "Beta"),
                _solver_axis_derivative_value(stab, roll_moment_coef, "p"),
                _solver_axis_derivative_value(stab, roll_moment_coef, "r"),
            ],
            [
                _solver_axis_derivative_value(stab, yaw_moment_coef, "Beta"),
                _solver_axis_derivative_value(stab, yaw_moment_coef, "p"),
                _solver_axis_derivative_value(stab, yaw_moment_coef, "r"),
            ],
        ],
        dtype=float,
    )
    rhs = -np.array(
        [
            _solver_axis_derivative_value(stab, side_force_coef, rudder_column),
            _solver_axis_derivative_value(stab, roll_moment_coef, rudder_column),
            _solver_axis_derivative_value(stab, yaw_moment_coef, rudder_column),
        ],
        dtype=float,
    )
    determinant = float(np.linalg.det(matrix))
    if not np.isfinite(determinant) or abs(determinant) < 1.0e-14:
        raise np.linalg.LinAlgError(f"Rudder coupled-gain matrix is singular or ill-conditioned. det={determinant}")

    beta_per_delta_r, phat_per_delta_r, rhat_per_delta_r = np.linalg.solve(matrix, rhs)
    p_per_delta_r = 2.0 * stab.V0 / stab.Bref * phat_per_delta_r
    return {
        "stab_path": str(stab_path),
        "control_columns": control_columns,
        "rudder_column": rudder_column,
        "side_force_coef": side_force_coef,
        "roll_moment_coef": roll_moment_coef,
        "yaw_moment_coef": yaw_moment_coef,
        "matrix_det": determinant,
        "beta_per_delta_r": float(beta_per_delta_r),
        "phat_per_delta_r": float(phat_per_delta_r),
        "rhat_per_delta_r": float(rhat_per_delta_r),
        "p_per_delta_r": float(p_per_delta_r),
        "V": float(stab.V0),
        "Sref": float(stab.Sref),
        "Cref": float(stab.Cref),
        "Bref": float(stab.Bref),
        "rho": None if stab.rho0 is None else float(stab.rho0),
    }

def calculate_linear_lateral_response_indices_from_stab(
    stab_path: str | Path,
    *,
    mass: float,
    Ixx: float,
    Izz: float,
    Ixz: float = 0.0,
    delta_r: float,
    t_final: float = 5.0,
    target_delta_phi: float = math.radians(5.0),
    taylor_tau_eval: float = 1.0,
    rho: float | None = None,
    g: float = 9.80665,
    control_map: Mapping[str, str] | None = None,
    side_force_coef: str = "CY",
    roll_moment_coef: str = "CMl",
    yaw_moment_coef: str = "CMn",
    max_step: float | None = None,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
) -> dict[str, Any]:
    """Calculate linear finite-time lateral response to a rudder step.

    The model uses the four small-disturbance lateral states

        x = [beta, phat, rhat, phi]

    with nondimensional time tau = 2 V / b * t.  It is the linear counterpart
    of the finite-time 6DOF signed-delta-phi index used in the Vv-Gamma
    workflow.  The same function also returns the initial derivative indices
    because they are simply the first Taylor coefficients of this linear model.
    """
    stab = read_vspaero_stab(stab_path)
    control_columns = resolve_control_columns_from_stab(stab, control_map)
    rudder_column = control_columns.get("delta_r")
    if rudder_column is None:
        raise KeyError("Could not resolve rudder control column. Pass control_map={'delta_r': 'ConGrp_N'} explicitly.")

    rho_used = float(stab.rho0 if rho is None else rho)
    if not math.isfinite(rho_used) or rho_used <= 0.0:
        raise ValueError("rho must be positive, or the .stab file must contain a positive Rho value.")

    V = float(stab.V0)
    if V <= 0.0:
        raise ValueError("The .stab file must contain a positive Vinf for linear lateral response indices.")
    mass = float(mass)
    if mass <= 0.0:
        raise ValueError("mass must be positive.")
    t_final = float(t_final)
    if t_final <= 0.0:
        raise ValueError("t_final must be positive.")
    delta_r = float(delta_r)
    if abs(delta_r) < 1.0e-14:
        raise ValueError("delta_r must be non-zero.")
    target_delta_phi = float(target_delta_phi)
    if abs(target_delta_phi) < 1.0e-14:
        raise ValueError("target_delta_phi must be non-zero.")

    Ixx = float(Ixx)
    Izz = float(Izz)
    Ixz = float(Ixz)
    inertia_det = Ixx * Izz - Ixz * Ixz
    if inertia_det <= 0.0:
        raise ValueError("Ixx * Izz - Ixz**2 must be positive.")

    Sref = float(stab.Sref)
    Bref = float(stab.Bref)
    qbar = 0.5 * rho_used * V * V
    tau_per_second = 2.0 * V / Bref
    second_per_tau = 1.0 / tau_per_second
    tau_final = tau_per_second * t_final
    mu_y = qbar * Sref * Bref / (2.0 * mass * V * V)
    mu_inertia = qbar * Sref * Bref**3 / (4.0 * V * V * inertia_det)
    lambda_g = float(g) * Bref / (2.0 * V * V)
    alpha0 = float(stab.alpha0)

    CY_beta = _solver_axis_derivative_value(stab, side_force_coef, "Beta")
    CY_phat = _solver_axis_derivative_value(stab, side_force_coef, "p")
    CY_rhat = _solver_axis_derivative_value(stab, side_force_coef, "r")
    CY_delta_r = _solver_axis_derivative_value(stab, side_force_coef, rudder_column)

    Cl_beta = _solver_axis_derivative_value(stab, roll_moment_coef, "Beta")
    Cl_phat = _solver_axis_derivative_value(stab, roll_moment_coef, "p")
    Cl_rhat = _solver_axis_derivative_value(stab, roll_moment_coef, "r")
    Cl_delta_r = _solver_axis_derivative_value(stab, roll_moment_coef, rudder_column)

    Cn_beta = _solver_axis_derivative_value(stab, yaw_moment_coef, "Beta")
    Cn_phat = _solver_axis_derivative_value(stab, yaw_moment_coef, "p")
    Cn_rhat = _solver_axis_derivative_value(stab, yaw_moment_coef, "r")
    Cn_delta_r = _solver_axis_derivative_value(stab, yaw_moment_coef, rudder_column)

    A = np.array(
        [
            [
                mu_y * CY_beta,
                mu_y * CY_phat + alpha0,
                mu_y * CY_rhat - 1.0,
                lambda_g,
            ],
            [
                mu_inertia * (Izz * Cl_beta + Ixz * Cn_beta),
                mu_inertia * (Izz * Cl_phat + Ixz * Cn_phat),
                mu_inertia * (Izz * Cl_rhat + Ixz * Cn_rhat),
                0.0,
            ],
            [
                mu_inertia * (Ixz * Cl_beta + Ixx * Cn_beta),
                mu_inertia * (Ixz * Cl_phat + Ixx * Cn_phat),
                mu_inertia * (Ixz * Cl_rhat + Ixx * Cn_rhat),
                0.0,
            ],
            [0.0, 1.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    B = np.array(
        [
            mu_y * CY_delta_r,
            mu_inertia * (Izz * Cl_delta_r + Ixz * Cn_delta_r),
            mu_inertia * (Ixz * Cl_delta_r + Ixx * Cn_delta_r),
            0.0,
        ],
        dtype=float,
    )

    beta_prime = float(B[0])
    phat_prime = float(B[1])
    rhat_prime = float(B[2])
    phat_prime_direct_roll_part = float(mu_inertia * Izz * Cl_delta_r)
    phat_prime_inertia_cross_yaw_part = float(mu_inertia * Ixz * Cn_delta_r)
    rhat_prime_direct_roll_cross_part = float(mu_inertia * Ixz * Cl_delta_r)
    rhat_prime_yaw_part = float(mu_inertia * Ixx * Cn_delta_r)

    phat_double_prime_beta_part = float(A[1, 0] * B[0])
    phat_double_prime_phat_part = float(A[1, 1] * B[1])
    phat_double_prime_rhat_part = float(A[1, 2] * B[2])
    phat_double_prime = float(A[1, :] @ B)

    # Exact Taylor coefficients of the same 4-state linear model.
    #
    #   Delta phi / delta_r
    #       = tau^2/2  * e_phi^T A B
    #       + tau^3/6  * e_phi^T A^2 B
    #       + tau^4/24 * e_phi^T A^3 B + ...
    #
    # These columns are short-time algebraic explanation terms.  They must not
    # be evaluated at the long 6DOF/linear response t_final, because the Taylor
    # series is a local expansion about tau = 0.
    AB = A @ B
    A2B = A @ AB
    A3B = A @ A2B
    linear_taylor_K1 = float(AB[3])
    linear_taylor_K2 = float(A2B[3])
    linear_taylor_K3 = float(A3B[3])
    # Simplified algebraic explanation terms.
    #
    # Assumptions used only in these simple_* columns:
    #   Ixz = 0
    #   Cl_delta_r = 0
    #   CY_phat = 0, CY_rhat = 0
    #   Cn_phat = 0
    #   alpha0 = 0
    # The kinematic beta' term -rhat is retained.
    #
    # Under these assumptions the dominant third-order bank-angle build-up is
    # split into two readable paths:
    #   rudder -> side force -> beta -> dihedral roll
    #   rudder -> yawing moment -> rhat -> yaw-rate roll
    simple_taylor_assumptions = (
        "Ixz=0; Cl_delta_r=0 for K2/K3 reduced paths; "
        "CY_phat=0; CY_rhat=0; Cn_phat=0; alpha0=0; "
        "retain beta_prime_kinematic_-rhat"
    )
    simple_K1_direct_roll = float(mu_inertia * Izz * Cl_delta_r)
    simple_K2_sideforce_dihedral = float(mu_inertia * mu_y * Izz * Cl_beta * CY_delta_r)
    simple_K2_yawrate_roll = float(mu_inertia**2 * Ixx * Izz * Cl_rhat * Cn_delta_r)
    simple_K2_total = simple_K2_sideforce_dihedral + simple_K2_yawrate_roll

    simple_K3_yawrate_to_beta_dihedral = float(
        -mu_inertia**2 * Ixx * Izz * Cl_beta * Cn_delta_r
    )

    # Minimal yaw-rate-roll / yaw-rate-to-beta-dihedral explanation terms.
    # These columns deliberately do not use an evaluation time such as tau_eff.
    # They are pure Taylor coefficients of the reduced path model:
    #   phi/delta_r = K2*tau**3/6 + K3*tau**4/24 + O(tau**5)
    simple_K2_plus_K3_yawrate_roll_beta_dihedral = float(
        simple_K2_yawrate_roll + simple_K3_yawrate_to_beta_dihedral
    )
    simple_roll_damping_rate = float(mu_inertia * Izz * Cl_phat)
    simple_yaw_damping_rate = float(mu_inertia * Ixx * Cn_rhat)
    simple_K3_yawrate_roll_damping_correction = float(
        simple_K2_yawrate_roll * (simple_roll_damping_rate + simple_yaw_damping_rate)
    )
    simple_K3_yawrate_beta_dihedral_damped = float(
        simple_K3_yawrate_to_beta_dihedral + simple_K3_yawrate_roll_damping_correction
    )
    simple_K3_reduced = float(
        simple_K3_yawrate_to_beta_dihedral
        + simple_roll_damping_rate * simple_K2_total
        + simple_yaw_damping_rate * simple_K2_yawrate_roll
    )

    # Normalized rudder-only turn tendency proxies.
    #
    # These are not turn-rate solutions.  They are dimensionless signed
    # coupling indicators obtained by normalizing the simple Taylor path
    # coefficients by mu_inertia**2 * Ixx * Izz:
    #
    #   simple_turn_rate_full
    #       = Cl_delta_r/(mu_inertia*Ixx)
    #       + mu_y*Cl_beta*CY_delta_r/(mu_inertia*Ixx)
    #       + (Cl_rhat - Cl_beta)*Cn_delta_r
    #
    #   simple_turn_rate
    #       = (Cl_rhat - Cl_beta)*Cn_delta_r
    #
    # The latter keeps only the yawing-moment -> yaw-rate -> roll and
    # yawing-moment -> yaw-rate -> beta -> dihedral-roll paths.
    simple_turn_rate_full = float(
        (
            simple_K1_direct_roll
            + simple_K2_sideforce_dihedral
            + simple_K2_yawrate_roll
            + simple_K3_yawrate_to_beta_dihedral
        )
        / (mu_inertia**2 * Ixx * Izz)
    )
    simple_turn_rate = float(
        simple_K2_plus_K3_yawrate_roll_beta_dihedral
        / (mu_inertia**2 * Ixx * Izz)
    )

    def rhs(tau: float, state: np.ndarray) -> np.ndarray:
        return A @ state + B * delta_r

    sign = 1.0 if target_delta_phi > 0.0 else -1.0

    def target_delta_phi_event(tau: float, state: np.ndarray) -> float:
        return sign * (float(state[3]) - target_delta_phi)

    target_delta_phi_event.terminal = False
    target_delta_phi_event.direction = 1.0

    kwargs: dict[str, Any] = {"rtol": rtol, "atol": atol, "events": target_delta_phi_event}
    if max_step is not None:
        kwargs["max_step"] = float(max_step) * tau_per_second

    solution = solve_ivp(rhs, (0.0, tau_final), np.zeros(4, dtype=float), **kwargs)
    if not solution.success:
        raise RuntimeError(solution.message)

    beta_final, phat_final, rhat_final, phi_final = [float(value) for value in solution.y[:, -1]]
    tau_actual = float(solution.t[-1])
    t_actual = tau_actual * second_per_tau
    reference_index = phi_final / tau_actual / delta_r
    reference_phi_rate = phi_final / t_actual

    reached = bool(solution.t_events and len(solution.t_events[0]) > 0)
    tau_reach = math.nan
    t_reach = math.nan
    target_index = math.nan
    target_phi_rate = math.nan
    target_phi_rate_per_delta_r = math.nan
    error = "target_delta_phi_not_reached"
    if reached:
        tau_reach = float(solution.t_events[0][0])
        t_reach = tau_reach * second_per_tau
        if tau_reach <= 0.0:
            error = "target_reached_at_or_before_start"
        else:
            target_index = target_delta_phi / tau_reach / delta_r
            target_phi_rate = target_delta_phi / t_reach
            target_phi_rate_per_delta_r = target_phi_rate / delta_r
            error = ""

    result: dict[str, Any] = {
        "linear_response_source": "small_disturbance_lateral_finite_time_rudder_step",
        "linear_success": True,
        "linear_message": "",
        "linear_rudder_column": rudder_column,
        "linear_delta_r": delta_r,
        "linear_target_delta_phi": target_delta_phi,
        "linear_target_delta_phi_deg": math.degrees(target_delta_phi),
        "linear_t_start": 0.0,
        "linear_tau_start": 0.0,
        "linear_t_final_requested": t_final,
        "linear_tau_final_requested": tau_final,
        "linear_t_final": t_actual,
        "linear_tau_final": tau_actual,
        "linear_tau_per_second": tau_per_second,
        "linear_second_per_tau": second_per_tau,
        "linear_index_Vinf": V,
        "linear_index_Bref": Bref,
        "linear_phi0": 0.0,
        "linear_phi_final": phi_final,
        "linear_phi_delta_final": phi_final,
        "linear_beta_final": beta_final,
        "linear_phat_final": phat_final,
        "linear_rhat_final": rhat_final,
        "linear_roll_response_reached": reached,
        "linear_t_reach": t_reach,
        "linear_tau_reach": tau_reach,
        "linear_dt_reach": t_reach,
        "linear_roll_response_phi_rate": target_phi_rate,
        "linear_roll_response_phi_rate_per_delta_r": target_phi_rate_per_delta_r,
        "linear_finite_time_roll_index": target_index,
        "linear_roll_response_error": error,
        "linear_roll_response_reference_phi_rate_to_t_final": reference_phi_rate,
        "linear_roll_response_reference_phi_rate_per_delta_r_to_t_final": reference_phi_rate / delta_r,
        "linear_roll_response_index_reference": reference_index,
        "linear_roll_response_fraction_of_target": phi_final / target_delta_phi,
        "linear_taylor_source": "exact_taylor_coefficients_from_full_linear_A_B",
        "linear_taylor_K1": linear_taylor_K1,
        "linear_taylor_K2": linear_taylor_K2,
        "linear_taylor_K3": linear_taylor_K3,
        "simple_taylor_source": "simplified_lateral_response_path_coefficients",
        "simple_taylor_assumptions": simple_taylor_assumptions,
        "simple_K1_direct_roll": simple_K1_direct_roll,
        "simple_K2_sideforce_dihedral": simple_K2_sideforce_dihedral,
        "simple_K2_yawrate_roll": simple_K2_yawrate_roll,
        "simple_K2_total": simple_K2_total,
        "simple_K3_yawrate_to_beta_dihedral": simple_K3_yawrate_to_beta_dihedral,
        "simple_K2_plus_K3_yawrate_roll_beta_dihedral": simple_K2_plus_K3_yawrate_roll_beta_dihedral,
        "simple_roll_damping_rate": simple_roll_damping_rate,
        "simple_yaw_damping_rate": simple_yaw_damping_rate,
        "simple_K3_yawrate_roll_damping_correction": simple_K3_yawrate_roll_damping_correction,
        "simple_K3_yawrate_beta_dihedral_damped": simple_K3_yawrate_beta_dihedral_damped,
        "simple_K3_reduced": simple_K3_reduced,
        "simple_turn_rate_full": simple_turn_rate_full,
        "simple_turn_rate": simple_turn_rate,
        "initial_response_source": "small_disturbance_lateral_taylor_terms_from_linear_model",
        "initial_tau_per_second": tau_per_second,
        "initial_second_per_tau": second_per_tau,
        "initial_mu_y": float(mu_y),
        "initial_mu_inertia": float(mu_inertia),
        "initial_inertia_det": float(inertia_det),
        "initial_rho": rho_used,
        "initial_V": V,
        "initial_Sref": Sref,
        "initial_Bref": Bref,
        "initial_rudder_column": rudder_column,
        "initial_beta_prime_per_delta_r": beta_prime,
        "initial_phat_prime_per_delta_r": phat_prime,
        "initial_rhat_prime_per_delta_r": rhat_prime,
        "initial_phat_double_prime_per_delta_r": phat_double_prime,
        "initial_phi_double_prime_per_delta_r": phat_prime,
        "initial_phi_triple_prime_per_delta_r": phat_double_prime,
        "initial_phat_prime_direct_roll_part": phat_prime_direct_roll_part,
        "initial_phat_prime_inertia_cross_yaw_part": phat_prime_inertia_cross_yaw_part,
        "initial_rhat_prime_direct_roll_cross_part": rhat_prime_direct_roll_cross_part,
        "initial_rhat_prime_yaw_part": rhat_prime_yaw_part,
        "initial_phat_double_prime_beta_part": phat_double_prime_beta_part,
        "initial_phat_double_prime_phat_part": phat_double_prime_phat_part,
        "initial_phat_double_prime_rhat_part": phat_double_prime_rhat_part,
        "initial_CY_beta": float(CY_beta),
        "initial_CY_phat": float(CY_phat),
        "initial_CY_rhat": float(CY_rhat),
        "initial_CY_delta_r": float(CY_delta_r),
        "initial_Cl_beta": float(Cl_beta),
        "initial_Cl_phat": float(Cl_phat),
        "initial_Cl_rhat": float(Cl_rhat),
        "initial_Cl_delta_r": float(Cl_delta_r),
        "initial_Cn_beta": float(Cn_beta),
        "initial_Cn_phat": float(Cn_phat),
        "initial_Cn_rhat": float(Cn_rhat),
        "initial_Cn_delta_r": float(Cn_delta_r),
    }

    state_names = ["beta", "phat", "rhat", "phi"]
    for row_index, row_name in enumerate(state_names):
        for col_index, col_name in enumerate(state_names):
            result[f"linear_A_{row_name}_{col_name}"] = float(A[row_index, col_index])
        result[f"linear_B_{row_name}"] = float(B[row_index])
    return result

def _initial_state_from_stab(
    stab,
    *,
    phi0: float = 0.0,
    theta0: float | None = None,
    psi0: float = 0.0,
    position0: Sequence[float] = (0.0, 0.0, 0.0),
) -> np.ndarray:
    V = float(stab.V0)
    alpha = float(stab.alpha0)
    beta = float(stab.beta0)
    theta = alpha if theta0 is None else float(theta0)
    u = V * math.cos(alpha) * math.cos(beta)
    v = V * math.sin(beta)
    w = V * math.sin(alpha) * math.cos(beta)
    x_e, y_e, z_e = [float(value) for value in position0]
    return np.array([u, v, w, 0.0, 0.0, 0.0, x_e, y_e, z_e, float(phi0), theta, float(psi0)], dtype=float)

def _linear_aero_from_stab(
    stab,
    state: Sequence[float],
    controls: Mapping[str, float],
    control_columns: Mapping[str, str],
    *,
    gust_body_velocity: Sequence[float] = (0.0, 0.0, 0.0),
) -> dict[str, float]:
    """Evaluate solver-axis aero coefficients from .stab linear derivatives.

    state[:3] is the aircraft velocity in solver body axes.
    gust_body_velocity is the local atmospheric velocity vector in the same
    axes.  Aerodynamic alpha, beta, reduced rates, and dynamic pressure are
    evaluated from the relative air velocity

        [u_air, v_air, w_air] = [u, v, w] - gust_body_velocity.

    Crosswind-facing sign conventions are handled by the caller.  In
    simulate_6dof_crosswind_gust_from_stab(), a positive crosswind input means
    wind from the aircraft right side and is converted to a negative
    atmospheric y-body velocity before this function is called.
    """
    u, v, w, p, q, r = [float(value) for value in state[:6]]
    u_g, v_g, w_g = [float(value) for value in gust_body_velocity]
    u_air = u - u_g
    v_air = v - v_g
    w_air = w - w_g

    V = max(math.sqrt(u_air * u_air + v_air * v_air + w_air * w_air), 1.0e-12)
    alpha = math.atan2(w_air, u_air)
    beta = math.asin(max(-1.0, min(1.0, v_air / V)))
    increments: dict[str, float] = {
        "Alpha": alpha - stab.alpha0,
        "Beta": beta - stab.beta0,
        "p": p * stab.Bref / (2.0 * V),
        "q": q * stab.Cref / (2.0 * V),
        "r": r * stab.Bref / (2.0 * V),
        "Mach": 0.0,
        "U": 0.0 if abs(stab.V0) < 1.0e-12 else (V - stab.V0) / stab.V0,
    }
    for delta_name, column_name in control_columns.items():
        if column_name in stab.derivatives.columns:
            increments[column_name] = float(controls.get(delta_name, 0.0))
    coefficients = _evaluate_stab_linear_coefficients(stab, increments, alpha=alpha, beta=beta)
    coefficients.update(
        {
            "V": V,
            "alpha": alpha,
            "beta": beta,
            "u_air": u_air,
            "v_air": v_air,
            "w_air": w_air,
            "u_gust": u_g,
            "v_gust": v_g,
            "w_gust": w_g,
        }
    )
    return coefficients

def one_cosine_gust_velocity(
    time: float,
    *,
    Uds: float,
    H: float,
    V: float,
    start_time: float = 0.0,
) -> float:
    """Return a signed 1-cosine gust velocity at time.

    For crosswind use, positive Uds means a gust blowing from the aircraft
    right side toward the left side.  The returned scalar uses the same
    positive-from-right convention.  Uds has the same length/time units as V.
    If a Part 25 design gust is specified in EAS, convert it to the physical/TAS
    gust velocity at the analysis density before calling this function.

    H is the gust-gradient distance.  The spatial gust coordinate is
    s = V * (time - start_time).  Outside 0 <= s <= 2H the gust velocity is zero.
    """
    H = float(H)
    V = float(V)
    if H <= 0.0:
        raise ValueError("H must be positive for a 1-cosine gust.")
    if V <= 0.0:
        raise ValueError("V must be positive for a 1-cosine gust.")
    s = V * (float(time) - float(start_time))
    if s < 0.0 or s > 2.0 * H:
        return 0.0
    return 0.5 * float(Uds) * (1.0 - math.cos(math.pi * s / H))

def _resolve_initial_controls(
    stab,
    state0: np.ndarray,
    control_columns: Mapping[str, str],
    *,
    delta_a: float,
    delta_e: float | None,
    delta_r: float,
    trim_elevator: bool,
) -> dict[str, float]:
    """Resolve delta_e/delta_a/delta_r, optionally trimming initial CMm with elevator."""
    controls = {"delta_a": float(delta_a), "delta_e": 0.0, "delta_r": float(delta_r)}
    if delta_e is not None:
        controls["delta_e"] = float(delta_e)
        return controls
    if not trim_elevator:
        controls["delta_e"] = 0.0
        return controls

    elevator_column = control_columns.get("delta_e")
    if elevator_column is None or elevator_column not in stab.derivatives.columns:
        raise KeyError("Could not resolve elevator control column for trim_elevator=True.")

    # Evaluate CMm with delta_e = 0 while retaining delta_a and delta_r.
    cm_without_elevator = _linear_aero_from_stab(stab, state0, controls, control_columns)["CMm"]
    cm_delta_e = _stab_value(stab, "CMm", elevator_column)
    if abs(cm_delta_e) < 1.0e-14:
        raise ZeroDivisionError(f"Elevator pitch derivative is too small: CMm/{elevator_column}={cm_delta_e}")
    controls["delta_e"] = -cm_without_elevator / cm_delta_e
    return controls

def _resolve_initial_thrust(
    stab,
    state0: np.ndarray,
    controls: Mapping[str, float],
    control_columns: Mapping[str, str],
    *,
    mass: float,
    rho: float,
    g: float,
    thrust: float | None,
    trim_thrust: bool,
) -> float:
    """Resolve thrust, optionally trimming initial body-x force balance."""
    if thrust is not None:
        return float(thrust)
    if not trim_thrust:
        return 0.0

    u, v, w, p, q, r, x_e, y_e, z_e, phi, theta, psi = [float(value) for value in state0]
    coeff = _linear_aero_from_stab(stab, state0, controls, control_columns)
    qbar_s = 0.5 * rho * coeff["V"] * coeff["V"] * stab.Sref
    required_x = mass * (q * w - r * v + g * math.sin(theta))
    return required_x - qbar_s * coeff["CX"]

def simulate_6dof_rudder_step_from_stab(
    stab_path: str | Path,
    *,
    mass: float,
    Ixx: float,
    Iyy: float,
    Izz: float,
    Ixz: float = 0.0,
    delta_r: float,
    t_final: float = 5.0,
    control_map: Mapping[str, str] | None = None,
    delta_a: float = 0.0,
    delta_e: float | None = None,
    trim_elevator: bool = True,
    theta_hold: bool = False,
    theta_ref: float | None = None,
    theta_hold_kp: float = 0.0,
    theta_hold_kq: float = 0.0,
    delta_e_min: float | None = None,
    delta_e_max: float | None = None,
    thrust: float | None = None,
    trim_thrust: bool = True,
    g: float = 9.80665,
    rho: float | None = None,
    phi0: float = 0.0,
    theta0: float | None = None,
    psi0: float = 0.0,
    position0: Sequence[float] = (0.0, 0.0, 0.0),
    max_step: float | None = None,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
    stop_at_target_delta_phi: bool = False,
    target_delta_phi: float = math.radians(5.0),
) -> pd.DataFrame:
    """
    Integrate a 6DOF nonlinear rigid-body model using .stab linear aero.

    If delta_e is None and trim_elevator=True, the initial elevator deflection is
    computed from the .stab CMm row so that initial CMm is zero. If theta_hold is
    False, this elevator deflection is fixed after initialization. If theta_hold
    is True, the initial value is used as delta_e0 and the elevator is updated
    inside the ODE right-hand side to reduce theta and pitch-rate drift. If
    thrust is None and trim_thrust=True, the initial body-x force balance is used
    to set thrust. These trims improve the initial condition but do not guarantee
    a sustained free-flight trim over time.

    If stop_at_target_delta_phi=True, integration stops when phi(t) - phi0
    reaches the signed target_delta_phi. This is intended for response-time
    indices; leave it False when a full time history to t_final is needed.
    """
    stab = read_vspaero_stab(stab_path)
    control_columns = resolve_control_columns_from_stab(stab, control_map)
    rho_used = float(stab.rho0 if rho is None else rho)
    state0 = _initial_state_from_stab(stab, phi0=phi0, theta0=theta0, psi0=psi0, position0=position0)
    controls = _resolve_initial_controls(
        stab,
        state0,
        control_columns,
        delta_a=delta_a,
        delta_e=delta_e,
        delta_r=delta_r,
        trim_elevator=trim_elevator,
    )
    delta_e0 = float(controls["delta_e"])
    theta_reference = float(state0[10] if theta_ref is None else theta_ref)
    elevator_column = control_columns.get("delta_e")
    cm_delta_e = None
    if theta_hold:
        if elevator_column is None or elevator_column not in stab.derivatives.columns:
            raise KeyError("Could not resolve elevator control column for theta_hold=True.")
        cm_delta_e = _stab_value(stab, "CMm", elevator_column)
        if abs(cm_delta_e) < 1.0e-14:
            raise ZeroDivisionError(f"Elevator pitch derivative is too small: CMm/{elevator_column}={cm_delta_e}")

    thrust_used = _resolve_initial_thrust(
        stab,
        state0,
        controls,
        control_columns,
        mass=mass,
        rho=rho_used,
        g=g,
        thrust=thrust,
        trim_thrust=trim_thrust,
    )

    inertia_denominator = Ixx * Izz - Ixz * Ixz
    if abs(inertia_denominator) < 1.0e-14:
        raise ValueError("Ixx * Izz - Ixz**2 must not be near zero.")

    def rhs(t: float, state: np.ndarray) -> np.ndarray:
        u, v, w, p, q, r, x_e, y_e, z_e, phi, theta, psi = state
        delta_e_current = delta_e0
        if theta_hold:
            theta_rate = q * math.cos(phi) - r * math.sin(phi)
            pitch_hold_moment = -theta_hold_kp * (theta - theta_reference) - theta_hold_kq * theta_rate
            delta_e_current = delta_e0 + pitch_hold_moment / cm_delta_e
            if delta_e_min is not None:
                delta_e_current = max(float(delta_e_min), delta_e_current)
            if delta_e_max is not None:
                delta_e_current = min(float(delta_e_max), delta_e_current)

        controls_current = {"delta_a": controls["delta_a"], "delta_e": delta_e_current, "delta_r": controls["delta_r"]}
        coeff = _linear_aero_from_stab(stab, state, controls_current, control_columns)
        qbar_s = 0.5 * rho_used * coeff["V"] * coeff["V"] * stab.Sref
        X = thrust_used + qbar_s * coeff["CX"]
        Y = qbar_s * coeff["CY"]
        Z = qbar_s * coeff["CZ"]
        L = qbar_s * stab.Bref * coeff["Cl"]
        M = qbar_s * stab.Cref * coeff["Cm"]
        N = qbar_s * stab.Bref * coeff["Cn"]

        u_dot = r * v - q * w + X / mass - g * math.sin(theta)
        v_dot = p * w - r * u + Y / mass + g * math.sin(phi) * math.cos(theta)
        w_dot = q * u - p * v + Z / mass + g * math.cos(phi) * math.cos(theta)
        p_dot = (Izz * (L + (Iyy - Izz) * q * r + Ixz * p * q) + Ixz * (N + (Ixx - Iyy) * p * q - Ixz * q * r)) / inertia_denominator
        q_dot = (M + (Izz - Ixx) * p * r - Ixz * (p * p - r * r)) / Iyy
        r_dot = (Ixz * (L + (Iyy - Izz) * q * r + Ixz * p * q) + Ixx * (N + (Ixx - Iyy) * p * q - Ixz * q * r)) / inertia_denominator

        cos_theta = math.cos(theta)
        if abs(cos_theta) < 1.0e-10:
            cos_theta = math.copysign(1.0e-10, cos_theta)

        x_dot = u * math.cos(theta) * math.cos(psi) + v * (math.sin(phi) * math.sin(theta) * math.cos(psi) - math.cos(phi) * math.sin(psi)) + w * (math.cos(phi) * math.sin(theta) * math.cos(psi) + math.sin(phi) * math.sin(psi))
        y_dot = u * math.cos(theta) * math.sin(psi) + v * (math.sin(phi) * math.sin(theta) * math.sin(psi) + math.cos(phi) * math.cos(psi)) + w * (math.cos(phi) * math.sin(theta) * math.sin(psi) - math.sin(phi) * math.cos(psi))
        z_dot = -u * math.sin(theta) + v * math.sin(phi) * math.cos(theta) + w * math.cos(phi) * math.cos(theta)
        phi_dot = p + q * math.sin(phi) * math.tan(theta) + r * math.cos(phi) * math.tan(theta)
        theta_dot = q * math.cos(phi) - r * math.sin(phi)
        psi_dot = (q * math.sin(phi) + r * math.cos(phi)) / cos_theta
        return np.array([u_dot, v_dot, w_dot, p_dot, q_dot, r_dot, x_dot, y_dot, z_dot, phi_dot, theta_dot, psi_dot], dtype=float)

    kwargs: dict[str, Any] = {"rtol": rtol, "atol": atol}
    if max_step is not None:
        kwargs["max_step"] = float(max_step)

    target_reached = False
    if stop_at_target_delta_phi:
        target = float(target_delta_phi)
        if abs(target) < 1.0e-14:
            raise ValueError("target_delta_phi must be non-zero when stop_at_target_delta_phi=True.")
        sign = 1.0 if target > 0.0 else -1.0
        phi_reference = float(phi0)

        def target_delta_phi_event(t: float, state: np.ndarray) -> float:
            return sign * ((float(state[9]) - phi_reference) - target)

        target_delta_phi_event.terminal = True
        target_delta_phi_event.direction = 1.0
        kwargs["events"] = target_delta_phi_event

    solution = solve_ivp(rhs, (0.0, float(t_final)), state0, **kwargs)
    if not solution.success:
        raise RuntimeError(solution.message)
    if stop_at_target_delta_phi and solution.t_events:
        target_reached = bool(len(solution.t_events[0]) > 0)

    columns = ["u", "v", "w", "p", "q", "r", "x_e", "y_e", "z_e", "phi", "theta", "psi"]
    history = pd.DataFrame(solution.y.T, columns=columns)
    history.insert(0, "time", solution.t)

    phidot_values, phat_values, qhat_values, rhat_values = [], [], [], []
    V_values, alpha_values, beta_values = [], [], []
    delta_e_values = []
    CX_values, CY_values, CZ_values = [], [], []
    CL_values, CD_values, CS_values = [], [], []

    for row in history.itertuples(index=False):
        state = np.array([row.u, row.v, row.w, row.p, row.q, row.r, row.x_e, row.y_e, row.z_e, row.phi, row.theta, row.psi])
        delta_e_current = delta_e0
        if theta_hold:
            theta_rate = row.q * math.cos(row.phi) - row.r * math.sin(row.phi)
            pitch_hold_moment = -theta_hold_kp * (row.theta - theta_reference) - theta_hold_kq * theta_rate
            delta_e_current = delta_e0 + pitch_hold_moment / cm_delta_e
            if delta_e_min is not None:
                delta_e_current = max(float(delta_e_min), delta_e_current)
            if delta_e_max is not None:
                delta_e_current = min(float(delta_e_max), delta_e_current)

        controls_current = {"delta_a": controls["delta_a"], "delta_e": delta_e_current, "delta_r": controls["delta_r"]}
        coeff = _linear_aero_from_stab(stab, state, controls_current, control_columns)
        V = coeff["V"]
        delta_e_values.append(delta_e_current)
        phidot_values.append(row.p + row.q * math.sin(row.phi) * math.tan(row.theta) + row.r * math.cos(row.phi) * math.tan(row.theta))
        phat_values.append(row.p * stab.Bref / (2.0 * V))
        qhat_values.append(row.q * stab.Cref / (2.0 * V))
        rhat_values.append(row.r * stab.Bref / (2.0 * V))
        V_values.append(V)
        alpha_values.append(coeff["alpha"])
        beta_values.append(coeff["beta"])
        CX_values.append(coeff["CX"])
        CY_values.append(coeff["CY"])
        CZ_values.append(coeff["CZ"])
        CL_values.append(coeff["CL"])
        CD_values.append(coeff["CD"])
        CS_values.append(coeff["CS"])

    history["V"] = V_values
    history["alpha"] = alpha_values
    history["beta"] = beta_values
    history["CX"] = CX_values
    history["CY"] = CY_values
    history["CZ"] = CZ_values
    history["CL"] = CL_values
    history["CD"] = CD_values
    history["CS"] = CS_values
    history["phi_dot"] = phidot_values
    history["phat"] = phat_values
    history["qhat"] = qhat_values
    history["rhat"] = rhat_values
    history["delta_e"] = delta_e_values
    history["delta_a"] = controls["delta_a"]
    history["delta_r"] = controls["delta_r"]
    history["thrust"] = thrust_used
    history.attrs.update({
        "stab_path": str(stab_path),
        "mass": float(mass),
        "Ixx": float(Ixx),
        "Iyy": float(Iyy),
        "Izz": float(Izz),
        "Ixz": float(Ixz),
        "rho": rho_used,
        "Sref": float(stab.Sref),
        "Cref": float(stab.Cref),
        "Bref": float(stab.Bref),
        "control_columns": control_columns,
        "delta_e0": delta_e0,
        "delta_e": float(delta_e_values[0]) if delta_e_values else delta_e0,
        "delta_a": controls["delta_a"],
        "delta_r": controls["delta_r"],
        "trim_elevator": bool(delta_e is None and trim_elevator),
        "theta_hold": bool(theta_hold),
        "theta_ref": theta_reference,
        "theta_hold_kp": float(theta_hold_kp),
        "theta_hold_kq": float(theta_hold_kq),
        "delta_e_min": None if delta_e_min is None else float(delta_e_min),
        "delta_e_max": None if delta_e_max is None else float(delta_e_max),
        "thrust": thrust_used,
        "trim_thrust": bool(thrust is None and trim_thrust),
        "stop_at_target_delta_phi": bool(stop_at_target_delta_phi),
        "target_delta_phi": float(target_delta_phi),
        "target_delta_phi_reached": bool(target_reached),
    })
    return history

def simulate_6dof_crosswind_gust_from_stab(
    stab_path: str | Path,
    *,
    mass: float,
    Ixx: float,
    Iyy: float,
    Izz: float,
    Ixz: float = 0.0,
    Uds: float,
    H: float,
    t_final: float | None = None,
    control_map: Mapping[str, str] | None = None,
    delta_a: float = 0.0,
    delta_e: float | None = None,
    delta_r: float = 0.0,
    trim_elevator: bool = True,
    theta_hold: bool = False,
    theta_ref: float | None = None,
    theta_hold_kp: float = 0.0,
    theta_hold_kq: float = 0.0,
    delta_e_min: float | None = None,
    delta_e_max: float | None = None,
    thrust: float | None = None,
    trim_thrust: bool = True,
    g: float = 9.80665,
    rho: float | None = None,
    phi0: float = 0.0,
    theta0: float | None = None,
    psi0: float = 0.0,
    position0: Sequence[float] = (0.0, 0.0, 0.0),
    gust_start_time: float = 0.0,
    max_step: float | None = None,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
) -> pd.DataFrame:
    """Integrate a 6DOF nonlinear rigid-body model with a lateral 1-cosine gust.

    The public crosswind sign convention is positive from the aircraft right
    side.  Therefore Uds > 0 and gust_velocity_y > 0 mean wind blowing from
    right to left, and beta_g > 0 is the corresponding positive gust sideslip.

    Solver body axes remain +x forward, +y right, +z down.  The actual local
    atmospheric velocity used by the relative-air calculation is therefore

        [0, -gust_velocity_y, 0].

    The aircraft states remain the rigid-body velocities and attitudes.  The
    gust is applied only when evaluating aerodynamic coefficients, by
    subtracting this atmospheric velocity vector from the aircraft body-axis
    velocity.  Thus v_air = v + gust_velocity_y and, for small disturbances,
    beta_air ~= beta_body + beta_g.

    Uds is the signed peak physical gust speed in the same length/time units as
    V0.  If a Part 25 design gust is specified in EAS, convert it to the
    physical/TAS gust speed at the analysis density before calling this
    function.  H is the gust-gradient distance.  The gust runs over
    0 <= s <= 2H, where s = V0 * (t - gust_start_time).
    """
    stab = read_vspaero_stab(stab_path)
    control_columns = resolve_control_columns_from_stab(stab, control_map)
    rho_used = float(stab.rho0 if rho is None else rho)
    Vref = float(stab.V0)
    if Vref <= 0.0:
        raise ValueError("The .stab file must contain a positive Vinf for gust simulation.")
    if float(H) <= 0.0:
        raise ValueError("H must be positive.")
    gust_duration = 2.0 * float(H) / Vref
    if t_final is None:
        t_final = float(gust_start_time) + gust_duration
    if float(t_final) <= 0.0:
        raise ValueError("t_final must be positive.")

    state0 = _initial_state_from_stab(stab, phi0=phi0, theta0=theta0, psi0=psi0, position0=position0)
    controls = _resolve_initial_controls(
        stab,
        state0,
        control_columns,
        delta_a=delta_a,
        delta_e=delta_e,
        delta_r=delta_r,
        trim_elevator=trim_elevator,
    )
    delta_e0 = float(controls["delta_e"])
    theta_reference = float(state0[10] if theta_ref is None else theta_ref)
    elevator_column = control_columns.get("delta_e")
    cm_delta_e = None
    if theta_hold:
        if elevator_column is None or elevator_column not in stab.derivatives.columns:
            raise KeyError("Could not resolve elevator control column for theta_hold=True.")
        cm_delta_e = _stab_value(stab, "CMm", elevator_column)
        if abs(cm_delta_e) < 1.0e-14:
            raise ZeroDivisionError(f"Elevator pitch derivative is too small: CMm/{elevator_column}={cm_delta_e}")

    thrust_used = _resolve_initial_thrust(
        stab,
        state0,
        controls,
        control_columns,
        mass=mass,
        rho=rho_used,
        g=g,
        thrust=thrust,
        trim_thrust=trim_thrust,
    )

    inertia_denominator = Ixx * Izz - Ixz * Ixz
    if abs(inertia_denominator) < 1.0e-14:
        raise ValueError("Ixx * Izz - Ixz**2 must not be near zero.")

    def gust_velocity_y_at(t: float) -> float:
        """Return the signed crosswind speed, positive when blowing from right."""
        return one_cosine_gust_velocity(
            t,
            Uds=Uds,
            H=H,
            V=Vref,
            start_time=gust_start_time,
        )

    def airmass_body_velocity_at(t: float) -> tuple[float, float, float]:
        """Convert positive-from-right crosswind speed to body-axis air velocity."""
        gust_velocity_y = gust_velocity_y_at(t)
        return 0.0, -gust_velocity_y, 0.0

    def elevator_deflection_for_state(state: Sequence[float]) -> float:
        if not theta_hold:
            return delta_e0
        _, _, _, _, q, r, _, _, _, phi, theta, _ = [float(value) for value in state]
        theta_rate = q * math.cos(phi) - r * math.sin(phi)
        pitch_hold_moment = -theta_hold_kp * (theta - theta_reference) - theta_hold_kq * theta_rate
        value = delta_e0 + pitch_hold_moment / cm_delta_e
        if delta_e_min is not None:
            value = max(float(delta_e_min), value)
        if delta_e_max is not None:
            value = min(float(delta_e_max), value)
        return float(value)

    def rhs(t: float, state: np.ndarray) -> np.ndarray:
        u, v, w, p, q, r, x_e, y_e, z_e, phi, theta, psi = state
        delta_e_current = elevator_deflection_for_state(state)
        controls_current = {"delta_a": controls["delta_a"], "delta_e": delta_e_current, "delta_r": controls["delta_r"]}
        coeff = _linear_aero_from_stab(
            stab,
            state,
            controls_current,
            control_columns,
            gust_body_velocity=airmass_body_velocity_at(t),
        )
        qbar_s = 0.5 * rho_used * coeff["V"] * coeff["V"] * stab.Sref
        X = thrust_used + qbar_s * coeff["CX"]
        Y = qbar_s * coeff["CY"]
        Z = qbar_s * coeff["CZ"]
        L = qbar_s * stab.Bref * coeff["Cl"]
        M = qbar_s * stab.Cref * coeff["Cm"]
        N = qbar_s * stab.Bref * coeff["Cn"]

        u_dot = r * v - q * w + X / mass - g * math.sin(theta)
        v_dot = p * w - r * u + Y / mass + g * math.sin(phi) * math.cos(theta)
        w_dot = q * u - p * v + Z / mass + g * math.cos(phi) * math.cos(theta)
        p_dot = (Izz * (L + (Iyy - Izz) * q * r + Ixz * p * q) + Ixz * (N + (Ixx - Iyy) * p * q - Ixz * q * r)) / inertia_denominator
        q_dot = (M + (Izz - Ixx) * p * r - Ixz * (p * p - r * r)) / Iyy
        r_dot = (Ixz * (L + (Iyy - Izz) * q * r + Ixz * p * q) + Ixx * (N + (Ixx - Iyy) * p * q - Ixz * q * r)) / inertia_denominator

        cos_theta = math.cos(theta)
        if abs(cos_theta) < 1.0e-10:
            cos_theta = math.copysign(1.0e-10, cos_theta)

        x_dot = u * math.cos(theta) * math.cos(psi) + v * (math.sin(phi) * math.sin(theta) * math.cos(psi) - math.cos(phi) * math.sin(psi)) + w * (math.cos(phi) * math.sin(theta) * math.cos(psi) + math.sin(phi) * math.sin(psi))
        y_dot = u * math.cos(theta) * math.sin(psi) + v * (math.sin(phi) * math.sin(theta) * math.sin(psi) + math.cos(phi) * math.cos(psi)) + w * (math.cos(phi) * math.sin(theta) * math.sin(psi) - math.sin(phi) * math.cos(psi))
        z_dot = -u * math.sin(theta) + v * math.sin(phi) * math.cos(theta) + w * math.cos(phi) * math.cos(theta)
        phi_dot = p + q * math.sin(phi) * math.tan(theta) + r * math.cos(phi) * math.tan(theta)
        theta_dot = q * math.cos(phi) - r * math.sin(phi)
        psi_dot = (q * math.sin(phi) + r * math.cos(phi)) / cos_theta
        return np.array([u_dot, v_dot, w_dot, p_dot, q_dot, r_dot, x_dot, y_dot, z_dot, phi_dot, theta_dot, psi_dot], dtype=float)

    kwargs: dict[str, Any] = {"rtol": rtol, "atol": atol}
    if max_step is not None:
        kwargs["max_step"] = float(max_step)

    solution = solve_ivp(rhs, (0.0, float(t_final)), state0, **kwargs)
    if not solution.success:
        raise RuntimeError(solution.message)

    columns = ["u", "v", "w", "p", "q", "r", "x_e", "y_e", "z_e", "phi", "theta", "psi"]
    history = pd.DataFrame(solution.y.T, columns=columns)
    history.insert(0, "time", solution.t)

    phidot_values, phat_values, qhat_values, rhat_values = [], [], [], []
    V_values, alpha_values, beta_values = [], [], []
    V_body_values, alpha_body_values, beta_body_values = [], [], []
    delta_e_values = []
    CX_values, CY_values, CZ_values = [], [], []
    CL_values, CD_values, CS_values = [], [], []
    gust_y_values, airmass_y_values, beta_g_values = [], [], []
    u_air_values, v_air_values, w_air_values = [], [], []

    for row in history.itertuples(index=False):
        state = np.array([row.u, row.v, row.w, row.p, row.q, row.r, row.x_e, row.y_e, row.z_e, row.phi, row.theta, row.psi])
        delta_e_current = elevator_deflection_for_state(state)
        controls_current = {"delta_a": controls["delta_a"], "delta_e": delta_e_current, "delta_r": controls["delta_r"]}
        gust_velocity_y = gust_velocity_y_at(row.time)
        airmass_body_velocity = airmass_body_velocity_at(row.time)
        coeff = _linear_aero_from_stab(
            stab,
            state,
            controls_current,
            control_columns,
            gust_body_velocity=airmass_body_velocity,
        )
        V = coeff["V"]
        V_body = max(math.sqrt(row.u * row.u + row.v * row.v + row.w * row.w), 1.0e-12)
        alpha_body = math.atan2(row.w, row.u)
        beta_body = math.asin(max(-1.0, min(1.0, row.v / V_body)))

        delta_e_values.append(delta_e_current)
        phidot_values.append(row.p + row.q * math.sin(row.phi) * math.tan(row.theta) + row.r * math.cos(row.phi) * math.tan(row.theta))
        phat_values.append(row.p * stab.Bref / (2.0 * V))
        qhat_values.append(row.q * stab.Cref / (2.0 * V))
        rhat_values.append(row.r * stab.Bref / (2.0 * V))
        V_values.append(V)
        alpha_values.append(coeff["alpha"])
        beta_values.append(coeff["beta"])
        V_body_values.append(V_body)
        alpha_body_values.append(alpha_body)
        beta_body_values.append(beta_body)
        CX_values.append(coeff["CX"])
        CY_values.append(coeff["CY"])
        CZ_values.append(coeff["CZ"])
        CL_values.append(coeff["CL"])
        CD_values.append(coeff["CD"])
        CS_values.append(coeff["CS"])
        gust_y_values.append(gust_velocity_y)
        airmass_y_values.append(airmass_body_velocity[1])
        beta_g_values.append(math.atan2(gust_velocity_y, Vref))
        u_air_values.append(coeff["u_air"])
        v_air_values.append(coeff["v_air"])
        w_air_values.append(coeff["w_air"])

    history["V"] = V_values
    history["alpha"] = alpha_values
    history["beta"] = beta_values
    history["V_body"] = V_body_values
    history["alpha_body"] = alpha_body_values
    history["beta_body"] = beta_body_values
    history["u_air"] = u_air_values
    history["v_air"] = v_air_values
    history["w_air"] = w_air_values
    history["gust_velocity_y"] = gust_y_values
    history["airmass_velocity_y"] = airmass_y_values
    history["beta_g"] = beta_g_values
    history["crosswind_sign_convention"] = "positive_from_right"
    history["CX"] = CX_values
    history["CY"] = CY_values
    history["CZ"] = CZ_values
    history["CL"] = CL_values
    history["CD"] = CD_values
    history["CS"] = CS_values
    history["phi_dot"] = phidot_values
    history["phat"] = phat_values
    history["qhat"] = qhat_values
    history["rhat"] = rhat_values
    history["delta_e"] = delta_e_values
    history["delta_a"] = controls["delta_a"]
    history["delta_r"] = controls["delta_r"]
    history["thrust"] = thrust_used
    history["phi_delta"] = history["phi"] - float(phi0)
    history.attrs.update({
        "simulation_source": "6dof_crosswind_1_cosine_gust_from_stab",
        "stab_path": str(stab_path),
        "mass": float(mass),
        "Ixx": float(Ixx),
        "Iyy": float(Iyy),
        "Izz": float(Izz),
        "Ixz": float(Ixz),
        "rho": rho_used,
        "Sref": float(stab.Sref),
        "Cref": float(stab.Cref),
        "Bref": float(stab.Bref),
        "control_columns": control_columns,
        "delta_e0": delta_e0,
        "delta_e": float(delta_e_values[0]) if delta_e_values else delta_e0,
        "delta_a": controls["delta_a"],
        "delta_r": controls["delta_r"],
        "trim_elevator": bool(delta_e is None and trim_elevator),
        "theta_hold": bool(theta_hold),
        "theta_ref": theta_reference,
        "theta_hold_kp": float(theta_hold_kp),
        "theta_hold_kq": float(theta_hold_kq),
        "delta_e_min": None if delta_e_min is None else float(delta_e_min),
        "delta_e_max": None if delta_e_max is None else float(delta_e_max),
        "thrust": thrust_used,
        "trim_thrust": bool(thrust is None and trim_thrust),
        "gust_Uds": float(Uds),
        "gust_H": float(H),
        "gust_start_time": float(gust_start_time),
        "gust_end_time": float(gust_start_time) + gust_duration,
        "gust_duration": gust_duration,
        "gust_reference_V": Vref,
        "gust_beta_peak": math.atan2(float(Uds), Vref),
    })
    return history

def _history_peak(history: pd.DataFrame, column: str) -> tuple[float, float, float]:
    values = history[column].to_numpy(dtype=float)
    if len(values) == 0 or np.all(np.isnan(values)):
        return math.nan, math.nan, math.nan
    index = int(np.nanargmax(np.abs(values)))
    return float(values[index]), float(abs(values[index])), float(history["time"].iloc[index])

def calculate_crosswind_gust_response_indices(
    history: pd.DataFrame,
    *,
    phi0: float | None = None,
    Vref: float | None = None,
    Bref: float | None = None,
) -> dict[str, Any]:
    """Calculate compact response indices from a positive-from-right gust history."""
    sign_convention = history.attrs.get("crosswind_sign_convention")
    if sign_convention is None and "crosswind_sign_convention" in history.columns:
        values = history["crosswind_sign_convention"].dropna().astype(str).unique()
        if len(values) == 1:
            sign_convention = values[0]

    for column in ("time", "phi", "beta", "p", "r", "phat", "rhat"):
        if column not in history.columns:
            raise KeyError(f"history is missing required column: {column!r}")

    time_values = history["time"].to_numpy(dtype=float)
    if len(time_values) < 2 or float(time_values[-1]) <= float(time_values[0]):
        raise ValueError("history must contain at least two increasing time samples.")

    phi_reference = float(history["phi"].iloc[0] if phi0 is None else phi0)
    V_used = float(history.attrs.get("gust_reference_V", Vref if Vref is not None else history["V"].iloc[0]))
    Bref_used = float(Bref if Bref is not None else history.attrs.get("Bref", history.attrs.get("gust_Bref", np.nan)))
    if not math.isfinite(Bref_used) and "phat" in history.columns and "p" in history.columns:
        Bref_used = math.nan

    Uds = float(history.attrs.get("gust_Uds", np.nan))
    H = float(history.attrs.get("gust_H", np.nan))
    gust_start = float(history.attrs.get("gust_start_time", time_values[0]))
    gust_end = float(history.attrs.get("gust_end_time", time_values[-1]))
    beta_g_peak = float(
        history.attrs.get(
            "gust_beta_peak",
            math.atan2(Uds, V_used)
            if math.isfinite(Uds) and V_used > 0.0
            else np.nan,
        )
    )
    beta_scale = abs(beta_g_peak) if math.isfinite(beta_g_peak) and abs(beta_g_peak) > 1.0e-14 else math.nan

    phi_delta = history["phi"].to_numpy(dtype=float) - phi_reference
    phi_delta_final = float(phi_delta[-1])
    phi_delta_at_gust_end = float(np.interp(gust_end, time_values, phi_delta))
    max_phi_index = int(np.nanargmax(np.abs(phi_delta)))
    max_phi_delta = float(phi_delta[max_phi_index])
    max_abs_phi_delta = float(abs(max_phi_delta))
    max_phi_time = float(time_values[max_phi_index])

    beta_peak, beta_peak_abs, beta_peak_time = _history_peak(history, "beta")
    p_peak, p_peak_abs, p_peak_time = _history_peak(history, "p")
    r_peak, r_peak_abs, r_peak_time = _history_peak(history, "r")
    phat_peak, phat_peak_abs, phat_peak_time = _history_peak(history, "phat")
    rhat_peak, rhat_peak_abs, rhat_peak_time = _history_peak(history, "rhat")

    result: dict[str, Any] = {
        "crosswind_gust_success": True,
        "crosswind_gust_message": "",
        "crosswind_gust_Uds": Uds,
        "crosswind_gust_H": H,
        "crosswind_gust_start_time": gust_start,
        "crosswind_gust_end_time": gust_end,
        "crosswind_gust_duration": float(history.attrs.get("gust_duration", gust_end - gust_start)),
        "crosswind_gust_reference_V": V_used,
        "crosswind_gust_beta_peak_input": beta_g_peak,
        "crosswind_gust_phi0": phi_reference,
        "crosswind_gust_t_start": float(time_values[0]),
        "crosswind_gust_t_final": float(time_values[-1]),
        "crosswind_gust_phi_final": float(history["phi"].iloc[-1]),
        "crosswind_gust_phi_delta_final": phi_delta_final,
        "crosswind_gust_phi_delta_at_gust_end": phi_delta_at_gust_end,
        "crosswind_gust_max_phi_delta": max_phi_delta,
        "crosswind_gust_max_abs_phi_delta": max_abs_phi_delta,
        "crosswind_gust_max_abs_phi_delta_time": max_phi_time,
        "crosswind_gust_peak_beta": beta_peak,
        "crosswind_gust_max_abs_beta": beta_peak_abs,
        "crosswind_gust_max_abs_beta_time": beta_peak_time,
        "crosswind_gust_peak_p": p_peak,
        "crosswind_gust_max_abs_p": p_peak_abs,
        "crosswind_gust_max_abs_p_time": p_peak_time,
        "crosswind_gust_peak_r": r_peak,
        "crosswind_gust_max_abs_r": r_peak_abs,
        "crosswind_gust_max_abs_r_time": r_peak_time,
        "crosswind_gust_peak_phat": phat_peak,
        "crosswind_gust_max_abs_phat": phat_peak_abs,
        "crosswind_gust_max_abs_phat_time": phat_peak_time,
        "crosswind_gust_peak_rhat": rhat_peak,
        "crosswind_gust_max_abs_rhat": rhat_peak_abs,
        "crosswind_gust_max_abs_rhat_time": rhat_peak_time,
        "crosswind_gust_bank_index_abs_per_beta_g": math.nan,
        "crosswind_gust_bank_index_final_per_beta_g": math.nan,
        "crosswind_gust_bank_index_at_gust_end_per_beta_g": math.nan,
        "crosswind_gust_phat_index_abs_per_beta_g": math.nan,
        "crosswind_gust_rhat_index_abs_per_beta_g": math.nan,
    }
    if math.isfinite(beta_scale):
        result["crosswind_gust_bank_index_abs_per_beta_g"] = max_abs_phi_delta / beta_scale
        result["crosswind_gust_phat_index_abs_per_beta_g"] = phat_peak_abs / beta_scale
        result["crosswind_gust_rhat_index_abs_per_beta_g"] = rhat_peak_abs / beta_scale
        result["crosswind_gust_bank_index_final_per_beta_g"] = phi_delta_final / beta_g_peak
        result["crosswind_gust_bank_index_at_gust_end_per_beta_g"] = phi_delta_at_gust_end / beta_g_peak
    if math.isfinite(Bref_used):
        result["crosswind_gust_index_Bref"] = Bref_used
    return result

def write_6dof_history_csv(history: pd.DataFrame, csv_path: str | Path) -> Path:
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    history.to_csv(csv_path, index=False)
    return csv_path

def plot_6dof_history(history: pd.DataFrame, *, plot_path: str | Path | None = None, show: bool = False, degrees: bool = True):
    import matplotlib.pyplot as plt

    t = history["time"]
    angle_factor = 180.0 / math.pi if degrees else 1.0
    angle_unit = "deg" if degrees else "rad"
    fig, axes = plt.subplots(5, 1, sharex=True, figsize=(10, 12))
    axes[0].plot(t, history["V"], label="V")
    axes[0].set_ylabel("V")
    axes[1].plot(t, history["alpha"] * angle_factor, label=f"alpha [{angle_unit}]")
    axes[1].plot(t, history["beta"] * angle_factor, label=f"beta [{angle_unit}]")
    axes[1].set_ylabel(f"alpha, beta [{angle_unit}]")
    axes[2].plot(t, history["p"], label="p [rad/s]")
    axes[2].plot(t, history["q"], label="q [rad/s]")
    axes[2].plot(t, history["r"], label="r [rad/s]")
    axes[2].set_ylabel("p, q, r [rad/s]")
    axes[3].plot(t, history["phi"] * angle_factor, label=f"phi [{angle_unit}]")
    axes[3].plot(t, history["theta"] * angle_factor, label=f"theta [{angle_unit}]")
    axes[3].plot(t, history["psi"] * angle_factor, label=f"psi [{angle_unit}]")
    axes[3].set_ylabel(f"Euler angles [{angle_unit}]")
    axes[4].plot(t, history["delta_e"] * angle_factor, label=f"delta_e [{angle_unit}]")
    axes[4].plot(t, history["delta_a"] * angle_factor, label=f"delta_a [{angle_unit}]")
    axes[4].plot(t, history["delta_r"] * angle_factor, label=f"delta_r [{angle_unit}]")
    axes[4].set_ylabel(f"controls [{angle_unit}]")
    axes[4].set_xlabel("t")
    for ax in axes:
        ax.legend(loc="best")
        ax.grid(True)
    fig.tight_layout()
    if plot_path is not None:
        plot_path = Path(plot_path)
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig

def simulate_reduced_lateral_response_from_stab(
    stab_path: str | Path,
    *,
    Ixx: float,
    Izz: float,
    delta_r: float,
    time: Sequence[float],
    Ixz: float = 0.0,
    rho: float | None = None,
    control_map: Mapping[str, str] | None = None,
    roll_moment_coef: str = "CMl",
    yaw_moment_coef: str = "CMn",
    include_roll_damping: bool = False,
    include_yaw_damping: bool = False,
    rtol: float = 1.0e-9,
    atol: float = 1.0e-11,
) -> pd.DataFrame:
    """Return the reduced lateral response on a given time grid.

    The base reduced model is

        beta' = -rhat
        phat' = mu_I Izz Cl_beta beta + mu_I Izz Cl_r rhat
        rhat' = mu_I Ixx Cn_delta_r delta_r
        phi' = phat

    with tau = 2 V / Bref * t.  Optional roll and yaw damping add

        + mu_I Izz Cl_p phat
        + mu_I Ixx Cn_r rhat

    respectively.  Side force, direct rudder roll, gravity, and Ixz coupling
    are intentionally excluded so this remains a readable explanatory model
    against the 6DOF history, not the main linear lateral model.
    """
    stab = read_vspaero_stab(stab_path)
    control_columns = resolve_control_columns_from_stab(stab, control_map)
    rudder_column = control_columns.get("delta_r")
    if rudder_column is None:
        raise KeyError("Could not resolve rudder control column. Pass control_map={'delta_r': 'ConGrp_N'} explicitly.")

    rho_used = float(stab.rho0 if rho is None else rho)
    if not math.isfinite(rho_used) or rho_used <= 0.0:
        raise ValueError("rho must be positive, or the .stab file must contain a positive Rho value.")

    Ixx = float(Ixx)
    Izz = float(Izz)
    Ixz = float(Ixz)
    inertia_det = Ixx * Izz - Ixz * Ixz
    if inertia_det <= 0.0:
        raise ValueError("Ixx * Izz - Ixz**2 must be positive.")

    time_values = np.asarray(time, dtype=float)
    if time_values.ndim != 1 or len(time_values) < 2:
        raise ValueError("time must be a one-dimensional sequence with at least two samples.")
    if np.any(~np.isfinite(time_values)) or np.any(np.diff(time_values) <= 0.0):
        raise ValueError("time must contain finite, strictly increasing values.")

    V = float(stab.V0)
    Bref = float(stab.Bref)
    Sref = float(stab.Sref)
    if V <= 0.0 or Bref <= 0.0:
        raise ValueError("The .stab file must contain positive Vinf and Bref values.")

    qbar = 0.5 * rho_used * V * V
    tau_per_second = 2.0 * V / Bref
    mu_inertia = qbar * Sref * Bref**3 / (4.0 * V * V * inertia_det)

    Cl_beta = _solver_axis_derivative_value(stab, roll_moment_coef, "Beta")
    Cl_phat = _solver_axis_derivative_value(stab, roll_moment_coef, "p")
    Cl_rhat = _solver_axis_derivative_value(stab, roll_moment_coef, "r")
    Cn_rhat = _solver_axis_derivative_value(stab, yaw_moment_coef, "r")
    Cn_delta_r = _solver_axis_derivative_value(stab, yaw_moment_coef, rudder_column)

    tau = tau_per_second * (time_values - time_values[0])
    delta_r = float(delta_r)
    rate_scale = 2.0 * V / Bref

    roll_damping = mu_inertia * Izz * Cl_phat if include_roll_damping else 0.0
    yaw_damping = mu_inertia * Ixx * Cn_rhat if include_yaw_damping else 0.0
    A = np.array(
        [
            [0.0, 0.0, -1.0, 0.0],
            [mu_inertia * Izz * Cl_beta, roll_damping, mu_inertia * Izz * Cl_rhat, 0.0],
            [0.0, 0.0, yaw_damping, 0.0],
            [0.0, 1.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    B = np.array(
        [
            0.0,
            0.0,
            mu_inertia * Ixx * Cn_delta_r,
            0.0,
        ],
        dtype=float,
    )

    if not include_roll_damping and not include_yaw_damping:
        # Keep the default minimal model algebraic and exactly identical to the
        # K2/K3 explanatory formula.  The solve_ivp path below is used only when
        # damping terms are requested.
        a = mu_inertia * Izz * Cl_beta
        b = mu_inertia * Izz * Cl_rhat
        c = mu_inertia * Ixx * Cn_delta_r
        rhat = c * delta_r * tau
        beta = -0.5 * c * delta_r * tau**2
        phat = c * delta_r * (0.5 * b * tau**2 - (a / 6.0) * tau**3)
        phi = c * delta_r * ((b / 6.0) * tau**3 - (a / 24.0) * tau**4)
    else:
        def rhs(tau_value: float, state: np.ndarray) -> np.ndarray:
            return A @ state + B * delta_r

        solution = solve_ivp(
            rhs,
            (float(tau[0]), float(tau[-1])),
            np.zeros(4, dtype=float),
            t_eval=tau,
            rtol=rtol,
            atol=atol,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        beta, phat, rhat, phi = solution.y

    history = pd.DataFrame(
        {
            "time": time_values,
            "tau": tau,
            "beta": beta,
            "phat": phat,
            "rhat": rhat,
            "phi": phi,
            "p": rate_scale * phat,
            "r": rate_scale * rhat,
        }
    )
    history.attrs.update(
        {
            "stab_path": str(stab_path),
            "rho": rho_used,
            "V": V,
            "Sref": Sref,
            "Bref": Bref,
            "Ixx": Ixx,
            "Izz": Izz,
            "Ixz": Ixz,
            "inertia_det": inertia_det,
            "mu_inertia": mu_inertia,
            "tau_per_second": tau_per_second,
            "delta_r": delta_r,
            "rudder_column": rudder_column,
            "Cl_beta": float(Cl_beta),
            "Cl_phat": float(Cl_phat),
            "Cl_rhat": float(Cl_rhat),
            "Cn_rhat": float(Cn_rhat),
            "Cn_delta_r": float(Cn_delta_r),
            "include_roll_damping": bool(include_roll_damping),
            "include_yaw_damping": bool(include_yaw_damping),
            "roll_damping": float(roll_damping),
            "yaw_damping": float(yaw_damping),
            "A_reduced": A.tolist(),
            "B_reduced": B.tolist(),
            "model": "reduced_beta_phat_rhat_phi",
        }
    )
    return history

def plot_vv_gamma_row_6dof_vs_reduced_response(
    indices: pd.DataFrame | str | Path,
    row_index: int,
    *,
    mass: float,
    inertia: Mapping[str, float],
    rho: float | None = None,
    control_map: Mapping[str, str] | None = None,
    history_path_column: str = "sixdof_history_csv_path",
    recompute_6dof_if_missing: bool = True,
    plot_path: str | Path | None = None,
    show: bool = False,
    degrees: bool = True,
    max_step: float | None = 0.01,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
    include_roll_damping: bool = False,
    include_yaw_damping: bool = False,
):
    """Plot one vv_gamma_indices row: 6DOF history versus reduced response.

    The selected row supplies the .stab path, rudder step, target-delta-phi
    settings, and optionally an existing sixdof_history_csv_path.  The reduced
    response is evaluated on the same time grid as the 6DOF history, then beta,
    p, r, and phi are plotted in a four-row shared-x figure.  Optional roll
    damping and yaw damping affect only the reduced model.  Each y-axis range
    is taken from the 6DOF history so the simplified response is judged against
    the nonlinear response scale.
    """
    table = indices.copy() if isinstance(indices, pd.DataFrame) else pd.read_csv(indices)
    row = table.iloc[int(row_index)]

    stab_path = Path(str(row["stab_path"]))
    delta_r = float(row.get("sixdof_delta_r", row.get("linear_delta_r", math.radians(5.0))))
    t_final = float(row.get("sixdof_t_final", row.get("linear_t_final_requested", 10.0)))
    target_delta_phi = float(row.get("sixdof_target_delta_phi", row.get("linear_target_delta_phi", math.radians(5.0))))
    stop_at_target = bool(row.get("sixdof_stop_at_target_delta_phi", False))

    Ixx = float(inertia["Ixx"])
    Izz = float(inertia["Izz"])
    Ixz = float(inertia.get("Ixz", 0.0))

    history_path = None
    if history_path_column in row.index and not pd.isna(row[history_path_column]):
        candidate = Path(str(row[history_path_column]))
        if candidate.exists():
            history_path = candidate

    if history_path is not None:
        sixdof = pd.read_csv(history_path)
    else:
        if not recompute_6dof_if_missing:
            raise FileNotFoundError(
                f"No existing 6DOF history was found in column {history_path_column!r}. "
                "Set recompute_6dof_if_missing=True or write histories during postprocess."
            )
        if "Iyy" not in inertia:
            raise KeyError("inertia must contain Iyy when recomputing the 6DOF history.")
        sixdof = simulate_6dof_rudder_step_from_stab(
            stab_path,
            mass=float(mass),
            Ixx=Ixx,
            Iyy=float(inertia["Iyy"]),
            Izz=Izz,
            Ixz=Ixz,
            delta_r=delta_r,
            t_final=t_final,
            control_map=control_map,
            theta_hold=True,
            theta_hold_kp=0.3,
            theta_hold_kq=0.8,
            rho=rho,
            max_step=max_step,
            rtol=rtol,
            atol=atol,
            stop_at_target_delta_phi=stop_at_target,
            target_delta_phi=target_delta_phi,
        )

    reduced = simulate_reduced_lateral_response_from_stab(
        stab_path,
        Ixx=Ixx,
        Izz=Izz,
        Ixz=Ixz,
        delta_r=delta_r,
        time=sixdof["time"].to_numpy(dtype=float),
        rho=rho,
        control_map=control_map,
        include_roll_damping=include_roll_damping,
        include_yaw_damping=include_yaw_damping,
        rtol=rtol,
        atol=atol,
    )

    import matplotlib.pyplot as plt

    angle_factor = 180.0 / math.pi if degrees else 1.0
    angle_unit = "deg" if degrees else "rad"
    plot_specs = [
        ("beta", angle_factor, rf"$\beta$ [{angle_unit}]"),
        ("p", 1.0, r"$p$ [rad/s]"),
        ("r", 1.0, r"$r$ [rad/s]"),
        ("phi", angle_factor, rf"$\phi$ [{angle_unit}]"),
    ]

    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(10, 9))
    for ax, (column, factor, ylabel) in zip(axes, plot_specs):
        y_6dof = sixdof[column].to_numpy(dtype=float) * factor
        y_reduced = reduced[column].to_numpy(dtype=float) * factor
        ax.plot(sixdof["time"], y_6dof, label="6DOF")
        reduced_label = "reduced"
        if include_roll_damping or include_yaw_damping:
            enabled = []
            if include_roll_damping:
                enabled.append("roll damping")
            if include_yaw_damping:
                enabled.append("yaw damping")
            reduced_label += " (" + ", ".join(enabled) + ")"
        ax.plot(reduced["time"], y_reduced, linestyle="--", label=reduced_label)
        ymin = float(np.nanmin(y_6dof))
        ymax = float(np.nanmax(y_6dof))
        if math.isclose(ymin, ymax, rel_tol=0.0, abs_tol=1.0e-12):
            margin = max(abs(ymin), 1.0) * 0.05
        else:
            margin = 0.05 * (ymax - ymin)
        ax.set_ylim(ymin - margin, ymax + margin)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.35)
        ax.legend(loc="best")

    case_name = row.get("case", row_index)
    axes[-1].set_xlabel("t [s]")
    fig.suptitle(f"6DOF vs reduced lateral response: {case_name}")
    fig.tight_layout()
    if plot_path is not None:
        plot_path = Path(plot_path)
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)

    return {
        "fig": fig,
        "axes": axes,
        "sixdof_history": sixdof,
        "reduced_history": reduced,
        "row": row,
        "history_path": None if history_path is None else str(history_path),
        "include_roll_damping": bool(include_roll_damping),
        "include_yaw_damping": bool(include_yaw_damping),
    }

def estimate_roll_rate_gain_from_history(history: pd.DataFrame, delta_r: float, *, evaluation_window: tuple[float, float] | None = None) -> dict[str, float]:
    if abs(delta_r) < 1.0e-14:
        raise ValueError("delta_r must be non-zero for gain estimation.")
    if evaluation_window is None:
        t0 = float(history["time"].iloc[0])
        t1 = float(history["time"].iloc[-1])
    else:
        t0, t1 = [float(value) for value in evaluation_window]

    window = history[(history["time"] >= t0) & (history["time"] <= t1)].copy()
    if len(window) < 2:
        raise ValueError("evaluation_window must contain at least two time samples.")

    time = window["time"].to_numpy(dtype=float)
    phi = window["phi"].to_numpy(dtype=float)
    phi_slope = float(np.polyfit(time, phi, 1)[0])
    p_mean = float(window["p"].mean())
    phidot_mean = float(window["phi_dot"].mean())
    phat_mean = float(window["phat"].mean())
    return {
        "evaluation_t0": t0,
        "evaluation_t1": t1,
        "sample_count": int(len(window)),
        "V_mean": float(window["V"].mean()),
        "V_min": float(window["V"].min()),
        "V_max": float(window["V"].max()),
        "alpha_mean": float(window["alpha"].mean()),
        "beta_mean": float(window["beta"].mean()),
        "theta_mean": float(window["theta"].mean()),
        "phi_mean": float(window["phi"].mean()),
        "delta_e": float(window["delta_e"].iloc[0]),
        "delta_a": float(window["delta_a"].iloc[0]),
        "delta_r": float(window["delta_r"].iloc[0]),
        "thrust": float(window["thrust"].iloc[0]),
        "p_mean": p_mean,
        "phi_dot_mean": phidot_mean,
        "phi_slope": phi_slope,
        "phat_mean": phat_mean,
        "p_per_delta_r_sim": p_mean / float(delta_r),
        "phidot_per_delta_r_sim": phidot_mean / float(delta_r),
        "phi_slope_per_delta_r_sim": phi_slope / float(delta_r),
        "phat_per_delta_r_sim": phat_mean / float(delta_r),
    }

def calculate_roll_response_index_by_delta_phi(
    history: pd.DataFrame,
    *,
    delta_r: float,
    target_delta_phi: float = math.radians(5.0),
    V: float,
    Bref: float,
    phi0: float | None = None,
) -> dict[str, Any]:
    """Calculate a signed delta-phi roll-response gain from a 6DOF history.

    The target is reached when phi(t) - phi0 first equals the signed
    target_delta_phi.  The reach time is linearly interpolated between the two
    surrounding samples.  If the target is not reached, the primary index is
    NaN and a reference index based on (phi_final - phi0) / t_final is still
    returned.  The V used in b/(2V) is supplied by the caller, normally the
    .stab Vinf, not the time-history average.
    """
    for column in ("time", "phi"):
        if column not in history.columns:
            raise KeyError(f"history is missing required column: {column!r}")
    if abs(float(delta_r)) < 1.0e-14:
        raise ValueError("delta_r must be non-zero.")
    if abs(float(target_delta_phi)) < 1.0e-14:
        raise ValueError("target_delta_phi must be non-zero.")
    if float(V) <= 0.0:
        raise ValueError("V must be positive.")
    if float(Bref) <= 0.0:
        raise ValueError("Bref must be positive.")

    time_values = history["time"].to_numpy(dtype=float)
    phi_values = history["phi"].to_numpy(dtype=float)
    if len(time_values) < 2:
        raise ValueError("history must contain at least two time samples.")

    start_time = float(time_values[0])
    end_time = float(time_values[-1])
    total_dt = end_time - start_time
    if total_dt <= 0.0:
        raise ValueError("history time must increase.")

    phi_start = float(phi_values[0] if phi0 is None else phi0)
    phi_delta = phi_values - phi_start
    target = float(target_delta_phi)
    sign = 1.0 if target > 0.0 else -1.0

    reached_indices = np.flatnonzero(sign * phi_delta >= sign * target)
    reached = bool(len(reached_indices) > 0)
    t_reach = math.nan
    if reached:
        idx = int(reached_indices[0])
        if idx == 0:
            t_reach = float(time_values[0])
        else:
            t0 = float(time_values[idx - 1])
            t1 = float(time_values[idx])
            y0 = float(phi_delta[idx - 1])
            y1 = float(phi_delta[idx])
            if abs(y1 - y0) < 1.0e-14:
                t_reach = t1
            else:
                t_reach = t0 + (target - y0) * (t1 - t0) / (y1 - y0)

    phi_final = float(phi_values[-1])
    phi_delta_final = phi_final - phi_start
    reference_phi_rate = phi_delta_final / total_dt
    scale = float(Bref) / (2.0 * float(V))
    result: dict[str, Any] = {
        "sixdof_roll_response_reached": reached,
        "sixdof_target_delta_phi": target,
        "sixdof_target_delta_phi_deg": math.degrees(target),
        "sixdof_phi0": phi_start,
        "sixdof_phi_final": phi_final,
        "sixdof_phi_delta_final": phi_delta_final,
        "sixdof_t_start": start_time,
        "sixdof_t_final": end_time,
        "sixdof_t_reach": t_reach,
        "sixdof_dt_reach": math.nan,
        "sixdof_delta_r": float(delta_r),
        "sixdof_index_Vinf": float(V),
        "sixdof_index_Bref": float(Bref),
        "sixdof_roll_response_phi_rate": math.nan,
        "sixdof_roll_response_phi_rate_per_delta_r": math.nan,
        "sixdof_finite_time_roll_index": math.nan,
        "sixdof_roll_response_error": "" if reached else "target_delta_phi_not_reached",
        "sixdof_roll_response_reference_phi_rate_to_t_final": reference_phi_rate,
        "sixdof_roll_response_reference_phi_rate_per_delta_r_to_t_final": reference_phi_rate / float(delta_r),
        "sixdof_roll_response_index_reference": scale * reference_phi_rate / float(delta_r),
        "sixdof_roll_response_fraction_of_target": phi_delta_final / target,
    }
    if reached:
        dt_reach = t_reach - start_time
        result["sixdof_dt_reach"] = dt_reach
        if dt_reach <= 0.0:
            result["sixdof_roll_response_error"] = "target_reached_at_or_before_start"
        else:
            phi_rate = target / dt_reach
            result["sixdof_roll_response_phi_rate"] = phi_rate
            result["sixdof_roll_response_phi_rate_per_delta_r"] = phi_rate / float(delta_r)
            result["sixdof_finite_time_roll_index"] = scale * phi_rate / float(delta_r)
    return result

def compare_quasi_steady_and_6dof_rudder_roll_gain(
    stab_path: str | Path,
    *,
    mass: float,
    Ixx: float,
    Iyy: float,
    Izz: float,
    Ixz: float = 0.0,
    delta_r: float,
    t_final: float = 5.0,
    evaluation_mode: str = "delta_phi",
    target_delta_phi: float = math.radians(5.0),
    stop_at_target_delta_phi: bool = False,
    evaluation_window: tuple[float, float] | None = None,
    index_V: float | None = None,
    index_Bref: float | None = None,
    control_map: Mapping[str, str] | None = None,
    delta_a: float = 0.0,
    delta_e: float | None = None,
    trim_elevator: bool = True,
    theta_hold: bool = False,
    theta_ref: float | None = None,
    theta_hold_kp: float = 0.0,
    theta_hold_kq: float = 0.0,
    delta_e_min: float | None = None,
    delta_e_max: float | None = None,
    thrust: float | None = None,
    trim_thrust: bool = True,
    g: float = 9.80665,
    rho: float | None = None,
    phi0: float = 0.0,
    theta0: float | None = None,
    psi0: float = 0.0,
    max_step: float | None = None,
    write_history_csv: bool = False,
    history_csv_path: str | Path | None = None,
    plot_history: bool = False,
    history_plot_path: str | Path | None = None,
    show_plot: bool = False,
    plot_degrees: bool = True,
) -> dict[str, Any]:
    """Compare quasi-steady and finite-time 6DOF rudder roll response.

    evaluation_mode='delta_phi' uses the signed target-delta-phi definition
    used by the Vv-Gamma post-processing workflow.  evaluation_mode='window'
    preserves the older window-average estimate for legacy checks.
    """
    if evaluation_mode not in {"delta_phi", "window"}:
        raise ValueError("evaluation_mode must be 'delta_phi' or 'window'.")

    quasi = calculate_quasi_steady_rudder_roll_gain_from_stab(stab_path, control_map=control_map)
    history = simulate_6dof_rudder_step_from_stab(
        stab_path,
        mass=mass,
        Ixx=Ixx,
        Iyy=Iyy,
        Izz=Izz,
        Ixz=Ixz,
        delta_r=delta_r,
        t_final=t_final,
        control_map=control_map,
        delta_a=delta_a,
        delta_e=delta_e,
        trim_elevator=trim_elevator,
        theta_hold=theta_hold,
        theta_ref=theta_ref,
        theta_hold_kp=theta_hold_kp,
        theta_hold_kq=theta_hold_kq,
        delta_e_min=delta_e_min,
        delta_e_max=delta_e_max,
        thrust=thrust,
        trim_thrust=trim_thrust,
        g=g,
        rho=rho,
        phi0=phi0,
        theta0=theta0,
        psi0=psi0,
        max_step=max_step,
        stop_at_target_delta_phi=bool(stop_at_target_delta_phi and evaluation_mode == "delta_phi"),
        target_delta_phi=target_delta_phi,
    )

    if evaluation_mode == "window":
        estimated = estimate_roll_rate_gain_from_history(history, delta_r, evaluation_window=evaluation_window)
    else:
        estimated = calculate_roll_response_index_by_delta_phi(
            history,
            delta_r=delta_r,
            target_delta_phi=target_delta_phi,
            V=float(quasi["V"] if index_V is None else index_V),
            Bref=float(quasi["Bref"] if index_Bref is None else index_Bref),
            phi0=phi0,
        )

    output_paths: dict[str, str] = {}
    if write_history_csv:
        if history_csv_path is None:
            history_csv_path = Path(stab_path).with_suffix(".6dof_history.csv")
        output_paths["history_csv_path"] = str(write_6dof_history_csv(history, history_csv_path))
    if plot_history or history_plot_path is not None or show_plot:
        if history_plot_path is None:
            history_plot_path = Path(stab_path).with_suffix(".6dof_history.png")
        plot_6dof_history(history, plot_path=history_plot_path, show=show_plot, degrees=plot_degrees)
        output_paths["history_plot_path"] = str(history_plot_path)

    def relative_error(reference: float, value: float) -> float:
        return math.nan if abs(reference) < 1.0e-14 else (value - reference) / reference

    if evaluation_mode == "window":
        differences = {
            "p_per_delta_r_sim_minus_stab": estimated["p_per_delta_r_sim"] - quasi["p_per_delta_r"],
            "p_per_delta_r_relative_error": relative_error(quasi["p_per_delta_r"], estimated["p_per_delta_r_sim"]),
            "phat_per_delta_r_sim_minus_stab": estimated["phat_per_delta_r_sim"] - quasi["phat_per_delta_r"],
            "phat_per_delta_r_relative_error": relative_error(quasi["phat_per_delta_r"], estimated["phat_per_delta_r_sim"]),
            "phidot_per_delta_r_sim_minus_p_stab": estimated["phidot_per_delta_r_sim"] - quasi["p_per_delta_r"],
            "phi_slope_per_delta_r_sim_minus_p_stab": estimated["phi_slope_per_delta_r_sim"] - quasi["p_per_delta_r"],
        }
    else:
        index = float(estimated["sixdof_finite_time_roll_index"])
        index_reference = float(estimated["sixdof_roll_response_index_reference"])
        p_scale = 2.0 * float(estimated["sixdof_index_Vinf"]) / float(estimated["sixdof_index_Bref"])
        differences = {
            "phat_per_delta_r_sim_minus_stab": index - quasi["phat_per_delta_r"],
            "phat_per_delta_r_relative_error": relative_error(quasi["phat_per_delta_r"], index),
            "p_per_delta_r_from_delta_phi_minus_stab": p_scale * index - quasi["p_per_delta_r"],
            "p_per_delta_r_from_delta_phi_relative_error": relative_error(quasi["p_per_delta_r"], p_scale * index),
            "phat_reference_per_delta_r_minus_stab": index_reference - quasi["phat_per_delta_r"],
            "phat_reference_per_delta_r_relative_error": relative_error(quasi["phat_per_delta_r"], index_reference),
        }

    return {
        "evaluation_mode": evaluation_mode,
        "quasi_steady": quasi,
        "simulation_gain": estimated,
        "differences": differences,
        "output_paths": output_paths,
        "history": history,
    }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Calculate quasi-steady rudder roll-rate gain from a VSPAERO .stab file.")
    parser.add_argument("stab_path")
    args = parser.parse_args()
    result = calculate_quasi_steady_rudder_roll_gain_from_stab(args.stab_path)
    for key in ["beta_per_delta_r", "phat_per_delta_r", "rhat_per_delta_r", "p_per_delta_r"]:
        print(f"{key}: {result[key]:.10g}")
