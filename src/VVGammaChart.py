"""
Vv-Gamma chart utilities for OpenVSP / VSPAERO rudder-only turn studies.

The code keeps the main flow deliberately direct:

    VTail Vv update -> wing elastic dihedral update -> VSPAERO stability run
    -> .stab parsing -> Vv-Gamma metrics -> optional turn-trim metric

Angles are radians in the numerical/flight-mechanics functions and degrees only
when values are written into OpenVSP wing-section parameters.
"""

from __future__ import annotations

import importlib
import math
import os
import shutil
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Small import helpers
# ---------------------------------------------------------------------------

def _import_openvsp():
    try:
        return importlib.import_module("openvsp")
    except ImportError as exc:
        raise ImportError(
            "openvsp is required for .vsp3 geometry editing and VSPAERO runs. "
            "Install/use this module inside an OpenVSP Python environment."
        ) from exc

def _import_turntrim():
    try:
        return importlib.import_module(".TurnTrim", package=__package__)
    except Exception:
        return importlib.import_module("TurnTrim")

def _import_analysis_vspaero():
    try:
        return importlib.import_module(".AnalysisVSPAERO", package=__package__)
    except Exception:
        return importlib.import_module("AnalysisVSPAERO")

@contextmanager
def _workdir(path: str | os.PathLike):
    previous = Path.cwd()
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)

# ---------------------------------------------------------------------------
# OpenVSP low-level helpers
# ---------------------------------------------------------------------------

def _find_one_geom(vsp, geom_name: str) -> str:
    geom_ids = list(vsp.FindGeomsWithName(geom_name))
    if len(geom_ids) != 1:
        raise ValueError(f"Geom '{geom_name}' must exist exactly once. found={geom_ids}")
    return geom_ids[0]

def _get_xsec_surf_id(vsp, geom_id: str, surf_index: int = 0) -> str:
    return vsp.GetXSecSurf(geom_id, surf_index)

def _xsec_indices(vsp, xsec_surf_id: str) -> range:
    # OpenVSP wing section parameters live on XSec 1..N-1.
    return range(1, int(vsp.GetNumXSec(xsec_surf_id)))

def _get_xsec_parm_id(vsp, xsec_id: str, parm_names: Sequence[str]) -> tuple[str, str]:
    for name in parm_names:
        try:
            parm_id = vsp.GetXSecParm(xsec_id, name)
        except Exception:
            parm_id = ""
        if parm_id:
            return parm_id, name
    raise KeyError(f"None of these XSec parameters were found: {list(parm_names)}")

def _get_xsec_value(vsp, xsec_id: str, parm_names: Sequence[str], default: float | None = None) -> float:
    try:
        parm_id, _ = _get_xsec_parm_id(vsp, xsec_id, parm_names)
    except KeyError:
        if default is None:
            raise
        return float(default)
    return float(vsp.GetParmVal(parm_id))

def _set_xsec_value(vsp, xsec_id: str, parm_names: Sequence[str], value: float) -> None:
    parm_id, _ = _get_xsec_parm_id(vsp, xsec_id, parm_names)
    vsp.SetParmVal(parm_id, float(value))

def _set_wing_section_driver_for_scaling(vsp, geom_id: str, section_index: int) -> None:
    # Keep this helper narrow and explicit: VTail scaling is easier to read when
    # every section is driven by Span, Root_Chord, and Tip_Chord directly.
    try:
        vsp.SetDriverGroup(
            geom_id,
            section_index,
            vsp.SPAN_WSECT_DRIVER,
            vsp.ROOTC_WSECT_DRIVER,
            vsp.TIPC_WSECT_DRIVER,
        )
    except Exception:
        # Some OpenVSP builds may already have the desired driver group and may
        # object to setting it before all parameters are valid.  The subsequent
        # SetParmVal calls will still fail loudly if the section cannot be edited.
        pass

def _geom_x_location(vsp, geom_id: str) -> float:
    parm_id = vsp.FindParm(geom_id, "X_Rel_Location", "XForm")
    if not parm_id:
        parm_id = vsp.FindParm(geom_id, "X_Location", "XForm")
    if not parm_id:
        raise KeyError("Could not find X_Rel_Location/X_Location on Geom XForm.")
    return float(vsp.GetParmVal(parm_id))

def _set_geom_x_location(vsp, geom_id: str, value: float) -> None:
    parm_id = vsp.FindParm(geom_id, "X_Rel_Location", "XForm")
    if not parm_id:
        parm_id = vsp.FindParm(geom_id, "X_Location", "XForm")
    if not parm_id:
        raise KeyError("Could not find X_Rel_Location/X_Location on Geom XForm.")
    vsp.SetParmVal(parm_id, float(value))

# Common OpenVSP wing-section parameter names.  The first name is preferred;
# aliases keep the code usable across small API/XML naming differences.
SPAN_NAMES = ("Span",)
ROOT_CHORD_NAMES = ("Root_Chord", "RootC")
TIP_CHORD_NAMES = ("Tip_Chord", "TipC")
SWEEP_NAMES = ("Sweep",)
SWEEP_LOCATION_NAMES = ("Sweep_Location", "SweepLoc")
DIHEDRAL_NAMES = ("Dihedral",)
AREA_NAMES = ("Area",)

# ---------------------------------------------------------------------------
# Geometry summaries used by both VTail scaling and wing deflection
# ---------------------------------------------------------------------------

def _read_wing_section_table_from_loaded(vsp, geom_id: str, *, surf_index: int = 0) -> pd.DataFrame:
    xsec_surf_id = _get_xsec_surf_id(vsp, geom_id, surf_index)

    rows = []
    root_le_x_local = 0.0
    root_y = 0.0
    for section_index in _xsec_indices(vsp, xsec_surf_id):
        xsec_id = vsp.GetXSec(xsec_surf_id, section_index)
        span = _get_xsec_value(vsp, xsec_id, SPAN_NAMES)
        root_chord = _get_xsec_value(vsp, xsec_id, ROOT_CHORD_NAMES)
        tip_chord = _get_xsec_value(vsp, xsec_id, TIP_CHORD_NAMES)
        sweep_deg = _get_xsec_value(vsp, xsec_id, SWEEP_NAMES, default=0.0)
        sweep_location = _get_xsec_value(vsp, xsec_id, SWEEP_LOCATION_NAMES, default=0.25)
        if sweep_location > 1.0:
            sweep_location = sweep_location / 100.0
        dihedral_deg = _get_xsec_value(vsp, xsec_id, DIHEDRAL_NAMES, default=0.0)

        tip_le_x_local = (
            root_le_x_local
            + span * math.tan(math.radians(sweep_deg))
            - sweep_location * (tip_chord - root_chord)
        )
        area = 0.5 * (root_chord + tip_chord) * span
        root_qc_x = root_le_x_local + 0.25 * root_chord
        tip_qc_x = tip_le_x_local + 0.25 * tip_chord
        qc_x = 0.5 * (root_qc_x + tip_qc_x)

        rows.append(
            {
                "section_index": section_index,
                "span": span,
                "root_chord": root_chord,
                "tip_chord": tip_chord,
                "taper": tip_chord / root_chord if abs(root_chord) > 1e-12 else np.nan,
                "sweep_deg": sweep_deg,
                "sweep_location": sweep_location,
                "dihedral_deg": dihedral_deg,
                "area": area,
                "root_y": root_y,
                "tip_y": root_y + span,
                "center_y": root_y + 0.5 * span,
                "root_le_x_local": root_le_x_local,
                "tip_le_x_local": tip_le_x_local,
                "quarter_chord_x_local": qc_x,
            }
        )

        root_le_x_local = tip_le_x_local
        root_y += span

    return pd.DataFrame(rows)

def read_wing_section_table(vsp3_path: str | os.PathLike, geom_name: str, *, surf_index: int = 0) -> pd.DataFrame:
    """Read editable wing-section geometry from an OpenVSP WING Geom.

    The returned rows correspond to OpenVSP section indices 1..GetNumXSec()-1.
    This function needs OpenVSP at runtime but does not run VSPAERO.
    """

    vsp = _import_openvsp()
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(vsp3_path))
    vsp.Update()

    geom_id = _find_one_geom(vsp, geom_name)
    return _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)

def wing_planform_summary(section_table: pd.DataFrame) -> dict[str, float]:
    area = float(section_table["area"].sum())
    span = float(section_table["span"].sum())
    quarter_chord_x_local = float(
        (section_table["area"] * section_table["quarter_chord_x_local"]).sum() / area
    )
    aspect_ratio = span * span / area if area > 0.0 else math.nan
    return {
        "area": area,
        "span": span,
        "aspect_ratio": aspect_ratio,
        "quarter_chord_x_local": quarter_chord_x_local,
    }

# ---------------------------------------------------------------------------
# VTail Vv update
# ---------------------------------------------------------------------------

def update_vtail_volume(
    input_vsp3_path: str | os.PathLike,
    output_vsp3_path: str | os.PathLike,
    vv: float,
    *,
    lv: float,
    wing_area: float,
    wing_span: float,
    xcg: float,
    vtail_name: str = "VTailGeom",
    surf_index: int = 0,
) -> dict[str, Any]:
    """Scale VTailGeom to a target vertical-tail volume and keep lv fixed.

    The VTail planform is scaled about the root leading edge.  Every section's
    Span, Root_Chord, and Tip_Chord is multiplied by the same factor, preserving
    section taper ratios and total aspect ratio.  After scaling, the whole Geom
    is translated in X so that the area-weighted quarter-chord representative
    aerodynamic center satisfies x_ac - xcg = lv.
    """

    if lv <= 0.0:
        raise ValueError("lv must be positive.")
    if wing_area <= 0.0 or wing_span <= 0.0:
        raise ValueError("wing_area and wing_span must be positive.")
    if vv <= 0.0:
        raise ValueError("vv must be positive.")

    vsp = _import_openvsp()
    input_vsp3_path = Path(input_vsp3_path)
    output_vsp3_path = Path(output_vsp3_path)
    output_vsp3_path.parent.mkdir(parents=True, exist_ok=True)

    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(input_vsp3_path))
    vsp.Update()

    geom_id = _find_one_geom(vsp, vtail_name)
    xsec_surf_id = _get_xsec_surf_id(vsp, geom_id, surf_index)

    before_sections = _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)
    before_summary = wing_planform_summary(before_sections)
    target_area = float(vv) * float(wing_area) * float(wing_span) / float(lv)
    scale = math.sqrt(target_area / before_summary["area"])

    for section_index in _xsec_indices(vsp, xsec_surf_id):
        xsec_id = vsp.GetXSec(xsec_surf_id, section_index)
        _set_wing_section_driver_for_scaling(vsp, geom_id, section_index)
        _set_xsec_value(vsp, xsec_id, SPAN_NAMES, _get_xsec_value(vsp, xsec_id, SPAN_NAMES) * scale)
        _set_xsec_value(vsp, xsec_id, ROOT_CHORD_NAMES, _get_xsec_value(vsp, xsec_id, ROOT_CHORD_NAMES) * scale)
        _set_xsec_value(vsp, xsec_id, TIP_CHORD_NAMES, _get_xsec_value(vsp, xsec_id, TIP_CHORD_NAMES) * scale)

    vsp.Update()
    after_sections = _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)
    after_summary = wing_planform_summary(after_sections)

    x_location_before_shift = _geom_x_location(vsp, geom_id)
    x_ac_before_shift = x_location_before_shift + after_summary["quarter_chord_x_local"]
    x_ac_target = float(xcg) + float(lv)
    delta_x = x_ac_target - x_ac_before_shift
    _set_geom_x_location(vsp, geom_id, x_location_before_shift + delta_x)
    vsp.Update()
    vsp.WriteVSPFile(os.fspath(output_vsp3_path))

    return {
        "input_vsp3_path": os.fspath(input_vsp3_path),
        "output_vsp3_path": os.fspath(output_vsp3_path),
        "vtail_name": vtail_name,
        "vv_target": float(vv),
        "lv_target": float(lv),
        "wing_area": float(wing_area),
        "wing_span": float(wing_span),
        "scale": scale,
        "area_before": before_summary["area"],
        "area_target": target_area,
        "area_after": after_summary["area"],
        "aspect_ratio_before": before_summary["aspect_ratio"],
        "aspect_ratio_after": after_summary["aspect_ratio"],
        "x_ac_target": x_ac_target,
        "x_ac_before_shift": x_ac_before_shift,
        "x_shift": delta_x,
        "x_location_after": x_location_before_shift + delta_x,
    }

# ---------------------------------------------------------------------------
# Elastic wing deflection and Dihedral approximation
# ---------------------------------------------------------------------------

def elastic_wing_deflection_distribution(
    tip_deflection: float,
    semispan: float,
    *,
    n_span: int = 101,
) -> pd.DataFrame:
    """Return elliptical-load, constant-EI cantilever deflection distribution.

    tip_deflection is the dimensional wing-tip deflection.  semispan uses the
    same length unit.  Positive tip_deflection produces positive elastic
    dihedral.
    """

    if semispan <= 0.0:
        raise ValueError("semispan must be positive.")
    if n_span < 2:
        raise ValueError("n_span must be at least 2.")

    eta = np.linspace(0.0, 1.0, int(n_span))
    eta_clipped = np.clip(eta, 0.0, 1.0)
    sqrt_term = np.sqrt(np.maximum(1.0 - eta_clipped * eta_clipped, 0.0))
    asin_term = np.arcsin(eta_clipped)

    # Dimensionless shape terms from the article derivation.  The common
    # load/EI factor is eliminated by normalizing z(1) to tip_deflection.
    theta_shape = (
        ((4.0 * eta_clipped**2 + 1.0) / 16.0) * asin_term
        - (math.pi * eta_clipped**2) / 8.0
        + eta_clipped * sqrt_term * (13.0 + 2.0 * eta_clipped**2) / 48.0
    )
    z_shape = (
        eta_clipped * (4.0 * eta_clipped**2 + 3.0) * asin_term / 48.0
        - math.pi * eta_clipped**3 / 24.0
        + sqrt_term * (6.0 * eta_clipped**4 + 83.0 * eta_clipped**2 + 16.0) / 720.0
        - 1.0 / 45.0
    )

    z_tip_shape = float(z_shape[-1])
    if abs(z_tip_shape) < 1e-15:
        raise ZeroDivisionError("The normalized deflection shape has zero tip value.")

    z = float(tip_deflection) * z_shape / z_tip_shape
    theta_rad = float(tip_deflection) * theta_shape / (float(semispan) * z_tip_shape)

    return pd.DataFrame(
        {
            "eta": eta,
            "y": eta * float(semispan),
            "z": z,
            "theta_rad": theta_rad,
            "theta_deg": np.degrees(theta_rad),
        }
    )

def apply_wing_deflection_as_dihedral(
    input_vsp3_path: str | os.PathLike,
    output_vsp3_path: str | os.PathLike,
    deflection_distribution: pd.DataFrame,
    *,
    wing_name: str = "WingGeom",
    surf_index: int = 0,
) -> dict[str, Any]:
    """Approximate elastic deflection by adding local theta_b to Dihedral.

    Each editable wing section receives the elastic angle evaluated at the
    section center.  The existing rigid Dihedral is preserved and incremented.
    """

    required_columns = {"y", "theta_deg"}
    missing = required_columns - set(deflection_distribution.columns)
    if missing:
        raise ValueError(f"deflection_distribution is missing columns: {sorted(missing)}")

    vsp = _import_openvsp()
    input_vsp3_path = Path(input_vsp3_path)
    output_vsp3_path = Path(output_vsp3_path)
    output_vsp3_path.parent.mkdir(parents=True, exist_ok=True)

    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(input_vsp3_path))
    vsp.Update()

    geom_id = _find_one_geom(vsp, wing_name)
    xsec_surf_id = _get_xsec_surf_id(vsp, geom_id, surf_index)
    original_sections = _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)

    y_values = deflection_distribution["y"].to_numpy(dtype=float)
    theta_values = deflection_distribution["theta_deg"].to_numpy(dtype=float)

    updates = []
    for row in original_sections.itertuples(index=False):
        section_index = int(row.section_index)
        xsec_id = vsp.GetXSec(xsec_surf_id, section_index)
        elastic_dihedral_deg = float(np.interp(float(row.center_y), y_values, theta_values))
        new_dihedral_deg = float(row.dihedral_deg) + elastic_dihedral_deg
        _set_xsec_value(vsp, xsec_id, DIHEDRAL_NAMES, new_dihedral_deg)
        updates.append(
            {
                "section_index": section_index,
                "center_y": float(row.center_y),
                "rigid_dihedral_deg": float(row.dihedral_deg),
                "elastic_dihedral_deg": elastic_dihedral_deg,
                "new_dihedral_deg": new_dihedral_deg,
            }
        )

    vsp.Update()
    vsp.WriteVSPFile(os.fspath(output_vsp3_path))

    return {
        "input_vsp3_path": os.fspath(input_vsp3_path),
        "output_vsp3_path": os.fspath(output_vsp3_path),
        "wing_name": wing_name,
        "updates": updates,
    }

# ---------------------------------------------------------------------------
# .stab based metrics
# ---------------------------------------------------------------------------

def _stab_derivative(stab, coef: str, column: str) -> float:
    if coef not in stab.derivatives.index:
        raise KeyError(f"Coefficient row '{coef}' was not found in .stab derivatives.")
    if column not in stab.derivatives.columns:
        raise KeyError(f"Derivative column '{column}' was not found in .stab derivatives.")
    return float(stab.derivatives.loc[coef, column])

def _safe_divide(numerator: float, denominator: float) -> float:
    return math.nan if abs(float(denominator)) < 1e-14 else float(numerator) / float(denominator)

def calculate_vv_gamma_metrics_from_stab(
    stab_path: str | os.PathLike,
    *,
    vv: float | None = None,
    tip_deflection: float | None = None,
    control_map: Mapping[str, str] | None = None,
    turn_trim_mode: str = "none",
    turn_trim_fixed: Mapping[str, float] | None = None,
    mass: float | None = None,
    rho: float | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    residual_tol: float = 1e-6,
    g: float = 9.80665,
) -> dict[str, Any]:
    """Calculate Vv-Gamma chart metrics from a VSPAERO .stab file."""

    turntrim = _import_turntrim()
    stab = turntrim.read_vspaero_stab(stab_path)
    control_columns = turntrim._control_columns_from_stab(stab, control_map)
    rudder_column = control_columns.get("delta_r")
    aileron_column = control_columns.get("delta_a")
    elevator_column = control_columns.get("delta_e")

    cl_alpha = _stab_derivative(stab, "CL", "Alpha")
    cl_beta = _stab_derivative(stab, "CMl", "Beta")
    cn_beta = _stab_derivative(stab, "CMn", "Beta")
    cl_p = _stab_derivative(stab, "CMl", "p")
    cl_r = _stab_derivative(stab, "CMl", "r")
    cn_p = _stab_derivative(stab, "CMn", "p")
    cn_r = _stab_derivative(stab, "CMn", "r")
    cy_beta = _stab_derivative(stab, "CS", "Beta")

    cl_delta_r = _stab_derivative(stab, "CMl", rudder_column) if rudder_column else math.nan
    cn_delta_r = _stab_derivative(stab, "CMn", rudder_column) if rudder_column else math.nan
    cy_delta_r = _stab_derivative(stab, "CS", rudder_column) if rudder_column else math.nan

    gamma_eff_rad = -3.0 * math.pi * cl_beta / (2.0 * cl_alpha)
    spiral_margin = cl_beta * cn_r - cn_beta * cl_r
    beta_per_delta_r_yaw_static = -_safe_divide(cn_delta_r, cn_beta)
    delta_r_per_beta_yaw_static = _safe_divide(1.0, beta_per_delta_r_yaw_static)
    p_hat_per_delta_r_beta_only = -_safe_divide(
        cl_delta_r + cl_beta * beta_per_delta_r_yaw_static,
        cl_p,
    )

    metrics: dict[str, Any] = {
        "stab_path": os.fspath(stab_path),
        "Vv": np.nan if vv is None else float(vv),
        "tip_deflection": np.nan if tip_deflection is None else float(tip_deflection),
        "Gamma_eff_rad": gamma_eff_rad,
        "Gamma_eff_deg": math.degrees(gamma_eff_rad),
        "CL_alpha": cl_alpha,
        "Cl_beta": cl_beta,
        "Cn_beta": cn_beta,
        "CY_beta": cy_beta,
        "Cl_p": cl_p,
        "Cl_r": cl_r,
        "Cn_p": cn_p,
        "Cn_r": cn_r,
        "Cl_delta_r": cl_delta_r,
        "Cn_delta_r": cn_delta_r,
        "CY_delta_r": cy_delta_r,
        "spiral_margin": spiral_margin,
        "beta_per_delta_r_yaw_static": beta_per_delta_r_yaw_static,
        "delta_r_per_beta_yaw_static": delta_r_per_beta_yaw_static,
        "p_hat_per_delta_r_beta_only": p_hat_per_delta_r_beta_only,
        "control_column_delta_a": aileron_column,
        "control_column_delta_e": elevator_column,
        "control_column_delta_r": rudder_column,
        "control_groups": dict(stab.control_groups),
        "Sref": stab.Sref,
        "Bref": stab.Bref,
        "Cref": stab.Cref,
        "Xcg": stab.references.get("Xcg"),
        "Ycg": stab.references.get("Ycg"),
        "Zcg": stab.references.get("Zcg"),
        "AoA_deg": stab.base_condition.get("AoA"),
        "Beta_deg": stab.base_condition.get("Beta"),
        "Mach": stab.base_condition.get("Mach"),
        "Rho": stab.base_condition.get("Rho"),
        "Vinf": stab.base_condition.get("Vinf"),
        "turn_trim_mode": turn_trim_mode,
        "turn_trim_passed": np.nan,
    }

    if turn_trim_mode not in {"none", "level", "gliding"}:
        raise ValueError("turn_trim_mode must be 'none', 'level', or 'gliding'.")

    if turn_trim_mode != "none":
        if mass is None:
            raise ValueError("mass is required when turn_trim_mode is 'level' or 'gliding'.")
        if turn_trim_fixed is None:
            raise ValueError("turn_trim_fixed is required when turn_trim_mode is 'level' or 'gliding'.")
        try:
            if turn_trim_mode == "level":
                trim = turntrim.solve_steady_level_turn(
                    turn_trim_fixed,
                    stab_path,
                    mass,
                    rho=rho,
                    initial_guess=initial_guess,
                    g=g,
                    bounds=bounds,
                    control_map=control_map,
                    residual_tol=residual_tol,
                    verbose=0,
                )
            else:
                trim = turntrim.solve_steady_gliding_turn(
                    turn_trim_fixed,
                    stab_path,
                    mass,
                    rho=rho,
                    initial_guess=initial_guess,
                    g=g,
                    bounds=bounds,
                    control_map=control_map,
                    residual_tol=residual_tol,
                    verbose=0,
                )
            metrics.update(_flatten_turn_trim(trim))
        except Exception as exc:
            metrics.update(
                {
                    "turn_trim_passed": False,
                    "turn_trim_error": repr(exc),
                }
            )

    return metrics

def _flatten_turn_trim(trim: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "turn_trim_passed": bool(trim.get("passed", False)),
        "turn_trim_message": trim.get("message", ""),
        "turn_trim_max_abs_residual": trim.get("max_abs_residual", np.nan),
        "turn_trim_cost": trim.get("cost", np.nan),
        "turn_trim_nfev": trim.get("nfev", np.nan),
    }
    for source_name in ("solution", "derived", "coefficients", "residuals"):
        values = trim.get(source_name, {}) or {}
        for key, value in values.items():
            if isinstance(value, (int, float, np.floating, np.integer, bool, str)) or value is None:
                result[f"turn_trim_{source_name}_{key}"] = value
    if "height_residual" in trim:
        result["turn_trim_height_residual"] = trim["height_residual"]
    return result

# ---------------------------------------------------------------------------
# VSPAERO run wrapper and chart loop
# ---------------------------------------------------------------------------

def _find_case_stab_file(case_dir: Path, fallback_stem: str | None = None) -> Path:
    candidates = sorted(case_dir.glob("*.stab"), key=lambda path: path.stat().st_mtime, reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No .stab file was found in {case_dir}.")
    if fallback_stem:
        for path in candidates:
            if path.stem == fallback_stem:
                return path
    return candidates[0]

def run_vspaero_stability_case(
    vsp3_path: str | os.PathLike,
    case_dir: str | os.PathLike,
    *,
    alpha_deg: float,
    mach: float,
    reynolds: float,
    verbose: int | bool = 0,
) -> dict[str, Any]:
    """Run the existing VSPAERO stability function in a case directory."""

    analysis_vspaero = _import_analysis_vspaero()
    case_dir = Path(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    vsp3_path = Path(vsp3_path).resolve()

    with _workdir(case_dir):
        report = analysis_vspaero.vsp_stability_derivatives(
            os.fspath(vsp3_path),
            alpha=float(alpha_deg),
            mach=float(mach),
            reynolds=float(reynolds),
            verbose=verbose,
        )

    report = dict(report)
    report["case_dir"] = os.fspath(case_dir)
    if report.get("passed"):
        stab_path = _find_case_stab_file(case_dir)
        report["stab_path"] = os.fspath(stab_path)
    else:
        report["stab_path"] = ""
    return report

def run_vv_gamma_chart(
    base_vsp3_path: str | os.PathLike,
    vv_values: Sequence[float],
    tip_deflections: Sequence[float],
    flight_condition: Mapping[str, float],
    geometry_config: Mapping[str, Any],
    output_dir: str | os.PathLike,
    *,
    turn_trim_mode: str = "none",
    turn_trim_fixed: Mapping[str, float] | None = None,
    mass: float | None = None,
    rho: float | None = None,
    initial_guess: Mapping[str, float] | None = None,
    bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    residual_tol: float = 1e-6,
    control_map: Mapping[str, str] | None = None,
    validate_base_model: bool = True,
    verbose: int | bool = 0,
) -> pd.DataFrame:
    """Generate a Vv-Gamma result table by sweeping Vv and tip deflection.

    geometry_config requires:
        lv, wing_area, wing_span, xcg

    Optional geometry_config keys:
        wing_name, vtail_name, n_span
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    base_vsp3_path = Path(base_vsp3_path).resolve()

    if turn_trim_mode not in {"none", "level", "gliding"}:
        raise ValueError("turn_trim_mode must be 'none', 'level', or 'gliding'.")

    if validate_base_model:
        analysis_vspaero = _import_analysis_vspaero()
        validation = analysis_vspaero.validate_vsp3_for_stability_derivatives(
            os.fspath(base_vsp3_path),
            verbose=verbose,
        )
        if not validation.get("passed", False):
            raise RuntimeError(f"Base .vsp3 validation failed: {validation.get('errors')}")

    rows = []
    for vv in vv_values:
        for tip_deflection in tip_deflections:
            case_name = f"vv_{float(vv):.5f}_wtip_{float(tip_deflection):.5f}".replace("-", "m").replace(".", "p")
            case_dir = output_dir / case_name
            case_dir.mkdir(parents=True, exist_ok=True)

            vtail_vsp3_path = case_dir / f"{case_name}_vtail.vsp3"
            deformed_vsp3_path = case_dir / f"{case_name}.vsp3"

            row: dict[str, Any] = {
                "case": case_name,
                "Vv": float(vv),
                "tip_deflection": float(tip_deflection),
                "case_dir": os.fspath(case_dir),
                "passed": False,
            }

            try:
                vtail_report = update_vtail_volume(
                    base_vsp3_path,
                    vtail_vsp3_path,
                    float(vv),
                    lv=float(geometry_config["lv"]),
                    wing_area=float(geometry_config["wing_area"]),
                    wing_span=float(geometry_config["wing_span"]),
                    xcg=float(geometry_config["xcg"]),
                    vtail_name=str(geometry_config.get("vtail_name", "VTailGeom")),
                )
                row.update({f"vtail_{key}": value for key, value in vtail_report.items() if isinstance(value, (int, float, str))})

                semispan = float(geometry_config["wing_span"]) / 2.0
                deflection = elastic_wing_deflection_distribution(
                    float(tip_deflection),
                    semispan,
                    n_span=int(geometry_config.get("n_span", 101)),
                )
                deflection.to_csv(case_dir / f"{case_name}_deflection.csv", index=False)

                wing_report = apply_wing_deflection_as_dihedral(
                    vtail_vsp3_path,
                    deformed_vsp3_path,
                    deflection,
                    wing_name=str(geometry_config.get("wing_name", "WingGeom")),
                )
                row["vsp3_path"] = os.fspath(deformed_vsp3_path)
                row["wing_dihedral_update_count"] = len(wing_report["updates"])

                stability_report = run_vspaero_stability_case(
                    deformed_vsp3_path,
                    case_dir,
                    alpha_deg=float(flight_condition["alpha_deg"]),
                    mach=float(flight_condition["mach"]),
                    reynolds=float(flight_condition["reynolds"]),
                    verbose=verbose,
                )
                row["stab_path"] = stability_report.get("stab_path", "")
                row["vspaero_passed"] = bool(stability_report.get("passed", False))
                if not stability_report.get("passed", False):
                    row["error"] = stability_report.get("errors", [])
                    rows.append(row)
                    continue

                metrics = calculate_vv_gamma_metrics_from_stab(
                    stability_report["stab_path"],
                    vv=float(vv),
                    tip_deflection=float(tip_deflection),
                    control_map=control_map,
                    turn_trim_mode=turn_trim_mode,
                    turn_trim_fixed=turn_trim_fixed,
                    mass=mass,
                    rho=rho,
                    initial_guess=initial_guess,
                    bounds=bounds,
                    residual_tol=residual_tol,
                )
                row.update(metrics)
                row["passed"] = True

            except Exception as exc:
                row["error"] = repr(exc)

            rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(output_dir / "vv_gamma_results.csv", index=False)
    return result

# ---------------------------------------------------------------------------
# Lightweight plotting helpers for notebooks
# ---------------------------------------------------------------------------

def plot_vv_gamma_contour(
    results: pd.DataFrame,
    value_column: str,
    *,
    ax=None,
    levels: int | Sequence[float] = 12,
):
    """Plot a Vv-Gamma contour from a completed results DataFrame."""

    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata

    data = results[results.get("passed", True).astype(bool)].copy()
    data = data.dropna(subset=["Vv", "Gamma_eff_deg", value_column])
    if data.empty:
        raise ValueError("No passed rows with Vv, Gamma_eff_deg, and requested value_column were found.")

    xi = np.linspace(data["Vv"].min(), data["Vv"].max(), 80)
    yi = np.linspace(data["Gamma_eff_deg"].min(), data["Gamma_eff_deg"].max(), 80)
    xx, yy = np.meshgrid(xi, yi)
    zz = griddata(
        data[["Vv", "Gamma_eff_deg"]].to_numpy(),
        data[value_column].to_numpy(),
        (xx, yy),
        method="linear",
    )

    if ax is None:
        _, ax = plt.subplots()
    contour = ax.contourf(xx, yy, zz, levels=levels)
    ax.scatter(data["Vv"], data["Gamma_eff_deg"], s=16)
    ax.set_xlabel("Vv")
    ax.set_ylabel("Gamma_eff [deg]")
    ax.set_title(value_column)
    return ax, contour
