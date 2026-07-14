"""
Vv-Gamma chart utilities for OpenVSP / VSPAERO rudder-only turn studies.

The main workflow is intentionally split into two readable stages:

    1. run_vv_wtip_stability_sweep:
       VTail Vv update -> wing elastic dihedral update -> VSPAERO stability run
       -> case table containing .vsp3 and .stab paths

    2. collect_vv_wtip_sweep_progress, optional while a sweep is running:
       expected Vv/wtip case table -> existing .vsp3/.stab paths

    3. postprocess_vv_gamma_cases:
       tip_deflection + semispan -> analytic elliptical-planform weighted
       Gamma_eff, then .stab parsing -> stability indices -> simple rudder-roll index
       -> 6DOF signed-delta-phi roll-response index ->
       crosswind-gust short-time roll-path index ->
       simple maximum-bank estimate -> optional fixed-bank and
       rudder-limit level/gliding turn-trim indices

Angles are radians in the numerical/flight-mechanics functions.  The final
postprocess_vv_gamma_cases() output DataFrame/CSV is a plotting and review
artifact: ordinary angles are written in degrees and ordinary angular rates / turn
rates are written in degrees per second.  Stability derivatives, reduced rates,
gains, and ratios remain in their original nondimensional or per-radian form.
"""

from __future__ import annotations

import importlib
import math
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from .util import (
    DIHEDRAL_NAMES,
    ROOT_CHORD_NAMES,
    SPAN_NAMES,
    SWEEP_LOCATION_NAMES,
    SWEEP_NAMES,
    TIP_CHORD_NAMES,
    XFORM_GROUP_NAMES,
    X_LOCATION_NAMES,
    find_one_geom,
    get_geom_parm_value,
    get_xsec_value,
    import_openvsp,
    set_geom_parm_value,
    set_wing_section_driver_for_scaling,
    set_xsec_value,
    workdir,
)
from .VSPAEROStab import (
    read_vspaero_stab,
    resolve_control_columns_from_stab,
    stab_coefficient_value,
    evaluate_stab_linear_aero, 
    solver_axis_derivative_value,
)
from .TrimTurnSolver import (
    solve_rudder_limit_turn,
    solve_steady_gliding_turn,
    solve_steady_level_turn,
)

# ---------------------------------------------------------------------------
# Small import and logging helpers
# ---------------------------------------------------------------------------

def _import_analysis_vspaero():
    # Prefer the revised local module when this file is used as a generated pair.
    # Fall back to the original module so the file remains usable after manual
    # renaming or in an existing package layout.
    try:
        return importlib.import_module(".AnalysisVSPAERO_revised", package=__package__)
    except Exception:
        try:
            return importlib.import_module("AnalysisVSPAERO_revised")
        except Exception:
            try:
                return importlib.import_module(".AnalysisVSPAERO", package=__package__)
            except Exception:
                return importlib.import_module("AnalysisVSPAERO")

def _vprint(verbose: int | bool, level: int, message: str) -> None:
    if verbose and int(verbose) >= level:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {message}", flush=True)

# ---------------------------------------------------------------------------
# Geometry summaries used by both VTail scaling and wing deflection
# ---------------------------------------------------------------------------

def _read_wing_section_table_from_loaded(vsp, geom_id: str, *, surf_index: int = 0) -> pd.DataFrame:
    xsec_surf_id = vsp.GetXSecSurf(geom_id, surf_index)

    rows = []
    root_le_x_local = 0.0
    root_y = 0.0
    # OpenVSP wing section parameters live on XSec 1..N-1.
    for section_index in range(1, int(vsp.GetNumXSec(xsec_surf_id))):
        xsec_id = vsp.GetXSec(xsec_surf_id, section_index)
        span = get_xsec_value(vsp, xsec_id, SPAN_NAMES)
        root_chord = get_xsec_value(vsp, xsec_id, ROOT_CHORD_NAMES)
        tip_chord = get_xsec_value(vsp, xsec_id, TIP_CHORD_NAMES)
        sweep_deg = get_xsec_value(vsp, xsec_id, SWEEP_NAMES, default=0.0)
        sweep_location = get_xsec_value(vsp, xsec_id, SWEEP_LOCATION_NAMES, default=0.25)
        if sweep_location > 1.0:
            sweep_location = sweep_location / 100.0
        dihedral_deg = get_xsec_value(vsp, xsec_id, DIHEDRAL_NAMES, default=0.0)
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
    """Read editable wing-section geometry from an OpenVSP WING Geom."""

    vsp = import_openvsp()
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(vsp3_path))
    vsp.Update()

    geom_id = find_one_geom(vsp, geom_name)
    return _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)

def wing_planform_summary(section_table: pd.DataFrame) -> dict[str, float]:
    if section_table.empty:
        raise ValueError("Wing section table is empty.")
    area = float(section_table["area"].sum())
    if not math.isfinite(area) or area <= 0.0:
        raise ValueError(f"Wing planform area must be positive. area={area}")
    span = float(section_table["span"].sum())
    quarter_chord_x_local = float(
        (section_table["area"] * section_table["quarter_chord_x_local"]).sum() / area
    )
    aspect_ratio = span * span / area
    return {
        "area": area,
        "span": span,
        "aspect_ratio": aspect_ratio,
        "quarter_chord_x_local": quarter_chord_x_local,
    }

def _reference_wing_summary_from_loaded(
    vsp,
    *,
    wing_name: str = "WingGeom",
    surf_index: int = 0,
) -> dict[str, Any]:
    """Return the reference wing area and span used for Vv scaling.

    The preferred source is OpenVSP's VSPAERO reference wing ID.  If it is not
    available in the loaded model, this function falls back to a named WING
    geometry.  Area and span are read from the WING Geom total parameters when
    available; the editable section table is used only as a fallback.
    """

    geom_id = ""
    source = ""
    try:
        geom_id = vsp.GetVSPAERORefWingID()
    except Exception:
        geom_id = ""

    if geom_id and str(geom_id).upper() != "NONE":
        try:
            vsp.GetXSecSurf(geom_id, surf_index)
        except Exception:
            geom_id = ""
    else:
        geom_id = ""

    if geom_id:
        source = "vspaero_reference_wing"
    else:
        geom_id = find_one_geom(vsp, wing_name)
        source = "wing_name"

    sections = _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)
    section_summary = wing_planform_summary(sections)

    area = get_geom_parm_value(
        vsp,
        geom_id,
        ("TotalArea", "Total_Area", "Area"),
        ("WingGeom",),
        default=section_summary["area"],
    )
    span = get_geom_parm_value(
        vsp,
        geom_id,
        ("TotalSpan", "Total_Span", "Span"),
        ("WingGeom",),
        default=section_summary["span"],
    )
    if not math.isfinite(area) or area <= 0.0:
        raise ValueError(f"Reference wing area must be positive. area={area}")
    if not math.isfinite(span) or span <= 0.0:
        raise ValueError(f"Reference wing span must be positive. span={span}")

    result = dict(section_summary)
    result.update(
        {
            "geom_id": geom_id,
            "source": source,
            "area": float(area),
            "span": float(span),
            "section_area": section_summary["area"],
            "section_span": section_summary["span"],
        }
    )
    return result

def read_reference_wing_summary(
    vsp3_path: str | os.PathLike,
    *,
    wing_name: str = "WingGeom",
    surf_index: int = 0,
) -> dict[str, Any]:
    """Read the reference wing area and span from a .vsp3 model."""

    vsp = import_openvsp()
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(vsp3_path))
    vsp.Update()
    return _reference_wing_summary_from_loaded(vsp, wing_name=wing_name, surf_index=surf_index)

def calculate_gamma_eff_from_tip_deflection(
    tip_deflection: float,
    semispan: float,
    *,
    n_span: int = 1001,
) -> dict[str, Any]:
    """Calculate Gamma_eff from an analytic elliptical-planform deflection model.

    This definition deliberately does not read a .vsp3 file.  The chart axis is
    tied directly to the design input tip_deflection and a chosen semispan.

    Assumptions:
    - chord distribution is elliptical: c(x) = c0 * sqrt(1 - x**2)
    - rigid dihedral is zero
    - local lift-curve slope is constant
    - elastic dihedral is the elliptical-load, constant-EI beam slope returned
      by elastic_wing_deflection_distribution()

    With x = y / semispan, c0 cancels from the yc(y)-weighted average:

        Gamma_eff = integral x sqrt(1-x^2) theta_b(x) dx
                  / integral x sqrt(1-x^2) dx
    """

    semispan = float(semispan)
    if semispan <= 0.0:
        raise ValueError("semispan must be positive for Gamma_eff calculation.")
    if int(n_span) < 11:
        raise ValueError("n_span must be at least 11 for Gamma_eff integration.")

    deflection = elastic_wing_deflection_distribution(
        float(tip_deflection),
        semispan,
        n_span=int(n_span),
    )
    x = deflection["eta"].to_numpy(dtype=float)
    theta_rad = deflection["theta_rad"].to_numpy(dtype=float)
    weight = x * np.sqrt(np.maximum(1.0 - x * x, 0.0))
    weight_sum = float(np.trapezoid(weight, x))
    if weight_sum <= 0.0:
        raise ValueError("Could not calculate Gamma_eff because the elliptical-planform weight is zero.")

    gamma_elastic_rad = float(np.trapezoid(weight * theta_rad, x) / weight_sum)
    gamma_rigid_rad = 0.0
    gamma_eff_rad = gamma_elastic_rad
    return {
        "Gamma_eff_rad": gamma_eff_rad,
        "Gamma_eff_deg": math.degrees(gamma_eff_rad),
        "Gamma_eff_rigid_rad": gamma_rigid_rad,
        "Gamma_eff_rigid_deg": 0.0,
        "Gamma_eff_elastic_rad": gamma_elastic_rad,
        "Gamma_eff_elastic_deg": math.degrees(gamma_elastic_rad),
        "Gamma_eff_weight_sum": weight_sum,
        "Gamma_eff_semispan": semispan,
        "Gamma_eff_n_span": int(n_span),
        "Gamma_eff_source": "analytic_elliptic_planform_zero_rigid_dihedral_constant_ei_tip_deflection_yc_weighted",
    }

# ---------------------------------------------------------------------------
# VTail Vv update
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

def update_vtail_volume(
    input_vsp3_path: str | os.PathLike,
    output_vsp3_path: str | os.PathLike,
    vv: float,
    *,
    xcg: float,
    wing_name: str = "WingGeom",
    vtail_name: str = "VTailGeom",
    vtail_area_scale_mode: str = "fixed_aspect_ratio",
    surf_index: int = 0,
    verbose: int | bool = 0,
) -> dict[str, Any]:
    """Scale VTailGeom to a target vertical-tail volume.

    The reference wing area/span are read from the loaded .vsp3 model.  The
    vertical-tail moment arm is not supplied by the caller; it is calculated
    from the user-specified Xcg and the original VTailGeom area-weighted
    quarter-chord position.  After scaling, the VTailGeom is shifted so this
    original approximate aerodynamic-center position is preserved.

    vtail_area_scale_mode selects how the requested area change is applied:
    "fixed_aspect_ratio" scales span and chords by the same factor, while
    "fixed_span" keeps section spans unchanged and scales root/tip chords.
    """

    if vv <= 0.0:
        raise ValueError("vv must be positive.")
    if vtail_area_scale_mode not in {"fixed_aspect_ratio", "fixed_span"}:
        raise ValueError("vtail_area_scale_mode must be 'fixed_aspect_ratio' or 'fixed_span'.")

    vsp = import_openvsp()
    input_vsp3_path = Path(input_vsp3_path)
    output_vsp3_path = Path(output_vsp3_path)
    output_vsp3_path.parent.mkdir(parents=True, exist_ok=True)

    _vprint(verbose, 2, f"VTail update start: Vv={vv:.6g}, Xcg={float(xcg):.6g}, output={output_vsp3_path}")

    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(input_vsp3_path))
    vsp.Update()

    wing_summary = _reference_wing_summary_from_loaded(vsp, wing_name=wing_name, surf_index=surf_index)
    wing_area = float(wing_summary["area"])
    wing_span = float(wing_summary["span"])

    geom_id = find_one_geom(vsp, vtail_name)
    xsec_surf_id = vsp.GetXSecSurf(geom_id, surf_index)

    before_sections = _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)
    before_summary = wing_planform_summary(before_sections)

    x_location_original = get_geom_parm_value(
        vsp, geom_id, X_LOCATION_NAMES, XFORM_GROUP_NAMES
    )
    x_ac_original = x_location_original + before_summary["quarter_chord_x_local"]
    lv_original = x_ac_original - float(xcg)
    if not math.isfinite(lv_original) or lv_original <= 0.0:
        raise ValueError(
            "The original VTailGeom approximate aerodynamic center must be aft of Xcg. "
            f"x_ac_original={x_ac_original}, xcg={float(xcg)}, lv_original={lv_original}"
        )

    target_area = float(vv) * wing_area * wing_span / lv_original
    area_ratio = target_area / before_summary["area"]
    if vtail_area_scale_mode == "fixed_aspect_ratio":
        span_scale = math.sqrt(area_ratio)
        chord_scale = span_scale
    else:
        span_scale = 1.0
        chord_scale = area_ratio

    for section_index in range(1, int(vsp.GetNumXSec(xsec_surf_id))):
        xsec_id = vsp.GetXSec(xsec_surf_id, section_index)
        set_wing_section_driver_for_scaling(vsp, geom_id, section_index)
        set_xsec_value(vsp, xsec_id, SPAN_NAMES, get_xsec_value(vsp, xsec_id, SPAN_NAMES) * span_scale)
        set_xsec_value(vsp, xsec_id, ROOT_CHORD_NAMES, get_xsec_value(vsp, xsec_id, ROOT_CHORD_NAMES) * chord_scale)
        set_xsec_value(vsp, xsec_id, TIP_CHORD_NAMES, get_xsec_value(vsp, xsec_id, TIP_CHORD_NAMES) * chord_scale)

    vsp.Update()
    after_sections = _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)
    after_summary = wing_planform_summary(after_sections)

    x_location_before_shift = get_geom_parm_value(
        vsp, geom_id, X_LOCATION_NAMES, XFORM_GROUP_NAMES
    )
    x_ac_before_shift = x_location_before_shift + after_summary["quarter_chord_x_local"]
    delta_x = x_ac_original - x_ac_before_shift
    x_location_after = x_location_before_shift + delta_x
    set_geom_parm_value(
        vsp, geom_id, X_LOCATION_NAMES, XFORM_GROUP_NAMES, x_location_after
    )
    vsp.Update()

    x_ac_after = x_location_after + after_summary["quarter_chord_x_local"]
    lv_actual = x_ac_after - float(xcg)
    vv_actual = after_summary["area"] * lv_actual / (wing_area * wing_span)

    vsp.WriteVSPFile(os.fspath(output_vsp3_path))
    _vprint(
        verbose,
        2,
        "VTail update done: "
        f"area_before={before_summary['area']:.6g}, area_target={target_area:.6g}, "
        f"area_after={after_summary['area']:.6g}, mode={vtail_area_scale_mode}, "
        f"span_scale={span_scale:.6g}, chord_scale={chord_scale:.6g}, "
        f"x_ac_original={x_ac_original:.6g}, x_shift={delta_x:.6g}, vv_actual={vv_actual:.6g}",
    )

    return {
        "input_vsp3_path": os.fspath(input_vsp3_path),
        "output_vsp3_path": os.fspath(output_vsp3_path),
        "vtail_name": vtail_name,
        "wing_reference_source": wing_summary["source"],
        "vv_target": float(vv),
        "vv_actual": float(vv_actual),
        "xcg": float(xcg),
        "lv_original": float(lv_original),
        "lv_actual": float(lv_actual),
        "wing_area": wing_area,
        "wing_span": wing_span,
        "wing_section_area": float(wing_summary["section_area"]),
        "wing_section_span": float(wing_summary["section_span"]),
        "vtail_area_scale_mode": vtail_area_scale_mode,
        "area_ratio": float(area_ratio),
        "span_scale": float(span_scale),
        "chord_scale": float(chord_scale),
        "scale": float(span_scale) if vtail_area_scale_mode == "fixed_aspect_ratio" else float(chord_scale),
        "area_before": before_summary["area"],
        "area_target": target_area,
        "area_after": after_summary["area"],
        "span_before": before_summary["span"],
        "span_after": after_summary["span"],
        "aspect_ratio_before": before_summary["aspect_ratio"],
        "aspect_ratio_after": after_summary["aspect_ratio"],
        "x_ac_original": x_ac_original,
        "x_ac_target": x_ac_original,
        "x_ac_before_shift": x_ac_before_shift,
        "x_ac_after": x_ac_after,
        "x_shift": delta_x,
        "x_location_original": x_location_original,
        "x_location_before_shift": x_location_before_shift,
        "x_location_after": x_location_after,
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
    """Return elliptical-load, constant-EI cantilever deflection distribution."""

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
    verbose: int | bool = 0,
) -> dict[str, Any]:
    """Approximate elastic deflection by adding local theta_b to Dihedral."""

    required_columns = {"y", "theta_deg"}
    missing = required_columns - set(deflection_distribution.columns)
    if missing:
        raise ValueError(f"deflection_distribution is missing columns: {sorted(missing)}")

    y_values = deflection_distribution["y"].to_numpy(dtype=float)
    theta_values = deflection_distribution["theta_deg"].to_numpy(dtype=float)
    if len(y_values) < 2:
        raise ValueError("deflection_distribution must contain at least two span stations.")
    if np.any(~np.isfinite(y_values)) or np.any(~np.isfinite(theta_values)):
        raise ValueError("deflection_distribution contains non-finite y or theta_deg values.")
    if np.any(np.diff(y_values) <= 0.0):
        raise ValueError("deflection_distribution['y'] must be strictly increasing.")

    vsp = import_openvsp()
    input_vsp3_path = Path(input_vsp3_path)
    output_vsp3_path = Path(output_vsp3_path)
    output_vsp3_path.parent.mkdir(parents=True, exist_ok=True)

    _vprint(verbose, 2, f"Wing deflection update start: input={input_vsp3_path}, output={output_vsp3_path}")

    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(input_vsp3_path))
    vsp.Update()

    geom_id = find_one_geom(vsp, wing_name)
    xsec_surf_id = vsp.GetXSecSurf(geom_id, surf_index)
    original_sections = _read_wing_section_table_from_loaded(vsp, geom_id, surf_index=surf_index)

    updates = []
    for row in original_sections.itertuples(index=False):
        section_index = int(row.section_index)
        xsec_id = vsp.GetXSec(xsec_surf_id, section_index)
        elastic_dihedral_deg = float(np.interp(float(row.center_y), y_values, theta_values))
        new_dihedral_deg = float(row.dihedral_deg) + elastic_dihedral_deg
        set_xsec_value(vsp, xsec_id, DIHEDRAL_NAMES, new_dihedral_deg)
        update = {
            "section_index": section_index,
            "center_y": float(row.center_y),
            "rigid_dihedral_deg": float(row.dihedral_deg),
            "elastic_dihedral_deg": elastic_dihedral_deg,
            "new_dihedral_deg": new_dihedral_deg,
        }
        updates.append(update)
        _vprint(
            verbose,
            3,
            "Wing section update: "
            f"section={section_index}, y={update['center_y']:.6g}, "
            f"rigid={update['rigid_dihedral_deg']:.6g} deg, "
            f"elastic={elastic_dihedral_deg:.6g} deg, new={new_dihedral_deg:.6g} deg",
        )

    vsp.Update()
    vsp.WriteVSPFile(os.fspath(output_vsp3_path))
    _vprint(verbose, 2, f"Wing deflection update done: sections={len(updates)}")

    return {
        "input_vsp3_path": os.fspath(input_vsp3_path),
        "output_vsp3_path": os.fspath(output_vsp3_path),
        "wing_name": wing_name,
        "updates": updates,
    }

# ---------------------------------------------------------------------------
# .stab based indices and post-processing
# ---------------------------------------------------------------------------

def _import_roll_rudder_gain():
    try:
        return importlib.import_module(".RollRudderGain", package=__package__)
    except Exception:
        return importlib.import_module("RollRudderGain")

def calculate_stab_basic_indices(
    stab_path: str | os.PathLike,
    *,
    vv: float | None = None,
    tip_deflection: float | None = None,
    control_map: Mapping[str, str] | None = None,
    verbose: int | bool = 0,
) -> dict[str, Any]:
    """Read one .stab file and return direct stability indices."""

    _vprint(verbose, 2, f"Basic .stab indices start: {stab_path}")
    stab = read_vspaero_stab(stab_path)
    control_columns = resolve_control_columns_from_stab(stab, control_map)
    rudder_column = control_columns.get("delta_r")
    aileron_column = control_columns.get("delta_a")
    elevator_column = control_columns.get("delta_e")

    cl_alpha = stab_coefficient_value(stab, "CL", "Alpha")
    cl_beta = stab_coefficient_value(stab, "CMl", "Beta")
    cn_beta = stab_coefficient_value(stab, "CMn", "Beta")
    cl_p = stab_coefficient_value(stab, "CMl", "p")
    cl_r = stab_coefficient_value(stab, "CMl", "r")
    cn_p = stab_coefficient_value(stab, "CMn", "p")
    cn_r = stab_coefficient_value(stab, "CMn", "r")
    cs_beta = stab_coefficient_value(stab, "CS", "Beta")
    cs_p = stab_coefficient_value(stab, "CS", "p")
    cs_r = stab_coefficient_value(stab, "CS", "r")
    cy_beta = solver_axis_derivative_value(stab, "CY", "Beta")
    cy_p = solver_axis_derivative_value(stab, "CY", "p")
    cy_r = solver_axis_derivative_value(stab, "CY", "r")
    cl_delta_r = stab_coefficient_value(stab, "CMl", rudder_column) if rudder_column else math.nan
    cn_delta_r = stab_coefficient_value(stab, "CMn", rudder_column) if rudder_column else math.nan
    cs_delta_r = stab_coefficient_value(stab, "CS", rudder_column) if rudder_column else math.nan
    cy_delta_r = solver_axis_derivative_value(stab, "CY", rudder_column) if rudder_column else math.nan

    spiral_margin = cl_beta * cn_r - cn_beta * cl_r
    simple_rudder_roll_index = (cl_r - cl_beta) * cn_delta_r
    simple_turn_trim_denominator = cn_delta_r * cl_r
    simple_turn_trim_delta_r_per_beta = math.nan
    if (
        math.isfinite(spiral_margin)
        and math.isfinite(simple_turn_trim_denominator)
        and abs(simple_turn_trim_denominator) >= 1.0e-14
    ):
        simple_turn_trim_delta_r_per_beta = spiral_margin / simple_turn_trim_denominator

    indices: dict[str, Any] = {
        "stab_path": os.fspath(stab_path),
        "Vv": np.nan if vv is None else float(vv),
        "tip_deflection": np.nan if tip_deflection is None else float(tip_deflection),
        "CL_alpha": cl_alpha,
        "Cl_beta": cl_beta,
        "Cn_beta": cn_beta,
        "CS_beta": cs_beta,
        "CS_p": cs_p,
        "CS_r": cs_r,
        "CY_beta": cy_beta,
        "CY_p": cy_p,
        "CY_r": cy_r,
        "Cl_p": cl_p,
        "Cl_r": cl_r,
        "Cn_p": cn_p,
        "Cn_r": cn_r,
        "Cl_delta_r": cl_delta_r,
        "Cn_delta_r": cn_delta_r,
        "CS_delta_r": cs_delta_r,
        "CY_delta_r": cy_delta_r,
        "spiral_margin": spiral_margin,
        "simple_rudder_roll_index": simple_rudder_roll_index,
        "simple_turn_trim_delta_r_per_beta": simple_turn_trim_delta_r_per_beta,
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
    }

    # Export the full .stab derivative table in the same column convention used
    # by the wake-iteration sensitivity notebook: CL_Base, CL_Alpha, CMl_Beta,
    # CMm_ConGrp_2, etc.  The compact aliases above are kept only for the
    # existing Vv-Gamma indices and labels.
    for coef_name in stab.derivatives.index:
        row = stab.derivatives.loc[coef_name]
        for derivative_name, value in row.items():
            column_name = f"{coef_name}_Base" if derivative_name == "Base" else f"{coef_name}_{derivative_name}"
            indices[column_name] = float(value)

    # Export solver-axis CY from the same CL/CD/CS conversion used by the trim
    # and 6DOF modules.  This keeps the chart table and simulations on one
    # coefficient convention.
    base_solver_aero = evaluate_stab_linear_aero(
        stab, {}, alpha=stab.alpha0, beta=stab.beta0
    )
    indices["CY_Base"] = base_solver_aero["CY"]
    for derivative_name in stab.derivatives.columns:
        if derivative_name != "Base":
            indices[f"CY_{derivative_name}"] = solver_axis_derivative_value(
                stab, "CY", derivative_name
            )

    return indices

def _flatten_turn_trim(trim: Mapping[str, Any]) -> dict[str, Any]:
    passed = bool(trim.get("passed", False))
    result: dict[str, Any] = {
        "turn_trim_passed": passed,
        "turn_trim_message": trim.get("message", ""),
        "turn_trim_max_abs_residual": trim.get("max_abs_residual", np.nan),
        "turn_trim_cost": trim.get("cost", np.nan),
        "turn_trim_nfev": trim.get("nfev", np.nan),
    }

    # Keep raw solver outputs for diagnosis even when the trim did not pass.
    # Direct turn_trim_* columns are populated only for accepted solutions.
    for source_name in ("solution", "derived", "coefficients", "residuals"):
        for key, value in (trim.get(source_name, {}) or {}).items():
            if isinstance(value, (int, float, np.floating, np.integer, bool, str)) or value is None:
                result[f"turn_trim_{source_name}_{key}"] = value
    if "height_residual" in trim:
        result["turn_trim_height_residual"] = trim["height_residual"]

    solution = trim.get("solution", {}) or {}
    for parameter_name in ("V", "alpha", "beta", "phi", "theta", "Omega", "delta_e", "delta_a", "delta_r", "T"):
        result[f"turn_trim_{parameter_name}"] = solution.get(parameter_name, np.nan) if passed else np.nan

    result["turn_trim_delta_r_per_beta"] = np.nan
    if passed:
        beta_value = float(solution.get("beta", np.nan))
        delta_r_value = float(solution.get("delta_r", np.nan))
        if math.isfinite(beta_value) and math.isfinite(delta_r_value) and abs(beta_value) >= 1.0e-14:
            result["turn_trim_delta_r_per_beta"] = delta_r_value / beta_value
    return result

def calculate_turn_trim_indices_from_stab(
    stab_path: str | os.PathLike,
    *,
    mode: str = "gliding",
    mass: float,
    rho: float | None = None,
    phi: float = math.radians(5.0),
    delta_a: float = 0.0,
    initial_guess: Mapping[str, float] | None = None,
    bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    residual_tol: float = 1.0e-6,
    control_map: Mapping[str, str] | None = None,
    g: float = 9.80665,
    verbose: int | bool = 0,
) -> dict[str, Any]:
    """Solve V=.stab Vinf, fixed phi, and fixed delta_a turn trim."""

    if mode not in {"none", "level", "gliding"}:
        raise ValueError("mode must be 'none', 'level', or 'gliding'.")

    result: dict[str, Any] = {
        "turn_trim_mode": mode,
        "turn_trim_fixed_phi": float(phi),
        "turn_trim_fixed_phi_deg": math.degrees(float(phi)),
        "turn_trim_fixed_delta_a": float(delta_a),
        "turn_trim_passed": np.nan,
    }
    if mode == "none":
        return result

    stab = read_vspaero_stab(stab_path)
    fixed = {"V": float(stab.V0), "phi": float(phi), "delta_a": float(delta_a)}
    result["turn_trim_fixed_V"] = fixed["V"]
    trim_verbose = 1 if verbose and int(verbose) >= 3 else 0

    if mode == "level":
        trim = solve_steady_level_turn(
            fixed, stab_path, mass, rho=rho, initial_guess=initial_guess, g=g,
            bounds=bounds, control_map=control_map, residual_tol=residual_tol,
            verbose=trim_verbose,
        )
    else:
        trim = solve_steady_gliding_turn(
            fixed, stab_path, mass, rho=rho, initial_guess=initial_guess, g=g,
            bounds=bounds, control_map=control_map, residual_tol=residual_tol,
            verbose=trim_verbose,
        )
    result.update(_flatten_turn_trim(trim))
    return result

def _flatten_rudder_limit_turn(trim: Mapping[str, Any]) -> dict[str, Any]:
    """Flatten the selected rudder-limit trim and both endpoint summaries."""

    prefix = "rudder_limit_turn"
    passed = bool(trim.get("passed", False))
    result: dict[str, Any] = {
        f"{prefix}_passed": passed,
        f"{prefix}_complete": bool(trim.get("rudder_limit_complete", False)),
        f"{prefix}_message": trim.get("message", ""),
        f"{prefix}_selected_side": trim.get("selected_side", ""),
        f"{prefix}_delta_r_max": trim.get("delta_r_max", np.nan),
        f"{prefix}_limiting_delta_r": trim.get("limiting_delta_r", np.nan),
        f"{prefix}_max_abs_phi": trim.get("max_abs_phi", np.nan),
        f"{prefix}_fixed_V": trim.get("fixed_V", np.nan),
        f"{prefix}_fixed_delta_a": trim.get("fixed_delta_a", np.nan),
        f"{prefix}_max_abs_residual": trim.get("max_abs_residual", np.nan),
        f"{prefix}_cost": trim.get("cost", np.nan),
        f"{prefix}_nfev": trim.get("nfev", np.nan),
    }

    # Keep raw selected-solution outputs for diagnosis.  Direct columns below
    # are populated only when at least one endpoint produced an accepted trim.
    for source_name in ("solution", "derived", "coefficients", "residuals"):
        for key, value in (trim.get(source_name, {}) or {}).items():
            if isinstance(value, (int, float, np.floating, np.integer, bool, str)) or value is None:
                result[f"{prefix}_{source_name}_{key}"] = value
    if "height_residual" in trim:
        result[f"{prefix}_height_residual"] = trim["height_residual"]

    solution = trim.get("solution", {}) or {}
    for parameter_name in ("V", "alpha", "beta", "phi", "theta", "Omega", "delta_e", "delta_a", "delta_r", "T"):
        result[f"{prefix}_{parameter_name}"] = solution.get(parameter_name, np.nan) if passed else np.nan

    result[f"{prefix}_delta_r_per_beta"] = np.nan
    if passed:
        beta_value = float(solution.get("beta", np.nan))
        delta_r_value = float(solution.get("delta_r", np.nan))
        if math.isfinite(beta_value) and math.isfinite(delta_r_value) and abs(beta_value) >= 1.0e-14:
            result[f"{prefix}_delta_r_per_beta"] = delta_r_value / beta_value

    for side in ("negative", "positive"):
        endpoint = trim.get(f"{side}_trim", {}) or {}
        endpoint_solution = endpoint.get("solution", {}) or {}
        result[f"{prefix}_{side}_passed"] = bool(endpoint.get("passed", False))
        result[f"{prefix}_{side}_message"] = endpoint.get("message", "")
        result[f"{prefix}_{side}_max_abs_residual"] = endpoint.get("max_abs_residual", np.nan)
        for parameter_name in ("beta", "phi", "Omega", "delta_r"):
            result[f"{prefix}_{side}_{parameter_name}"] = endpoint_solution.get(parameter_name, np.nan)

    return result

def calculate_rudder_limit_turn_indices_from_stab(
    stab_path: str | os.PathLike,
    *,
    mode: str = "gliding",
    mass: float,
    delta_r_max: float,
    rho: float | None = None,
    delta_a: float = 0.0,
    initial_guess: Mapping[str, float] | None = None,
    bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    residual_tol: float = 1.0e-6,
    control_map: Mapping[str, str] | None = None,
    g: float = 9.80665,
    verbose: int | bool = 0,
) -> dict[str, Any]:
    """Return maximum-bank indices from the two rudder-limit endpoint trims."""

    if mode not in {"none", "level", "gliding"}:
        raise ValueError("mode must be 'none', 'level', or 'gliding'.")

    result: dict[str, Any] = {
        "rudder_limit_turn_mode": mode,
        "rudder_limit_turn_delta_r_max": float(delta_r_max),
        "rudder_limit_turn_fixed_delta_a": float(delta_a),
        "rudder_limit_turn_passed": np.nan,
    }
    if mode == "none":
        return result

    trim = solve_rudder_limit_turn(
        stab_path,
        mass,
        delta_r_max,
        mode=mode,
        delta_a=delta_a,
        rho=rho,
        initial_guess=initial_guess,
        g=g,
        bounds=bounds,
        control_map=control_map,
        residual_tol=residual_tol,
        verbose=1 if verbose and int(verbose) >= 3 else 0,
    )
    result.update(_flatten_rudder_limit_turn(trim))
    return result

def _find_case_stab_file(case_dir: Path, fallback_stem: str | None = None) -> Path:
    candidates = sorted(case_dir.glob("*.stab"), key=lambda path: path.stat().st_mtime, reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No .stab file was found in {case_dir}.")
    if fallback_stem:
        for path in candidates:
            if path.stem == fallback_stem:
                return path
    return candidates[0]

def vv_wtip_case_name(vv: float, tip_deflection: float) -> str:
    """Return the case-directory name shared by sweep and progress collection."""

    return f"vv_{float(vv):.5f}_wtip_{float(tip_deflection):.5f}".replace("-", "m").replace(".", "p")

def collect_vv_wtip_sweep_progress(
    vv_values: Sequence[float],
    tip_deflections: Sequence[float],
    output_dir: str | os.PathLike,
    *,
    gamma_semispan: float | None = None,
    include_incomplete_cases: bool = False,
    output_csv_path: str | os.PathLike | None = None,
) -> pd.DataFrame:
    """Build a minimal sweep table from completed case folders.

    This is intended for a second notebook while run_vv_wtip_stability_sweep()
    is still running.  The user already knows the Vv and wtip grids from the
    notebook inputs, so this function does not try to parse design variables
    from file names.  It only reconstructs the expected case names and records
    which case folders already contain a .stab file.
    """

    vv_values = list(vv_values)
    tip_deflections = list(tip_deflections)
    if not vv_values:
        raise ValueError("vv_values must not be empty.")
    if not tip_deflections:
        raise ValueError("tip_deflections must not be empty.")

    output_dir = Path(output_dir)
    total_cases = len(vv_values) * len(tip_deflections)
    rows: list[dict[str, Any]] = []
    case_index = 0

    for vv in vv_values:
        for tip_deflection in tip_deflections:
            case_index += 1
            case_name = vv_wtip_case_name(vv, tip_deflection)
            case_dir = output_dir / case_name
            deformed_vsp3_path = case_dir / f"{case_name}.vsp3"
            deflection_csv_path = case_dir / f"{case_name}_deflection.csv"

            stab_path = ""
            passed = False
            error = ""
            if not case_dir.exists():
                error = "case_dir_not_found"
            else:
                try:
                    stab_path = os.fspath(_find_case_stab_file(case_dir, fallback_stem=case_name))
                    passed = True
                except FileNotFoundError:
                    error = "stab_file_not_found"

            if passed or include_incomplete_cases:
                rows.append(
                    {
                        "case": case_name,
                        "case_index": case_index,
                        "case_count": total_cases,
                        "Vv": float(vv),
                        "tip_deflection": float(tip_deflection),
                        "case_dir": os.fspath(case_dir),
                        "gamma_semispan": np.nan if gamma_semispan is None else float(gamma_semispan),
                        "vsp3_path": os.fspath(deformed_vsp3_path),
                        "deflection_csv_path": os.fspath(deflection_csv_path),
                        "stab_path": stab_path,
                        "passed": passed,
                        "vspaero_passed": passed,
                        "error": error,
                    }
                )

    result = pd.DataFrame(rows)
    if output_csv_path is not None:
        output_csv_path = Path(output_csv_path)
        output_csv_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_csv_path, index=False)
    return result

def run_vspaero_stability_case(
    vsp3_path: str | os.PathLike,
    case_dir: str | os.PathLike,
    *,
    alpha_deg: float,
    mach: float,
    reynolds: float,
    ncpu: int | None = None,
    wake_num_iter: int | None = None,
    fixed_wake_flag: bool | None = None,
    redirect_file: str | None = "",
    stop_before_run: bool = False,
    vspaero_verbose: int | bool = 0,
    verbose: int | bool = 0,
) -> dict[str, Any]:
    """Run one VSPAERO stability analysis in a case directory."""

    analysis_vspaero = _import_analysis_vspaero()
    case_dir = Path(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    vsp3_path = Path(vsp3_path).resolve()
    start = time.perf_counter()
    with workdir(case_dir):
        report = analysis_vspaero.vsp_stability_derivatives(
            os.fspath(vsp3_path), alpha=float(alpha_deg), mach=float(mach), reynolds=float(reynolds),
            ncpu=ncpu, wake_num_iter=wake_num_iter, fixed_wake_flag=fixed_wake_flag,
            redirect_file=redirect_file, stop_before_run=stop_before_run,
            vspaero_verbose=vspaero_verbose, verbose=verbose,
        )
    report = dict(report)
    report["case_dir"] = os.fspath(case_dir)
    report["wall_elapsed_s"] = time.perf_counter() - start
    report["stab_path"] = os.fspath(_find_case_stab_file(case_dir)) if report.get("passed") else ""
    return report

def run_vv_wtip_stability_sweep(
    base_vsp3_path: str | os.PathLike,
    vv_values: Sequence[float],
    tip_deflections: Sequence[float],
    flight_condition: Mapping[str, float],
    geometry_config: Mapping[str, Any],
    output_dir: str | os.PathLike,
    *,
    validate_base_model: bool = True,
    ncpu: int | None = None,
    wake_num_iter: int | None = None,
    fixed_wake_flag: bool | None = None,
    redirect_file: str | None = "",
    stop_before_run: bool = False,
    vspaero_verbose: int | bool = 0,
    output_csv_name: str = "vv_wtip_stability_sweep.csv",
    verbose: int | bool = 0,
) -> pd.DataFrame:
    """Create Vv/wtip models and run VSPAERO stability derivatives only."""

    vv_values = list(vv_values)
    tip_deflections = list(tip_deflections)
    if not vv_values:
        raise ValueError("vv_values must not be empty.")
    if not tip_deflections:
        raise ValueError("tip_deflections must not be empty.")
    for key in ("alpha_deg", "mach", "reynolds"):
        if key not in flight_condition:
            raise KeyError(f"flight_condition is missing required key: {key}")
    if "xcg" not in geometry_config:
        raise KeyError("geometry_config is missing required key: 'xcg'")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    base_vsp3_path = Path(base_vsp3_path).resolve()
    total_cases = len(vv_values) * len(tip_deflections)
    sweep_start = time.perf_counter()

    if validate_base_model:
        analysis_vspaero = _import_analysis_vspaero()
        validation = analysis_vspaero.validate_vsp3_for_stability_derivatives(os.fspath(base_vsp3_path), verbose=verbose)
        if not validation.get("passed", False):
            raise RuntimeError(f"Base .vsp3 validation failed: {validation.get('errors')}")

    reference_wing_summary = read_reference_wing_summary(
        base_vsp3_path,
        wing_name=str(geometry_config.get("wing_name", "WingGeom")),
        surf_index=int(geometry_config.get("surf_index", 0)),
    )
    reference_wing_span = float(reference_wing_summary["span"])
    rows: list[dict[str, Any]] = []
    case_index = 0

    for vv in vv_values:
        for tip_deflection in tip_deflections:
            case_index += 1
            case_start = time.perf_counter()
            case_name = vv_wtip_case_name(vv, tip_deflection)
            case_dir = output_dir / case_name
            case_dir.mkdir(parents=True, exist_ok=True)
            vtail_vsp3_path = case_dir / f"{case_name}_vtail.vsp3"
            deformed_vsp3_path = case_dir / f"{case_name}.vsp3"
            row: dict[str, Any] = {
                "case": case_name,
                "case_index": case_index,
                "case_count": total_cases,
                "Vv": float(vv),
                "tip_deflection": float(tip_deflection),
                "gamma_semispan": reference_wing_span / 2.0,
                "case_dir": os.fspath(case_dir),
                "base_vsp3_path": os.fspath(base_vsp3_path),
                "vsp3_path": os.fspath(deformed_vsp3_path),
                "stab_path": "",
                "passed": False,
                "vspaero_passed": False,
            }
            try:
                vtail_report = update_vtail_volume(
                    base_vsp3_path, vtail_vsp3_path, float(vv), xcg=float(geometry_config["xcg"]),
                    wing_name=str(geometry_config.get("wing_name", "WingGeom")),
                    vtail_name=str(geometry_config.get("vtail_name", "VTailGeom")),
                    vtail_area_scale_mode=str(geometry_config.get("vtail_area_scale_mode", "fixed_aspect_ratio")),
                    surf_index=int(geometry_config.get("surf_index", 0)), verbose=verbose,
                )
                row.update({f"vtail_{key}": value for key, value in vtail_report.items() if isinstance(value, (int, float, str))})
                deflection = elastic_wing_deflection_distribution(
                    float(tip_deflection), reference_wing_span / 2.0, n_span=int(geometry_config.get("n_span", 101))
                )
                deflection_csv_path = case_dir / f"{case_name}_deflection.csv"
                deflection.to_csv(deflection_csv_path, index=False)
                row["deflection_csv_path"] = os.fspath(deflection_csv_path)
                wing_report = apply_wing_deflection_as_dihedral(
                    vtail_vsp3_path, deformed_vsp3_path, deflection,
                    wing_name=str(geometry_config.get("wing_name", "WingGeom")),
                    surf_index=int(geometry_config.get("surf_index", 0)), verbose=verbose,
                )
                row["wing_dihedral_update_count"] = len(wing_report["updates"])
                stability_report = run_vspaero_stability_case(
                    deformed_vsp3_path, case_dir,
                    alpha_deg=float(flight_condition["alpha_deg"]), mach=float(flight_condition["mach"]),
                    reynolds=float(flight_condition["reynolds"]), ncpu=ncpu, wake_num_iter=wake_num_iter,
                    fixed_wake_flag=fixed_wake_flag, redirect_file=redirect_file,
                    stop_before_run=stop_before_run, vspaero_verbose=vspaero_verbose, verbose=verbose,
                )
                row["stab_path"] = stability_report.get("stab_path", "")
                row["vspaero_passed"] = bool(stability_report.get("passed", False))
                row["vspaero_wall_elapsed_s"] = stability_report.get("wall_elapsed_s", math.nan)
                for timing_key, timing_value in (stability_report.get("timing", {}) or {}).items():
                    row[f"vspaero_{timing_key}"] = timing_value
                row["passed"] = bool(stability_report.get("passed", False))
                if not row["passed"]:
                    row["error"] = "VSPAERO stopped before solver execution by request." if stability_report.get("stopped_before_run") else repr(stability_report.get("errors", []))
            except Exception as exc:
                row["error"] = repr(exc)
                _vprint(verbose, 1, f"Case {case_index}/{total_cases} failed: {repr(exc)}")
            row["elapsed_s"] = time.perf_counter() - case_start
            rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(output_dir / output_csv_name, index=False)
    _vprint(verbose, 1, f"Vv-wtip stability sweep done: elapsed={time.perf_counter() - sweep_start:.1f} s")
    return result

def run_vv_gamma_chart(*args, **kwargs) -> pd.DataFrame:
    """Deprecated name. Runs only the separated Vv/wtip stability sweep."""

    removed_keywords = {"turn_trim_mode", "turn_trim_fixed", "rudder_limit_turn_mode", "rudder_limit_turn_delta_r_max", "mass", "rho", "initial_guess", "bounds", "residual_tol", "control_map"}
    used_removed = sorted(key for key in removed_keywords if key in kwargs)
    if used_removed:
        raise TypeError(
            "run_vv_gamma_chart no longer accepts post-processing arguments: "
            f"{used_removed}. Use run_vv_wtip_stability_sweep(...) followed by postprocess_vv_gamma_cases(...)."
        )
    return run_vv_wtip_stability_sweep(*args, **kwargs)

def convert_vv_gamma_postprocess_output_units(result: pd.DataFrame) -> pd.DataFrame:
    """Convert the final Vv-Gamma postprocess table to display units.

    All solvers and response calculations return angles in rad and angular rates
    in rad/s.  The postprocess table is the human-facing CSV/plotting artifact,
    so this final step converts only ordinary angles and ordinary angular rates.

    Deliberately unchanged:
    - stability derivatives such as Cl_beta, Cn_r, CL_Alpha, CMl_Beta
    - reduced rates such as p_hat, q_hat, r_hat, phat, qhat, rhat
    - gains and ratios such as *_per_delta_r, *_per_beta, finite-time indices
    - nondimensional response indices and coefficients
    """

    result = result.copy()

    # Explicit *_deg columns are already in degrees.  Explicit *_rad columns are
    # removed from the final table so the CSV contains only the degree version of
    # the same Gamma_eff quantities.
    result = result.drop(
        columns=["Gamma_eff_rad", "Gamma_eff_rigid_rad", "Gamma_eff_elastic_rad"],
        errors="ignore",
    )

    angle_columns_from_rad = [
        "alpha",
        "beta",
        "phi",
        "theta",
        "psi",
        "delta_e",
        "delta_a",
        "delta_r",
        "sixdof_delta_r",
        "sixdof_target_delta_phi",
        "sixdof_phi0",
        "sixdof_phi_final",
        "sixdof_phi_delta_final",
        "sixdof_delta_e_initial",
        "crosswind_gust_beta_peak_input",
        "crosswind_gust_phi0",
        "crosswind_gust_phi_final",
        "crosswind_gust_phi_delta_final",
        "crosswind_gust_phi_delta_at_gust_end",
        "crosswind_gust_max_phi_delta",
        "crosswind_gust_max_abs_phi_delta",
        "crosswind_gust_peak_beta",
        "crosswind_gust_max_abs_beta",
        "crosswind_gust_delta_e_initial",
        "crosswind_gust_delta_a",
        "crosswind_gust_delta_r",
        "turn_trim_fixed_phi",
        "turn_trim_fixed_delta_a",
        "turn_trim_alpha",
        "turn_trim_beta",
        "turn_trim_phi",
        "turn_trim_theta",
        "turn_trim_delta_e",
        "turn_trim_delta_a",
        "turn_trim_delta_r",
        "turn_trim_solution_alpha",
        "turn_trim_solution_beta",
        "turn_trim_solution_phi",
        "turn_trim_solution_theta",
        "turn_trim_solution_delta_e",
        "turn_trim_solution_delta_a",
        "turn_trim_solution_delta_r",
        "simple_rudder_limit_turn_phi_at_positive_delta_r_max",
        "simple_rudder_limit_turn_max_abs_phi",
        "rudder_limit_turn_delta_r_max",
        "rudder_limit_turn_limiting_delta_r",
        "rudder_limit_turn_max_abs_phi",
        "rudder_limit_turn_fixed_delta_a",
        "rudder_limit_turn_alpha",
        "rudder_limit_turn_beta",
        "rudder_limit_turn_phi",
        "rudder_limit_turn_theta",
        "rudder_limit_turn_delta_e",
        "rudder_limit_turn_delta_a",
        "rudder_limit_turn_delta_r",
        "rudder_limit_turn_solution_alpha",
        "rudder_limit_turn_solution_beta",
        "rudder_limit_turn_solution_phi",
        "rudder_limit_turn_solution_theta",
        "rudder_limit_turn_solution_delta_e",
        "rudder_limit_turn_solution_delta_a",
        "rudder_limit_turn_solution_delta_r",
        "rudder_limit_turn_negative_beta",
        "rudder_limit_turn_negative_phi",
        "rudder_limit_turn_negative_delta_r",
        "rudder_limit_turn_positive_beta",
        "rudder_limit_turn_positive_phi",
        "rudder_limit_turn_positive_delta_r",
    ]
    angular_rate_columns_from_rad_s = [
        "p",
        "q",
        "r",
        "Omega",
        "phi_dot",
        "sixdof_roll_response_phi_rate",
        "sixdof_roll_response_reference_phi_rate_to_t_final",
        "crosswind_gust_peak_p",
        "crosswind_gust_max_abs_p",
        "crosswind_gust_peak_r",
        "crosswind_gust_max_abs_r",
        "turn_trim_Omega",
        "turn_trim_solution_Omega",
        "turn_trim_derived_p",
        "turn_trim_derived_q",
        "turn_trim_derived_r",
        "rudder_limit_turn_Omega",
        "rudder_limit_turn_solution_Omega",
        "rudder_limit_turn_derived_p",
        "rudder_limit_turn_derived_q",
        "rudder_limit_turn_derived_r",
        "rudder_limit_turn_negative_Omega",
        "rudder_limit_turn_positive_Omega",
    ]

    for column in angle_columns_from_rad:
        if column in result.columns:
            result[column] = np.degrees(result[column].astype(float))

    for column in angular_rate_columns_from_rad_s:
        if column in result.columns:
            result[column] = np.degrees(result[column].astype(float))

    return result

def postprocess_vv_gamma_cases(
    results: pd.DataFrame | str | os.PathLike,
    *,
    stab_path_column: str = "stab_path",
    gamma_semispan: float | None = None,
    gamma_semispan_column: str = "gamma_semispan",
    gamma_n_span: int = 1001,
    control_map: Mapping[str, str] | None = None,
    mass: float | None = None,
    inertia: Mapping[str, float] | None = None,
    rho: float | None = None,
    g: float = 9.80665,
    run_6dof: bool = True,
    delta_r: float = math.radians(5.0),
    target_delta_phi: float = math.radians(5.0),
    stop_6dof_at_target_delta_phi: bool = True,
    t_final: float = 10.0,
    delta_a: float = 0.0,
    delta_e: float | None = None,
    trim_elevator: bool = True,
    theta_hold: bool = True,
    theta_ref: float | None = None,
    theta_hold_kp: float = 0.3,
    theta_hold_kq: float = 0.8,
    delta_e_min: float | None = None,
    delta_e_max: float | None = None,
    thrust: float | None = None,
    trim_thrust: bool = True,
    phi0: float = 0.0,
    theta0: float | None = None,
    psi0: float = 0.0,
    max_step: float | None = 0.01,
    rtol: float = 1.0e-8,
    atol: float = 1.0e-10,
    turn_trim_mode: str = "gliding",
    turn_trim_phi: float = math.radians(5.0),
    turn_trim_delta_a: float = 0.0,
    turn_trim_initial_guess: Mapping[str, float] | None = None,
    turn_trim_bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    turn_trim_residual_tol: float = 1.0e-6,
    rudder_limit_turn_mode: str = "none",
    rudder_limit_turn_delta_r_max: float | None = None,
    rudder_limit_turn_delta_a: float = 0.0,
    rudder_limit_turn_initial_guess: Mapping[str, float] | None = None,
    rudder_limit_turn_bounds: tuple[Mapping[str, float], Mapping[str, float]] | None = None,
    rudder_limit_turn_residual_tol: float = 1.0e-6,
    write_6dof_history: bool = False,
    plot_6dof_history: bool = False,
    run_crosswind_gust_6dof: bool = False,
    crosswind_gust_Uds: float | None = None,
    crosswind_gust_H: float | None = None,
    crosswind_gust_t_final: float | None = None,
    crosswind_gust_start_time: float = 0.0,
    crosswind_gust_delta_a: float = 0.0,
    crosswind_gust_delta_e: float | None = None,
    crosswind_gust_delta_r: float = 0.0,
    write_crosswind_gust_history: bool = False,
    plot_crosswind_gust_history: bool = False,
    history_output_dir: str | os.PathLike | None = None,
    output_csv_path: str | os.PathLike | None = None,
    verbose: int | bool = 0,
) -> pd.DataFrame:
    """Add Gamma_eff, .stab indices, response indices, and turn-trim indices.

    Internal calculations stay in rad and rad/s. The returned DataFrame and
    optional CSV are converted at the end: angles are deg, angular rates and
    turn rates are deg/s, while stability derivatives, reduced rates, gains,
    and ratios are left unchanged.

    crosswind_gust_Uds uses the positive-from-right convention.  A positive
    value means a lateral gust blowing from the aircraft right side toward the
    left side.
    """

    if turn_trim_mode not in {"none", "level", "gliding"}:
        raise ValueError("turn_trim_mode must be 'none', 'level', or 'gliding'.")
    if rudder_limit_turn_mode not in {"none", "level", "gliding"}:
        raise ValueError("rudder_limit_turn_mode must be 'none', 'level', or 'gliding'.")
    if rudder_limit_turn_mode != "none":
        if rudder_limit_turn_delta_r_max is None:
            raise ValueError(
                "rudder_limit_turn_delta_r_max is required when rudder-limit turn indices are enabled."
            )
        if not math.isfinite(float(rudder_limit_turn_delta_r_max)) or float(rudder_limit_turn_delta_r_max) <= 0.0:
            raise ValueError("rudder_limit_turn_delta_r_max must be positive and finite.")
    table = results.copy() if isinstance(results, pd.DataFrame) else pd.read_csv(results)
    if stab_path_column not in table.columns:
        raise KeyError(f"results is missing required .stab path column: {stab_path_column!r}")

    gamma_semispan_value = gamma_semispan
    if gamma_semispan_value is None and gamma_semispan_column in table.columns:
        candidates = [
            value
            for value in table[gamma_semispan_column].tolist()
            if not pd.isna(value) and str(value).strip()
        ]
        if candidates:
            gamma_semispan_value = float(candidates[0])
    if gamma_semispan_value is None:
        raise ValueError(
            "gamma_semispan is required unless results contains a non-empty "
            f"{gamma_semispan_column!r} column."
        )
    gamma_semispan_value = float(gamma_semispan_value)
    if gamma_semispan_value <= 0.0:
        raise ValueError("gamma_semispan must be positive.")

    needs_mass = bool(
        run_6dof
        or run_crosswind_gust_6dof
        or turn_trim_mode != "none"
        or rudder_limit_turn_mode != "none"
    )
    if needs_mass and mass is None:
        raise ValueError(
            "mass is required when 6DOF, gust 6DOF, fixed-bank turn, "
            "or rudder-limit turn indices are enabled."
        )
    if inertia is not None:
        if "Izz" not in inertia:
            raise KeyError(
                "inertia is missing required key 'Izz' for crosswind_gust_roll_index."
            )
        if not math.isfinite(float(inertia["Izz"])) or float(inertia["Izz"]) <= 0.0:
            raise ValueError("inertia['Izz'] must be positive and finite.")
    if run_6dof or run_crosswind_gust_6dof:
        if inertia is None:
            raise ValueError("inertia is required when 6DOF or gust 6DOF indices are enabled.")
        for key in ("Ixx", "Iyy", "Izz"):
            if key not in inertia:
                raise KeyError(f"inertia is missing required key: {key}")
    if run_crosswind_gust_6dof:
        if crosswind_gust_Uds is None:
            raise ValueError("crosswind_gust_Uds is required when run_crosswind_gust_6dof=True.")
        if crosswind_gust_H is None:
            raise ValueError("crosswind_gust_H is required when run_crosswind_gust_6dof=True.")

    rudder_gain = _import_roll_rudder_gain() if (run_6dof or run_crosswind_gust_6dof) else None
    Ixx = float(inertia["Ixx"]) if inertia is not None and "Ixx" in inertia else math.nan
    Iyy = float(inertia["Iyy"]) if inertia is not None and "Iyy" in inertia else math.nan
    Izz = float(inertia["Izz"]) if inertia is not None else math.nan
    Ixz = float(inertia.get("Ixz", 0.0)) if inertia is not None else 0.0
    history_output_dir_path = None if history_output_dir is None else Path(history_output_dir)
    if history_output_dir_path is not None:
        history_output_dir_path.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for row_index, row in table.iterrows():
        case_name = str(row.get("case", f"case_{row_index}"))
        output_row = dict(row)
        output_row["postprocess_passed"] = False
        output_row["postprocess_error"] = ""
        output_row["Gamma_eff_rad"] = np.nan
        output_row["Gamma_eff_deg"] = np.nan
        output_row["Gamma_eff_source"] = "analytic_elliptic_planform_zero_rigid_dihedral_constant_ei_tip_deflection_yc_weighted"
        output_row["Gamma_eff_rigid_rad"] = np.nan
        output_row["Gamma_eff_rigid_deg"] = np.nan
        output_row["Gamma_eff_elastic_rad"] = np.nan
        output_row["Gamma_eff_elastic_deg"] = np.nan
        output_row["Gamma_eff_weight_sum"] = np.nan
        output_row["Gamma_eff_semispan"] = gamma_semispan_value
        output_row["Gamma_eff_n_span"] = int(gamma_n_span)
        output_row["simple_rudder_roll_index"] = np.nan
        output_row["simple_turn_trim_delta_r_per_beta"] = np.nan
        output_row["simple_rudder_turn_K_phi"] = np.nan
        output_row["simple_rudder_turn_phi_per_delta_r"] = np.nan
        output_row["simple_rudder_limit_turn_phi_at_positive_delta_r_max"] = np.nan
        output_row["simple_rudder_limit_turn_max_abs_phi"] = np.nan
        output_row["crosswind_gust_roll_index"] = np.nan
        output_row["crosswind_gust_roll_index_abs"] = np.nan
        # Keep the output schema stable even when every completed case fails
        # before the 6DOF or turn-trim step.  This makes downstream plotting
        # fail with "no valid rows" instead of "missing column".
        if run_6dof:
            output_row.update({
                "sixdof_success": False,
                "sixdof_message": "not_evaluated",
                "sixdof_roll_response_reached": False,
                "sixdof_target_delta_phi": float(target_delta_phi),
                "sixdof_target_delta_phi_deg": math.degrees(float(target_delta_phi)),
                "sixdof_stop_at_target_delta_phi": bool(stop_6dof_at_target_delta_phi),
                "sixdof_stopped_at_target_delta_phi": False,
                "sixdof_phi0": float(phi0),
                "sixdof_phi_final": np.nan,
                "sixdof_phi_delta_final": np.nan,
                "sixdof_t_start": np.nan,
                "sixdof_t_final": float(t_final),
                "sixdof_t_reach": np.nan,
                "sixdof_dt_reach": np.nan,
                "sixdof_delta_r": float(delta_r),
                "sixdof_index_Vinf": np.nan,
                "sixdof_index_Bref": np.nan,
                "sixdof_roll_response_phi_rate": np.nan,
                "sixdof_roll_response_phi_rate_per_delta_r": np.nan,
                "sixdof_finite_time_roll_index": np.nan,
                "sixdof_roll_response_error": "not_evaluated",
                "sixdof_roll_response_reference_phi_rate_to_t_final": np.nan,
                "sixdof_roll_response_reference_phi_rate_per_delta_r_to_t_final": np.nan,
                "sixdof_roll_response_index_reference": np.nan,
                "sixdof_roll_response_fraction_of_target": np.nan,
                "sixdof_delta_e_initial": np.nan,
                "sixdof_thrust": np.nan,
            })
        if run_crosswind_gust_6dof:
            output_row.update({
                "crosswind_gust_success": False,
                "crosswind_gust_message": "not_evaluated",
                "crosswind_gust_Uds": np.nan if crosswind_gust_Uds is None else float(crosswind_gust_Uds),
                "crosswind_gust_H": np.nan if crosswind_gust_H is None else float(crosswind_gust_H),
                "crosswind_gust_start_time": float(crosswind_gust_start_time),
                "crosswind_gust_end_time": np.nan,
                "crosswind_gust_duration": np.nan,
                "crosswind_gust_reference_V": np.nan,
                "crosswind_gust_beta_peak_input": np.nan,
                "crosswind_gust_phi0": float(phi0),
                "crosswind_gust_t_start": np.nan,
                "crosswind_gust_t_final": np.nan if crosswind_gust_t_final is None else float(crosswind_gust_t_final),
                "crosswind_gust_phi_final": np.nan,
                "crosswind_gust_phi_delta_final": np.nan,
                "crosswind_gust_phi_delta_at_gust_end": np.nan,
                "crosswind_gust_max_phi_delta": np.nan,
                "crosswind_gust_max_abs_phi_delta": np.nan,
                "crosswind_gust_max_abs_phi_delta_time": np.nan,
                "crosswind_gust_max_abs_beta": np.nan,
                "crosswind_gust_max_abs_beta_time": np.nan,
                "crosswind_gust_max_abs_p": np.nan,
                "crosswind_gust_max_abs_p_time": np.nan,
                "crosswind_gust_max_abs_r": np.nan,
                "crosswind_gust_max_abs_r_time": np.nan,
                "crosswind_gust_max_abs_phat": np.nan,
                "crosswind_gust_max_abs_phat_time": np.nan,
                "crosswind_gust_max_abs_rhat": np.nan,
                "crosswind_gust_max_abs_rhat_time": np.nan,
                "crosswind_gust_bank_index_abs_per_beta_g": np.nan,
                "crosswind_gust_bank_index_final_per_beta_g": np.nan,
                "crosswind_gust_bank_index_at_gust_end_per_beta_g": np.nan,
                "crosswind_gust_phat_index_abs_per_beta_g": np.nan,
                "crosswind_gust_rhat_index_abs_per_beta_g": np.nan,
                "crosswind_gust_delta_e_initial": np.nan,
                "crosswind_gust_delta_a": float(crosswind_gust_delta_a),
                "crosswind_gust_delta_r": float(crosswind_gust_delta_r),
                "crosswind_gust_thrust": np.nan,
            })
        if turn_trim_mode != "none":
            output_row.update({
                "turn_trim_mode": turn_trim_mode,
                "turn_trim_passed": False,
                "turn_trim_V": np.nan,
                "turn_trim_alpha": np.nan,
                "turn_trim_beta": np.nan,
                "turn_trim_phi": np.nan,
                "turn_trim_theta": np.nan,
                "turn_trim_Omega": np.nan,
                "turn_trim_delta_e": np.nan,
                "turn_trim_delta_a": np.nan,
                "turn_trim_delta_r": np.nan,
                "turn_trim_T": np.nan,
                "turn_trim_delta_r_per_beta": np.nan,
            })
        if rudder_limit_turn_mode != "none":
            output_row.update({
                "rudder_limit_turn_mode": rudder_limit_turn_mode,
                "rudder_limit_turn_passed": False,
                "rudder_limit_turn_complete": False,
                "rudder_limit_turn_delta_r_max": float(rudder_limit_turn_delta_r_max),
                "rudder_limit_turn_fixed_delta_a": float(rudder_limit_turn_delta_a),
                "rudder_limit_turn_selected_side": "",
                "rudder_limit_turn_max_abs_phi": np.nan,
                "rudder_limit_turn_limiting_delta_r": np.nan,
                "rudder_limit_turn_V": np.nan,
                "rudder_limit_turn_alpha": np.nan,
                "rudder_limit_turn_beta": np.nan,
                "rudder_limit_turn_phi": np.nan,
                "rudder_limit_turn_theta": np.nan,
                "rudder_limit_turn_Omega": np.nan,
                "rudder_limit_turn_delta_e": np.nan,
                "rudder_limit_turn_delta_a": np.nan,
                "rudder_limit_turn_delta_r": np.nan,
                "rudder_limit_turn_T": np.nan,
                "rudder_limit_turn_delta_r_per_beta": np.nan,
                "rudder_limit_turn_negative_passed": False,
                "rudder_limit_turn_positive_passed": False,
            })

        stab_path_value = row.get(stab_path_column, "")
        if pd.isna(stab_path_value) or str(stab_path_value).strip() == "":
            output_row["postprocess_error"] = "missing stab_path"
            rows.append(output_row)
            continue
        stab_path = Path(str(stab_path_value))
        if not stab_path.exists():
            output_row["postprocess_error"] = f"stab file not found: {stab_path}"
            rows.append(output_row)
            continue

        try:
            gamma_eff = calculate_gamma_eff_from_tip_deflection(
                float(row.get("tip_deflection")),
                gamma_semispan_value,
                n_span=int(gamma_n_span),
            )
            output_row.update(gamma_eff)
            basic = calculate_stab_basic_indices(stab_path, vv=row.get("Vv", None), tip_deflection=row.get("tip_deflection", None), control_map=control_map, verbose=verbose)
            output_row.update(basic)

            # Small-angle steady rudder-only turn approximation from the article:
            #
            #   phi / delta_r = q S / (m g) * K_phi
            #
            #   K_phi = -Cn_delta_r / spiral_margin
            #           * [4m/(rho S b) Cl_beta + CY_beta Cl_r]
            #
            # The signed value uses +delta_r_max.  The absolute value is the
            # simple estimate of the maximum trimmable bank angle.  Near a zero
            # spiral denominator the approximation is intentionally left NaN.
            if rudder_limit_turn_mode != "none":
                rho_used = float(basic["Rho"] if rho is None else rho)
                spiral_margin = float(basic["spiral_margin"])
                if (
                    math.isfinite(rho_used)
                    and rho_used > 0.0
                    and math.isfinite(spiral_margin)
                    and abs(spiral_margin) >= 1.0e-14
                ):
                    simple_K_phi = (
                        -float(basic["Cn_delta_r"])
                        / spiral_margin
                        * (
                            4.0
                            * float(mass)
                            / (rho_used * float(basic["Sref"]) * float(basic["Bref"]))
                            * float(basic["Cl_beta"])
                            + float(basic["CY_beta"]) * float(basic["Cl_r"])
                        )
                    )
                    simple_phi_per_delta_r = (
                        0.5
                        * rho_used
                        * float(basic["Vinf"]) ** 2
                        * float(basic["Sref"])
                        / (float(mass) * float(g))
                        * simple_K_phi
                    )
                    simple_phi_at_positive_limit = (
                        simple_phi_per_delta_r
                        * float(rudder_limit_turn_delta_r_max)
                    )
                    output_row["simple_rudder_turn_K_phi"] = np.abs(simple_K_phi)
                    output_row["simple_rudder_turn_phi_per_delta_r"] = simple_phi_per_delta_r
                    output_row["simple_rudder_limit_turn_phi_at_positive_delta_r_max"] = simple_phi_at_positive_limit
                    output_row["simple_rudder_limit_turn_max_abs_phi"] = abs(simple_phi_at_positive_limit)

            # Inertia-retaining crosswind-gust roll-path index:
            #
            #   K_gust,roll =
            #       8 Izz / (rho S b^3) * Cl_beta
            #       + Cl_rhat * Cn_beta
            #
            # Positive gust input means wind from the aircraft right side.
            # This reduced index assumes Ixz ~= 0 and removes the common positive
            # factor mu_x * mu_z.  It is a static .stab-derived signed path index,
            # so it does not require a time-marching gust simulation.
            if inertia is not None:
                rho_used = float(basic["Rho"] if rho is None else rho)
                Sref = float(basic["Sref"])
                Bref = float(basic["Bref"])
                if (
                    math.isfinite(rho_used)
                    and rho_used > 0.0
                    and math.isfinite(Sref)
                    and Sref > 0.0
                    and math.isfinite(Bref)
                    and Bref > 0.0
                ):
                    output_row["crosswind_gust_roll_index"] = (
                        8.0
                        * Izz
                        / (rho_used * Sref * Bref**3)
                        * float(basic["Cl_beta"])
                        + float(basic["Cl_r"])
                        * float(basic["Cn_beta"])
                    )
                    output_row["crosswind_gust_roll_index_abs"] = abs(
                        output_row["crosswind_gust_roll_index"]
                    )

            if run_6dof:
                history = rudder_gain.simulate_6dof_rudder_step_from_stab(
                    stab_path, mass=float(mass), Ixx=Ixx, Iyy=Iyy, Izz=Izz, Ixz=Ixz,
                    delta_r=float(delta_r), t_final=float(t_final), control_map=control_map,
                    delta_a=float(delta_a), delta_e=delta_e, trim_elevator=trim_elevator,
                    theta_hold=theta_hold, theta_ref=theta_ref, theta_hold_kp=float(theta_hold_kp),
                    theta_hold_kq=float(theta_hold_kq), delta_e_min=delta_e_min, delta_e_max=delta_e_max,
                    thrust=thrust, trim_thrust=trim_thrust, g=float(g), rho=rho,
                    phi0=float(phi0), theta0=theta0, psi0=float(psi0), max_step=max_step,
                    rtol=rtol, atol=atol,
                    stop_at_target_delta_phi=bool(stop_6dof_at_target_delta_phi),
                    target_delta_phi=float(target_delta_phi),
                )
                response = rudder_gain.calculate_roll_response_index_by_delta_phi(
                    history, delta_r=float(delta_r), target_delta_phi=float(target_delta_phi),
                    V=float(basic["Vinf"]), Bref=float(basic["Bref"]), phi0=float(phi0),
                )
                output_row.update(response)
                output_row["sixdof_stop_at_target_delta_phi"] = bool(history.attrs.get("stop_at_target_delta_phi", False))
                output_row["sixdof_stopped_at_target_delta_phi"] = bool(history.attrs.get("target_delta_phi_reached", False))
                output_row["sixdof_success"] = bool(response["sixdof_roll_response_reached"])
                output_row["sixdof_message"] = response.get("sixdof_roll_response_error", "")
                output_row["sixdof_delta_e_initial"] = float(history["delta_e"].iloc[0]) if "delta_e" in history.columns else math.nan
                output_row["sixdof_thrust"] = float(history["thrust"].iloc[0]) if "thrust" in history.columns else math.nan
                if history_output_dir_path is not None:
                    history_csv_path = history_output_dir_path / f"{case_name}_6dof_history.csv"
                    history_plot_path = history_output_dir_path / f"{case_name}_6dof_history.png"
                    if write_6dof_history:
                        output_row["sixdof_history_csv_path"] = str(rudder_gain.write_6dof_history_csv(history, history_csv_path))
                    if plot_6dof_history:
                        rudder_gain.plot_6dof_history(history, plot_path=history_plot_path, show=False, degrees=True)
                        output_row["sixdof_history_plot_path"] = str(history_plot_path)
            if run_crosswind_gust_6dof:
                gust_history = rudder_gain.simulate_6dof_crosswind_gust_from_stab(
                    stab_path, mass=float(mass), Ixx=Ixx, Iyy=Iyy, Izz=Izz, Ixz=Ixz,
                    Uds=float(crosswind_gust_Uds), H=float(crosswind_gust_H),
                    t_final=None if crosswind_gust_t_final is None else float(crosswind_gust_t_final),
                    control_map=control_map,
                    delta_a=float(crosswind_gust_delta_a),
                    delta_e=crosswind_gust_delta_e,
                    delta_r=float(crosswind_gust_delta_r),
                    trim_elevator=trim_elevator,
                    theta_hold=theta_hold, theta_ref=theta_ref,
                    theta_hold_kp=float(theta_hold_kp), theta_hold_kq=float(theta_hold_kq),
                    delta_e_min=delta_e_min, delta_e_max=delta_e_max,
                    thrust=thrust, trim_thrust=trim_thrust, g=float(g), rho=rho,
                    phi0=float(phi0), theta0=theta0, psi0=float(psi0),
                    gust_start_time=float(crosswind_gust_start_time),
                    max_step=max_step, rtol=rtol, atol=atol,
                )
                gust_response = rudder_gain.calculate_crosswind_gust_response_indices(
                    gust_history, phi0=float(phi0), Vref=float(basic["Vinf"]), Bref=float(basic["Bref"])
                )
                output_row.update(gust_response)
                output_row["crosswind_gust_delta_e_initial"] = float(gust_history["delta_e"].iloc[0]) if "delta_e" in gust_history.columns else math.nan
                output_row["crosswind_gust_delta_a"] = float(gust_history["delta_a"].iloc[0]) if "delta_a" in gust_history.columns else float(crosswind_gust_delta_a)
                output_row["crosswind_gust_delta_r"] = float(gust_history["delta_r"].iloc[0]) if "delta_r" in gust_history.columns else float(crosswind_gust_delta_r)
                output_row["crosswind_gust_thrust"] = float(gust_history["thrust"].iloc[0]) if "thrust" in gust_history.columns else math.nan
                if history_output_dir_path is not None:
                    gust_csv_path = history_output_dir_path / f"{case_name}_crosswind_gust_6dof_history.csv"
                    gust_plot_path = history_output_dir_path / f"{case_name}_crosswind_gust_6dof_history.png"
                    if write_crosswind_gust_history:
                        output_row["crosswind_gust_history_csv_path"] = str(rudder_gain.write_6dof_history_csv(gust_history, gust_csv_path))
                    if plot_crosswind_gust_history:
                        rudder_gain.plot_6dof_history(gust_history, plot_path=gust_plot_path, show=False, degrees=True)
                        output_row["crosswind_gust_history_plot_path"] = str(gust_plot_path)
            if turn_trim_mode != "none":
                try:
                    output_row.update(calculate_turn_trim_indices_from_stab(
                        stab_path,
                        mode=turn_trim_mode,
                        mass=float(mass),
                        rho=rho,
                        phi=float(turn_trim_phi),
                        delta_a=float(turn_trim_delta_a),
                        initial_guess=turn_trim_initial_guess,
                        bounds=turn_trim_bounds,
                        residual_tol=float(turn_trim_residual_tol),
                        control_map=control_map,
                        g=float(g),
                        verbose=verbose,
                    ))
                except Exception as exc:
                    output_row.update({
                        "turn_trim_mode": turn_trim_mode,
                        "turn_trim_passed": False,
                        "turn_trim_V": np.nan,
                        "turn_trim_alpha": np.nan,
                        "turn_trim_beta": np.nan,
                        "turn_trim_phi": np.nan,
                        "turn_trim_theta": np.nan,
                        "turn_trim_Omega": np.nan,
                        "turn_trim_delta_e": np.nan,
                        "turn_trim_delta_a": np.nan,
                        "turn_trim_delta_r": np.nan,
                        "turn_trim_T": np.nan,
                        "turn_trim_delta_r_per_beta": np.nan,
                        "turn_trim_error": repr(exc),
                    })
            if rudder_limit_turn_mode != "none":
                try:
                    output_row.update(calculate_rudder_limit_turn_indices_from_stab(
                        stab_path,
                        mode=rudder_limit_turn_mode,
                        mass=float(mass),
                        delta_r_max=float(rudder_limit_turn_delta_r_max),
                        rho=rho,
                        delta_a=float(rudder_limit_turn_delta_a),
                        initial_guess=rudder_limit_turn_initial_guess,
                        bounds=rudder_limit_turn_bounds,
                        residual_tol=float(rudder_limit_turn_residual_tol),
                        control_map=control_map,
                        g=float(g),
                        verbose=verbose,
                    ))
                except Exception as exc:
                    output_row.update({
                        "rudder_limit_turn_mode": rudder_limit_turn_mode,
                        "rudder_limit_turn_passed": False,
                        "rudder_limit_turn_complete": False,
                        "rudder_limit_turn_V": np.nan,
                        "rudder_limit_turn_alpha": np.nan,
                        "rudder_limit_turn_beta": np.nan,
                        "rudder_limit_turn_phi": np.nan,
                        "rudder_limit_turn_theta": np.nan,
                        "rudder_limit_turn_Omega": np.nan,
                        "rudder_limit_turn_delta_e": np.nan,
                        "rudder_limit_turn_delta_a": np.nan,
                        "rudder_limit_turn_delta_r": np.nan,
                        "rudder_limit_turn_T": np.nan,
                        "rudder_limit_turn_delta_r_per_beta": np.nan,
                        "rudder_limit_turn_error": repr(exc),
                    })
            output_row["postprocess_passed"] = True
        except Exception as exc:
            output_row["postprocess_error"] = repr(exc)
            _vprint(verbose, 1, f"Postprocess case failed: {case_name}: {repr(exc)}")
        rows.append(output_row)

    result = pd.DataFrame(rows)
    result = convert_vv_gamma_postprocess_output_units(result)
    if output_csv_path is not None:
        output_csv_path = Path(output_csv_path)
        output_csv_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(output_csv_path, index=False)
    return result

# ---------------------------------------------------------------------------
# Lightweight plotting helpers for notebooks
# ---------------------------------------------------------------------------

VV_GAMMA_LATEX_LABELS = {
    "Gamma_eff_deg": r"$\Gamma_{\mathrm{eff}}\ [\mathrm{deg}]$",
    "CL_alpha": r"$C_{L\alpha}$",
    "Cl_beta": r"$C_{l\beta}$",
    "Cn_beta": r"$C_{n\beta}$",
    "CS_beta": r"$C_{S\beta}$",
    "CS_p": r"$C_{S\hat{p}}$",
    "CS_r": r"$C_{S\hat{r}}$",
    "CY_beta": r"$C_{Y\beta}$",
    "CY_p": r"$C_{Y\hat{p}}$",
    "CY_r": r"$C_{Y\hat{r}}$",
    "Cl_p": r"$C_{l\hat{p}}$",
    "Cl_r": r"$C_{l\hat{r}}$",
    "Cn_p": r"$C_{n\hat{p}}$",
    "Cn_r": r"$C_{n\hat{r}}$",
    "Cl_delta_r": r"$C_{l\delta_r}$",
    "Cn_delta_r": r"$C_{n\delta_r}$",
    "CS_delta_r": r"$C_{S\delta_r}$",
    "CY_delta_r": r"$C_{Y\delta_r}$",
    "spiral_margin": r"$C_{l\beta}C_{n\hat{r}}-C_{n\beta}C_{l\hat{r}}$",
    "simple_turn_trim_delta_r_per_beta": r"$\left(\delta_r/\beta\right)=\frac{C_{l\beta}C_{n\hat{r}}-C_{n\beta}C_{l\hat{r}}}{C_{n\delta_r}C_{l\hat{r}}}$",
    "simple_rudder_turn_K_phi": r"$\left|K_\phi\right|=\left|-\frac{C_{n\delta_r}}{C_{l\beta}C_{n\hat{r}}-C_{n\beta}C_{l\hat{r}}}\left(\frac{4m}{\rho S b}C_{l\beta}+C_{Y\beta}C_{l\hat{r}}\right)\right|$",
    "simple_rudder_turn_phi_per_delta_r": r"$(\phi/\delta_r)_{simple}$",
    "simple_rudder_limit_turn_phi_at_positive_delta_r_max": r"$\phi_{simple}(+\delta_{r,\max})\ [\mathrm{deg}]$",
    "simple_rudder_limit_turn_max_abs_phi": r"$|\phi|_{\max,simple}\ [\mathrm{deg}]$",
    "simple_rudder_roll_index": r"$K_{rudder,roll}=(C_{l\hat{r}}-C_{l\beta})C_{n\delta_r}$",
    "sixdof_finite_time_roll_index": r"$\frac{b}{2V}\frac{\Delta\phi/\Delta t}{\delta_r}$",
    "sixdof_roll_response_index_reference": r"$\frac{b}{2V}\frac{(\phi_f-\phi_0)/t_f}{\delta_r}$",
    "crosswind_gust_roll_index": r"$K_{gust,roll}=\frac{8I_z}{\rho S b^3}C_{l\beta}+C_{l\hat{r}}C_{n\beta}$",
    "crosswind_gust_roll_index_abs": r"$|K_{gust,roll}|$",
    "sixdof_roll_response_phi_rate_per_delta_r": r"$\{(\Delta\phi/\Delta t)/\delta_r\}_{6DOF}$",
    "sixdof_roll_response_fraction_of_target": r"$(\phi_f-\phi_0)/\Delta\phi_{target}$",
    "turn_trim_V": r"$V_{trim}$",
    "turn_trim_alpha": r"$\alpha_{trim}\ [\mathrm{deg}]$",
    "turn_trim_beta": r"$\beta_{trim}\ [\mathrm{deg}]$",
    "turn_trim_phi": r"$\phi_{trim}\ [\mathrm{deg}]$",
    "turn_trim_theta": r"$\theta_{trim}\ [\mathrm{deg}]$",
    "turn_trim_Omega": r"$\Omega_{trim}\ [\mathrm{deg/s}]$",
    "turn_trim_delta_e": r"$\delta_{e,trim}\ [\mathrm{deg}]$",
    "turn_trim_delta_a": r"$\delta_{a,trim}\ [\mathrm{deg}]$",
    "turn_trim_delta_r": r"$\delta_{r,trim}\ [\mathrm{deg}]$",
    "turn_trim_T": r"$T_{trim}$",
    "turn_trim_delta_r_per_beta": r"$(\delta_r/\beta)_{trim}$",
    "turn_trim_solution_beta": r"$\beta\ [\mathrm{deg}]$",
    "turn_trim_solution_delta_r": r"$\delta_r\ [\mathrm{deg}]$",
    "turn_trim_solution_Omega": r"$\Omega\ [\mathrm{deg/s}]$",
    "rudder_limit_turn_V": r"$V_{\delta_r,lim}$",
    "rudder_limit_turn_alpha": r"$\alpha_{\delta_r,lim}\ [\mathrm{deg}]$",
    "rudder_limit_turn_beta": r"$\beta_{\delta_r,lim}\ [\mathrm{deg}]$",
    "rudder_limit_turn_phi": r"$\phi_{\delta_r,lim}\ [\mathrm{deg}]$",
    "rudder_limit_turn_max_abs_phi": r"$|\phi|_{\max,\delta_r}\ [\mathrm{deg}]$",
    "rudder_limit_turn_theta": r"$\theta_{\delta_r,lim}\ [\mathrm{deg}]$",
    "rudder_limit_turn_Omega": r"$\Omega_{\delta_r,lim}\ [\mathrm{deg/s}]$",
    "rudder_limit_turn_delta_e": r"$\delta_{e,\delta_r,lim}\ [\mathrm{deg}]$",
    "rudder_limit_turn_delta_a": r"$\delta_{a,\delta_r,lim}\ [\mathrm{deg}]$",
    "rudder_limit_turn_delta_r": r"$\delta_{r,lim}\ [\mathrm{deg}]$",
    "rudder_limit_turn_delta_r_max": r"$\delta_{r,\max}\ [\mathrm{deg}]$",
    "rudder_limit_turn_limiting_delta_r": r"$\delta_{r,selected}\ [\mathrm{deg}]$",
    "rudder_limit_turn_T": r"$T_{\delta_r,lim}$",
    "rudder_limit_turn_delta_r_per_beta": r"$(\delta_r/\beta)_{\delta_r,lim}$",
    "rudder_limit_turn_solution_beta": r"$\beta\ [\mathrm{deg}]$",
    "rudder_limit_turn_solution_delta_r": r"$\delta_r\ [\mathrm{deg}]$",
    "rudder_limit_turn_solution_Omega": r"$\Omega\ [\mathrm{deg/s}]$",
    "sixdof_delta_r": r"$\delta_r\ [\mathrm{deg}]$",
    "sixdof_target_delta_phi": r"$\Delta\phi_{target}\ [\mathrm{deg}]$",
    "sixdof_roll_response_phi_rate": r"$(\Delta\phi/\Delta t)_{6DOF}\ [\mathrm{deg/s}]$",
    "turn_trim_derived_sink_rate": r"$w_{\mathrm{sink}}$",
    "rudder_limit_turn_derived_sink_rate": r"$w_{\mathrm{sink}}$",
    "rudder_limit_turn_max_abs_residual": r"$\max|r_i|$",
}

def vv_gamma_latex_label(column_name: str) -> str:
    """Return a LaTeX label for a Vv-Gamma result column."""

    if column_name in VV_GAMMA_LATEX_LABELS:
        return VV_GAMMA_LATEX_LABELS[column_name]

    coefficient_labels = {
        "CL": r"C_L",
        "CD": r"C_D",
        "CS": r"C_S",
        "CY": r"C_Y",
        "CMl": r"C_l",
        "CMm": r"C_m",
        "CMn": r"C_n",
        "CFx": r"C_X",
        "CFy": r"C_Y",
        "CFz": r"C_Z",
        "CMx": r"C_x",
        "CMy": r"C_y",
        "CMz": r"C_z",
    }
    variable_labels = {
        "U": r"U",
        "Alpha": r"\alpha",
        "Beta": r"\beta",
        "p": r"\hat{p}",
        "q": r"\hat{q}",
        "r": r"\hat{r}",
        "Mach": r"M",
        "ConGrp_1": r"\delta_{1}",
        "ConGrp_2": r"\delta_{2}",
        "ConGrp_3": r"\delta_{3}",
    }

    for coef_name, coef_label in coefficient_labels.items():
        prefix = f"{coef_name}_"
        if not column_name.startswith(prefix):
            continue
        suffix = column_name[len(prefix):]
        if suffix == "Base":
            return f"${coef_label}$"
        if suffix in variable_labels:
            return rf"$\partial {coef_label}/\partial {variable_labels[suffix]}$"

    return "$" + str(column_name).replace("_", r"\_") + "$"

def build_vv_gamma_response_surface(
    results: pd.DataFrame,
    value_column: str,
    *,
    x_column: str = "Vv",
    y_column: str = "Gamma_eff_deg",
    grid_size: int = 80,
    method: str = "linear",
    fallback_method: str = "nearest",
) -> dict[str, Any]:
    """Build a scattered-interpolation response surface for a Vv-Gamma chart.

    The returned surface is a plotting-oriented response surface, not a
    validated surrogate model.  It is built from completed Vv-Gamma result
    points and can be passed directly to contour plotting.
    """

    required_columns = {x_column, y_column, value_column}
    missing_columns = required_columns - set(results.columns)
    if missing_columns:
        available = ", ".join(map(str, results.columns[:30]))
        raise KeyError(
            f"results is missing required column(s): {sorted(missing_columns)}. "
            "Pass the postprocess_vv_gamma_cases(...) result, not the sweep/progress table. "
            f"Available columns start with: {available}"
        )

    if "passed" in results.columns:
        data = results[results["passed"].astype(bool)].copy()
    else:
        data = results.copy()

    data = data.dropna(subset=[x_column, y_column, value_column])
    if data.empty:
        raise ValueError(f"No rows with {x_column}, {y_column}, and {value_column} were found.")

    x = data[x_column].to_numpy(dtype=float)
    y = data[y_column].to_numpy(dtype=float)
    z = data[value_column].to_numpy(dtype=float)

    if len(data) < 3 or data[x_column].nunique() < 2 or data[y_column].nunique() < 2:
        return {
            "data": data,
            "x": x,
            "y": y,
            "z": z,
            "xx": None,
            "yy": None,
            "zz": None,
            "x_column": x_column,
            "y_column": y_column,
            "value_column": value_column,
            "method": "points_only",
        }

    from scipy.interpolate import griddata

    xi = np.linspace(float(np.min(x)), float(np.max(x)), int(grid_size))
    yi = np.linspace(float(np.min(y)), float(np.max(y)), int(grid_size))
    xx, yy = np.meshgrid(xi, yi)

    selected_method = method
    try:
        zz = griddata(
            data[[x_column, y_column]].to_numpy(dtype=float),
            z,
            (xx, yy),
            method=method,
        )
        if np.all(np.isnan(zz)):
            raise ValueError(f"{method} griddata returned only NaN values.")
    except Exception:
        selected_method = fallback_method
        zz = griddata(
            data[[x_column, y_column]].to_numpy(dtype=float),
            z,
            (xx, yy),
            method=fallback_method,
        )

    return {
        "data": data,
        "x": x,
        "y": y,
        "z": z,
        "xx": xx,
        "yy": yy,
        "zz": zz,
        "x_column": x_column,
        "y_column": y_column,
        "value_column": value_column,
        "method": selected_method,
    }

def _format_vv_gamma_axes(
    ax,
    *,
    label: str,
    x_ticks: Sequence[float] | None = None,
    y_ticks: Sequence[float] | None = None,
) -> None:
    """Apply the common Vv-Gamma axis labels, title, and optional ticks."""

    ax.set_xlabel(r"$V_v$")
    ax.set_ylabel(r"$\Gamma_{\mathrm{eff}}\ [\mathrm{deg}]$")
    ax.set_title(label)
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
    if y_ticks is not None:
        ax.set_yticks(y_ticks)

def plot_vv_gamma_contour(
    results: pd.DataFrame,
    value_column: str,
    *,
    ax=None,
    levels: int | Sequence[float] = 12,
    x_column: str = "Vv",
    y_column: str = "Gamma_eff_deg",
    grid_size: int = 80,
    method: str = "linear",
    fallback_method: str = "nearest",
    value_label: str | None = None,
    colorbar_ticks: Sequence[float] | None = None,
    x_ticks: Sequence[float] | None = None,
    y_ticks: Sequence[float] | None = None,
    show_points: bool = True,
    show_colorbar: bool = True,
):
    """Plot a Vv-Gamma response-surface contour from completed results.

    When colorbar_ticks is specified and levels is an integer, the first and
    last ticks define the contour color range.  levels controls the number of
    contour intervals within that range.  When levels is an explicit sequence,
    that sequence continues to define the contour boundaries directly.

    x_ticks and y_ticks optionally override the Vv and Gamma_eff axis ticks.
    """

    import matplotlib.pyplot as plt

    surface = build_vv_gamma_response_surface(
        results,
        value_column,
        x_column=x_column,
        y_column=y_column,
        grid_size=grid_size,
        method=method,
        fallback_method=fallback_method,
    )

    if ax is None:
        _, ax = plt.subplots()

    label = vv_gamma_latex_label(value_column) if value_label is None else value_label

    contour = None
    colorbar = None
    if surface["zz"] is not None:
        contour_levels = levels
        extend = "neither"

        if colorbar_ticks is not None and isinstance(levels, (int, np.integer)):
            colorbar_ticks = np.asarray(colorbar_ticks, dtype=float)
            if (
                colorbar_ticks.size < 2
                or not np.all(np.isfinite(colorbar_ticks))
                or np.any(np.diff(colorbar_ticks) <= 0.0)
            ):
                raise ValueError(
                    "colorbar_ticks must contain at least two finite, "
                    "strictly increasing values."
                )

            contour_min = float(colorbar_ticks[0])
            contour_max = float(colorbar_ticks[-1])
            contour_levels = np.linspace(
                contour_min,
                contour_max,
                int(levels) + 1,
            )

            data_min = float(np.nanmin(surface["zz"]))
            data_max = float(np.nanmax(surface["zz"]))
            if data_min < contour_min and data_max > contour_max:
                extend = "both"
            elif data_min < contour_min:
                extend = "min"
            elif data_max > contour_max:
                extend = "max"

        contour = ax.contourf(
            surface["xx"],
            surface["yy"],
            surface["zz"],
            levels=contour_levels,
            cmap="jet",
            extend=extend,
        )
        if show_colorbar:
            colorbar = ax.figure.colorbar(contour, ax=ax, ticks=colorbar_ticks)
            colorbar.set_label(label)

    if show_points:
        ax.scatter(surface["x"], surface["y"], s=16)

    _format_vv_gamma_axes(ax, label=label, x_ticks=x_ticks, y_ticks=y_ticks)

    return ax, contour, colorbar, surface

def plot_vv_gamma_contour_panel(
    results: pd.DataFrame,
    column: str,
    *,
    ax=None,
    label: str | None = None,
    levels: int | Sequence[float] = 12,
    x_column: str = "Vv",
    y_column: str = "Gamma_eff_deg",
    grid_size: int = 80,
    method: str = "linear",
    fallback_method: str = "nearest",
    colorbar_ticks_by_column: Mapping[str, Sequence[float]] | None = None,
    x_ticks: Sequence[float] | None = None,
    y_ticks: Sequence[float] | None = None,
    show_points: bool = True,
    show_colorbar: bool = True,
    print_errors: bool = True,
):
    """Plot one Vv-Gamma contour panel and handle missing/failed columns.

    colorbar_ticks_by_column is keyed by the result column name, not the LaTeX
    label.  This keeps notebook-side plot settings stable even when labels are
    renamed.
    """

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()

    panel_label = vv_gamma_latex_label(column) if label is None else label
    if column not in results.columns:
        ax.text(
            0.5,
            0.5,
            f"missing\n{column}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        _format_vv_gamma_axes(ax, label=panel_label, x_ticks=x_ticks, y_ticks=y_ticks)
        return ax, None, None, None

    colorbar_ticks = None
    if colorbar_ticks_by_column is not None:
        colorbar_ticks = colorbar_ticks_by_column.get(column)

    try:
        return plot_vv_gamma_contour(
            results,
            column,
            ax=ax,
            levels=levels,
            x_column=x_column,
            y_column=y_column,
            grid_size=grid_size,
            method=method,
            fallback_method=fallback_method,
            value_label=panel_label,
            colorbar_ticks=colorbar_ticks,
            x_ticks=x_ticks,
            y_ticks=y_ticks,
            show_points=show_points,
            show_colorbar=show_colorbar,
        )
    except Exception as exc:
        ax.text(
            0.5,
            0.5,
            f"failed\n{column}\n{type(exc).__name__}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        _format_vv_gamma_axes(ax, label=panel_label, x_ticks=x_ticks, y_ticks=y_ticks)
        if print_errors:
            print(f"{column}: {exc!r}")
        return ax, None, None, None

def plot_vv_gamma_contour_lines(
    results: pd.DataFrame,
    column_settings: Mapping[str, Mapping[str, Any]],
    *,
    ax=None,
    x_column: str = "Vv",
    y_column: str = "Gamma_eff_deg",
    grid_size: int = 80,
    method: str = "linear",
    fallback_method: str = "nearest",
    x_ticks: Sequence[float] | None = None,
    y_ticks: Sequence[float] | None = None,
    show_points: bool = False,
    title: str = "",
):
    """Overlay line contours for multiple Vv-Gamma result columns.

    ``column_settings`` is keyed by result-column name. Each value may specify
    ``levels``, ``color``, ``linestyle``, ``linewidth_min``,
    ``linewidth_max``, ``linewidth_scale``, ``show_labels``,
    ``label_format``, and ``label``.
    Color identifies the result column. Line width represents the contour
    value within that column. ``linewidth_scale`` accepts ``value``,
    ``absolute``, ``value-``, and ``absolute-``; the default is ``value``.
    """
    import warnings
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    if not column_settings:
        raise ValueError("column_settings must contain at least one result column.")
    if ax is None:
        _, ax = plt.subplots()
    default_colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["C0"])
    contours: dict[str, Any] = {}
    surfaces: dict[str, dict[str, Any]] = {}
    legend_handles = []
    first_surface = None
    for column_index, (column, settings) in enumerate(column_settings.items()):
        settings = dict(settings)
        surface = build_vv_gamma_response_surface(
            results,
            column,
            x_column=x_column,
            y_column=y_column,
            grid_size=grid_size,
            method=method,
            fallback_method=fallback_method,
        )
        surfaces[column] = surface
        if first_surface is None:
            first_surface = surface
        if surface["zz"] is None:
            warnings.warn(
                f"{column}: a two-dimensional response surface could not be built; contour lines were skipped.",
                RuntimeWarning,
                stacklevel=2,
            )
            contours[column] = None
            continue
        finite_values = surface["zz"][np.isfinite(surface["zz"])]
        if finite_values.size == 0:
            warnings.warn(
                f"{column}: the response surface has no finite values; contour lines were skipped.",
                RuntimeWarning,
                stacklevel=2,
            )
            contours[column] = None
            continue
        levels = settings.get("levels", 5)
        if isinstance(levels, (int, np.integer)):
            if int(levels) < 1:
                raise ValueError(f"{column}: integer levels must be at least 1.")
            value_min = float(np.min(finite_values))
            value_max = float(np.max(finite_values))
            if math.isclose(value_min, value_max, rel_tol=0.0, abs_tol=1.0e-14):
                warnings.warn(
                    f"{column}: the response surface is constant; contour lines were skipped.",
                    RuntimeWarning,
                    stacklevel=2,
                )
                contours[column] = None
                continue
            contour_levels = np.linspace(value_min, value_max, int(levels))
        else:
            contour_levels = np.asarray(levels, dtype=float)
            if (
                contour_levels.ndim != 1
                or contour_levels.size < 1
                or not np.all(np.isfinite(contour_levels))
                or np.any(np.diff(contour_levels) <= 0.0)
            ):
                raise ValueError(
                    f"{column}: explicit levels must be a finite, strictly increasing one-dimensional sequence."
                )
        linewidth_min = float(settings.get("linewidth_min", 1.0))
        linewidth_max = float(settings.get("linewidth_max", 4.0))
        if linewidth_min <= 0.0 or linewidth_max < linewidth_min:
            raise ValueError(f"{column}: require 0 < linewidth_min <= linewidth_max.")
        linewidth_scale = settings.get("linewidth_scale", "value")
        if linewidth_scale not in {"value", "absolute", "value-", "absolute-"}:
            raise ValueError(
                f"{column}: linewidth_scale must be 'value', 'absolute', 'value-', or 'absolute-'."
            )
        width_values = (
            np.abs(contour_levels)
            if linewidth_scale.startswith("absolute")
            else contour_levels.copy()
        )
        width_value_min = float(np.min(width_values))
        width_value_max = float(np.max(width_values))
        if math.isclose(width_value_min, width_value_max, rel_tol=0.0, abs_tol=1.0e-14):
            normalized_widths = np.full(contour_levels.shape, 0.5, dtype=float)
        else:
            normalized_widths = (
                (width_values - width_value_min)
                / (width_value_max - width_value_min)
            )
        if linewidth_scale.endswith("-"):
            normalized_widths = 1.0 - normalized_widths
        linewidths = linewidth_min + normalized_widths * (linewidth_max - linewidth_min)
        color = settings.get("color", default_colors[column_index % len(default_colors)])
        linestyle = settings.get("linestyle", "-")
        contour = ax.contour(
            surface["xx"],
            surface["yy"],
            surface["zz"],
            levels=contour_levels,
            colors=[color],
            linestyles=linestyle,
            linewidths=linewidths,
        )
        contours[column] = contour
        if settings.get("show_labels", False):
            ax.clabel(
                contour,
                contour.levels,
                inline=True,
                fmt=settings.get("label_format", "%.2g"),
            )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                linestyle=linestyle,
                linewidth=0.5 * (linewidth_min + linewidth_max),
                label=settings.get("label", vv_gamma_latex_label(column)),
            )
        )
    if show_points and first_surface is not None:
        ax.scatter(first_surface["x"], first_surface["y"], s=16)
    _format_vv_gamma_axes(ax, label=title, x_ticks=x_ticks, y_ticks=y_ticks)
    if legend_handles:
        ax.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left', )
    return ax, contours, surfaces
