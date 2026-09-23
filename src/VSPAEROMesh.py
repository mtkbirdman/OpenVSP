"""OpenVSP/VSPAERO surface-tessellation initialization and convergence tools.

``equalize_vspaero_tessellation`` edits a .vsp3 model so that physical
surface-mesh edge lengths are close to the chordwise edge length of a reference
wing. ``optimize_vspaero_tessellation`` is model-name agnostic: the user picks
ThinWing / ThickWing / Hybrid / ThickAll, the active Geoms are obtained from
that representation's existing Geom Sets, and W then U tessellation are varied
in separate convergence stages while geometry, clustering and wake settings
remain fixed.

The equalization workflow uses a reference wing as follows.

The reference wing's Tess_W is supplied by the caller.  That setting is applied
first and OpenVSP's actual 3D W-edge lengths define the target mesh size.  U
resolution is then derived from geometry rather than from the model's existing
SectTess_U proportions:

1. Measure each XSec section's 3D surface length at several W stations.
2. Choose SectTess_U from section_length / target_edge_size.
3. Update OpenVSP and measure the actual tessellated U edges in each section.
4. Apply a small number of section-by-section corrections when needed.

For Geoms without SectTess_U, the same idea is applied to global Tess_U.  Tess_W
for non-reference Geoms is initialized from representative 3D W-curve length
and then corrected from OpenVSP's actual tessellation.  Active Wing end caps
are included through CapUMinTess and use the same physical-edge feedback.

Clustering parameters are deliberately left unchanged.  This module controls
Tess_W, U-direction tessellation counts, and Wing end-cap tessellation.
Post-intersection VSPAERO NGon quality remains the responsibility of
VSPAEROMeshQuality.py.

OpenVSP implementation/API basis:
- Tess_W is the Geom W-direction tessellation control.
- SectTess_U is the number of tessellated U curves for one XSec section.
- Wing/Fuselage SectTess_U values feed the U-direction tessellation vector.
- GetUWTess01 returns OpenVSP's actual tessellated U/W stations.
- CompVecPnt01 maps normalized surface coordinates to 3D model coordinates.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import time
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from .util import (
    find_container_parm,
    find_geom_parm,
    find_one_geom,
    get_xsec_parm_id,
    import_openvsp,
    results_dataframe,
    workdir,
)

_TESS_GROUP_NAMES = ("Shape",)
_TESS_W_NAMES = ("Tess_W",)
_TESS_U_NAMES = ("Tess_U",)
_SECT_TESS_U_NAMES = ("SectTess_U",)

# Interior samples avoid LE/TE coincidence and degenerate nose/tip locations.
_LENGTH_W_SAMPLES = (0.125, 0.25, 0.5, 0.75, 0.875)
_LENGTH_U_SAMPLES = (0.15, 0.30, 0.50, 0.70, 0.85)
_SECTION_POLYLINE_POINTS = 17
_W_POLYLINE_POINTS = 65

def _surface_tessellation_metrics(vsp, geom_id: str) -> dict[str, float | int]:
    """Measure actual 3D U/W edge lengths on every main surface of a Geom."""

    u_edges: list[np.ndarray] = []
    w_edges: list[np.ndarray] = []
    surface_count = int(vsp.GetNumMainSurfs(geom_id))

    for surf_index in range(surface_count):
        u_tess, w_tess = vsp.GetUWTess01(geom_id, surf_index)
        u_tess = np.asarray(u_tess, dtype=float)
        w_tess = np.asarray(w_tess, dtype=float)
        if u_tess.size < 2 or w_tess.size < 2:
            continue

        u_grid = np.repeat(u_tess, w_tess.size)
        w_grid = np.tile(w_tess, u_tess.size)
        points = vsp.CompVecPnt01(geom_id, surf_index, u_grid.tolist(), w_grid.tolist())
        xyz = np.asarray([[p.x(), p.y(), p.z()] for p in points], dtype=float)
        xyz = xyz.reshape(u_tess.size, w_tess.size, 3)

        u_edges.append(np.linalg.norm(xyz[1:, :, :] - xyz[:-1, :, :], axis=2).ravel())
        w_edges.append(np.linalg.norm(xyz[:, 1:, :] - xyz[:, :-1, :], axis=2).ravel())

    if not u_edges or not w_edges:
        raise ValueError(f"Geom '{vsp.GetGeomName(geom_id)}' has no measurable tessellated surface.")

    u_values = np.concatenate(u_edges)
    w_values = np.concatenate(w_edges)
    u_values = u_values[np.isfinite(u_values) & (u_values > 1.0e-12)]
    w_values = w_values[np.isfinite(w_values) & (w_values > 1.0e-12)]
    if u_values.size == 0 or w_values.size == 0:
        raise ValueError(f"Geom '{vsp.GetGeomName(geom_id)}' has only degenerate tessellation edges.")

    u_median = float(np.median(u_values))
    w_median = float(np.median(w_values))
    return {
        "surface_count": surface_count,
        "u_edge_count": int(u_values.size),
        "w_edge_count": int(w_values.size),
        "u_edge_median": u_median,
        "w_edge_median": w_median,
        "edge_aspect_ratio": float(max(u_median, w_median) / min(u_median, w_median)),
    }

def _u_tessellation_parms(vsp, geom_id: str) -> tuple[str, list[str]]:
    """Return the editable U-tessellation mode and Parm IDs for one Geom."""

    try:
        xsec_surf_count = int(vsp.GetNumXSecSurfs(geom_id))
    except Exception:
        xsec_surf_count = 0

    # Wing, Fuselage, and Stack use one XSecSurf.  Avoid pretending that one
    # section sequence maps cleanly onto a more exotic multi-XSecSurf Geom.
    if xsec_surf_count == 1:
        xsec_surf_id = vsp.GetXSecSurf(geom_id, 0)
        section_parm_ids: list[str] = []
        for xsec_index in range(1, int(vsp.GetNumXSec(xsec_surf_id))):
            xsec_id = vsp.GetXSec(xsec_surf_id, xsec_index)
            try:
                parm_id, _ = get_xsec_parm_id(vsp, xsec_id, _SECT_TESS_U_NAMES)
            except KeyError:
                section_parm_ids = []
                break
            section_parm_ids.append(parm_id)
        if section_parm_ids:
            return "SectTess_U", section_parm_ids

    tess_u_id = find_geom_parm(vsp, geom_id, _TESS_U_NAMES, _TESS_GROUP_NAMES)
    if tess_u_id:
        return "Tess_U", [tess_u_id]
    return "", []

def _u_intervals(vsp, geom_id: str, section_count: int) -> list[tuple[str, int | None, float, float]]:
    """Return normalized U intervals for active end caps and XSec sections."""

    cap_min_id = find_container_parm(vsp, geom_id, "CapUMinOption")
    cap_max_id = find_container_parm(vsp, geom_id, "CapUMaxOption")
    cap_min = bool(cap_min_id and int(round(vsp.GetParmVal(cap_min_id))) != 0)
    cap_max = bool(cap_max_id and int(round(vsp.GetParmVal(cap_max_id))) != 0)

    interval_count = section_count + int(cap_min) + int(cap_max)
    if interval_count <= 0:
        return []

    intervals: list[tuple[str, int | None, float, float]] = []
    offset = 0
    if cap_min:
        intervals.append(("cap_min", None, 0.0, 1.0 / interval_count))
        offset = 1

    for section_index in range(section_count):
        u0 = (offset + section_index) / interval_count
        u1 = (offset + section_index + 1) / interval_count
        intervals.append(("section", section_index, u0, u1))

    if cap_max:
        intervals.append(
            (
                "cap_max",
                None,
                (interval_count - 1) / interval_count,
                1.0,
            )
        )

    return intervals

def _u_interval_surface_lengths(
    vsp,
    geom_id: str,
    intervals: Sequence[tuple[str, int | None, float, float]],
) -> list[float]:
    """Measure representative 3D U length for each requested surface interval."""

    surface_count = int(vsp.GetNumMainSurfs(geom_id))
    interval_lengths: list[float] = []

    for _kind, _index, u0, u1 in intervals:
        measured: list[float] = []
        u_line = np.linspace(u0, u1, _SECTION_POLYLINE_POINTS)
        for surf_index in range(surface_count):
            for w in _LENGTH_W_SAMPLES:
                points = vsp.CompVecPnt01(
                    geom_id,
                    surf_index,
                    u_line.tolist(),
                    [float(w)] * u_line.size,
                )
                xyz = np.asarray([[p.x(), p.y(), p.z()] for p in points], dtype=float)
                length = float(np.linalg.norm(xyz[1:] - xyz[:-1], axis=1).sum())
                if np.isfinite(length) and length > 1.0e-12:
                    measured.append(length)
        if not measured:
            raise ValueError(
                f"Geom '{vsp.GetGeomName(geom_id)}' has an unmeasurable U interval."
            )
        interval_lengths.append(float(np.median(measured)))

    return interval_lengths

def _global_u_surface_length(vsp, geom_id: str) -> float:
    """Measure representative full-surface U length for a global Tess_U Geom."""

    measured: list[float] = []
    u_line = np.linspace(0.0, 1.0, _SECTION_POLYLINE_POINTS * 2 - 1)
    for surf_index in range(int(vsp.GetNumMainSurfs(geom_id))):
        for w in _LENGTH_W_SAMPLES:
            points = vsp.CompVecPnt01(
                geom_id,
                surf_index,
                u_line.tolist(),
                [float(w)] * u_line.size,
            )
            xyz = np.asarray([[p.x(), p.y(), p.z()] for p in points], dtype=float)
            length = float(np.linalg.norm(xyz[1:] - xyz[:-1], axis=1).sum())
            if np.isfinite(length) and length > 1.0e-12:
                measured.append(length)
    if not measured:
        raise ValueError(f"Geom '{vsp.GetGeomName(geom_id)}' has no measurable U surface length.")
    return float(np.median(measured))

def _representative_w_surface_length(vsp, geom_id: str) -> float:
    """Measure representative 3D W-curve length for choosing Tess_W."""

    measured: list[float] = []
    w_line = np.linspace(0.0, 1.0, _W_POLYLINE_POINTS)
    for surf_index in range(int(vsp.GetNumMainSurfs(geom_id))):
        for u in _LENGTH_U_SAMPLES:
            points = vsp.CompVecPnt01(
                geom_id,
                surf_index,
                [float(u)] * w_line.size,
                w_line.tolist(),
            )
            xyz = np.asarray([[p.x(), p.y(), p.z()] for p in points], dtype=float)
            length = float(np.linalg.norm(xyz[1:] - xyz[:-1], axis=1).sum())
            if np.isfinite(length) and length > 1.0e-12:
                measured.append(length)
    if not measured:
        raise ValueError(f"Geom '{vsp.GetGeomName(geom_id)}' has no measurable W surface length.")
    return float(np.median(measured))

def _u_interval_edge_medians(
    vsp,
    geom_id: str,
    intervals: Sequence[tuple[str, int | None, float, float]],
) -> list[float]:
    """Measure actual tessellated U-edge median for every requested interval."""

    values: list[list[float]] = [[] for _ in intervals]

    for surf_index in range(int(vsp.GetNumMainSurfs(geom_id))):
        u_tess, w_tess = vsp.GetUWTess01(geom_id, surf_index)
        u_tess = np.asarray(u_tess, dtype=float)
        w_tess = np.asarray(w_tess, dtype=float)
        if u_tess.size < 2 or w_tess.size < 1:
            continue

        u_grid = np.repeat(u_tess, w_tess.size)
        w_grid = np.tile(w_tess, u_tess.size)
        points = vsp.CompVecPnt01(geom_id, surf_index, u_grid.tolist(), w_grid.tolist())
        xyz = np.asarray([[p.x(), p.y(), p.z()] for p in points], dtype=float)
        xyz = xyz.reshape(u_tess.size, w_tess.size, 3)
        edge_lengths = np.linalg.norm(xyz[1:, :, :] - xyz[:-1, :, :], axis=2)
        edge_mid_u = 0.5 * (u_tess[:-1] + u_tess[1:])

        for interval_index, (_kind, _index, u0, u1) in enumerate(intervals):
            if interval_index == len(intervals) - 1:
                mask = (edge_mid_u >= u0) & (edge_mid_u <= u1)
            else:
                mask = (edge_mid_u >= u0) & (edge_mid_u < u1)
            if not np.any(mask):
                continue
            interval_values = edge_lengths[mask, :].ravel()
            interval_values = interval_values[
                np.isfinite(interval_values) & (interval_values > 1.0e-12)
            ]
            values[interval_index].extend(interval_values.tolist())

    medians: list[float] = []
    for interval, interval_values in zip(intervals, values):
        if not interval_values:
            kind, index, _u0, _u1 = interval
            label = kind if index is None else f"section {index + 1}"
            raise ValueError(
                f"Geom '{vsp.GetGeomName(geom_id)}' {label} has no measurable "
                "tessellated U edges."
            )
        medians.append(float(np.median(interval_values)))
    return medians

def _interval_values_by_kind(
    intervals: Sequence[tuple[str, int | None, float, float]],
    values: Sequence[float],
    kind: str,
) -> list[float]:
    """Return values whose interval kind matches ``kind`` in interval order."""

    return [
        float(value)
        for interval, value in zip(intervals, values)
        if interval[0] == kind
    ]

def _cap_value(
    intervals: Sequence[tuple[str, int | None, float, float]],
    values: Sequence[float],
    kind: str,
) -> float | None:
    """Return one cap interval value, or None when that cap is inactive."""

    for interval, value in zip(intervals, values):
        if interval[0] == kind:
            return float(value)
    return None

def equalize_vspaero_tessellation(
    input_vsp3_path: str | os.PathLike,
    output_vsp3_path: str | os.PathLike,
    reference_wing_name: str,
    reference_tess_w: int,
    *,
    geom_names: Sequence[str] | None = None,
    tolerance: float = 0.10,
    max_iterations: int = 4,
) -> pd.DataFrame:
    """Equalize OpenVSP surface tessellation around a reference-wing mesh size.

    The reference Wing Tess_W is fixed by the caller.  Its actual median 3D
    W-edge length becomes the physical target size.  Every XSec section then
    receives a SectTess_U derived from its own measured 3D surface length;
    existing SectTess_U proportions are not preserved.  Clustering parameters
    are preserved.

    Geoms with a global Tess_U instead of SectTess_U use their representative
    full-surface U length.  Non-reference Tess_W is initialized from a measured
    representative W-curve length.  Active Wing end caps also participate via
    ``CapUMinTess``.  OpenVSP's resulting tessellation is then re-measured and
    corrected for at most ``max_iterations`` passes.
    """

    input_path = Path(input_vsp3_path)
    output_path = Path(output_vsp3_path)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if reference_tess_w < 1:
        raise ValueError("reference_tess_w must be a positive integer.")
    if not 0.0 < tolerance < 1.0:
        raise ValueError("tolerance must be between 0 and 1.")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least 1.")

    vsp = import_openvsp()
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(input_path))
    vsp.Update()

    reference_id = find_one_geom(vsp, reference_wing_name)
    if str(vsp.GetGeomTypeName(reference_id)).lower() != "wing":
        raise ValueError(f"Reference Geom '{reference_wing_name}' is not a Wing.")

    reference_tess_w_id = find_geom_parm(vsp, reference_id, _TESS_W_NAMES, _TESS_GROUP_NAMES)
    if not reference_tess_w_id:
        raise KeyError(f"Reference Wing '{reference_wing_name}' has no Tess_W Parm.")

    requested_names = None if geom_names is None else set(geom_names)
    if requested_names is not None:
        requested_names.add(reference_wing_name)

    rows: list[dict] = []
    selected: list[dict[str, object]] = []
    for geom_id in vsp.FindGeoms():
        geom_name = str(vsp.GetGeomName(geom_id))
        if requested_names is not None and geom_name not in requested_names:
            continue

        geom_type = str(vsp.GetGeomTypeName(geom_id))
        try:
            surface_count = int(vsp.GetNumMainSurfs(geom_id))
        except Exception:
            surface_count = 0
        tess_w_id = find_geom_parm(vsp, geom_id, _TESS_W_NAMES, _TESS_GROUP_NAMES)
        u_mode, u_parm_ids = _u_tessellation_parms(vsp, geom_id)

        if surface_count < 1 or not tess_w_id:
            rows.append(
                {
                    "geom_name": geom_name,
                    "geom_type": geom_type,
                    "status": "skipped",
                    "reason": "no main surface" if surface_count < 1 else "no Tess_W Parm",
                }
            )
            continue

        section_count = len(u_parm_ids) if u_mode == "SectTess_U" else 0
        if str(geom_type).lower() == "wing" and int(vsp.GetNumXSecSurfs(geom_id)) == 1:
            xsec_surf_id = vsp.GetXSecSurf(geom_id, 0)
            section_count = max(0, int(vsp.GetNumXSec(xsec_surf_id)) - 1)

        u_intervals = _u_intervals(vsp, geom_id, section_count) if section_count else []
        has_active_cap = any(interval[0].startswith("cap_") for interval in u_intervals)
        cap_tess_id = (
            find_container_parm(vsp, geom_id, "CapUMinTess")
            if str(geom_type).lower() == "wing" and has_active_cap
            else ""
        )

        selected.append(
            {
                "geom_id": geom_id,
                "geom_name": geom_name,
                "geom_type": geom_type,
                "tess_w_id": tess_w_id,
                "u_mode": u_mode,
                "u_parm_ids": u_parm_ids,
                "u_intervals": u_intervals,
                "cap_tess_id": cap_tess_id,
            }
        )

    selected_ids = {item["geom_id"] for item in selected}
    if reference_id not in selected_ids:
        raise ValueError(f"Reference Wing '{reference_wing_name}' is not an adjustable surface Geom.")
    if requested_names is not None:
        found_names = {str(item["geom_name"]) for item in selected} | {row["geom_name"] for row in rows}
        missing_names = sorted(requested_names - found_names)
        if missing_names:
            raise ValueError(f"Geom name(s) not found: {missing_names}")

    before: dict[str, dict] = {}
    for item in selected:
        geom_id = str(item["geom_id"])
        u_parm_ids = list(item["u_parm_ids"])
        u_intervals = list(item["u_intervals"])
        interval_edge_medians = (
            _u_interval_edge_medians(vsp, geom_id, u_intervals)
            if u_intervals
            else []
        )
        cap_tess_id = str(item["cap_tess_id"])

        before[geom_id] = {
            "geom_name": item["geom_name"],
            "geom_type": item["geom_type"],
            "u_mode": item["u_mode"],
            "tess_w": int(round(vsp.GetParmVal(str(item["tess_w_id"])))),
            "u_values": [int(round(vsp.GetParmVal(pid))) for pid in u_parm_ids],
            "cap_tess": int(round(vsp.GetParmVal(cap_tess_id))) if cap_tess_id else None,
            "cap_min_edge_median": _cap_value(
                u_intervals, interval_edge_medians, "cap_min"
            ),
            "cap_max_edge_median": _cap_value(
                u_intervals, interval_edge_medians, "cap_max"
            ),
            **_surface_tessellation_metrics(vsp, geom_id),
        }

    # The caller defines reference-wing chordwise resolution.  OpenVSP's actual
    # allowed value and actual 3D spacing, not the requested integer alone, are
    # authoritative from here onward.
    vsp.SetParmVal(reference_tess_w_id, float(reference_tess_w))
    vsp.Update()
    target_edge_size = float(_surface_tessellation_metrics(vsp, reference_id)["w_edge_median"])

    geometry_targets: dict[str, dict] = {}
    for item in selected:
        geom_id = str(item["geom_id"])
        u_mode = str(item["u_mode"])
        u_parm_ids = list(item["u_parm_ids"])
        u_intervals = list(item["u_intervals"])
        tess_w_id = str(item["tess_w_id"])
        cap_tess_id = str(item["cap_tess_id"])

        w_length = _representative_w_surface_length(vsp, geom_id)
        if geom_id != reference_id:
            panel_count = max(1, int(round(w_length / target_edge_size)))
            vsp.SetParmVal(tess_w_id, float(panel_count + 1))

        interval_lengths = (
            _u_interval_surface_lengths(vsp, geom_id, u_intervals)
            if u_intervals
            else []
        )
        section_lengths = _interval_values_by_kind(
            u_intervals, interval_lengths, "section"
        )

        if u_mode == "SectTess_U":
            for parm_id, length in zip(u_parm_ids, section_lengths):
                panel_count = max(1, int(round(length / target_edge_size)))
                vsp.SetParmVal(parm_id, float(panel_count + 1))
            u_surface_lengths = section_lengths
        elif u_mode == "Tess_U":
            u_length = _global_u_surface_length(vsp, geom_id)
            panel_count = max(1, int(round(u_length / target_edge_size)))
            vsp.SetParmVal(u_parm_ids[0], float(panel_count + 1))
            u_surface_lengths = [u_length]
        else:
            u_surface_lengths = []

        cap_lengths = [
            float(length)
            for interval, length in zip(u_intervals, interval_lengths)
            if interval[0].startswith("cap_")
        ]
        if cap_tess_id and cap_lengths:
            representative_cap_length = float(np.median(cap_lengths))
            panel_count = max(1, int(round(representative_cap_length / target_edge_size)))
            vsp.SetParmVal(cap_tess_id, float(panel_count + 1))

        geometry_targets[geom_id] = {
            "w_surface_length": w_length,
            "u_surface_lengths": u_surface_lengths,
            "cap_min_surface_length": _cap_value(
                u_intervals, interval_lengths, "cap_min"
            ),
            "cap_max_surface_length": _cap_value(
                u_intervals, interval_lengths, "cap_max"
            ),
        }

    vsp.Update()

    # Feedback uses the mesh OpenVSP actually produced.  W is corrected per Geom;
    # SectTess_U is corrected per section.  Active Wing end caps share one
    # CapUMinTess control, so their measured edge medians are pooled into one
    # representative correction while both cap medians remain visible in output.
    for item in selected:
        geom_id = str(item["geom_id"])
        u_mode = str(item["u_mode"])
        u_parm_ids = list(item["u_parm_ids"])
        u_intervals = list(item["u_intervals"])
        tess_w_id = str(item["tess_w_id"])
        cap_tess_id = str(item["cap_tess_id"])
        is_reference = geom_id == reference_id

        for _ in range(max_iterations):
            metrics = _surface_tessellation_metrics(vsp, geom_id)
            old_w = int(round(vsp.GetParmVal(tess_w_id)))
            old_u = [int(round(vsp.GetParmVal(pid))) for pid in u_parm_ids]
            old_cap = int(round(vsp.GetParmVal(cap_tess_id))) if cap_tess_id else None
            interval_edge_medians = (
                _u_interval_edge_medians(vsp, geom_id, u_intervals)
                if u_intervals
                else []
            )
            changed = False

            if not is_reference:
                w_ratio = float(metrics["w_edge_median"]) / target_edge_size
                if abs(w_ratio - 1.0) > tolerance:
                    new_w = max(2, 1 + int(round(max(1, old_w - 1) * w_ratio)))
                    vsp.SetParmVal(tess_w_id, float(new_w))
                    changed = True

            if u_mode == "SectTess_U":
                section_edge_medians = _interval_values_by_kind(
                    u_intervals, interval_edge_medians, "section"
                )
                for parm_id, old_value, edge_median in zip(
                    u_parm_ids, old_u, section_edge_medians
                ):
                    ratio = edge_median / target_edge_size
                    if abs(ratio - 1.0) <= tolerance:
                        continue
                    new_value = max(
                        2,
                        1 + int(round(max(1, old_value - 1) * ratio)),
                    )
                    vsp.SetParmVal(parm_id, float(new_value))
                    changed = True
            elif u_mode == "Tess_U":
                u_ratio = float(metrics["u_edge_median"]) / target_edge_size
                if abs(u_ratio - 1.0) > tolerance:
                    old_value = old_u[0]
                    new_value = max(
                        2,
                        1 + int(round(max(1, old_value - 1) * u_ratio)),
                    )
                    vsp.SetParmVal(u_parm_ids[0], float(new_value))
                    changed = True

            if cap_tess_id:
                cap_edge_medians = [
                    float(edge_median)
                    for interval, edge_median in zip(u_intervals, interval_edge_medians)
                    if interval[0].startswith("cap_")
                ]
                if cap_edge_medians:
                    cap_ratio = float(np.median(cap_edge_medians)) / target_edge_size
                    if abs(cap_ratio - 1.0) > tolerance:
                        new_cap = max(
                            2,
                            1 + int(round(max(1, old_cap - 1) * cap_ratio)),
                        )
                        vsp.SetParmVal(cap_tess_id, float(new_cap))
                        changed = True

            if not changed:
                break

            vsp.Update()
            new_w = int(round(vsp.GetParmVal(tess_w_id)))
            new_u = [int(round(vsp.GetParmVal(pid))) for pid in u_parm_ids]
            new_cap = int(round(vsp.GetParmVal(cap_tess_id))) if cap_tess_id else None
            if new_w == old_w and new_u == old_u and new_cap == old_cap:
                break

    vsp.Update()

    for item in selected:
        geom_id = str(item["geom_id"])
        geom_name = str(item["geom_name"])
        u_mode = str(item["u_mode"])
        u_parm_ids = list(item["u_parm_ids"])
        u_intervals = list(item["u_intervals"])
        tess_w_id = str(item["tess_w_id"])
        cap_tess_id = str(item["cap_tess_id"])

        after_metrics = _surface_tessellation_metrics(vsp, geom_id)
        before_row = before[geom_id]
        after_interval_edge_medians = (
            _u_interval_edge_medians(vsp, geom_id, u_intervals)
            if u_intervals
            else []
        )
        if u_mode == "SectTess_U":
            after_section_edge_medians = _interval_values_by_kind(
                u_intervals, after_interval_edge_medians, "section"
            )
        elif u_mode == "Tess_U":
            after_section_edge_medians = [float(after_metrics["u_edge_median"])]
        else:
            after_section_edge_medians = []

        after_cap_min = _cap_value(
            u_intervals, after_interval_edge_medians, "cap_min"
        )
        after_cap_max = _cap_value(
            u_intervals, after_interval_edge_medians, "cap_max"
        )
        after_cap_values = [
            value for value in (after_cap_min, after_cap_max) if value is not None
        ]

        rows.append(
            {
                "geom_name": geom_name,
                "geom_type": item["geom_type"],
                "status": "adjusted",
                "reason": (
                    "reference wing"
                    if geom_id == reference_id
                    else ("no U tessellation Parm" if not u_parm_ids else "")
                ),
                "u_tessellation_mode": u_mode,
                "target_edge_size": target_edge_size,
                "geometry_w_surface_length": geometry_targets[geom_id]["w_surface_length"],
                "geometry_u_surface_lengths": geometry_targets[geom_id]["u_surface_lengths"],
                "geometry_cap_min_surface_length": geometry_targets[geom_id]["cap_min_surface_length"],
                "geometry_cap_max_surface_length": geometry_targets[geom_id]["cap_max_surface_length"],
                "before_tess_w": before_row["tess_w"],
                "after_tess_w": int(round(vsp.GetParmVal(tess_w_id))),
                "before_u_tess": before_row["u_values"],
                "after_u_tess": [int(round(vsp.GetParmVal(pid))) for pid in u_parm_ids],
                "before_cap_tess": before_row["cap_tess"],
                "after_cap_tess": (
                    int(round(vsp.GetParmVal(cap_tess_id)))
                    if cap_tess_id
                    else None
                ),
                "after_section_u_edge_medians": after_section_edge_medians,
                "before_cap_min_u_edge_median": before_row["cap_min_edge_median"],
                "after_cap_min_u_edge_median": after_cap_min,
                "before_cap_max_u_edge_median": before_row["cap_max_edge_median"],
                "after_cap_max_u_edge_median": after_cap_max,
                "before_u_edge_median": before_row["u_edge_median"],
                "after_u_edge_median": after_metrics["u_edge_median"],
                "before_w_edge_median": before_row["w_edge_median"],
                "after_w_edge_median": after_metrics["w_edge_median"],
                "before_edge_aspect_ratio": before_row["edge_aspect_ratio"],
                "after_edge_aspect_ratio": after_metrics["edge_aspect_ratio"],
                "u_target_ratio": float(after_metrics["u_edge_median"]) / target_edge_size,
                "w_target_ratio": float(after_metrics["w_edge_median"]) / target_edge_size,
                "cap_target_ratio": (
                    float(np.median(after_cap_values)) / target_edge_size
                    if after_cap_values
                    else None
                ),
                "surface_count": after_metrics["surface_count"],
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    vsp.WriteVSPFile(os.fspath(output_path), vsp.SET_ALL)
    vsp.Update()

    return pd.DataFrame(rows)

# VSPAERO mesh-convergence optimization ---------------------------------------

_REPRESENTATIONS = ("ThinWing", "ThickWing", "Hybrid", "ThickAll")

def _canonical_representation_name(representation: str) -> str:
    """Return one supported representation name without model-specific aliases."""

    lookup = {name.lower(): name for name in _REPRESENTATIONS}
    try:
        return lookup[str(representation).strip().lower()]
    except KeyError as exc:
        raise ValueError(
            f"representation must be one of {list(_REPRESENTATIONS)}; "
            f"got {representation!r}."
        ) from exc

def resolve_vspaero_representation(
    vsp,
    representation: str,
    *,
    lifting_set_name: str = "ThinGeom",
    body_set_name: str = "ThickGeom",
    thick_all_set_name: str = "ThickAll",
) -> dict[str, object]:
    """Resolve a user-selected VSPAERO representation to thick/thin Geom Sets.

    The representation names describe how existing model sets are assigned to
    VSPAERO.  They do not encode any Boeing/G103A geometry names:

    ``ThinWing``
        lifting set -> ThinGeomSet, no thick geometry.
    ``ThickWing``
        lifting set -> GeomSet, no thin geometry.
    ``Hybrid``
        lifting set -> ThinGeomSet, body set -> GeomSet.
    ``ThickAll``
        all-thick set -> GeomSet, no thin geometry.

    The three set names are configurable so models do not need to use the
    repository's default ``ThinGeom / ThickGeom / ThickAll`` labels.
    """

    name = _canonical_representation_name(representation)
    set_none = int(getattr(vsp, "SET_NONE", -1))
    definitions = {
        "ThinWing": {"thin": lifting_set_name, "thick": None},
        "ThickWing": {"thin": None, "thick": lifting_set_name},
        "Hybrid": {"thin": lifting_set_name, "thick": body_set_name},
        "ThickAll": {"thin": None, "thick": thick_all_set_name},
    }

    result: dict[str, object] = {"name": name}
    active_geom_ids: list[str] = []
    for side in ("thin", "thick"):
        set_name = definitions[name][side]
        if set_name is None:
            result[f"{side}_set_name"] = None
            result[f"{side}_set_index"] = set_none
            result[f"{side}_geom_ids"] = []
            result[f"{side}_geom_names"] = []
            continue

        set_index = int(vsp.GetSetIndex(set_name))
        if set_index == set_none:
            raise ValueError(
                f"representation={name!r} requires Geom Set {set_name!r}, "
                "but that set does not exist in the loaded model."
            )
        geom_ids = list(vsp.GetGeomSetAtIndex(set_index))
        if not geom_ids:
            raise ValueError(
                f"representation={name!r} requires Geom Set {set_name!r}, "
                "but that set is empty."
            )
        geom_names = [str(vsp.GetGeomName(geom_id)) for geom_id in geom_ids]
        result[f"{side}_set_name"] = set_name
        result[f"{side}_set_index"] = set_index
        result[f"{side}_geom_ids"] = geom_ids
        result[f"{side}_geom_names"] = geom_names
        for geom_id in geom_ids:
            if geom_id not in active_geom_ids:
                active_geom_ids.append(geom_id)

    result["active_geom_ids"] = active_geom_ids
    result["active_geom_names"] = [str(vsp.GetGeomName(geom_id)) for geom_id in active_geom_ids]
    return result

def _set_saved_vspaero_representation(vsp, representation: Mapping[str, object]) -> None:
    """Persist the selected thick/thin Geom Sets in VSPAEROSettings."""

    settings_id = vsp.FindContainer("VSPAEROSettings", 0)
    if not settings_id:
        raise RuntimeError("The loaded model has no VSPAEROSettings container.")

    for parm_name, key in (
        ("GeomSet", "thick_set_index"),
        ("ThinGeomSet", "thin_set_index"),
    ):
        parm_id = find_container_parm(vsp, settings_id, parm_name)
        if not parm_id:
            if parm_name == "ThinGeomSet" and int(representation[key]) == int(
                getattr(vsp, "SET_NONE", -1)
            ):
                continue
            raise RuntimeError(
                f"VSPAEROSettings does not expose {parm_name!r} in this OpenVSP build."
            )
        vsp.SetParmVal(parm_id, float(representation[key]))
    vsp.Update()

def _scaled_tessellation_count(value: int | float, scale: float) -> int:
    """Scale an OpenVSP tessellation count while preserving its interval meaning."""

    value = int(round(float(value)))
    scale = float(scale)
    if scale <= 0.0:
        raise ValueError("Tessellation scale factors must be positive.")
    return max(2, 1 + int(round(max(1, value - 1) * scale)))

def _mesh_controls(vsp, geom_id: str) -> dict[str, object] | None:
    """Return editable tessellation controls for one active VSPAERO Geom."""

    tess_w_id = find_geom_parm(vsp, geom_id, _TESS_W_NAMES, _TESS_GROUP_NAMES)
    if not tess_w_id:
        return None
    try:
        if int(vsp.GetNumMainSurfs(geom_id)) < 1:
            return None
    except Exception:
        return None

    u_mode, u_parm_ids = _u_tessellation_parms(vsp, geom_id)
    cap_tess_id = find_container_parm(vsp, geom_id, "CapUMinTess")
    return {
        "geom_id": geom_id,
        "geom_name": str(vsp.GetGeomName(geom_id)),
        "geom_type": str(vsp.GetGeomTypeName(geom_id)),
        "tess_w_id": tess_w_id,
        "u_mode": u_mode,
        "u_parm_ids": list(u_parm_ids),
        "cap_tess_id": cap_tess_id,
    }

def _scale_mesh_direction(vsp, controls: Sequence[dict[str, object]], direction: str, scale: float) -> None:
    """Scale one mesh direction for all selected Geoms in the loaded model."""

    if direction not in {"W", "U"}:
        raise ValueError("direction must be 'W' or 'U'.")

    for control in controls:
        if direction == "W":
            parm_id = str(control["tess_w_id"])
            current = int(round(vsp.GetParmVal(parm_id)))
            vsp.SetParmVal(parm_id, float(_scaled_tessellation_count(current, scale)))
            continue

        for parm_id in control["u_parm_ids"]:
            current = int(round(vsp.GetParmVal(str(parm_id))))
            vsp.SetParmVal(
                str(parm_id),
                float(_scaled_tessellation_count(current, scale)),
            )

        # Wing end-cap tessellation participates in U-direction resolution.
        cap_tess_id = str(control["cap_tess_id"] or "")
        if cap_tess_id:
            current = int(round(vsp.GetParmVal(cap_tess_id)))
            vsp.SetParmVal(
                cap_tess_id,
                float(_scaled_tessellation_count(current, scale)),
            )

    vsp.Update()

def _topology_signature(diagnostics: dict) -> dict[str, object]:
    """Build a refinement-tolerant topology fingerprint from one diagnosis."""

    surfaces = sorted(
        {
            (
                row.get("component_id"),
                row.get("surface_name", ""),
            )
            for row in diagnostics.get("surface_metadata", [])
        }
    )

    junction_pairs = set()
    for row in diagnostics.get("junction_quality", []):
        pair = tuple(sorted(tuple(item) for item in row.get("component_surface_pairs", [])))
        if pair:
            junction_pairs.add(pair)

    kutta_surfaces = []
    kutta = diagnostics.get("kutta") or {}
    for row in kutta.get("surfaces", []):
        if int(row.get("n_kutta_lists", 0) or 0) <= 0:
            continue
        kutta_surfaces.append(
            (
                row.get("component_id"),
                row.get("surface_name", ""),
                bool(row.get("kutta_coverage_passed")),
            )
        )

    summary = diagnostics["summary"]
    return {
        "surfaces": surfaces,
        "junction_pairs": sorted(junction_pairs),
        "kutta_surfaces": sorted(kutta_surfaces),
        "same_surface_non_manifold_edge_count": summary["topology"].get(
            "same_surface_non_manifold_edge_count"
        ),
    }

def _same_topology_family(first: Mapping[str, object], second: Mapping[str, object]) -> bool:
    """Return True when two cases retain the same structural mesh topology."""

    return (
        first.get("surfaces") == second.get("surfaces")
        and first.get("junction_pairs") == second.get("junction_pairs")
        and first.get("kutta_surfaces") == second.get("kutta_surfaces")
        and first.get("same_surface_non_manifold_edge_count") == 0
        and second.get("same_surface_non_manifold_edge_count") == 0
    )

def _qoi_pair_convergence(
    coarse: Mapping[str, object],
    fine: Mapping[str, object],
    qoi_tolerances: Mapping[str, float],
    qoi_scales: Mapping[str, float] | None = None,
) -> tuple[bool, dict[str, float]]:
    """Compare two mesh cases using user-selected normalized QoI changes."""

    scales = qoi_scales or {}
    errors: dict[str, float] = {}
    for name, tolerance in qoi_tolerances.items():
        a = float(coarse[name])
        b = float(fine[name])
        if not math.isfinite(a) or not math.isfinite(b):
            return False, {name: math.inf}
        denominator = max(abs(a), abs(b), abs(float(scales.get(name, 1.0e-12))))
        error = abs(b - a) / denominator
        errors[name] = error
        if error > float(tolerance):
            return False, errors
    return True, errors

def _local_health_not_worse(coarse: Mapping[str, object], fine: Mapping[str, object]) -> bool:
    """Reject a coarse candidate when refinement removes an obvious pathology."""

    for name in (
        "strong_mesh_advisory_count",
        "n_local_cp_spikes",
        "lod_outlier_count",
    ):
        a = coarse.get(name)
        b = fine.get(name)
        if a is not None and b is not None and int(a) > int(b):
            return False
    return True

def _evaluate_candidate_pair(
    coarse: Mapping[str, object],
    fine: Mapping[str, object],
    qoi_tolerances: Mapping[str, float],
    qoi_scales: Mapping[str, float] | None,
) -> dict[str, object]:
    topology_ok = _same_topology_family(
        json.loads(str(coarse["topology_signature_json"])),
        json.loads(str(fine["topology_signature_json"])),
    )
    qoi_ok, errors = _qoi_pair_convergence(coarse, fine, qoi_tolerances, qoi_scales)
    local_ok = _local_health_not_worse(coarse, fine)
    return {
        "passed": bool(
            coarse.get("feasible")
            and fine.get("feasible")
            and topology_ok
            and qoi_ok
            and local_ok
        ),
        "topology_compatible": topology_ok,
        "qoi_converged": qoi_ok,
        "local_health_not_worse": local_ok,
        "qoi_errors": errors,
    }

def _select_converged_case(
    rows: Sequence[Mapping[str, object]],
    qoi_tolerances: Mapping[str, float],
    qoi_scales: Mapping[str, float] | None,
) -> dict[str, object]:
    """Select the coarsest case confirmed against every requested finer level."""

    candidates = sorted(rows, key=lambda row: float(row["scale"]))
    feasible = [
        row
        for row in candidates
        if row.get("status") == "completed" and row.get("feasible")
    ]
    if not feasible:
        raise RuntimeError("No feasible mesh candidate completed successfully.")

    for index, candidate in enumerate(candidates[:-1]):
        if not (
            candidate.get("status") == "completed"
            and candidate.get("feasible")
        ):
            continue

        finer_cases = candidates[index + 1 :]
        if any(
            case.get("status") != "completed" or not case.get("feasible")
            for case in finer_cases
        ):
            continue

        checks = [
            _evaluate_candidate_pair(
                candidate,
                fine_case,
                qoi_tolerances,
                qoi_scales,
            )
            for fine_case in finer_cases
        ]
        if not all(check["passed"] for check in checks):
            continue

        return {
            "selected": dict(candidate),
            "converged": True,
            "adjacent_check": checks[0],
            "confirmation_check": checks[-1] if len(checks) > 1 else None,
        }

    return {
        "selected": dict(feasible[-1]),
        "converged": False,
        "adjacent_check": None,
        "confirmation_check": None,
    }

def _run_saved_vspaero_mesh_case(
    vsp,
    case_vsp3_path: str | os.PathLike,
    case_dir: str | os.PathLike,
    representation: Mapping[str, object],
    *,
    alpha: float,
    mach: float,
    reynolds_number: float,
    ncpu: int | None,
    wake_num_iter: int | None,
    wake_num_nodes: int | None,
    fixed_wake_flag: bool | None,
    runtime_openvsp_version: str,
    verbose: int,
) -> dict[str, object]:
    """Run one saved VSPAERO case and return its polar and mesh diagnostics.

    Callers remain responsible for editing and saving the .vsp3.  This function
    owns the repeated execution contract shared by convergence studies and raw
    calibration datasets: one VSPAERO sweep, required artifact checks, and one
    VSPAEROMeshQuality diagnosis of exactly the saved model.
    """

    from .AnalysisVSPAERO import vsp_sweep
    from .VSPAEROMeshQuality import analyze_vspaero_mesh_quality

    case_vsp3 = Path(case_vsp3_path).resolve()
    case_path = Path(case_dir).resolve()
    start = time.perf_counter()
    with workdir(case_path):
        result_id = vsp_sweep(
            vsp=vsp,
            alpha=[float(alpha)],
            mach=[float(mach)],
            reynolds=[float(reynolds_number)],
            verbose=verbose,
            ncpu=ncpu,
            wake_num_iter=wake_num_iter,
            wake_num_nodes=wake_num_nodes,
            fixed_wake_flag=fixed_wake_flag,
            thick_geom_set=int(representation["thick_set_index"]),
            thin_geom_set=int(representation["thin_set_index"]),
        )
    elapsed_s = time.perf_counter() - start

    polar = results_dataframe(vsp, result_id, "VSPAERO_Polar")
    if polar.empty:
        polar = results_dataframe(vsp, result_id, "VSPAERO Polar")
    if polar.empty:
        raise RuntimeError(f"{case_vsp3.stem}: VSPAERO returned no polar result.")

    stem = case_vsp3.stem
    adb_path = case_path / f"{stem}.adb"
    vspgeom_path = case_path / f"{stem}.vspgeom"
    if not adb_path.is_file() or not vspgeom_path.is_file():
        raise FileNotFoundError(
            f"Expected {adb_path.name} and {vspgeom_path.name} after VSPAERO."
        )
    history_path = case_path / f"{stem}.history"
    polar_path = case_path / f"{stem}.polar"
    lod_path = case_path / f"{stem}.lod"
    vspaero_path = case_path / f"{stem}.vspaero"
    diagnostics = analyze_vspaero_mesh_quality(
        adb_path=adb_path,
        vspgeom_path=vspgeom_path,
        output_dir=case_path / "mesh_quality",
        history_path=history_path if history_path.is_file() else None,
        polar_path=polar_path if polar_path.is_file() else None,
        lod_path=lod_path if lod_path.is_file() else None,
        vspaero_path=vspaero_path if vspaero_path.is_file() else None,
        vsp3_path=case_vsp3,
        runtime_openvsp_version=runtime_openvsp_version,
    )
    return {
        "result_id": result_id,
        "polar": polar,
        "diagnostics": diagnostics,
        "elapsed_s": elapsed_s,
    }

def optimize_vspaero_tessellation(
    input_vsp3_path: str | os.PathLike,
    output_dir: str | os.PathLike,
    representation: str,
    *,
    qoi_tolerances: Mapping[str, float],
    mesh_targets: Sequence[str] | None = None,
    lifting_set_name: str = "ThinGeom",
    body_set_name: str = "ThickGeom",
    thick_all_set_name: str = "ThickAll",
    tess_w_scales: Sequence[float] = (0.6, 1.0, 1.5),
    u_scales: Sequence[float] = (0.6, 1.0, 1.5),
    qoi_scales: Mapping[str, float] | None = None,
    alpha: float = 2.0,
    mach: float = 0.1,
    reynolds_number: float = 4.4e6,
    ncpu: int | None = None,
    wake_num_iter: int | None = None,
    wake_num_nodes: int | None = None,
    fixed_wake_flag: bool | None = None,
    expected_openvsp_version: str | None = None,
    verbose: int = 1,
) -> dict[str, object]:
    """Select a low-cost VSPAERO tessellation by staged mesh convergence.

    The workflow is deliberately model-name agnostic.  The user selects one of
    ``ThinWing / ThickWing / Hybrid / ThickAll``; that representation is
    resolved from existing Geom Sets, and only active Geoms with editable
    surface tessellation are considered.  ``mesh_targets`` can optionally
    restrict the active Geoms by unique name or exact OpenVSP Geom ID.

    Search order is intentionally simple and traceable:

    1. Scale ``Tess_W`` for all targets while U controls and clustering remain
       fixed.
    2. Start from the selected W case and scale ``SectTess_U`` / ``Tess_U``
       (plus active Wing end-cap tessellation) while W remains fixed.
    3. Choose the coarsest case whose user-selected QoIs are converged against
       finer topology-compatible cases and whose obvious local diagnostics are
       not worse.  If no level satisfies the criterion, retain the finest
       feasible case and report that convergence was not demonstrated.

    Geometry, Thin/Thick set membership, clustering parameters, flight
    condition, and wake settings are never optimized here.  They remain fixed
    so this function measures mesh-discretization sensitivity rather than a
    mixture of mesh, geometry, model-form, and wake effects.
    """

    input_path = Path(input_vsp3_path).resolve()
    output_path = Path(output_dir).resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if not qoi_tolerances:
        raise ValueError("qoi_tolerances must contain at least one VSPAERO output column.")

    w_scales = sorted({float(value) for value in tess_w_scales})
    u_scale_values = sorted({float(value) for value in u_scales})
    if len(w_scales) < 2 or len(u_scale_values) < 2:
        raise ValueError("tess_w_scales and u_scales must each contain at least two values.")
    if any(not math.isfinite(value) or value <= 0.0 for value in [*w_scales, *u_scale_values]):
        raise ValueError("Mesh scale factors must be finite and positive.")
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in qoi_tolerances.values()):
        raise ValueError("qoi_tolerances must contain finite non-negative values.")
    if qoi_scales and any(
        not math.isfinite(float(value)) or float(value) <= 0.0
        for value in qoi_scales.values()
    ):
        raise ValueError("qoi_scales must contain finite positive values.")

    vsp = import_openvsp()
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(input_path))
    vsp.Update()

    actual_version = str(vsp.GetVSPVersion()) if hasattr(vsp, "GetVSPVersion") else "unknown"
    actual_version_number = actual_version.rsplit(" ", 1)[-1]
    if expected_openvsp_version and str(expected_openvsp_version) not in {
        actual_version,
        actual_version_number,
    }:
        raise RuntimeError(
            f"OpenVSP version mismatch: expected {expected_openvsp_version!r}, "
            f"actual {actual_version!r}."
        )

    representation_info = resolve_vspaero_representation(
        vsp,
        representation,
        lifting_set_name=lifting_set_name,
        body_set_name=body_set_name,
        thick_all_set_name=thick_all_set_name,
    )

    active_ids = list(representation_info["active_geom_ids"])
    active_id_set = set(active_ids)
    active_by_name: dict[str, list[str]] = {}
    for geom_id in active_ids:
        active_by_name.setdefault(str(vsp.GetGeomName(geom_id)), []).append(geom_id)

    if mesh_targets is None:
        target_ids = active_ids
    else:
        target_ids = []
        for target in mesh_targets:
            target = str(target)
            if target in active_id_set:
                target_ids.append(target)
                continue
            matches = active_by_name.get(target, [])
            if len(matches) == 1:
                target_ids.append(matches[0])
                continue
            if len(matches) > 1:
                raise ValueError(
                    f"mesh target name {target!r} is not unique; pass the OpenVSP Geom ID instead."
                )
            raise ValueError(
                f"mesh target {target!r} is not active in representation "
                f"{representation_info['name']!r}."
            )
        target_ids = list(dict.fromkeys(target_ids))

    target_controls = []
    skipped_active_geoms = []
    for geom_id in target_ids:
        control = _mesh_controls(vsp, geom_id)
        if control is None:
            skipped_active_geoms.append(str(vsp.GetGeomName(geom_id)))
        else:
            target_controls.append(control)
    if mesh_targets is not None and skipped_active_geoms:
        raise ValueError(
            "Explicit mesh target(s) have no editable surface Tess_W: "
            f"{skipped_active_geoms}"
        )
    if not target_controls:
        raise ValueError("The selected representation has no mesh targets with editable Tess_W.")

    output_path.mkdir(parents=True, exist_ok=True)
    for stage_dir_name in ("w", "u", "selected"):
        stage_dir = output_path / stage_dir_name
        if stage_dir.exists():
            shutil.rmtree(stage_dir)
    results_csv = output_path / "mesh_optimization_results.csv"
    rows: list[dict[str, object]] = []

    def run_case(base_vsp3_path: Path, stage: str, scale: float, case_index: int) -> dict[str, object]:
        case_name = f"{stage.lower()}_{case_index:02d}_{str(scale).replace('.', 'p')}"
        case_dir = output_path / stage.lower() / case_name
        if case_dir.exists():
            shutil.rmtree(case_dir)
        case_dir.mkdir(parents=True)
        case_vsp3 = case_dir / f"{input_path.stem}.{representation_info['name']}.{case_name}.vsp3"

        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(os.fspath(base_vsp3_path))
        vsp.Update()

        current_representation = resolve_vspaero_representation(
            vsp,
            str(representation_info["name"]),
            lifting_set_name=lifting_set_name,
            body_set_name=body_set_name,
            thick_all_set_name=thick_all_set_name,
        )
        _set_saved_vspaero_representation(vsp, current_representation)

        current_active_ids = set(current_representation["active_geom_ids"])
        current_controls = []
        for original in target_controls:
            geom_id = str(original["geom_id"])
            name = str(original["geom_name"])
            if geom_id not in current_active_ids:
                raise RuntimeError(
                    f"Geom {name!r} ({geom_id}) is no longer active after model reload."
                )
            control = _mesh_controls(vsp, geom_id)
            if control is None:
                raise RuntimeError(
                    f"Geom {name!r} lost editable tessellation controls after reload."
                )
            current_controls.append(control)

        _scale_mesh_direction(vsp, current_controls, stage, scale)
        vsp.WriteVSPFile(os.fspath(case_vsp3), vsp.SET_ALL)
        vsp.Update()

        # Reload exactly what was written before the solver run.  This prevents
        # unsaved in-memory parameter state from entering the comparison.
        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(os.fspath(case_vsp3))
        vsp.Update()

        effective_mesh_settings = {}
        for original in current_controls:
            geom_id = str(original["geom_id"])
            control = _mesh_controls(vsp, geom_id)
            if control is None:
                raise RuntimeError(
                    f"Geom {original['geom_name']!r} lost mesh controls after save/reload."
                )
            cap_id = str(control["cap_tess_id"] or "")
            effective_mesh_settings[geom_id] = {
                "geom_name": control["geom_name"],
                "Tess_W": int(round(vsp.GetParmVal(str(control["tess_w_id"])))),
                "U_mode": control["u_mode"],
                "U_values": [
                    int(round(vsp.GetParmVal(str(parm_id))))
                    for parm_id in control["u_parm_ids"]
                ],
                "CapUMinTess": int(round(vsp.GetParmVal(cap_id))) if cap_id else None,
            }

        case_result = _run_saved_vspaero_mesh_case(
            vsp,
            case_vsp3,
            case_dir,
            current_representation,
            alpha=float(alpha),
            mach=float(mach),
            reynolds_number=float(reynolds_number),
            ncpu=ncpu,
            wake_num_iter=wake_num_iter,
            wake_num_nodes=wake_num_nodes,
            fixed_wake_flag=fixed_wake_flag,
            runtime_openvsp_version=actual_version,
            verbose=verbose,
        )
        elapsed_s = float(case_result["elapsed_s"])
        polar_row = case_result["polar"].iloc[-1]
        diagnostics = case_result["diagnostics"]
        summary = diagnostics["summary"]
        mesh = summary["mesh"]
        cp = summary["cp"]
        lod = summary["lod"]
        checks = summary["checks"]

        provenance = summary["provenance"]
        row: dict[str, object] = {
            "stage": stage,
            "case_name": case_name,
            "scale": float(scale),
            "status": "completed",
            "feasible": bool(checks["geometry_topology_checks_passed"]),
            "representation": representation_info["name"],
            "thick_set_name": current_representation["thick_set_name"],
            "thick_set_index": current_representation["thick_set_index"],
            "thin_set_name": current_representation["thin_set_name"],
            "thin_set_index": current_representation["thin_set_index"],
            "mesh_targets_json": json.dumps(
                [control["geom_name"] for control in current_controls],
                ensure_ascii=False,
            ),
            "mesh_settings_json": json.dumps(
                effective_mesh_settings,
                ensure_ascii=False,
                sort_keys=True,
            ),
            "openvsp_version": actual_version,
            "elapsed_s": elapsed_s,
            "vsp3_path": str(case_vsp3),
            "vsp3_sha256": (provenance.get("vsp3_file") or {}).get("sha256"),
            "vspgeom_sha256": (provenance.get("vspgeom_file") or {}).get("sha256"),
            "adb_sha256": (provenance.get("adb_file") or {}).get("sha256"),
            "history_sha256": (provenance.get("history_file") or {}).get("sha256"),
            "lod_sha256": (provenance.get("lod_file") or {}).get("sha256"),
            "n_surface_triangles": mesh.get("surface_triangle_count"),
            "n_ngons": mesh.get("ngon_count"),
            "strong_mesh_advisory_count": mesh.get("strong_mesh_advisory_count"),
            "junction_edge_min": (mesh.get("junction_edge_length") or {}).get("min"),
            "junction_to_local_p50_ratio_min": (
                mesh.get("junction_to_local_p50_ratio") or {}
            ).get("min"),
            "n_local_cp_spikes": cp.get("local_spike_count"),
            "cp_min_actual": cp.get("min_actual"),
            "lod_outlier_count": lod.get("lod_outlier_count"),
            "geometry_topology_failure_categories_json": json.dumps(
                checks.get("geometry_topology_failure_categories", []),
                ensure_ascii=False,
            ),
            "topology_signature_json": json.dumps(
                _topology_signature(diagnostics),
                sort_keys=True,
                ensure_ascii=False,
            ),
        }

        missing_qoi = []
        for qoi_name in qoi_tolerances:
            value = polar_row.get(qoi_name)
            if value is None or not math.isfinite(float(value)):
                row[qoi_name] = math.nan
                missing_qoi.append(qoi_name)
            else:
                row[qoi_name] = float(value)
        if missing_qoi:
            row["feasible"] = False
            row["geometry_topology_failure_categories_json"] = json.dumps(
                [
                    *json.loads(str(row["geometry_topology_failure_categories_json"])),
                    f"missing_qoi:{','.join(missing_qoi)}",
                ],
                ensure_ascii=False,
            )

        return row

    # Stage 1: W-direction convergence from one immutable input baseline.
    w_rows = []
    for index, scale in enumerate(w_scales, start=1):
        if verbose:
            print(f"\n[mesh optimization] W {index}/{len(w_scales)} scale={scale:g}")
        try:
            row = run_case(input_path, "W", scale, index)
        except Exception as exc:
            row = {
                "stage": "W",
                "case_name": f"w_{index:02d}_{str(scale).replace('.', 'p')}",
                "scale": float(scale),
                "status": "failed",
                "feasible": False,
                "representation": representation_info["name"],
                "openvsp_version": actual_version,
                "error": repr(exc),
            }
            if verbose:
                print(f"  FAILED: {exc}")
        rows.append(row)
        w_rows.append(row)
        pd.DataFrame(rows).to_csv(results_csv, index=False)

    w_selection = _select_converged_case(w_rows, qoi_tolerances, qoi_scales)
    selected_w_path = Path(str(w_selection["selected"]["vsp3_path"]))

    # Stage 2: U-direction convergence begins from the selected W mesh.
    u_rows = []
    for index, scale in enumerate(u_scale_values, start=1):
        if verbose:
            print(f"\n[mesh optimization] U {index}/{len(u_scale_values)} scale={scale:g}")
        try:
            row = run_case(selected_w_path, "U", scale, index)
        except Exception as exc:
            row = {
                "stage": "U",
                "case_name": f"u_{index:02d}_{str(scale).replace('.', 'p')}",
                "scale": float(scale),
                "status": "failed",
                "feasible": False,
                "representation": representation_info["name"],
                "openvsp_version": actual_version,
                "error": repr(exc),
            }
            if verbose:
                print(f"  FAILED: {exc}")
        rows.append(row)
        u_rows.append(row)
        pd.DataFrame(rows).to_csv(results_csv, index=False)

    u_selection = _select_converged_case(u_rows, qoi_tolerances, qoi_scales)
    selected_row = u_selection["selected"]
    selected_source = Path(str(selected_row["vsp3_path"]))
    selected_dir = output_path / "selected"
    selected_dir.mkdir(parents=True, exist_ok=True)
    selected_vsp3 = selected_dir / f"{input_path.stem}.{representation_info['name']}.optimized.vsp3"
    shutil.copy2(selected_source, selected_vsp3)

    table = pd.DataFrame(rows)
    table["selected_w_stage"] = (
        (table["stage"] == "W")
        & (table["case_name"] == str(w_selection["selected"]["case_name"]))
    )
    table["selected_u_stage"] = (
        (table["stage"] == "U")
        & (table["case_name"] == str(u_selection["selected"]["case_name"]))
    )

    # Pairwise convergence columns make the CSV independently reviewable.
    for stage_name, stage_rows in (("W", w_rows), ("U", u_rows)):
        ordered = sorted(stage_rows, key=lambda row: float(row["scale"]))
        for first, second in zip(ordered[:-1], ordered[1:]):
            mask = (table["stage"] == stage_name) & (table["case_name"] == first["case_name"])
            if not (
                first.get("status") == "completed"
                and second.get("status") == "completed"
                and "topology_signature_json" in first
                and "topology_signature_json" in second
            ):
                table.loc[mask, "topology_compatible_to_next"] = False
                table.loc[mask, "qoi_converged_to_next"] = False
                table.loc[mask, "local_health_not_worse_than_next"] = False
                table.loc[mask, "pair_passed_to_next"] = False
                continue

            check = _evaluate_candidate_pair(first, second, qoi_tolerances, qoi_scales)
            table.loc[mask, "topology_compatible_to_next"] = bool(check["topology_compatible"])
            table.loc[mask, "qoi_converged_to_next"] = bool(check["qoi_converged"])
            table.loc[mask, "local_health_not_worse_than_next"] = bool(check["local_health_not_worse"])
            table.loc[mask, "pair_passed_to_next"] = bool(check["passed"])
            table.loc[mask, "qoi_error_to_next_json"] = json.dumps(
                check["qoi_errors"], sort_keys=True
            )

    table.to_csv(results_csv, index=False)

    with input_path.open("rb") as fp:
        input_sha256 = hashlib.file_digest(fp, "sha256").hexdigest()
    with selected_vsp3.open("rb") as fp:
        selected_sha256 = hashlib.file_digest(fp, "sha256").hexdigest()

    summary = {
        "input_vsp3": str(input_path),
        "input_vsp3_sha256": input_sha256,
        "representation": {
            key: value
            for key, value in representation_info.items()
            if key not in {"active_geom_ids", "thin_geom_ids", "thick_geom_ids"}
        },
        "mesh_targets": [control["geom_name"] for control in target_controls],
        "skipped_active_geoms": skipped_active_geoms,
        "openvsp_version": actual_version,
        "expected_openvsp_version": expected_openvsp_version,
        "flight_condition": {
            "alpha_deg": float(alpha),
            "mach": float(mach),
            "reynolds_number": float(reynolds_number),
        },
        "wake_settings": {
            "wake_num_iter": wake_num_iter,
            "wake_num_nodes": wake_num_nodes,
            "fixed_wake_flag": fixed_wake_flag,
        },
        "qoi_tolerances": dict(qoi_tolerances),
        "qoi_scales": dict(qoi_scales or {}),
        "tess_w_scales": w_scales,
        "u_scales": u_scale_values,
        "w_stage": {
            "converged": bool(w_selection["converged"]),
            "selected_case": w_selection["selected"]["case_name"],
            "selected_scale": w_selection["selected"]["scale"],
            "selected_case_vsp3": w_selection["selected"].get("vsp3_path"),
            "adjacent_check": w_selection["adjacent_check"],
            "confirmation_check": w_selection["confirmation_check"],
        },
        "u_stage": {
            "converged": bool(u_selection["converged"]),
            "selected_case": u_selection["selected"]["case_name"],
            "selected_scale": u_selection["selected"]["scale"],
            "selected_case_vsp3": u_selection["selected"].get("vsp3_path"),
            "adjacent_check": u_selection["adjacent_check"],
            "confirmation_check": u_selection["confirmation_check"],
        },
        "fully_converged": bool(w_selection["converged"] and u_selection["converged"]),
        "selected_vsp3": str(selected_vsp3),
        "selected_vsp3_sha256": selected_sha256,
        "results_csv": str(results_csv),
        "note": (
            "Selection is based on user-selected QoI convergence inside a structurally "
            "compatible mesh family. Geometry, representation, clustering and wake settings "
            "are fixed and are not optimization variables."
        ),
    }
    summary_json = output_path / "mesh_optimization_summary.json"
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    return {
        "summary": summary,
        "cases": table,
        "selected_vsp3_path": selected_vsp3,
        "results_csv_path": results_csv,
        "summary_json_path": summary_json,
    }
