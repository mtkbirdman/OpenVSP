"""OpenVSP/VSPAERO surface-tessellation equalization.

``equalize_vspaero_tessellation`` edits a .vsp3 model so that physical
surface-mesh edge lengths are close to the chordwise edge length of a reference
wing.

The caller sets the reference wing's ``Tess_W``. OpenVSP's actual 3D W-edge
median then becomes the target physical mesh size. U-direction tessellation is
derived from measured geometry instead of preserving the model's existing
``SectTess_U`` proportions:

1. Measure each XSec section's 3D surface length.
2. Choose ``SectTess_U`` or global ``Tess_U`` from length / target mesh size.
3. Initialize non-reference ``Tess_W`` from representative 3D W-curve length.
4. Include active Wing end caps through ``CapUMinTess``.
5. Re-measure OpenVSP's actual tessellation and apply a small number of
   corrections when the physical edge size is outside the requested tolerance.

Clustering parameters are deliberately left unchanged. This module only
initializes OpenVSP surface tessellation; VSPAERO convergence studies and
post-intersection mesh diagnostics are separate concerns.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .util import (
    find_container_parm,
    find_geom_parm,
    find_one_geom,
    get_xsec_parm_id,
    import_openvsp,
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
    """Measure the U/W edge medians used by the equalization feedback loop."""

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
        xyz = np.asarray([[point.x(), point.y(), point.z()] for point in points], dtype=float)
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
        "u_edge_median": u_median,
        "w_edge_median": w_median,
        "edge_aspect_ratio": max(u_median, w_median) / min(u_median, w_median),
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
