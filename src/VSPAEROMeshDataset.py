"""Build an append-only VSPAERO mesh-calibration dataset.

This module is intentionally a data generator, not an optimizer.  It varies all
editable mesh controls that are actually present on the active VSPAERO geometry,
runs reproducible cases, and preserves raw geometry, mesh, topology, solution,
and failure observations.  No case is labelled "good" or "bad" here; convergence
families and production sizing rules are derived later from the saved data.

The workflow is linear:

1. Load one source .vsp3 and resolve ThinWing / ThickWing / Hybrid / ThickAll.
2. Discover the supported mesh Parm whitelist on every active Geom/XSec.
3. Build a baseline, one-factor sweeps for every discovered Parm, and an
   all-parameter Latin-hypercube sample.  Optional explicit cases can be added.
4. For every case, reload the immutable source model, apply absolute requested
   Parm values, save/reload the case .vsp3, measure geometry and actual surface
   tessellation, run VSPAERO, and run VSPAEROMeshQuality diagnostics.
5. Commit every case independently inside its attempt directory.  Aggregate CSV
   tables are derived only after a normal run or by an explicit rebuild.  Completed
   deterministic case IDs are skipped on a later call, so a long campaign can be
   resumed or extended without depending on an aggregate checkpoint.

The saved parameter values are OpenVSP input controls.  Physical edge sizes,
spacing-growth metrics, curvature features, junction metrics, Cp/LOD diagnostics,
and aerodynamic results are observations.  In particular, SmallPanelW/MaxGrowth
are recomputed from the generated physical mesh instead of treated as DOE inputs.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from .VSPAEROMesh import (
    _global_u_surface_length,
    _representative_w_surface_length,
    _run_saved_vspaero_mesh_case,
    _scaled_tessellation_count,
    _set_saved_vspaero_representation,
    _surface_tessellation_metrics,
    _u_interval_edge_medians,
    _u_interval_surface_lengths,
    _u_intervals,
    _u_tessellation_parms,
    resolve_vspaero_representation,
)
from .util import import_openvsp

DATASET_SCHEMA_VERSION = 2

# Only parameters that change surface tessellation are eligible.  Geometry
# design parameters and cap-option geometry controls are deliberately excluded.
_GEOM_MESH_PARMS = {
    "Tess_W": "count",
    "Tess_U": "count",
    "LECluster": "continuous",
    "TECluster": "continuous",
    "CapUMinTess": "count",
}
_XSEC_MESH_PARMS = {
    "SectTess_U": "count",
    "InCluster": "continuous",
    "OutCluster": "continuous",
    "FwdCluster": "continuous",
    "AftCluster": "continuous",
}

_CASE_COLUMNS = [
    "case_id",
    "attempt",
    "experiment_type",
    "source_model_name",
    "source_model_path",
    "source_model_sha256",
    "representation",
    "openvsp_version",
    "alpha_deg",
    "mach",
    "reynolds_number",
    "ncpu",
    "wake_num_iter",
    "wake_num_nodes",
    "fixed_wake_flag",
    "started_at",
    "finished_at",
    "status",
    "geometry_topology_checks_passed",
    "solution_review_required",
    "elapsed_s",
    "surface_triangle_count",
    "ngon_count",
    "surface_edge_count",
    "strong_mesh_advisory_count",
    "junction_edge_count",
    "junction_edge_min",
    "junction_to_local_p50_ratio_min",
    "local_cp_spike_count",
    "cp_min_actual",
    "lod_outlier_count",
    "geometry_topology_failure_categories_json",
    "solution_advisories_json",
    "requested_parameters_json",
    "vsp3_path",
    "vsp3_sha256",
    "vspgeom_sha256",
    "adb_sha256",
    "history_sha256",
    "lod_sha256",
    "error",
]
_PARAMETER_COLUMNS = [
    "case_id",
    "attempt",
    "parameter_key",
    "scope",
    "parameter_kind",
    "geom_id",
    "geom_name",
    "geom_type",
    "xsec_surf_index",
    "xsec_index",
    "parameter_name",
    "parameter_group",
    "baseline_value",
    "requested_value",
    "effective_value",
    "lower_limit",
    "upper_limit",
]
_CATALOG_COLUMNS = [column for column in _PARAMETER_COLUMNS if column not in {
    "case_id", "attempt", "requested_value", "effective_value"
}]
_GEOM_COLUMNS = [
    "case_id",
    "attempt",
    "geom_id",
    "geom_name",
    "geom_type",
    "representation_role",
    "surface_count",
    "bbox_dx",
    "bbox_dy",
    "bbox_dz",
    "u_surface_length_median",
    "w_surface_length_median",
    "curvature_radius_min",
    "curvature_radius_p10",
    "curvature_radius_median",
    "u_edge_count",
    "w_edge_count",
    "u_edge_min",
    "u_edge_p10",
    "u_edge_median",
    "u_edge_p90",
    "u_edge_max",
    "w_edge_min",
    "w_edge_p10",
    "w_edge_median",
    "w_edge_p90",
    "w_edge_max",
    "edge_aspect_ratio_median",
    "u_adjacent_growth_p95",
    "u_adjacent_growth_max",
    "w_adjacent_growth_p95",
    "w_adjacent_growth_max",
    "small_panel_w",
    "max_growth_w",
]
_SECTION_COLUMNS = [
    "case_id",
    "attempt",
    "geom_id",
    "geom_name",
    "geom_type",
    "interval_kind",
    "section_index",
    "xsec_index",
    "u0",
    "u1",
    "surface_length",
    "actual_u_edge_median",
    "xsec_width",
    "xsec_height",
    "u_mode",
    "u_parameter_value",
    "in_cluster",
    "out_cluster",
    "fwd_cluster",
    "aft_cluster",
    "cap_u_min_tess",
]
_JUNCTION_COLUMNS = [
    "case_id",
    "attempt",
    "surface_a",
    "surface_b",
    "component_a",
    "component_b",
    "edge_length",
    "local_p50",
    "junction_to_local_p50_ratio",
    "raw_json",
]
_POLAR_COLUMNS = ["case_id", "attempt", "name", "numeric_value", "text_value"]

def _sha256_file(path: Path) -> str:
    with path.open("rb") as fp:
        return hashlib.file_digest(fp, "sha256").hexdigest()

def _write_json_atomic(path: Path, data: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    os.replace(temporary, path)

def _write_csv_atomic(path: Path, rows: Sequence[Mapping[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    for column in columns:
        if column not in frame.columns:
            frame[column] = None
    frame = frame[list(columns)]
    temporary = path.with_name(f"{path.name}.tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)

def _load_case_states(dataset_dir: Path) -> dict[str, dict[str, object]]:
    """Read resume state only from committed per-attempt case directories."""

    states: dict[str, dict[str, object]] = {}
    cases_dir = dataset_dir / "cases"
    if not cases_dir.is_dir():
        return states

    for case_dir in sorted(path for path in cases_dir.iterdir() if path.is_dir()):
        max_attempt = 0
        committed: list[tuple[int, dict[str, object]]] = []
        for attempt_dir in sorted(path for path in case_dir.iterdir() if path.is_dir() and path.name.startswith("attempt_")):
            try:
                attempt = int(attempt_dir.name.removeprefix("attempt_"))
            except ValueError:
                continue
            max_attempt = max(max_attempt, attempt)
            case_path = attempt_dir / "case.json"
            if case_path.is_file():
                committed.append((attempt, json.loads(case_path.read_text(encoding="utf-8"))))

        latest_status = None
        latest_attempt = None
        latest_elapsed_s = None
        if committed:
            latest_attempt, latest_case = max(committed, key=lambda item: item[0])
            latest_status = str(latest_case.get("status", ""))
            latest_elapsed_s = latest_case.get("elapsed_s")
        states[case_dir.name] = {
            "max_attempt": max_attempt,
            "latest_committed_attempt": latest_attempt,
            "latest_status": latest_status,
            "latest_elapsed_s": latest_elapsed_s,
        }
    return states

def rebuild_mesh_dataset_tables(output_dir: str | os.PathLike) -> dict[str, object]:
    """Rebuild aggregate CSV tables from committed case-attempt directories.

    ``cases/<case_id>/attempt_xx/case.json`` is the commit marker and source of
    truth.  Attempt directories without ``case.json`` are incomplete and are
    intentionally ignored.  This function can therefore be run after a normal
    campaign, after an interrupted campaign, or at any later time.
    """

    dataset_dir = Path(output_dir).resolve()
    cases_dir = dataset_dir / "cases"
    table_specs = {
        "cases": ("case.json", _CASE_COLUMNS),
        "parameters": ("parameters.csv", _PARAMETER_COLUMNS),
        "geoms": ("geoms.csv", _GEOM_COLUMNS),
        "sections": ("sections.csv", _SECTION_COLUMNS),
        "junctions": ("junctions.csv", _JUNCTION_COLUMNS),
        "polar": ("polar.csv", _POLAR_COLUMNS),
    }
    collected: dict[str, list[dict[str, object]]] = {name: [] for name in table_specs}

    if cases_dir.is_dir():
        for case_dir in sorted(path for path in cases_dir.iterdir() if path.is_dir()):
            for attempt_dir in sorted(path for path in case_dir.iterdir() if path.is_dir() and path.name.startswith("attempt_")):
                case_path = attempt_dir / "case.json"
                if not case_path.is_file():
                    continue
                collected["cases"].append(json.loads(case_path.read_text(encoding="utf-8")))
                for name, (filename, _columns) in table_specs.items():
                    if name == "cases":
                        continue
                    path = attempt_dir / filename
                    if path.is_file():
                        collected[name].extend(pd.read_csv(path).to_dict(orient="records"))

    paths: dict[str, Path] = {}
    for name, (_filename, columns) in table_specs.items():
        rows = collected[name]
        if rows:
            rows.sort(key=lambda row: (str(row.get("case_id", "")), int(row.get("attempt", 0))))
        path = dataset_dir / f"{name}.csv"
        _write_csv_atomic(path, rows, columns)
        paths[name] = path

    return {
        "committed_attempt_count": len(collected["cases"]),
        **{f"{name}_path": path for name, path in paths.items()},
    }

def _representation_role(geom_id: str, representation: Mapping[str, object]) -> str:
    thin = geom_id in set(representation.get("thin_geom_ids", []))
    thick = geom_id in set(representation.get("thick_geom_ids", []))
    if thin and thick:
        return "thin+thick"
    if thin:
        return "thin"
    if thick:
        return "thick"
    return "inactive"

def discover_vspaero_mesh_parameters(vsp, representation: Mapping[str, object]) -> list[dict[str, object]]:
    """Discover supported editable mesh controls on every active Geom/XSec.

    Discovery is intentionally whitelist-based.  This makes the function work
    across Wing/Fuselage/Stack/other Geom types without assuming that every type
    exposes the same controls, while preventing geometry-design Parms from being
    swept accidentally.
    """

    parameters: list[dict[str, object]] = []
    for geom_id in representation["active_geom_ids"]:
        geom_name = str(vsp.GetGeomName(geom_id))
        geom_type = str(vsp.GetGeomTypeName(geom_id))

        for parm_id in vsp.FindContainerParmIDs(geom_id):
            name = str(vsp.GetParmName(parm_id))
            kind = _GEOM_MESH_PARMS.get(name)
            if kind is None:
                continue
            parameters.append(
                {
                    "parameter_key": f"geom:{geom_id}:{name}",
                    "scope": "geom",
                    "parameter_kind": kind,
                    "geom_id": geom_id,
                    "geom_name": geom_name,
                    "geom_type": geom_type,
                    "xsec_surf_index": None,
                    "xsec_index": None,
                    "parameter_name": name,
                    "parameter_group": str(vsp.GetParmGroupName(parm_id)) if hasattr(vsp, "GetParmGroupName") else "",
                    "baseline_value": float(vsp.GetParmVal(parm_id)),
                    "lower_limit": float(vsp.GetParmLowerLimit(parm_id)) if hasattr(vsp, "GetParmLowerLimit") else -math.inf,
                    "upper_limit": float(vsp.GetParmUpperLimit(parm_id)) if hasattr(vsp, "GetParmUpperLimit") else math.inf,
                }
            )

        try:
            xsec_surf_count = int(vsp.GetNumXSecSurfs(geom_id))
        except Exception:
            xsec_surf_count = 0
        for surf_index in range(xsec_surf_count):
            xsec_surf_id = vsp.GetXSecSurf(geom_id, surf_index)
            for xsec_index in range(int(vsp.GetNumXSec(xsec_surf_id))):
                xsec_id = vsp.GetXSec(xsec_surf_id, xsec_index)
                for name, kind in _XSEC_MESH_PARMS.items():
                    try:
                        parm_id = vsp.GetXSecParm(xsec_id, name)
                    except Exception:
                        parm_id = ""
                    if not parm_id or str(parm_id).upper() == "NONE":
                        continue
                    parameters.append(
                        {
                            "parameter_key": f"xsec:{geom_id}:{surf_index}:{xsec_index}:{name}",
                            "scope": "xsec",
                            "parameter_kind": kind,
                            "geom_id": geom_id,
                            "geom_name": geom_name,
                            "geom_type": geom_type,
                            "xsec_surf_index": surf_index,
                            "xsec_index": xsec_index,
                            "parameter_name": name,
                            "parameter_group": str(vsp.GetParmGroupName(parm_id)) if hasattr(vsp, "GetParmGroupName") else "",
                            "baseline_value": float(vsp.GetParmVal(parm_id)),
                            "lower_limit": float(vsp.GetParmLowerLimit(parm_id)) if hasattr(vsp, "GetParmLowerLimit") else -math.inf,
                            "upper_limit": float(vsp.GetParmUpperLimit(parm_id)) if hasattr(vsp, "GetParmUpperLimit") else math.inf,
                        }
                    )

    # Some OpenVSP containers expose the same Parm through more than one API
    # path.  The stable parameter key is the dataset identity, so keep one row.
    unique: dict[str, dict[str, object]] = {}
    for row in parameters:
        unique[str(row["parameter_key"])] = row
    return [unique[key] for key in sorted(unique)]

def _resolve_parameter_id(vsp, parameter: Mapping[str, object]) -> str:
    if parameter["scope"] == "geom":
        for parm_id in vsp.FindContainerParmIDs(str(parameter["geom_id"])):
            if str(vsp.GetParmName(parm_id)) == str(parameter["parameter_name"]):
                return str(parm_id)
        raise KeyError(f"Missing Geom Parm {parameter['parameter_key']!r} after reload.")

    geom_id = str(parameter["geom_id"])
    surf_index = int(parameter["xsec_surf_index"])
    xsec_index = int(parameter["xsec_index"])
    xsec_surf_id = vsp.GetXSecSurf(geom_id, surf_index)
    xsec_id = vsp.GetXSec(xsec_surf_id, xsec_index)
    parm_id = vsp.GetXSecParm(xsec_id, str(parameter["parameter_name"]))
    if not parm_id or str(parm_id).upper() == "NONE":
        raise KeyError(f"Missing XSec Parm {parameter['parameter_key']!r} after reload.")
    return str(parm_id)

def _bounded_value(parameter: Mapping[str, object], value: float) -> float:
    value = float(value)
    lower = float(parameter["lower_limit"])
    upper = float(parameter["upper_limit"])
    if math.isfinite(lower):
        value = max(lower, value)
    if math.isfinite(upper):
        value = min(upper, value)
    if parameter["parameter_kind"] == "count":
        value = float(int(round(value)))
    return value

def _parameter_levels(
    parameter: Mapping[str, object],
    count_scales: Sequence[float],
    continuous_scales: Sequence[float],
) -> list[float]:
    baseline = float(parameter["baseline_value"])
    if parameter["parameter_kind"] == "count":
        values = [
            _bounded_value(parameter, _scaled_tessellation_count(baseline, float(scale)))
            for scale in count_scales
        ]
    else:
        values = [_bounded_value(parameter, baseline * float(scale)) for scale in continuous_scales]
        # A zero baseline cannot be explored multiplicatively.  In that unusual
        # case use finite Parm limits instead of silently making it constant.
        if max(values) - min(values) <= 1.0e-12:
            lower = float(parameter["lower_limit"])
            upper = float(parameter["upper_limit"])
            if math.isfinite(lower) and math.isfinite(upper) and upper > lower:
                values.extend([lower + 0.25 * (upper - lower), lower + 0.5 * (upper - lower), lower + 0.75 * (upper - lower)])
    values.append(_bounded_value(parameter, baseline))
    return sorted({float(value) for value in values})

def _latin_hypercube(sample_count: int, dimensions: int, seed: int) -> np.ndarray:
    if sample_count <= 0 or dimensions <= 0:
        return np.empty((0, max(0, dimensions)), dtype=float)
    rng = np.random.default_rng(int(seed))
    points = np.empty((sample_count, dimensions), dtype=float)
    for dimension in range(dimensions):
        values = (np.arange(sample_count, dtype=float) + rng.random(sample_count)) / sample_count
        rng.shuffle(values)
        points[:, dimension] = values
    return points

def _build_experiment_cases(
    parameters: Sequence[Mapping[str, object]],
    *,
    count_scales: Sequence[float],
    continuous_scales: Sequence[float],
    lhs_samples: int,
    lhs_seed: int,
    include_single_parameter_sweeps: bool,
    include_pairwise_extremes: bool,
    explicit_cases: Sequence[Mapping[str, float]] | None,
) -> list[dict[str, object]]:
    baseline = {str(row["parameter_key"]): float(row["baseline_value"]) for row in parameters}
    levels = {
        str(row["parameter_key"]): _parameter_levels(row, count_scales, continuous_scales)
        for row in parameters
    }
    cases: list[dict[str, object]] = [{"experiment_type": "baseline", "values": dict(baseline)}]

    if include_single_parameter_sweeps:
        for parameter in parameters:
            key = str(parameter["parameter_key"])
            for value in levels[key]:
                if math.isclose(value, baseline[key], rel_tol=0.0, abs_tol=1.0e-12):
                    continue
                values = dict(baseline)
                values[key] = value
                cases.append({"experiment_type": "single", "values": values})

    if include_pairwise_extremes:
        for i, first in enumerate(parameters):
            first_key = str(first["parameter_key"])
            first_values = [levels[first_key][0], levels[first_key][-1]]
            for second in parameters[i + 1 :]:
                second_key = str(second["parameter_key"])
                second_values = [levels[second_key][0], levels[second_key][-1]]
                for first_value in first_values:
                    for second_value in second_values:
                        values = dict(baseline)
                        values[first_key] = first_value
                        values[second_key] = second_value
                        cases.append({"experiment_type": "pairwise_extreme", "values": values})

    lhs = _latin_hypercube(int(lhs_samples), len(parameters), int(lhs_seed))
    for sample_index, sample in enumerate(lhs):
        values = dict(baseline)
        for parameter, unit_value in zip(parameters, sample):
            key = str(parameter["parameter_key"])
            parameter_levels = levels[key]
            if len(parameter_levels) == 1:
                values[key] = parameter_levels[0]
                continue
            low = parameter_levels[0]
            high = parameter_levels[-1]
            value = low + float(unit_value) * (high - low)
            values[key] = _bounded_value(parameter, value)
        cases.append({"experiment_type": f"lhs_{sample_index + 1:04d}", "values": values})

    known_keys = set(baseline)
    for index, explicit in enumerate(explicit_cases or (), start=1):
        unknown = set(explicit) - known_keys
        if unknown:
            raise KeyError(f"explicit case {index} contains unknown parameter keys: {sorted(unknown)}")
        values = dict(baseline)
        for key, value in explicit.items():
            parameter = next(row for row in parameters if str(row["parameter_key"]) == str(key))
            values[str(key)] = _bounded_value(parameter, float(value))
        cases.append({"experiment_type": f"explicit_{index:04d}", "values": values})

    # Integer rounding, limit clamping, and repeated experiment blocks can make
    # multiple requests identical.  Preserve the first human-readable origin.
    unique: dict[str, dict[str, object]] = {}
    for case in cases:
        key = json.dumps(case["values"], sort_keys=True, separators=(",", ":"))
        unique.setdefault(key, case)
    return list(unique.values())

def _case_id(
    source_sha256: str,
    openvsp_version: str,
    representation: str,
    values: Mapping[str, float],
    flight_condition: Mapping[str, object],
    solver_settings: Mapping[str, object],
) -> str:
    payload = {
        "schema": DATASET_SCHEMA_VERSION,
        "source_sha256": source_sha256,
        "openvsp_version": openvsp_version,
        "representation": representation,
        "values": {key: float(values[key]) for key in sorted(values)},
        "flight_condition": dict(flight_condition),
        "solver_settings": dict(solver_settings),
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return digest[:20]

def _polyline_curvature_radii(points: np.ndarray) -> list[float]:
    radii: list[float] = []
    for first, middle, last in zip(points[:-2], points[1:-1], points[2:]):
        a = float(np.linalg.norm(middle - first))
        b = float(np.linalg.norm(last - middle))
        c = float(np.linalg.norm(last - first))
        twice_area = float(np.linalg.norm(np.cross(middle - first, last - first)))
        if min(a, b, c) <= 1.0e-12 or twice_area <= 1.0e-12:
            continue
        radius = a * b * c / (2.0 * twice_area)
        if math.isfinite(radius) and radius > 1.0e-12:
            radii.append(radius)
    return radii

def _surface_geometry_features(vsp, geom_id: str) -> dict[str, float | int | None]:
    """Measure geometry-only scale and curvature features on normalized surfaces."""

    all_points: list[np.ndarray] = []
    radii: list[float] = []
    u_values = np.linspace(0.05, 0.95, 13)
    w_values = np.linspace(0.02, 0.98, 33)

    for surf_index in range(int(vsp.GetNumMainSurfs(geom_id))):
        u_grid = np.repeat(u_values, w_values.size)
        w_grid = np.tile(w_values, u_values.size)
        points = vsp.CompVecPnt01(geom_id, surf_index, u_grid.tolist(), w_grid.tolist())
        xyz = np.asarray([[point.x(), point.y(), point.z()] for point in points], dtype=float)
        xyz = xyz.reshape(u_values.size, w_values.size, 3)
        all_points.append(xyz.reshape(-1, 3))
        for row in xyz:
            radii.extend(_polyline_curvature_radii(row))
        for column in np.swapaxes(xyz, 0, 1):
            radii.extend(_polyline_curvature_radii(column))

    if not all_points:
        return {
            "bbox_dx": None,
            "bbox_dy": None,
            "bbox_dz": None,
            "curvature_radius_min": None,
            "curvature_radius_p10": None,
            "curvature_radius_median": None,
        }

    points = np.concatenate(all_points, axis=0)
    extent = np.max(points, axis=0) - np.min(points, axis=0)
    radius_array = np.asarray(radii, dtype=float)
    return {
        "bbox_dx": float(extent[0]),
        "bbox_dy": float(extent[1]),
        "bbox_dz": float(extent[2]),
        "curvature_radius_min": float(np.min(radius_array)) if radius_array.size else None,
        "curvature_radius_p10": float(np.quantile(radius_array, 0.10)) if radius_array.size else None,
        "curvature_radius_median": float(np.median(radius_array)) if radius_array.size else None,
    }

def _surface_mesh_distribution(vsp, geom_id: str) -> dict[str, float | int | None]:
    """Measure physical edge distributions and adjacent edge-size growth."""

    u_edges: list[float] = []
    w_edges: list[float] = []
    u_growth: list[float] = []
    w_growth: list[float] = []

    for surf_index in range(int(vsp.GetNumMainSurfs(geom_id))):
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
        u = np.linalg.norm(xyz[1:, :, :] - xyz[:-1, :, :], axis=2)
        w = np.linalg.norm(xyz[:, 1:, :] - xyz[:, :-1, :], axis=2)
        u_edges.extend(u[np.isfinite(u) & (u > 1.0e-12)].tolist())
        w_edges.extend(w[np.isfinite(w) & (w > 1.0e-12)].tolist())

        if u.shape[0] > 1:
            first = u[:-1, :]
            second = u[1:, :]
            mask = np.isfinite(first) & np.isfinite(second) & (first > 1.0e-12) & (second > 1.0e-12)
            ratio = np.maximum(first[mask] / second[mask], second[mask] / first[mask])
            u_growth.extend(ratio.tolist())
        if w.shape[1] > 1:
            first = w[:, :-1]
            second = w[:, 1:]
            mask = np.isfinite(first) & np.isfinite(second) & (first > 1.0e-12) & (second > 1.0e-12)
            ratio = np.maximum(first[mask] / second[mask], second[mask] / first[mask])
            w_growth.extend(ratio.tolist())

    def distribution(values: Sequence[float], prefix: str) -> dict[str, float | int | None]:
        array = np.asarray(values, dtype=float)
        if array.size == 0:
            return {f"{prefix}_{name}": None for name in ("min", "p10", "median", "p90", "max")}
        return {
            f"{prefix}_min": float(np.min(array)),
            f"{prefix}_p10": float(np.quantile(array, 0.10)),
            f"{prefix}_median": float(np.median(array)),
            f"{prefix}_p90": float(np.quantile(array, 0.90)),
            f"{prefix}_max": float(np.max(array)),
        }

    result = {
        "u_edge_count": len(u_edges),
        "w_edge_count": len(w_edges),
        **distribution(u_edges, "u_edge"),
        **distribution(w_edges, "w_edge"),
    }
    for values, prefix in ((u_growth, "u_adjacent_growth"), (w_growth, "w_adjacent_growth")):
        array = np.asarray(values, dtype=float)
        result[f"{prefix}_p95"] = float(np.quantile(array, 0.95)) if array.size else None
        result[f"{prefix}_max"] = float(np.max(array)) if array.size else None
    result["small_panel_w"] = result["w_edge_min"]
    result["max_growth_w"] = result["w_adjacent_growth_max"]
    return result

def _collect_geom_and_section_rows(
    vsp,
    representation: Mapping[str, object],
    case_id: str,
    attempt: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    geom_rows: list[dict[str, object]] = []
    section_rows: list[dict[str, object]] = []

    for geom_id in representation["active_geom_ids"]:
        geom_name = str(vsp.GetGeomName(geom_id))
        geom_type = str(vsp.GetGeomTypeName(geom_id))
        try:
            tess = _surface_tessellation_metrics(vsp, geom_id)
        except Exception:
            tess = {"surface_count": None, "edge_aspect_ratio": None}
        try:
            mesh_distribution = _surface_mesh_distribution(vsp, geom_id)
        except Exception:
            mesh_distribution = {}
        try:
            geometry = _surface_geometry_features(vsp, geom_id)
        except Exception:
            geometry = {}
        try:
            u_length = _global_u_surface_length(vsp, geom_id)
        except Exception:
            u_length = None
        try:
            w_length = _representative_w_surface_length(vsp, geom_id)
        except Exception:
            w_length = None

        geom_rows.append(
            {
                "case_id": case_id,
                "attempt": attempt,
                "geom_id": geom_id,
                "geom_name": geom_name,
                "geom_type": geom_type,
                "representation_role": _representation_role(str(geom_id), representation),
                "surface_count": tess.get("surface_count"),
                **geometry,
                "u_surface_length_median": u_length,
                "w_surface_length_median": w_length,
                **mesh_distribution,
                "edge_aspect_ratio_median": tess.get("edge_aspect_ratio"),
            }
        )

        try:
            u_mode, u_parm_ids = _u_tessellation_parms(vsp, geom_id)
            section_count = len(u_parm_ids) if u_mode == "SectTess_U" else 0
            intervals = _u_intervals(vsp, geom_id, section_count)
            if intervals:
                lengths = _u_interval_surface_lengths(vsp, geom_id, intervals)
                edge_medians = _u_interval_edge_medians(vsp, geom_id, intervals)
            else:
                lengths = []
                edge_medians = []
        except Exception:
            u_mode, u_parm_ids, intervals, lengths, edge_medians = "", [], [], [], []

        if intervals:
            xsec_surf_id = vsp.GetXSecSurf(geom_id, 0) if section_count else None
            for interval, length, edge_median in zip(intervals, lengths, edge_medians):
                kind, section_index, u0, u1 = interval
                row = {
                    "case_id": case_id,
                    "attempt": attempt,
                    "geom_id": geom_id,
                    "geom_name": geom_name,
                    "geom_type": geom_type,
                    "interval_kind": kind,
                    "section_index": section_index,
                    "xsec_index": None,
                    "u0": u0,
                    "u1": u1,
                    "surface_length": length,
                    "actual_u_edge_median": edge_median,
                    "u_mode": u_mode,
                    "cap_u_min_tess": None,
                }
                if kind == "section" and section_index is not None and xsec_surf_id:
                    xsec_index = int(section_index) + 1
                    xsec_id = vsp.GetXSec(xsec_surf_id, xsec_index)
                    row["xsec_index"] = xsec_index
                    for output_name, api_name in (("xsec_width", "GetXSecWidth"), ("xsec_height", "GetXSecHeight")):
                        try:
                            row[output_name] = float(getattr(vsp, api_name)(xsec_id))
                        except Exception:
                            row[output_name] = None
                    row["u_parameter_value"] = float(vsp.GetParmVal(str(u_parm_ids[int(section_index)])))
                    for output_name, parm_name in (
                        ("in_cluster", "InCluster"),
                        ("out_cluster", "OutCluster"),
                        ("fwd_cluster", "FwdCluster"),
                        ("aft_cluster", "AftCluster"),
                    ):
                        try:
                            parm_id = vsp.GetXSecParm(xsec_id, parm_name)
                            row[output_name] = float(vsp.GetParmVal(parm_id)) if parm_id and str(parm_id).upper() != "NONE" else None
                        except Exception:
                            row[output_name] = None
                else:
                    row.update({
                        "xsec_width": None,
                        "xsec_height": None,
                        "u_parameter_value": None,
                        "in_cluster": None,
                        "out_cluster": None,
                        "fwd_cluster": None,
                        "aft_cluster": None,
                    })
                    if kind in {"cap_min", "cap_max"}:
                        try:
                            cap_id = next(
                                parm_id
                                for parm_id in vsp.FindContainerParmIDs(geom_id)
                                if str(vsp.GetParmName(parm_id)) == "CapUMinTess"
                            )
                            row["cap_u_min_tess"] = float(vsp.GetParmVal(cap_id))
                        except Exception:
                            row["cap_u_min_tess"] = None
                section_rows.append(row)
        elif u_mode == "Tess_U" and u_parm_ids:
            try:
                edge_median = _surface_tessellation_metrics(vsp, geom_id)["u_edge_median"]
            except Exception:
                edge_median = None
            section_rows.append(
                {
                    "case_id": case_id,
                    "attempt": attempt,
                    "geom_id": geom_id,
                    "geom_name": geom_name,
                    "geom_type": geom_type,
                    "interval_kind": "global",
                    "section_index": None,
                    "xsec_index": None,
                    "u0": 0.0,
                    "u1": 1.0,
                    "surface_length": u_length,
                    "actual_u_edge_median": edge_median,
                    "xsec_width": None,
                    "xsec_height": None,
                    "u_mode": u_mode,
                    "u_parameter_value": float(vsp.GetParmVal(str(u_parm_ids[0]))),
                    "in_cluster": None,
                    "out_cluster": None,
                    "fwd_cluster": None,
                    "aft_cluster": None,
                    "cap_u_min_tess": None,
                }
            )

    return geom_rows, section_rows

def _evaluate_case(
    *,
    vsp,
    source_path: Path,
    source_sha256: str,
    dataset_dir: Path,
    case_id: str,
    attempt: int,
    experiment_type: str,
    requested_values: Mapping[str, float],
    parameter_catalog: Sequence[Mapping[str, object]],
    representation_name: str,
    lifting_set_name: str,
    body_set_name: str,
    thick_all_set_name: str,
    actual_version: str,
    alpha: float,
    mach: float,
    reynolds_number: float,
    ncpu: int | None,
    wake_num_iter: int | None,
    wake_num_nodes: int | None,
    fixed_wake_flag: bool | None,
) -> dict[str, object]:
    case_dir = dataset_dir / "cases" / case_id / f"attempt_{attempt:02d}"
    case_dir.mkdir(parents=True, exist_ok=True)
    case_vsp3 = case_dir / f"{source_path.stem}.{representation_name}.{case_id}.vsp3"

    parameter_rows: list[dict[str, object]] = [
        {
            "case_id": case_id,
            "attempt": attempt,
            **parameter,
            "requested_value": float(requested_values[str(parameter["parameter_key"])]),
            "effective_value": None,
        }
        for parameter in parameter_catalog
    ]
    geom_rows: list[dict[str, object]] = []
    section_rows: list[dict[str, object]] = []
    junction_rows: list[dict[str, object]] = []
    polar_rows: list[dict[str, object]] = []
    start = time.perf_counter()
    started_at = datetime.now().astimezone().isoformat()
    row: dict[str, object] = {
        "case_id": case_id,
        "attempt": attempt,
        "experiment_type": experiment_type,
        "source_model_name": source_path.name,
        "source_model_path": str(source_path),
        "source_model_sha256": source_sha256,
        "representation": representation_name,
        "openvsp_version": actual_version,
        "alpha_deg": float(alpha),
        "mach": float(mach),
        "reynolds_number": float(reynolds_number),
        "ncpu": ncpu,
        "wake_num_iter": wake_num_iter,
        "wake_num_nodes": wake_num_nodes,
        "fixed_wake_flag": fixed_wake_flag,
        "started_at": started_at,
        "status": "failed",
        "requested_parameters_json": json.dumps(dict(sorted(requested_values.items())), sort_keys=True),
        "vsp3_path": str(case_vsp3),
    }
    _write_json_atomic(
        case_dir / "request.json",
        {
            "schema_version": DATASET_SCHEMA_VERSION,
            "case_id": case_id,
            "attempt": attempt,
            "experiment_type": experiment_type,
            "source_model_path": str(source_path),
            "source_model_sha256": source_sha256,
            "representation": representation_name,
            "openvsp_version": actual_version,
            "flight_condition": {
                "alpha_deg": float(alpha),
                "mach": float(mach),
                "reynolds_number": float(reynolds_number),
            },
            "solver_settings": {
                "ncpu": ncpu,
                "wake_num_iter": wake_num_iter,
                "wake_num_nodes": wake_num_nodes,
                "fixed_wake_flag": fixed_wake_flag,
            },
            "requested_parameters": dict(sorted(requested_values.items())),
            "started_at": started_at,
        },
    )

    try:
        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(os.fspath(source_path))
        vsp.Update()
        representation = resolve_vspaero_representation(
            vsp,
            representation_name,
            lifting_set_name=lifting_set_name,
            body_set_name=body_set_name,
            thick_all_set_name=thick_all_set_name,
        )
        _set_saved_vspaero_representation(vsp, representation)

        for parameter in parameter_catalog:
            parm_id = _resolve_parameter_id(vsp, parameter)
            requested = float(requested_values[str(parameter["parameter_key"])])
            vsp.SetParmVal(parm_id, requested)
        vsp.Update()
        vsp.WriteVSPFile(os.fspath(case_vsp3), vsp.SET_ALL)
        vsp.Update()

        # All reported effective values and geometry features are measured after
        # reloading exactly the file that will be given to VSPAERO.
        vsp.ClearVSPModel()
        vsp.Update()
        vsp.ReadVSPFile(os.fspath(case_vsp3))
        vsp.Update()
        representation = resolve_vspaero_representation(
            vsp,
            representation_name,
            lifting_set_name=lifting_set_name,
            body_set_name=body_set_name,
            thick_all_set_name=thick_all_set_name,
        )
        _set_saved_vspaero_representation(vsp, representation)

        for parameter_row, parameter in zip(parameter_rows, parameter_catalog):
            parm_id = _resolve_parameter_id(vsp, parameter)
            parameter_row["effective_value"] = float(vsp.GetParmVal(parm_id))

        geom_rows, section_rows = _collect_geom_and_section_rows(vsp, representation, case_id, attempt)

        case_result = _run_saved_vspaero_mesh_case(
            vsp,
            case_vsp3,
            case_dir,
            representation,
            alpha=float(alpha),
            mach=float(mach),
            reynolds_number=float(reynolds_number),
            ncpu=ncpu,
            wake_num_iter=wake_num_iter,
            wake_num_nodes=wake_num_nodes,
            fixed_wake_flag=fixed_wake_flag,
            runtime_openvsp_version=actual_version,
            verbose=0,
        )
        polar_row = case_result["polar"].iloc[-1]
        for name, value in polar_row.items():
            try:
                numeric = float(value)
                if not math.isfinite(numeric):
                    numeric = None
            except (TypeError, ValueError):
                numeric = None
            polar_rows.append(
                {
                    "case_id": case_id,
                    "attempt": attempt,
                    "name": str(name),
                    "numeric_value": numeric,
                    "text_value": None if numeric is not None else str(value),
                }
            )

        diagnostics = case_result["diagnostics"]
        summary = diagnostics["summary"]
        mesh = summary["mesh"]
        cp = summary["cp"]
        lod = summary["lod"]
        checks = summary["checks"]
        provenance = summary["provenance"]

        for junction in diagnostics.get("junction_quality", []):
            junction_rows.append(
                {
                    "case_id": case_id,
                    "attempt": attempt,
                    "surface_a": junction.get("surface_a", junction.get("surface_1")),
                    "surface_b": junction.get("surface_b", junction.get("surface_2")),
                    "component_a": junction.get("component_a", junction.get("component_1")),
                    "component_b": junction.get("component_b", junction.get("component_2")),
                    "edge_length": junction.get("edge_length"),
                    "local_p50": junction.get("local_p50"),
                    "junction_to_local_p50_ratio": junction.get("junction_to_local_p50_ratio"),
                    "raw_json": json.dumps(junction, sort_keys=True, ensure_ascii=False),
                }
            )

        row.update(
            {
                "status": "completed",
                "geometry_topology_checks_passed": bool(checks["geometry_topology_checks_passed"]),
                "solution_review_required": bool(checks["solution_review_required"]),
                "surface_triangle_count": mesh.get("surface_triangle_count"),
                "ngon_count": mesh.get("ngon_count"),
                "surface_edge_count": mesh.get("surface_edge_count"),
                "strong_mesh_advisory_count": mesh.get("strong_mesh_advisory_count"),
                "junction_edge_count": mesh.get("junction_edge_count"),
                "junction_edge_min": (mesh.get("junction_edge_length") or {}).get("min"),
                "junction_to_local_p50_ratio_min": (mesh.get("junction_to_local_p50_ratio") or {}).get("min"),
                "local_cp_spike_count": cp.get("local_spike_count"),
                "cp_min_actual": cp.get("min_actual"),
                "lod_outlier_count": lod.get("lod_outlier_count"),
                "geometry_topology_failure_categories_json": json.dumps(checks.get("geometry_topology_failure_categories", []), ensure_ascii=False),
                "solution_advisories_json": json.dumps(checks.get("solution_advisories", []), ensure_ascii=False),
                "vsp3_sha256": (provenance.get("vsp3_file") or {}).get("sha256"),
                "vspgeom_sha256": (provenance.get("vspgeom_file") or {}).get("sha256"),
                "adb_sha256": (provenance.get("adb_file") or {}).get("sha256"),
                "history_sha256": (provenance.get("history_file") or {}).get("sha256"),
                "lod_sha256": (provenance.get("lod_file") or {}).get("sha256"),
            }
        )
    except Exception as exc:
        row["error"] = repr(exc)
        if case_vsp3.is_file():
            row["vsp3_sha256"] = _sha256_file(case_vsp3)
    finally:
        row["elapsed_s"] = time.perf_counter() - start
        row["finished_at"] = datetime.now().astimezone().isoformat()

    result = {
        "case": row,
        "parameters": parameter_rows,
        "geoms": geom_rows,
        "sections": section_rows,
        "junctions": junction_rows,
        "polar": polar_rows,
    }
    _write_csv_atomic(case_dir / "parameters.csv", parameter_rows, _PARAMETER_COLUMNS)
    _write_csv_atomic(case_dir / "geoms.csv", geom_rows, _GEOM_COLUMNS)
    _write_csv_atomic(case_dir / "sections.csv", section_rows, _SECTION_COLUMNS)
    _write_csv_atomic(case_dir / "junctions.csv", junction_rows, _JUNCTION_COLUMNS)
    _write_csv_atomic(case_dir / "polar.csv", polar_rows, _POLAR_COLUMNS)
    # case.json is written last and is the atomic commit marker for this attempt.
    _write_json_atomic(case_dir / "case.json", row)
    return result

def build_vspaero_mesh_rule_dataset(
    input_vsp3_path: str | os.PathLike,
    output_dir: str | os.PathLike,
    representation: str,
    *,
    lifting_set_name: str = "ThinGeom",
    body_set_name: str = "ThickGeom",
    thick_all_set_name: str = "ThickAll",
    count_scales: Sequence[float] = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
    continuous_scales: Sequence[float] = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0),
    lhs_samples: int = 128,
    lhs_seed: int = 0,
    include_single_parameter_sweeps: bool = True,
    include_pairwise_extremes: bool = False,
    explicit_cases: Sequence[Mapping[str, float]] | None = None,
    alpha: float = 2.0,
    mach: float = 0.1,
    reynolds_number: float = 4.4e6,
    ncpu: int | None = None,
    wake_num_iter: int | None = None,
    wake_num_nodes: int | None = None,
    fixed_wake_flag: bool | None = None,
    expected_openvsp_version: str | None = None,
    rerun_failed: bool = False,
    verbose: int = 1,
) -> dict[str, object]:
    """Append calibration cases for one model/representation to a raw dataset.

    Repeated calls may point to the same ``output_dir`` with another model,
    representation, sample count, or explicit case list.  Case IDs are derived
    from source-model content, runtime OpenVSP version, representation, flight
    and solver settings, and every requested mesh Parm value.  A completed case
    is therefore skipped automatically on resume.

    ``include_pairwise_extremes`` can create many cases (four cases for every
    parameter pair).  It is intentionally explicit rather than silently enabled;
    the default one-factor sweeps plus Latin-hypercube block already include all
    discovered freedoms and can be extended in a later call.
    """

    source_path = Path(input_vsp3_path).resolve()
    dataset_dir = Path(output_dir).resolve()
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    if lhs_samples < 0:
        raise ValueError("lhs_samples must be non-negative.")
    if any(float(value) <= 0.0 for value in [*count_scales, *continuous_scales]):
        raise ValueError("count_scales and continuous_scales must be positive.")

    dataset_dir.mkdir(parents=True, exist_ok=True)
    (dataset_dir / "cases").mkdir(exist_ok=True)
    manifest_path = dataset_dir / "dataset_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if int(manifest.get("schema_version", -1)) != DATASET_SCHEMA_VERSION:
            raise RuntimeError(
                f"Dataset schema mismatch in {manifest_path}: "
                f"expected {DATASET_SCHEMA_VERSION}, found {manifest.get('schema_version')!r}."
            )
    else:
        manifest_path.write_text(
            json.dumps(
                {
                    "schema_version": DATASET_SCHEMA_VERSION,
                    "purpose": "Raw OpenVSP/VSPAERO mesh-rule calibration observations. Derived convergence labels are intentionally excluded.",
                    "case_storage": {
                        "request.json": "attempt input record written before analysis starts",
                        "case.json": "atomic attempt commit marker written after all raw tables",
                        "parameters.csv": "requested/effective value of every control in this attempt",
                        "geoms.csv": "geometry and generated physical mesh metrics per active Geom",
                        "sections.csv": "U-section/cap geometry and generated edge metrics",
                        "junctions.csv": "post-intersection junction diagnostics",
                        "polar.csv": "VSPAERO polar outputs in long format",
                    },
                    "aggregate_tables": "Derived from committed case directories at normal completion or by rebuild_mesh_dataset_tables().",
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

    vsp = import_openvsp()
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(os.fspath(source_path))
    vsp.Update()
    actual_version = str(vsp.GetVSPVersion()) if hasattr(vsp, "GetVSPVersion") else "unknown"
    actual_version_number = actual_version.rsplit(" ", 1)[-1]
    if expected_openvsp_version and str(expected_openvsp_version) not in {actual_version, actual_version_number}:
        raise RuntimeError(
            f"OpenVSP version mismatch: expected {expected_openvsp_version!r}, actual {actual_version!r}."
        )

    representation_info = resolve_vspaero_representation(
        vsp,
        representation,
        lifting_set_name=lifting_set_name,
        body_set_name=body_set_name,
        thick_all_set_name=thick_all_set_name,
    )
    parameters = discover_vspaero_mesh_parameters(vsp, representation_info)
    if not parameters:
        raise ValueError("The selected representation exposes no supported editable mesh parameters.")

    source_sha256 = _sha256_file(source_path)
    catalog_path = dataset_dir / "parameter_catalog.csv"
    existing_catalog_keys: set[tuple[str, str, str]] = set()
    if catalog_path.exists():
        existing_catalog = pd.read_csv(catalog_path)
        if not existing_catalog.empty:
            existing_catalog_keys = set(
                zip(
                    existing_catalog["source_model_sha256"].astype(str),
                    existing_catalog["representation"].astype(str),
                    existing_catalog["parameter_key"].astype(str),
                )
            )
    catalog_rows = []
    for parameter in parameters:
        key = (source_sha256, str(representation_info["name"]), str(parameter["parameter_key"]))
        if key in existing_catalog_keys:
            continue
        catalog_rows.append(
            {
                "source_model_sha256": source_sha256,
                "representation": representation_info["name"],
                **parameter,
            }
        )
    catalog_columns = ["source_model_sha256", "representation", *_CATALOG_COLUMNS]
    if catalog_rows:
        existing_rows = existing_catalog.to_dict(orient="records") if catalog_path.exists() and not existing_catalog.empty else []
        _write_csv_atomic(catalog_path, [*existing_rows, *catalog_rows], catalog_columns)

    experiment_cases = _build_experiment_cases(
        parameters,
        count_scales=count_scales,
        continuous_scales=continuous_scales,
        lhs_samples=int(lhs_samples),
        lhs_seed=int(lhs_seed),
        include_single_parameter_sweeps=bool(include_single_parameter_sweeps),
        include_pairwise_extremes=bool(include_pairwise_extremes),
        explicit_cases=explicit_cases,
    )
    flight_condition = {
        "alpha_deg": float(alpha),
        "mach": float(mach),
        "reynolds_number": float(reynolds_number),
    }
    solver_settings = {
        "ncpu": ncpu,
        "wake_num_iter": wake_num_iter,
        "wake_num_nodes": wake_num_nodes,
        "fixed_wake_flag": fixed_wake_flag,
    }
    for case in experiment_cases:
        case["case_id"] = _case_id(
            source_sha256,
            actual_version,
            str(representation_info["name"]),
            case["values"],
            flight_condition,
            solver_settings,
        )

    plan_rows = [
        {
            "case_id": case["case_id"],
            "experiment_type": case["experiment_type"],
            "requested_parameters_json": json.dumps(case["values"], sort_keys=True),
        }
        for case in experiment_cases
    ]
    plan_path = dataset_dir / f"plan_{source_sha256[:12]}_{representation_info['name']}.csv"
    _write_csv_atomic(
        plan_path,
        plan_rows,
        ["case_id", "experiment_type", "requested_parameters_json"],
    )

    case_states = _load_case_states(dataset_dir)
    planned_run_count = 0
    recent_durations: list[float] = []
    for planned in experiment_cases:
        state = case_states.get(str(planned["case_id"]), {})
        status = state.get("latest_status")
        if status == "completed" or (status == "failed" and not rerun_failed):
            if status == "completed":
                elapsed = state.get("latest_elapsed_s")
                if elapsed is not None:
                    recent_durations.append(float(elapsed))
        else:
            planned_run_count += 1

    campaign_start = time.perf_counter()
    executed = 0
    completed = 0
    failed = 0
    skipped = 0

    for index, planned in enumerate(experiment_cases, start=1):
        case_id = str(planned["case_id"])
        state = case_states.get(case_id, {})
        status = state.get("latest_status")
        if status == "completed" or (status == "failed" and not rerun_failed):
            skipped += 1
            if verbose:
                now = datetime.now().astimezone()
                print(
                    f"[mesh dataset] {index}/{len(experiment_cases)} {case_id} "
                    f"skipped ({status}) | now={now.strftime('%Y-%m-%d %H:%M:%S %Z')}"
                )
            continue

        attempt = int(state.get("max_attempt", 0)) + 1
        if verbose:
            now = datetime.now().astimezone()
            elapsed = time.perf_counter() - campaign_start
            remaining_runs = planned_run_count - executed
            if recent_durations:
                typical_runtime = float(np.median(recent_durations[-20:]))
                eta = now + timedelta(seconds=typical_runtime * remaining_runs)
                eta_text = eta.strftime("%Y-%m-%d %H:%M:%S %Z")
            else:
                eta_text = "estimating"
            print(
                f"[mesh dataset] {index}/{len(experiment_cases)} "
                f"{planned['experiment_type']} case={case_id} attempt={attempt} "
                f"| now={now.strftime('%Y-%m-%d %H:%M:%S %Z')} "
                f"| elapsed={str(timedelta(seconds=int(elapsed)))} "
                f"| ETA={eta_text}"
            )

        result = _evaluate_case(
            vsp=vsp,
            source_path=source_path,
            source_sha256=source_sha256,
            dataset_dir=dataset_dir,
            case_id=case_id,
            attempt=attempt,
            experiment_type=str(planned["experiment_type"]),
            requested_values=planned["values"],
            parameter_catalog=parameters,
            representation_name=str(representation_info["name"]),
            lifting_set_name=lifting_set_name,
            body_set_name=body_set_name,
            thick_all_set_name=thick_all_set_name,
            actual_version=actual_version,
            alpha=float(alpha),
            mach=float(mach),
            reynolds_number=float(reynolds_number),
            ncpu=ncpu,
            wake_num_iter=wake_num_iter,
            wake_num_nodes=wake_num_nodes,
            fixed_wake_flag=fixed_wake_flag,
        )

        executed += 1
        case_runtime = float(result["case"].get("elapsed_s") or 0.0)
        if result["case"]["status"] == "completed" and case_runtime > 0.0:
            recent_durations.append(case_runtime)
        case_states[case_id] = {
            "max_attempt": attempt,
            "latest_committed_attempt": attempt,
            "latest_status": str(result["case"]["status"]),
            "latest_elapsed_s": case_runtime,
        }
        if result["case"]["status"] == "completed":
            completed += 1
        else:
            failed += 1

        if verbose:
            now = datetime.now().astimezone()
            elapsed = time.perf_counter() - campaign_start
            remaining_runs = planned_run_count - executed
            if recent_durations and remaining_runs > 0:
                typical_runtime = float(np.median(recent_durations[-20:]))
                eta = now + timedelta(seconds=typical_runtime * remaining_runs)
                eta_text = eta.strftime("%Y-%m-%d %H:%M:%S %Z")
            else:
                eta_text = "complete" if remaining_runs == 0 else "estimating"
            print(
                f"[mesh dataset] finished {index}/{len(experiment_cases)} "
                f"status={result['case']['status']} case={case_id} "
                f"| now={now.strftime('%Y-%m-%d %H:%M:%S %Z')} "
                f"| case={str(timedelta(seconds=int(case_runtime)))} "
                f"| elapsed={str(timedelta(seconds=int(elapsed)))} "
                f"| ETA={eta_text}"
            )
            if result["case"]["status"] != "completed":
                print(f"  FAILED: {result['case'].get('error')}")

    aggregates = rebuild_mesh_dataset_tables(dataset_dir)
    return {
        "dataset_dir": dataset_dir,
        "source_model_sha256": source_sha256,
        "openvsp_version": actual_version,
        "representation": representation_info,
        "parameter_catalog": pd.DataFrame(parameters),
        "plan": pd.DataFrame(plan_rows),
        "planned_case_count": len(experiment_cases),
        "executed_this_run": executed,
        "completed_this_run": completed,
        "failed_this_run": failed,
        "skipped_this_run": skipped,
        "manifest_path": manifest_path,
        "parameter_catalog_path": catalog_path,
        "plan_path": plan_path,
        **aggregates,
    }
