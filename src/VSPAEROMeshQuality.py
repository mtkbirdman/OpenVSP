"""
VSPAEROMeshQuality.py -- v4

OpenVSP / VSPAERO の surface mesh, Cp, Kutta / wake topology と、
.history / .polar / .lod に保存された solution-level diagnostics を同じレポートに整理する。

v4 の設計方針
----------------
1. ADB の surface triangle と Cp を authoritative geometry / pressure source とする。
2. VSPGeom の original NGon, alternate triangulation, UV, Kutta lists を対応付ける。
3. mesh / Cp と Kutta / wake topology を独立した verification 軸として保持する。
4. .history / .polar / .lod は solution-level verification の補助情報として読む。
5. surface-derived force と wake-derived force の不一致は advisory とし、それだけで FAIL にしない。
6. LOD の極端値や左右非対称も heuristic advisory とし、閾値は設定可能にする。
7. <=15 deg の child-triangle small-angle 件数は backward-compatible advisory として残し、
   severe mesh warning は極端な角度 / area transition / Cp anomaly の組み合わせで別集計する。
8. mixed thick/thin では cross-surface junction 自体は正常に存在し得るため、件数だけで FAIL にしない。
9. mapping failure, same-surface non-manifold, incomplete Kutta coverage は structural failure として扱う。
10. Rotor / nozzle を含む ADB は、未対応データを黙って誤読せず明示的に停止する。

Reference implementation baseline
---------------------------------
OpenVSP/OpenVSP tag OpenVSP_3.51.3
  src/vsp_aero/Viewer/glviewer.C
  src/vsp_aero/Solver/VSP_Solver.C
  src/vsp_aero/Solver/VSP_Geom.C
  src/geom_core/MeshGeom.cpp
  src/geom_core/TMesh.cpp

v4 は OpenVSP の validity criterion を新たに定義するものではない。
出力される threshold-based flags は verification を進めるための diagnostic / advisory である。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import struct
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Iterable

import numpy as np

MESH_QUALITY_TOOL_VERSION = 4
OPENVSP_REFERENCE_VERSION = "3.51.3"
ADB_MAGIC_V2 = -123789456
ADB_MAGIC_V3 = ADB_MAGIC_V2 + 3
WAKE_EDGE_TOLERANCE = 1.0e-12
DEFAULT_TE_PARAM_NEAR_MISS_BAND = 1.0e-3

class BinaryReader:
    """Small binary helper. ADB uses 32-bit int/float and 64-bit double."""

    def __init__(self, fp, endian: str):
        self.fp = fp
        self.endian = endian

    def _read(self, fmt: str):
        size = struct.calcsize(self.endian + fmt)
        raw = self.fp.read(size)
        if len(raw) != size:
            raise EOFError(f"Unexpected EOF at byte {self.fp.tell() - len(raw)}")
        return struct.unpack(self.endian + fmt, raw)[0]

    def i32(self) -> int:
        return self._read("i")

    def f32(self) -> float:
        return self._read("f")

    def f64(self) -> float:
        return self._read("d")

    def bytes(self, n: int) -> bytes:
        raw = self.fp.read(n)
        if len(raw) != n:
            raise EOFError(f"Unexpected EOF at byte {self.fp.tell() - len(raw)}")
        return raw

def _sha256(path: str | Path | None) -> str | None:
    if path is None:
        return None
    path = Path(path)
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as fp:
        for chunk in iter(lambda: fp.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def _detect_adb_endian_and_version(fp) -> tuple[str, int]:
    raw = fp.read(4)
    if len(raw) != 4:
        raise ValueError("ADB file is too short.")
    for endian in ("<", ">"):
        value = struct.unpack(endian + "i", raw)[0]
        if value == ADB_MAGIC_V3:
            return endian, 3
        if value == ADB_MAGIC_V2:
            return endian, 2
    raise ValueError("ADB magic number is not recognized.")

def _read_adb_header(reader: BinaryReader, version: int) -> dict:
    header = {
        "version": version,
        "model_type": reader.i32(),
        "symmetry_flag": reader.i32(),
        "time_accurate": reader.i32(),
        "number_of_vortex_loops": reader.i32(),
        "number_of_nodes": reader.i32(),
        "number_of_tris": reader.i32(),
        "number_of_surface_vortex_edges": reader.i32(),
        "sref": reader.f32(),
        "cref": reader.f32(),
        "bref": reader.f32(),
        "xcg": reader.f32(),
        "ycg": reader.f32(),
        "zcg": reader.f32(),
    }

    n_cart3d = reader.i32()
    cart3d = []
    for _ in range(n_cart3d):
        tag = reader.i32()
        name = reader.bytes(100).split(b"\0", 1)[0].decode("utf-8", errors="replace")
        component = reader.i32()
        cart3d.append({"tag": tag, "name": name, "component": component})
    header["cart3d_surfaces"] = cart3d
    return header

def _read_adb_geometry(reader: BinaryReader, header: dict) -> dict:
    ntris = header["number_of_tris"]
    nnodes = header["number_of_nodes"]

    triangles = []
    for triangle_id in range(1, ntris + 1):
        triangles.append(
            {
                "triangle_id": triangle_id,
                "node1": reader.i32(),
                "node2": reader.i32(),
                "node3": reader.i32(),
                "component_id": reader.i32(),
                "surface_id": reader.i32(),
                "min_valid_timestep": reader.i32(),
                "stored_area": reader.f32(),
            }
        )

    nodes = np.empty((nnodes + 1, 3), dtype=float)
    nodes[0] = np.nan
    for node_id in range(1, nnodes + 1):
        nodes[node_id] = (reader.f32(), reader.f32(), reader.f32())

    number_of_rotors = reader.i32()
    number_of_nozzles = reader.i32() if header["version"] == 3 else 0
    if number_of_rotors or number_of_nozzles:
        raise NotImplementedError(
            "VSPAEROMeshQuality does not yet skip rotor/nozzle binary STP records. "
            f"Found rotors={number_of_rotors}, nozzles={number_of_nozzles}."
        )

    number_of_mesh_levels = reader.i32()
    coarse_levels = []
    for level in range(1, number_of_mesh_levels + 1):
        n_nodes = reader.i32()
        n_edges = reader.i32()
        coordinates = np.array(
            [[reader.f32(), reader.f32(), reader.f32()] for _ in range(n_nodes)],
            dtype=float,
        )
        edges = []
        for _ in range(n_edges):
            signed_surface_id = reader.i32()
            min_valid_timestep = reader.i32()
            node1 = reader.i32()
            node2 = reader.i32()
            edges.append(
                {
                    "surface_id": abs(signed_surface_id),
                    "is_boundary_edge": signed_surface_id < 0,
                    "min_valid_timestep": min_valid_timestep,
                    "node1": node1,
                    "node2": node2,
                }
            )
        coarse_levels.append({"level": level, "nodes": coordinates, "edges": edges})

    number_of_kutta_edges = reader.i32()
    kutta_edges = [reader.i32() for _ in range(number_of_kutta_edges)]

    number_of_kutta_nodes = reader.i32()
    kutta_nodes = [reader.i32() for _ in range(number_of_kutta_nodes)]

    number_of_control_surfaces = reader.i32()
    controls = []
    for control_id in range(1, number_of_control_surfaces + 1):
        n_control_nodes = reader.i32()
        control_nodes = np.array(
            [[reader.f32(), reader.f32(), reader.f32()] for _ in range(n_control_nodes)],
            dtype=float,
        )
        hinge_node1 = np.array([reader.f32(), reader.f32(), reader.f32()])
        hinge_node2 = np.array([reader.f32(), reader.f32(), reader.f32()])
        hinge_vec = np.array([reader.f32(), reader.f32(), reader.f32()])
        n_loops = reader.i32()
        loop_ids = [reader.i32() for _ in range(n_loops)]
        controls.append(
            {
                "control_id": control_id,
                "nodes": control_nodes,
                "hinge_node1": hinge_node1,
                "hinge_node2": hinge_node2,
                "hinge_vec": hinge_vec,
                "loop_ids": loop_ids,
            }
        )

    return {
        "triangles": triangles,
        "nodes": nodes,
        "coarse_levels": coarse_levels,
        "kutta_edges": kutta_edges,
        "kutta_nodes": kutta_nodes,
        "controls": controls,
    }

def _polyline_metrics(points: list[tuple[float, float, float]]) -> dict:
    if not points:
        return {
            "arc_length": 0.0,
            "start_x": math.nan,
            "start_y": math.nan,
            "start_z": math.nan,
            "end_x": math.nan,
            "end_y": math.nan,
            "end_z": math.nan,
            "x_min": math.nan,
            "x_max": math.nan,
            "y_min": math.nan,
            "y_max": math.nan,
            "z_min": math.nan,
            "z_max": math.nan,
        }
    arr = np.asarray(points, dtype=float)
    arc = float(np.linalg.norm(np.diff(arr, axis=0), axis=1).sum()) if len(arr) > 1 else 0.0
    return {
        "arc_length": arc,
        "start_x": float(arr[0, 0]),
        "start_y": float(arr[0, 1]),
        "start_z": float(arr[0, 2]),
        "end_x": float(arr[-1, 0]),
        "end_y": float(arr[-1, 1]),
        "end_z": float(arr[-1, 2]),
        "x_min": float(arr[:, 0].min()),
        "x_max": float(arr[:, 0].max()),
        "y_min": float(arr[:, 1].min()),
        "y_max": float(arr[:, 1].max()),
        "z_min": float(arr[:, 2].min()),
        "z_max": float(arr[:, 2].max()),
    }

def _read_adb_solution(reader: BinaryReader, header: dict, geometry: dict) -> dict:
    mach = reader.f32()
    alpha_rad = reader.f32()
    beta_rad = reader.f32()
    cp_min_solver = reader.f32()
    cp_max_solver = reader.f32()

    nloops = header["number_of_vortex_loops"]
    computational_gamma = np.empty(nloops, dtype=float)
    computational_dcp_unsteady = np.empty(nloops, dtype=float)
    for i in range(nloops):
        computational_gamma[i] = reader.f64()
        computational_dcp_unsteady[i] = reader.f64()

    nedges = header["number_of_surface_vortex_edges"]
    edge_force = np.empty((nedges, 3), dtype=float)
    for i in range(nedges):
        edge_force[i] = (reader.f64(), reader.f64(), reader.f64())

    computational_velocity = np.empty((nloops, 3), dtype=float)
    for i in range(nloops):
        computational_velocity[i] = (reader.f64(), reader.f64(), reader.f64())

    ntris = header["number_of_tris"]
    cp = np.empty(ntris, dtype=float)
    cp_unsteady = np.empty(ntris, dtype=float)
    gamma = np.empty(ntris, dtype=float)
    for i in range(ntris):
        cp[i] = reader.f32()
        cp_unsteady[i] = reader.f32()
        gamma[i] = reader.f32()

    number_of_trailing_vortex_edges = reader.i32()
    trailing_vortex_edges = []
    for edge_id in range(1, number_of_trailing_vortex_edges + 1):
        wing_wake_node = reader.i32()
        span_location = reader.f64()
        n_sub_vortex_nodes = reader.i32()
        points = [(reader.f64(), reader.f64(), reader.f64()) for _ in range(n_sub_vortex_nodes)]
        trailing_vortex_edges.append(
            {
                "trailing_vortex_edge_id": edge_id,
                "wing_wake_node": wing_wake_node,
                "span_location": span_location,
                "n_sub_vortex_nodes": n_sub_vortex_nodes,
                "points": points,
                **_polyline_metrics(points),
            }
        )

    control_deflections = [reader.f32() for _ in geometry["controls"]]

    return {
        "mach": mach,
        "alpha_deg": math.degrees(alpha_rad),
        "beta_deg": math.degrees(beta_rad),
        "cp_min_solver": cp_min_solver,
        "cp_max_solver": cp_max_solver,
        "computational_gamma": computational_gamma,
        "computational_dcp_unsteady": computational_dcp_unsteady,
        "edge_force": edge_force,
        "computational_velocity": computational_velocity,
        "cp": cp,
        "cp_unsteady": cp_unsteady,
        "cp_steady": cp - cp_unsteady,
        "gamma": gamma,
        "trailing_vortex_edges": trailing_vortex_edges,
        "control_deflections": control_deflections,
    }

def read_adb_v3(path: str | Path, solution_case: int = 1) -> dict:
    """Read one solution case from an OpenVSP/VSPAERO ADB v3 file."""
    path = Path(path)
    if solution_case < 1:
        raise ValueError("solution_case must be >= 1.")

    with path.open("rb") as fp:
        endian, version = _detect_adb_endian_and_version(fp)
        if version != 3:
            raise ValueError(f"Only ADB v3 is supported; found v{version}.")
        reader = BinaryReader(fp, endian)
        header = _read_adb_header(reader, version)

        selected_geometry = None
        selected_solution = None
        for case_id in range(1, solution_case + 1):
            try:
                geometry = _read_adb_geometry(reader, header)
                solution = _read_adb_solution(reader, header, geometry)
            except EOFError as exc:
                raise ValueError(
                    f"solution_case={solution_case} does not exist or the ADB is truncated."
                ) from exc
            if case_id == solution_case:
                selected_geometry = geometry
                selected_solution = solution

    return {
        "path": str(path),
        "endian": "little" if endian == "<" else "big",
        "header": header,
        "geometry": selected_geometry,
        "solution": selected_solution,
        "solution_case": solution_case,
    }

def read_vspgeom_v3(path: str | Path) -> dict:
    """Read VSPGeom v3 NGons, UV data, Kutta lists and alternate triangulation."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines or lines[0].strip() != "# vspgeom v3":
        raise ValueError("Only '# vspgeom v3' files are supported.")

    tokens = iter(" ".join(lines[1:]).split())

    def next_int():
        return int(next(tokens))

    def next_float():
        return float(next(tokens))

    number_of_refine_levels = next_int()
    num_nodes = next_int()
    num_loops_header = next_int()
    wake_nodes = next_int()

    nodes = np.empty((num_nodes + 1, 3), dtype=float)
    nodes[0] = np.nan
    for node_id in range(1, num_nodes + 1):
        nodes[node_id] = (next_float(), next_float(), next_float())

    num_loops = next_int()
    if num_loops != num_loops_header:
        raise ValueError(
            f"VSPGeom loop count mismatch: header={num_loops_header}, body={num_loops}."
        )

    ngons = []
    for ngon_id in range(1, num_loops + 1):
        n = next_int()
        node_ids = [next_int() for _ in range(n)]
        ngons.append({"ngon_id": ngon_id, "node_ids": node_ids, "n_vertices": n})

    node_surface_uv = defaultdict(lambda: defaultdict(list))
    for ngon in ngons:
        surface_id = next_int()
        subsurface_id = next_int()
        uvs = []
        for node_id in ngon["node_ids"]:
            u = next_float()
            w = next_float()
            uvs.append((u, w))
            node_surface_uv[surface_id][node_id].append((u, w))
        ngon["vspgeom_surface_id"] = surface_id
        ngon["vspgeom_subsurface_id"] = subsurface_id
        ngon["uvs"] = uvs

    for ngon in ngons:
        ngon["parent_a"] = next_int()
        ngon["parent_b"] = next_int()

    number_of_kutta_lists = next_int()
    kutta_lists = []
    for kutta_list_id in range(1, number_of_kutta_lists + 1):
        signed_count = next_int()
        wake_part_num = next_int()
        node_ids = [next_int() for _ in range(abs(signed_count))]
        kutta_lists.append(
            {
                "kutta_list_id": kutta_list_id,
                "body_wake": signed_count < 0,
                "wake_part_num": wake_part_num,
                "node_ids": node_ids,
            }
        )

    child_to_parent = {}
    triangle_connectivity = {}
    child_triangle_id = 0
    for expected_ngon_id, ngon in enumerate(ngons, start=1):
        record_ngon_id = next_int()
        ntris = next_int()
        ngon["alternate_record_id"] = record_ngon_id
        child_ids = []
        for _ in range(ntris):
            child_triangle_id += 1
            tri_nodes = (next_int(), next_int(), next_int())
            child_to_parent[child_triangle_id] = expected_ngon_id
            triangle_connectivity[child_triangle_id] = tri_nodes
            child_ids.append(child_triangle_id)
        ngon["child_triangle_ids"] = child_ids
        ngon["n_child_triangles"] = ntris

    return {
        "path": str(path),
        "number_of_refine_levels": number_of_refine_levels,
        "num_nodes": num_nodes,
        "num_loops": num_loops,
        "wake_nodes": wake_nodes,
        "nodes": nodes,
        "ngons": ngons,
        "kutta_lists": kutta_lists,
        "node_surface_uv": node_surface_uv,
        "child_to_parent": child_to_parent,
        "triangle_connectivity": triangle_connectivity,
        "num_alternate_triangles": child_triangle_id,
    }

def _triangle_geometry(node_ids: tuple[int, int, int], nodes: np.ndarray) -> dict:
    p1, p2, p3 = (nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]])
    e12 = float(np.linalg.norm(p2 - p1))
    e23 = float(np.linalg.norm(p3 - p2))
    e31 = float(np.linalg.norm(p1 - p3))
    edges = np.array([e12, e23, e31], dtype=float)
    edge_min = float(edges.min())
    edge_max = float(edges.max())
    edge_ratio = edge_max / edge_min if edge_min > 0 else math.inf

    cross = np.cross(p2 - p1, p3 - p1)
    cross_mag = float(np.linalg.norm(cross))
    area = 0.5 * cross_mag
    normal = cross / cross_mag if cross_mag > 0 else np.array([math.nan] * 3)

    angles = []
    opposite = ((e23, e12, e31), (e31, e12, e23), (e12, e23, e31))
    for a, b, c in opposite:
        if b <= 0 or c <= 0:
            angles.append(math.nan)
        else:
            cosine = max(-1.0, min(1.0, (b * b + c * c - a * a) / (2 * b * c)))
            angles.append(math.degrees(math.acos(cosine)))

    centroid = (p1 + p2 + p3) / 3.0
    finite_angles = [a for a in angles if math.isfinite(a)]
    return {
        "centroid_x": float(centroid[0]),
        "centroid_y": float(centroid[1]),
        "centroid_z": float(centroid[2]),
        "area": area,
        "edge_min": edge_min,
        "edge_max": edge_max,
        "edge_ratio": edge_ratio,
        "angle_min_deg": min(finite_angles) if finite_angles else math.nan,
        "angle_max_deg": max(finite_angles) if finite_angles else math.nan,
        "normal_x": float(normal[0]),
        "normal_y": float(normal[1]),
        "normal_z": float(normal[2]),
    }

def _surface_metadata_from_adb(adb: dict) -> tuple[dict[int, dict], list[dict]]:
    rows = []
    by_surface = {}
    for record in adb["header"].get("cart3d_surfaces", []):
        surface_id = int(record["tag"])
        row = {
            "surface_id": surface_id,
            "surface_name": record["name"],
            "component_id": int(record["component"]),
        }
        rows.append(row)
        by_surface[surface_id] = row
    return by_surface, rows

def _build_component_metadata(surface_rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in surface_rows:
        grouped[row["component_id"]].append(row)
    out = []
    for component_id in sorted(grouped):
        items = sorted(grouped[component_id], key=lambda x: x["surface_id"])
        out.append(
            {
                "component_id": component_id,
                "n_surfaces": len(items),
                "surface_ids": [x["surface_id"] for x in items],
                "surface_names": [x["surface_name"] for x in items],
            }
        )
    return out

def _build_triangle_rows(adb: dict, surface_meta: dict[int, dict]) -> tuple[list[dict], dict[int, set[int]], dict]:
    geometry = adb["geometry"]
    solution = adb["solution"]
    nodes = geometry["nodes"]

    rows = []
    edge_to_triangles = defaultdict(list)
    for tri, cp, cp_unsteady, gamma in zip(
        geometry["triangles"], solution["cp"], solution["cp_unsteady"], solution["gamma"]
    ):
        node_ids = (tri["node1"], tri["node2"], tri["node3"])
        metrics = _triangle_geometry(node_ids, nodes)
        area_error = metrics["area"] - tri["stored_area"]
        area_rel_error = area_error / tri["stored_area"] if abs(tri["stored_area"]) > 0 else math.nan
        meta = surface_meta.get(tri["surface_id"], {})
        row = {
            **tri,
            **metrics,
            "surface_name": meta.get("surface_name", ""),
            "surface_component_id_metadata": meta.get("component_id"),
            "area_error": area_error,
            "area_rel_error": area_rel_error,
            "cp": float(cp),
            "cp_unsteady": float(cp_unsteady),
            "cp_steady": float(cp - cp_unsteady),
            "gamma": float(gamma),
            "parent_ngon_id": None,
        }
        rows.append(row)
        for a, b in ((node_ids[0], node_ids[1]), (node_ids[1], node_ids[2]), (node_ids[2], node_ids[0])):
            edge_to_triangles[tuple(sorted((a, b)))].append(tri["triangle_id"])

    by_id = {row["triangle_id"]: row for row in rows}
    neighbors = {row["triangle_id"]: set() for row in rows}
    boundary_edges = []
    edge_multiplicity_gt2 = []
    same_surface_non_manifold_edges = []
    cross_surface_junction_edges = []

    for edge, tri_ids in edge_to_triangles.items():
        if len(tri_ids) == 1:
            boundary_edges.append(edge)
        elif len(tri_ids) > 2:
            surface_keys = sorted(
                {(by_id[t]["component_id"], by_id[t]["surface_id"]) for t in tri_ids}
            )
            edge_info = {
                "edge": edge,
                "triangle_ids": list(tri_ids),
                "component_surface_pairs": surface_keys,
            }
            edge_multiplicity_gt2.append(edge_info)
            if len(surface_keys) == 1:
                same_surface_non_manifold_edges.append(edge_info)
            else:
                cross_surface_junction_edges.append(edge_info)
        for tri_id in tri_ids:
            neighbors[tri_id].update(other for other in tri_ids if other != tri_id)

    for row in rows:
        ratios = []
        normal_angles = []
        n1 = np.array([row["normal_x"], row["normal_y"], row["normal_z"]])
        for neighbor_id in neighbors[row["triangle_id"]]:
            other = by_id[neighbor_id]
            if row["area"] > 0 and other["area"] > 0:
                ratios.append(max(row["area"] / other["area"], other["area"] / row["area"]))
            n2 = np.array([other["normal_x"], other["normal_y"], other["normal_z"]])
            if np.all(np.isfinite(n1)) and np.all(np.isfinite(n2)):
                normal_angles.append(math.degrees(math.acos(float(np.clip(np.dot(n1, n2), -1.0, 1.0)))))
        row["neighbor_area_ratio"] = max(ratios) if ratios else math.nan
        row["max_neighbor_normal_angle_deg"] = max(normal_angles) if normal_angles else math.nan
        row["neighbor_count"] = len(neighbors[row["triangle_id"]])

    topology = {
        "boundary_edge_count": len(boundary_edges),
        "edge_multiplicity_gt2": edge_multiplicity_gt2,
        "same_surface_non_manifold_edges": same_surface_non_manifold_edges,
        "cross_surface_junction_edges": cross_surface_junction_edges,
    }
    return rows, neighbors, topology

def _validate_surface_metadata(triangle_rows: list[dict], surface_meta: dict[int, dict]) -> dict:
    surface_ids = sorted({r["surface_id"] for r in triangle_rows if r["surface_id"] > 0})
    missing = [sid for sid in surface_ids if sid not in surface_meta]
    mismatches = []
    for row in triangle_rows:
        if row["surface_id"] <= 0 or row["surface_id"] not in surface_meta:
            continue
        expected = surface_meta[row["surface_id"]]["component_id"]
        if row["component_id"] != expected:
            item = (row["surface_id"], row["component_id"], expected)
            if item not in mismatches:
                mismatches.append(item)
    return {
        "missing_surface_ids": missing,
        "component_mismatches": mismatches,
        "passed": not missing and not mismatches,
    }

def _attach_parent_ngons(triangle_rows: list[dict], vspgeom: dict) -> dict:
    by_id = {row["triangle_id"]: row for row in triangle_rows}
    mismatches = []
    for child_id, parent_id in vspgeom["child_to_parent"].items():
        if child_id not in by_id:
            mismatches.append(f"VSPGeom child triangle {child_id} is missing from ADB.")
            continue
        expected_nodes = tuple(sorted(vspgeom["triangle_connectivity"][child_id]))
        actual_nodes = tuple(sorted((by_id[child_id]["node1"], by_id[child_id]["node2"], by_id[child_id]["node3"])))
        if expected_nodes != actual_nodes:
            mismatches.append(
                f"Triangle {child_id}: VSPGeom connectivity={expected_nodes}, ADB={actual_nodes}."
            )
            continue
        by_id[child_id]["parent_ngon_id"] = parent_id
    return {
        "mapped_triangles": sum(row["parent_ngon_id"] is not None for row in triangle_rows),
        "vspgeom_alternate_triangles": vspgeom["num_alternate_triangles"],
        "mismatches": mismatches,
        "passed": not mismatches,
    }

def _normal_angle(row1: dict, row2: dict) -> float:
    n1 = np.array([row1["normal_x"], row1["normal_y"], row1["normal_z"]])
    n2 = np.array([row2["normal_x"], row2["normal_y"], row2["normal_z"]])
    if not np.all(np.isfinite(n1)) or not np.all(np.isfinite(n2)):
        return math.inf
    return math.degrees(math.acos(float(np.clip(np.dot(n1, n2), -1.0, 1.0))))

def _cp_outlier_metrics(cp: float, neighbor_cp: list[float], spike_z_threshold: float, spike_delta_threshold: float) -> dict:
    if not neighbor_cp:
        return {
            "cp_neighbor_median": math.nan,
            "cp_local_delta": math.nan,
            "cp_local_mad": math.nan,
            "cp_robust_z": math.nan,
            "cp_robust_z_status": "no_neighbors",
            "cp_zero_mad_outlier": False,
            "is_negative_cp_spike": False,
        }
    cp_median = float(median(neighbor_cp))
    mad = float(median([abs(value - cp_median) for value in neighbor_cp]))
    delta = cp - cp_median
    scale = 1.4826 * mad
    if scale > 0:
        robust_z = delta / scale
        status = "finite"
        zero_mad_outlier = False
        is_negative_spike = robust_z <= spike_z_threshold and delta <= spike_delta_threshold
    elif delta == 0:
        robust_z = 0.0
        status = "zero_mad_same"
        zero_mad_outlier = False
        is_negative_spike = False
    else:
        robust_z = math.nan
        status = "zero_mad_negative_delta" if delta < 0 else "zero_mad_positive_delta"
        zero_mad_outlier = True
        is_negative_spike = delta <= spike_delta_threshold
    return {
        "cp_neighbor_median": cp_median,
        "cp_local_delta": delta,
        "cp_local_mad": mad,
        "cp_robust_z": robust_z,
        "cp_robust_z_status": status,
        "cp_zero_mad_outlier": zero_mad_outlier,
        "is_negative_cp_spike": is_negative_spike,
    }

def _add_local_cp_metrics(rows: list[dict], neighbors: dict[int, set[int]], normal_angle_limit_deg: float, spike_z_threshold: float, spike_delta_threshold: float):
    by_id = {row["triangle_id"]: row for row in rows}
    for row in rows:
        cp_neighbors = []
        for neighbor_id in neighbors.get(row["triangle_id"], set()):
            other = by_id.get(neighbor_id)
            if other is None:
                continue
            if other["component_id"] != row["component_id"] or other["surface_id"] != row["surface_id"]:
                continue
            if _normal_angle(row, other) > normal_angle_limit_deg:
                continue
            cp_neighbors.append(other["cp"])
        row.update(_cp_outlier_metrics(row["cp"], cp_neighbors, spike_z_threshold, spike_delta_threshold))

def _mark_cp_floor_plateau(rows: list[dict], tolerance: float, min_cells: int, max_cp: float) -> dict:
    for row in rows:
        row["cp_floor_plateau_size"] = 0
        row["cp_floor_plateau_fraction"] = 0.0
        row["is_cp_floor_plateau"] = False
        row["is_cp_limiter_candidate"] = False
        row["cpcrit_equivalent_local_mach"] = math.nan
    if not rows:
        return {
            "cp_floor_value": None,
            "cp_floor_plateau_size": 0,
            "cp_floor_plateau_fraction": 0.0,
            "cp_floor_plateau_detected": False,
            "cp_floor_cpcrit_equivalent_local_mach": None,
        }
    cp_floor = min(row["cp"] for row in rows)
    members = [row for row in rows if abs(row["cp"] - cp_floor) <= tolerance]
    fraction = len(members) / len(rows)
    detected = cp_floor <= max_cp and len(members) >= min_cells
    k = 5.0 / 2.4
    equivalent_mach = math.sqrt(k / (k - cp_floor)) if cp_floor < 0 and k - cp_floor > 0 else math.nan
    for row in members:
        row["cp_floor_plateau_size"] = len(members)
        row["cp_floor_plateau_fraction"] = fraction
        row["is_cp_floor_plateau"] = detected
        row["is_cp_limiter_candidate"] = detected
        row["cpcrit_equivalent_local_mach"] = equivalent_mach
    return {
        "cp_floor_value": cp_floor,
        "cp_floor_plateau_size": len(members),
        "cp_floor_plateau_fraction": fraction,
        "cp_floor_plateau_detected": detected,
        "cp_floor_cpcrit_equivalent_local_mach": equivalent_mach,
    }

def _build_ngon_rows(triangle_rows: list[dict], vspgeom: dict, normal_angle_limit_deg: float, spike_z_threshold: float, spike_delta_threshold: float) -> list[dict]:
    tri_by_id = {row["triangle_id"]: row for row in triangle_rows}
    rows = []
    edge_to_ngons = defaultdict(list)
    for ngon in vspgeom["ngons"]:
        node_ids = ngon["node_ids"]
        for i in range(len(node_ids)):
            edge_to_ngons[tuple(sorted((node_ids[i], node_ids[(i + 1) % len(node_ids)])))].append(ngon["ngon_id"])
    neighbors = {ngon["ngon_id"]: set() for ngon in vspgeom["ngons"]}
    for ngon_ids in edge_to_ngons.values():
        for ngon_id in ngon_ids:
            neighbors[ngon_id].update(other for other in ngon_ids if other != ngon_id)

    for ngon in vspgeom["ngons"]:
        child_rows = [tri_by_id[c] for c in ngon["child_triangle_ids"] if c in tri_by_id]
        if not child_rows:
            continue
        area = sum(row["area"] for row in child_rows)
        cp_values = [row["cp"] for row in child_rows]
        weighted_normal = np.zeros(3)
        weighted_centroid = np.zeros(3)
        for row in child_rows:
            weighted_normal += row["area"] * np.array([row["normal_x"], row["normal_y"], row["normal_z"]])
            weighted_centroid += row["area"] * np.array([row["centroid_x"], row["centroid_y"], row["centroid_z"]])
        norm = np.linalg.norm(weighted_normal)
        weighted_normal = weighted_normal / norm if norm > 0 else np.array([math.nan] * 3)
        centroid = weighted_centroid / area if area > 0 else np.array([math.nan] * 3)
        child_min_angle = min(row["angle_min_deg"] for row in child_rows)
        rows.append(
            {
                "ngon_id": ngon["ngon_id"],
                "component_id": child_rows[0]["component_id"],
                "surface_id": child_rows[0]["surface_id"],
                "surface_name": child_rows[0].get("surface_name", ""),
                "n_vertices": ngon["n_vertices"],
                "n_child_triangles": len(child_rows),
                "centroid_x": float(centroid[0]),
                "centroid_y": float(centroid[1]),
                "centroid_z": float(centroid[2]),
                "area": area,
                "cp": float(median(cp_values)),
                "child_cp_range": max(cp_values) - min(cp_values),
                "child_min_angle_deg": child_min_angle,
                "child_max_edge_ratio": max(row["edge_ratio"] for row in child_rows),
                "child_max_neighbor_area_ratio": max(
                    (row["neighbor_area_ratio"] for row in child_rows if math.isfinite(row["neighbor_area_ratio"])),
                    default=math.nan,
                ),
                "normal_x": float(weighted_normal[0]),
                "normal_y": float(weighted_normal[1]),
                "normal_z": float(weighted_normal[2]),
            }
        )

    by_id = {row["ngon_id"]: row for row in rows}
    for row in rows:
        cp_neighbors = []
        for neighbor_id in neighbors[row["ngon_id"]]:
            other = by_id.get(neighbor_id)
            if other is None:
                continue
            if other["component_id"] != row["component_id"] or other["surface_id"] != row["surface_id"]:
                continue
            if _normal_angle(row, other) > normal_angle_limit_deg:
                continue
            cp_neighbors.append(other["cp"])
        row.update(_cp_outlier_metrics(row["cp"], cp_neighbors, spike_z_threshold, spike_delta_threshold))
    return rows

def _matches_filters(row: dict, component_ids: set[int] | None, surface_ids: set[int] | None, bbox: tuple[float, float, float, float, float, float] | None) -> bool:
    if component_ids is not None and row.get("component_id") not in component_ids:
        return False
    if surface_ids is not None and row.get("surface_id") not in surface_ids:
        return False
    if bbox is not None and "centroid_x" in row:
        xmin, xmax, ymin, ymax, zmin, zmax = bbox
        if not (xmin <= row["centroid_x"] <= xmax and ymin <= row["centroid_y"] <= ymax and zmin <= row["centroid_z"] <= zmax):
            return False
    return True

def _surface_geometry_from_vspgeom(vspgeom: dict, surface_meta: dict[int, dict]) -> dict[int, dict]:
    out = {}
    by_surface_ngons = defaultdict(list)
    for ngon in vspgeom["ngons"]:
        by_surface_ngons[ngon["vspgeom_surface_id"]].append(ngon)

    for surface_id, ngons in by_surface_ngons.items():
        node_ids = sorted({node_id for ngon in ngons for node_id in ngon["node_ids"]})
        coords = vspgeom["nodes"][node_ids]
        u_values = []
        w_values = []
        edges = {}
        for ngon in ngons:
            n = len(ngon["node_ids"])
            for i, node_id in enumerate(ngon["node_ids"]):
                u, w = ngon["uvs"][i]
                u_values.append(u)
                w_values.append(w)
                j = (i + 1) % n
                node2 = ngon["node_ids"][j]
                u2, w2 = ngon["uvs"][j]
                edge = tuple(sorted((node_id, node2)))
                if edge not in edges:
                    edges[edge] = {
                        "node1": node_id,
                        "node2": node2,
                        "ngon_id": ngon["ngon_id"],
                        "u1": u,
                        "w1": w,
                        "u2": u2,
                        "w2": w2,
                    }
        y_range = float(coords[:, 1].max() - coords[:, 1].min()) if len(coords) else 0.0
        z_range = float(coords[:, 2].max() - coords[:, 2].min()) if len(coords) else 0.0
        span_axis = "y" if y_range >= z_range else "z"
        span_values = coords[:, 1] if span_axis == "y" else coords[:, 2]
        meta = surface_meta.get(surface_id, {})
        out[surface_id] = {
            "surface_id": surface_id,
            "surface_name": meta.get("surface_name", ""),
            "component_id": meta.get("component_id"),
            "node_ids": node_ids,
            "node_id_set": set(node_ids),
            "n_nodes": len(node_ids),
            "x_min": float(coords[:, 0].min()),
            "x_max": float(coords[:, 0].max()),
            "y_min": float(coords[:, 1].min()),
            "y_max": float(coords[:, 1].max()),
            "z_min": float(coords[:, 2].min()),
            "z_max": float(coords[:, 2].max()),
            "span_axis": span_axis,
            "span_min": float(span_values.min()),
            "span_max": float(span_values.max()),
            "span_outer_abs": float(np.max(np.abs(span_values))),
            "param_u_min": min(u_values),
            "param_u_max": max(u_values),
            "param_w_min": min(w_values),
            "param_w_max": max(w_values),
            "edges": list(edges.values()),
        }
    return out

def _node_uv_summary(vspgeom: dict, surface_id: int, node_id: int) -> tuple[float, float, float, float]:
    vals = vspgeom["node_surface_uv"].get(surface_id, {}).get(node_id, [])
    if not vals:
        return math.nan, math.nan, math.nan, math.nan
    u = [x[0] for x in vals]
    w = [x[1] for x in vals]
    return min(u), max(u), min(w), max(w)

def _infer_kutta_surface(node_ids: list[int], surfaces: dict[int, dict]) -> tuple[int | None, float]:
    best_id = None
    best_fraction = 0.0
    if not node_ids:
        return None, 0.0
    nodes = set(node_ids)
    for surface_id, surface in surfaces.items():
        fraction = len(nodes & surface["node_id_set"]) / len(nodes)
        if fraction > best_fraction:
            best_id = surface_id
            best_fraction = fraction
    return best_id, best_fraction

def _relative_difference(a: float | None, b: float | None) -> float | None:
    if a is None or b is None or not math.isfinite(a) or not math.isfinite(b):
        return None
    scale = max(abs(a), abs(b), 1e-15)
    return abs(a - b) / scale

def _analyze_kutta(vspgeom: dict, surface_meta: dict[int, dict], coverage_min: float, symmetry_tolerance_fraction: float, te_param_near_miss_band: float) -> dict:
    surfaces = _surface_geometry_from_vspgeom(vspgeom, surface_meta)
    kutta_lines = []
    kutta_nodes = []
    kutta_node_set = set()

    for kutta in vspgeom["kutta_lists"]:
        node_ids = kutta["node_ids"]
        kutta_node_set.update(node_ids)
        surface_id, inference_fraction = _infer_kutta_surface(node_ids, surfaces)
        surface = surfaces.get(surface_id) if surface_id is not None else None
        points = [tuple(vspgeom["nodes"][node_id]) for node_id in node_ids]
        metrics = _polyline_metrics(points)
        span_axis = surface["span_axis"] if surface else ""
        if span_axis == "y":
            span_values = [p[1] for p in points]
        elif span_axis == "z":
            span_values = [p[2] for p in points]
        else:
            span_values = []
        span_min = min(span_values) if span_values else math.nan
        span_max = max(span_values) if span_values else math.nan
        outer_abs_span = max((abs(v) for v in span_values), default=math.nan)
        surface_outer = surface["span_outer_abs"] if surface else math.nan
        coverage = outer_abs_span / surface_outer if surface and surface_outer > 0 else math.nan
        side = ""
        if span_values:
            mean_span = sum(span_values) / len(span_values)
            side = "positive" if mean_span > 0 else "negative" if mean_span < 0 else "center"

        max_w_offset = math.nan
        if surface:
            offsets = []
            for node_id in node_ids:
                _, _, wmin, wmax = _node_uv_summary(vspgeom, surface_id, node_id)
                if math.isfinite(wmin):
                    offsets.append(max(abs(wmin - surface["param_w_min"]), abs(wmax - surface["param_w_min"])))
            max_w_offset = max(offsets) if offsets else math.nan

        kutta_lines.append(
            {
                "kutta_list_id": kutta["kutta_list_id"],
                "body_wake": kutta["body_wake"],
                "wake_part_num": kutta["wake_part_num"],
                "n_nodes": len(node_ids),
                "surface_id": surface_id,
                "surface_name": surface["surface_name"] if surface else "",
                "component_id": surface["component_id"] if surface else None,
                "surface_inference_fraction": inference_fraction,
                "span_axis": span_axis,
                "side": side,
                **metrics,
                "span_min": span_min,
                "span_max": span_max,
                "outer_abs_span": outer_abs_span,
                "surface_outer_abs_span": surface_outer,
                "outer_coverage_fraction": coverage,
                "kutta_param_w_offset_max": max_w_offset,
            }
        )

        for sequence_index, node_id in enumerate(node_ids, start=1):
            xyz = vspgeom["nodes"][node_id]
            umin = umax = wmin = wmax = math.nan
            surface_wmin = math.nan
            offset = math.nan
            if surface:
                umin, umax, wmin, wmax = _node_uv_summary(vspgeom, surface_id, node_id)
                surface_wmin = surface["param_w_min"]
                if math.isfinite(wmin):
                    offset = min(abs(wmin - surface_wmin), abs(wmax - surface_wmin))
            kutta_nodes.append(
                {
                    "kutta_list_id": kutta["kutta_list_id"],
                    "sequence_index": sequence_index,
                    "node_id": node_id,
                    "wake_part_num": kutta["wake_part_num"],
                    "surface_id": surface_id,
                    "surface_name": surface["surface_name"] if surface else "",
                    "component_id": surface["component_id"] if surface else None,
                    "span_axis": span_axis,
                    "side": side,
                    "x": float(xyz[0]),
                    "y": float(xyz[1]),
                    "z": float(xyz[2]),
                    "param_u_min": umin,
                    "param_u_max": umax,
                    "param_w_min": wmin,
                    "param_w_max": wmax,
                    "surface_param_w_min": surface_wmin,
                    "param_w_offset_min": offset,
                }
            )

    by_surface_lines = defaultdict(list)
    for row in kutta_lines:
        if row["surface_id"] is not None and not row["body_wake"]:
            by_surface_lines[row["surface_id"]].append(row)

    surface_rows = []
    near_miss_edges = []
    for surface_id, surface in surfaces.items():
        strict_count = 0
        near_count = 0
        near_miss_count = 0
        for edge in surface["edges"]:
            wmin = surface["param_w_min"]
            off1 = edge["w1"] - wmin
            off2 = edge["w2"] - wmin
            strict = off1 <= WAKE_EDGE_TOLERANCE and off2 <= WAKE_EDGE_TOLERANCE
            near = off1 <= te_param_near_miss_band and off2 <= te_param_near_miss_band
            if strict:
                strict_count += 1
            if near:
                near_count += 1
            if near and not strict:
                near_miss_count += 1
                p1 = vspgeom["nodes"][edge["node1"]]
                p2 = vspgeom["nodes"][edge["node2"]]
                mid = (p1 + p2) / 2
                near_miss_edges.append(
                    {
                        "surface_id": surface_id,
                        "surface_name": surface["surface_name"],
                        "node1": edge["node1"],
                        "node2": edge["node2"],
                        "ngon_id": edge["ngon_id"],
                        "midpoint_x": float(mid[0]),
                        "midpoint_y": float(mid[1]),
                        "midpoint_z": float(mid[2]),
                        "edge_length": float(np.linalg.norm(p2 - p1)),
                        "surface_param_w_min": wmin,
                        "param_w1": edge["w1"],
                        "param_w2": edge["w2"],
                        "param_w_offset1": off1,
                        "param_w_offset2": off2,
                        "param_w_offset_max": max(off1, off2),
                        "wake_edge_tolerance": WAKE_EDGE_TOLERANCE,
                        "te_param_near_miss_band": te_param_near_miss_band,
                        "node1_is_kutta_node": edge["node1"] in kutta_node_set,
                        "node2_is_kutta_node": edge["node2"] in kutta_node_set,
                        "is_kutta_node_pair": edge["node1"] in kutta_node_set and edge["node2"] in kutta_node_set,
                    }
                )

        lines = by_surface_lines.get(surface_id, [])
        positive = [x for x in lines if x["side"] == "positive"]
        negative = [x for x in lines if x["side"] == "negative"]
        n_nodes = sum(x["n_nodes"] for x in lines)
        outer = max((x["outer_abs_span"] for x in lines), default=math.nan)
        coverage = outer / surface["span_outer_abs"] if surface["span_outer_abs"] > 0 and math.isfinite(outer) else math.nan
        pos_outer = max((x["outer_abs_span"] for x in positive), default=math.nan)
        neg_outer = max((x["outer_abs_span"] for x in negative), default=math.nan)
        pos_cov = pos_outer / surface["span_outer_abs"] if math.isfinite(pos_outer) else math.nan
        neg_cov = neg_outer / surface["span_outer_abs"] if math.isfinite(neg_outer) else math.nan
        coverage_passed = True if not lines else (math.isfinite(coverage) and coverage >= coverage_min)
        surface_rows.append(
            {
                **{k: v for k, v in surface.items() if k not in {"node_ids", "node_id_set", "edges"}},
                "strict_parametric_te_edge_count": strict_count,
                "near_parametric_te_edge_count": near_count,
                "te_param_near_miss_edge_count": near_miss_count,
                "n_kutta_lists": len(lines),
                "n_kutta_nodes": n_nodes,
                "kutta_outer_abs_span": outer,
                "kutta_outer_coverage_fraction": coverage,
                "positive_kutta_outer_span": pos_outer,
                "negative_kutta_outer_span_abs": neg_outer,
                "positive_kutta_coverage_fraction": pos_cov,
                "negative_kutta_coverage_fraction": neg_cov,
                "is_y_symmetric_geometry": False,
                "geometry_outer_extent_relative_difference": None,
                "kutta_coverage_passed": coverage_passed,
            }
        )

    surface_row_by_id = {r["surface_id"]: r for r in surface_rows}
    for row in near_miss_edges:
        srow = surface_row_by_id.get(row["surface_id"])
        row["surface_kutta_coverage_passed"] = srow["kutta_coverage_passed"] if srow else None

    symmetry_rows = []
    grouped = defaultdict(list)
    for row in surface_rows:
        grouped[(row["component_id"], row["surface_name"])].append(row)
    for (component_id, surface_name), items in grouped.items():
        if len(items) != 2:
            continue
        a, b = items
        if a["span_axis"] != "y" or b["span_axis"] != "y":
            continue
        positive = a if a["y_min"] >= 0 else b if b["y_min"] >= 0 else None
        negative = a if a["y_max"] <= 0 else b if b["y_max"] <= 0 else None
        if positive is None or negative is None:
            continue
        outer_diff = _relative_difference(positive["kutta_outer_abs_span"], negative["kutta_outer_abs_span"])
        node_diff = _relative_difference(float(positive["n_kutta_nodes"]), float(negative["n_kutta_nodes"]))
        geometry_diff = _relative_difference(positive["span_outer_abs"], negative["span_outer_abs"])
        cov_diff = _relative_difference(positive["kutta_outer_coverage_fraction"], negative["kutta_outer_coverage_fraction"])
        passed = all(
            d is not None and d <= symmetry_tolerance_fraction
            for d in (outer_diff, node_diff, geometry_diff, cov_diff)
        )
        symmetry_rows.append(
            {
                "component_id": component_id,
                "surface_name": surface_name,
                "symmetry_mode": "mirrored_surface_pair",
                "positive_surface_id": positive["surface_id"],
                "negative_surface_id": negative["surface_id"],
                "positive_kutta_node_count": positive["n_kutta_nodes"],
                "negative_kutta_node_count": negative["n_kutta_nodes"],
                "positive_outer_span": positive["kutta_outer_abs_span"],
                "negative_outer_span_abs": negative["kutta_outer_abs_span"],
                "positive_coverage_fraction": positive["kutta_outer_coverage_fraction"],
                "negative_coverage_fraction": negative["kutta_outer_coverage_fraction"],
                "geometry_outer_extent_relative_difference": geometry_diff,
                "outer_extent_relative_difference": outer_diff,
                "node_count_relative_difference": node_diff,
                "coverage_relative_difference": cov_diff,
                "tolerance_fraction": symmetry_tolerance_fraction,
                "passed": passed,
            }
        )

    unresolved = [x for x in kutta_lines if not x["body_wake"] and (x["surface_id"] is None or x["surface_inference_fraction"] < 0.5)]
    coverage_failures = [x for x in surface_rows if x["n_kutta_lists"] > 0 and not x["kutta_coverage_passed"]]
    symmetry_failures = [x for x in symmetry_rows if not x["passed"]]
    return {
        "lines": kutta_lines,
        "nodes": kutta_nodes,
        "surfaces": surface_rows,
        "symmetry": symmetry_rows,
        "near_miss_edges": near_miss_edges,
        "coverage_failures": coverage_failures,
        "symmetry_failures": symmetry_failures,
        "unresolved": unresolved,
        "passed": not coverage_failures and not symmetry_failures and not unresolved,
    }

def _point_segment_distance(point: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    ab = b - a
    denom = float(np.dot(ab, ab))
    if denom <= 0:
        return float(np.linalg.norm(point - a))
    t = float(np.clip(np.dot(point - a, ab) / denom, 0.0, 1.0))
    return float(np.linalg.norm(point - (a + t * ab)))

def _attach_nearest_kutta_distance(rows: list[dict], kutta_result: dict, vspgeom: dict):
    segments = []
    endpoints = []
    for line in kutta_result["lines"]:
        if line["body_wake"]:
            continue
        node_ids = next((k["node_ids"] for k in vspgeom["kutta_lists"] if k["kutta_list_id"] == line["kutta_list_id"]), [])
        pts = [vspgeom["nodes"][n] for n in node_ids]
        for a, b in zip(pts, pts[1:]):
            segments.append((line["kutta_list_id"], np.asarray(a), np.asarray(b)))
        for p in pts[:1] + pts[-1:]:
            endpoints.append((line["kutta_list_id"], np.asarray(p)))
    for row in rows:
        p = np.array([row["centroid_x"], row["centroid_y"], row["centroid_z"]], dtype=float)
        nearest_distance = math.inf
        nearest_list = None
        for list_id, a, b in segments:
            d = _point_segment_distance(p, a, b)
            if d < nearest_distance:
                nearest_distance = d
                nearest_list = list_id
        endpoint_distance = min((float(np.linalg.norm(p - q)) for _, q in endpoints), default=math.inf)
        row["nearest_kutta_distance"] = nearest_distance if math.isfinite(nearest_distance) else math.nan
        row["nearest_kutta_endpoint_distance"] = endpoint_distance if math.isfinite(endpoint_distance) else math.nan
        row["nearest_kutta_list_id"] = nearest_list

def _build_junction_quality_rows(topology: dict, triangle_rows: list[dict], nodes: np.ndarray) -> list[dict]:
    by_id = {r["triangle_id"]: r for r in triangle_rows}
    rows = []
    for junction_id, item in enumerate(topology["cross_surface_junction_edges"], start=1):
        node1, node2 = item["edge"]
        p1, p2 = nodes[node1], nodes[node2]
        incident = [by_id[t] for t in item["triangle_ids"] if t in by_id]
        rows.append(
            {
                "junction_edge_id": junction_id,
                "node1": node1,
                "node2": node2,
                "midpoint_x": float((p1[0] + p2[0]) / 2),
                "midpoint_y": float((p1[1] + p2[1]) / 2),
                "midpoint_z": float((p1[2] + p2[2]) / 2),
                "edge_length": float(np.linalg.norm(p2 - p1)),
                "triangle_ids": item["triangle_ids"],
                "component_surface_pairs": item["component_surface_pairs"],
                "incident_triangle_count": len(incident),
                "incident_min_angle_deg": min((r["angle_min_deg"] for r in incident), default=math.nan),
                "incident_max_edge_ratio": max((r["edge_ratio"] for r in incident), default=math.nan),
                "incident_max_neighbor_area_ratio": max((r["neighbor_area_ratio"] for r in incident if math.isfinite(r["neighbor_area_ratio"])), default=math.nan),
                "incident_cp_min": min((r["cp"] for r in incident), default=math.nan),
                "incident_cp_max": max((r["cp"] for r in incident), default=math.nan),
            }
        )
    return rows

def _cp_force_integration(triangle_rows: list[dict], header: dict, alpha_deg: float) -> tuple[dict, list[dict]]:
    sref = float(header["sref"])
    cref = float(header["cref"])
    bref = float(header["bref"])
    ref = np.array([header["xcg"], header["ycg"], header["zcg"]], dtype=float)
    alpha = math.radians(alpha_deg)
    lift_dir = np.array([-math.sin(alpha), 0.0, math.cos(alpha)], dtype=float)

    accum = defaultdict(lambda: {"force": np.zeros(3), "moment": np.zeros(3), "area": 0.0})
    total_force = np.zeros(3)
    total_moment = np.zeros(3)
    total_area = 0.0
    for row in triangle_rows:
        if row["surface_id"] <= 0 or not math.isfinite(row["cp"]):
            continue
        normal = np.array([row["normal_x"], row["normal_y"], row["normal_z"]], dtype=float)
        if not np.all(np.isfinite(normal)):
            continue
        dcf = -row["cp"] * normal * row["area"] / sref
        centroid = np.array([row["centroid_x"], row["centroid_y"], row["centroid_z"]])
        dm_dimless = np.cross(centroid - ref, dcf)
        total_force += dcf
        total_moment += dm_dimless
        total_area += row["area"]
        bucket = accum[(row["component_id"], row.get("surface_name", ""))]
        bucket["force"] += dcf
        bucket["moment"] += dm_dimless
        bucket["area"] += row["area"]

    def make_row(component_id, name, data):
        f = data["force"]
        m = data["moment"]
        return {
            "component_id": component_id,
            "component_surface_name": name,
            "area": data["area"],
            "CFx_cp": float(f[0]),
            "CFy_cp": float(f[1]),
            "CFz_cp": float(f[2]),
            "CL_cp": float(np.dot(f, lift_dir)),
            "CMx_cp": float(m[0] / bref),
            "CMy_cp": float(m[1] / cref),
            "CMz_cp": float(m[2] / bref),
        }

    rows = [make_row(cid, name, data) for (cid, name), data in sorted(accum.items())]
    total = make_row(0, "TOTAL", {"force": total_force, "moment": total_moment, "area": total_area})
    return total, rows

def _parse_named_value_metadata(lines: list[str]) -> dict:
    out = {}
    pattern = re.compile(r"^\s*([A-Za-z0-9_]+)\s+([-+0-9.eE]+)\s*(.*)$")
    for line in lines:
        m = pattern.match(line)
        if not m:
            continue
        key, value, units = m.groups()
        try:
            out[key] = {"value": float(value), "units": units.strip()}
        except ValueError:
            continue
    return out

def _parse_whitespace_table(path: str | Path, header_starts: tuple[str, ...]) -> dict:
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    metadata = _parse_named_value_metadata(lines)
    header = None
    rows = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if any(stripped.startswith(prefix) for prefix in header_starts):
            candidate = stripped.split()
            if len(candidate) >= 2:
                header = candidate
            continue
        if header is None:
            continue
        parts = stripped.split()
        if len(parts) != len(header):
            continue
        try:
            values = [float(x) for x in parts]
        except ValueError:
            continue
        rows.append(dict(zip(header, values)))
    return {"path": str(path), "metadata": metadata, "header": header or [], "rows": rows}

def read_history(path: str | Path) -> dict:
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    metadata = _parse_named_value_metadata(lines)
    header = None
    rows = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("Iter "):
            normalized = stripped.replace("L2 Residual", "L2_Residual").replace("Max Residual", "Max_Residual")
            header = normalized.split()
            continue
        if header is None:
            continue
        parts = stripped.split()
        if len(parts) != len(header):
            continue
        try:
            values = [float(x) for x in parts]
        except ValueError:
            continue
        rows.append(dict(zip(header, values)))
    return {"path": str(path), "metadata": metadata, "header": header or [], "rows": rows}

def read_polar(path: str | Path) -> dict:
    return _parse_whitespace_table(path, ("Beta ", "Mach ", "AoA "))

def read_lod(path: str | Path) -> dict:
    return _parse_whitespace_table(path, ("Iter ",))

def read_vspaero_config(path: str | Path) -> dict:
    path = Path(path)
    data = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()
        try:
            number = float(value)
            data[key] = int(number) if number.is_integer() else number
        except ValueError:
            data[key] = value
    return {"path": str(path), "values": data}

def read_vsp3_metadata(path: str | Path) -> dict:
    path = Path(path)
    try:
        root = ET.parse(path).getroot()
    except Exception as exc:
        return {"path": str(path), "read_error": repr(exc)}
    version = root.findtext("Version")
    vehicle_name = root.findtext("./Vehicle/ParmContainer/Name")
    return {"path": str(path), "vsp3_schema_version": version, "vehicle_name": vehicle_name}

def _lod_diagnostics(lod: dict | None, cl_abs_warning: float, cmy_abs_warning: float, symmetry_relative_warning: float) -> dict:
    if not lod or not lod.get("rows"):
        return {"outliers": [], "symmetry": [], "summary": {"lod_available": False}}
    rows = lod["rows"]
    outliers = []
    for row in rows:
        cl = row.get("Cl", math.nan)
        cmy = row.get("Cmy", math.nan)
        cl_flag = math.isfinite(cl) and abs(cl) >= cl_abs_warning
        cmy_flag = math.isfinite(cmy) and abs(cmy) >= cmy_abs_warning
        if cl_flag or cmy_flag:
            outliers.append({**row, "warn_abs_cl": cl_flag, "warn_abs_cmy": cmy_flag})

    grouped = defaultdict(dict)
    for i, row in enumerate(rows):
        y = row.get("Yavg", math.nan)
        if not math.isfinite(y) or abs(y) < 1e-8:
            continue
        key = (
            round(row.get("Xavg", 0.0), 4),
            round(abs(y), 4),
            round(row.get("Zavg", 0.0), 4),
            round(row.get("Chord", 0.0), 4),
            round(row.get("dSpan", 0.0), 4),
        )
        grouped[key]["positive" if y > 0 else "negative"] = (i, row)

    symmetry = []
    for key, pair in grouped.items():
        if "positive" not in pair or "negative" not in pair:
            continue
        ip, p = pair["positive"]
        ineg, n = pair["negative"]
        cl_diff = _relative_difference(p.get("Cl"), n.get("Cl"))
        cdi_diff = _relative_difference(p.get("Cdi"), n.get("Cdi"))
        passed = cl_diff is not None and cl_diff <= symmetry_relative_warning
        symmetry.append(
            {
                "positive_row_index": ip,
                "negative_row_index": ineg,
                "Xavg": key[0],
                "abs_Yavg": key[1],
                "Zavg": key[2],
                "Chord": key[3],
                "dSpan": key[4],
                "positive_Cl": p.get("Cl"),
                "negative_Cl": n.get("Cl"),
                "Cl_relative_difference": cl_diff,
                "Cdi_relative_difference": cdi_diff,
                "warning_threshold": symmetry_relative_warning,
                "passed": passed,
            }
        )

    max_abs_cl = max((abs(r.get("Cl", 0.0)) for r in rows if math.isfinite(r.get("Cl", math.nan))), default=None)
    max_abs_cmy = max((abs(r.get("Cmy", 0.0)) for r in rows if math.isfinite(r.get("Cmy", math.nan))), default=None)
    return {
        "outliers": outliers,
        "symmetry": symmetry,
        "summary": {
            "lod_available": True,
            "lod_row_count": len(rows),
            "lod_outlier_count": len(outliers),
            "lod_symmetry_pair_count": len(symmetry),
            "lod_symmetry_warning_count": sum(not r["passed"] for r in symmetry),
            "lod_max_abs_cl": max_abs_cl,
            "lod_max_abs_cmy": max_abs_cmy,
            "lod_cl_abs_warning_threshold": cl_abs_warning,
            "lod_cmy_abs_warning_threshold": cmy_abs_warning,
            "lod_symmetry_relative_warning_threshold": symmetry_relative_warning,
        },
    }

def _solution_diagnostics(history: dict | None, polar: dict | None, vspaero: dict | None, surface_wake_cl_relative_warning: float) -> dict:
    source = None
    final = None
    if history and history.get("rows"):
        source = "history"
        final = history["rows"][-1]
    elif polar and polar.get("rows"):
        source = "polar"
        final = polar["rows"][-1]

    summary = {
        "solution_data_available": final is not None,
        "solution_primary_source": source,
        "surface_wake_cl_relative_warning_threshold": surface_wake_cl_relative_warning,
    }
    if final is None:
        return {"summary": summary, "final": None}

    for key in (
        "Mach", "AoA", "Beta", "CLo", "CLi", "CLtot", "CDo", "CDi", "CDtot",
        "CSo", "CSi", "CStot", "CMxtot", "CMytot", "CMztot", "CLwtot", "CDwtot",
        "CSwtot", "CLiw", "CDiw", "CSiw", "L/D", "E", "LoDw", "Ew", "L2_Residual", "Max_Residual",
    ):
        if key in final:
            summary[f"solution_{key.replace('/', '_over_').replace(' ', '_')}"] = final[key]

    cltot = final.get("CLtot")
    clw = final.get("CLwtot")
    cli = final.get("CLi")
    cliw = final.get("CLiw")
    cstot = final.get("CStot")
    csw = final.get("CSwtot")
    summary["surface_wake_cl_relative_difference"] = _relative_difference(cltot, clw)
    summary["surface_wake_cli_relative_difference"] = _relative_difference(cli, cliw)
    summary["surface_wake_cs_absolute_difference"] = abs(cstot - csw) if cstot is not None and csw is not None else None
    d = summary["surface_wake_cl_relative_difference"]
    summary["surface_wake_cl_difference_advisory"] = d is not None and d >= surface_wake_cl_relative_warning

    wake_iters = None
    implicit_wake = None
    if vspaero:
        wake_iters = vspaero["values"].get("WakeIters")
        implicit_wake = vspaero["values"].get("ImplicitWake")
    summary["wake_iterations"] = wake_iters
    summary["implicit_wake"] = implicit_wake
    if wake_iters == 0 and history and len(history.get("rows", [])) == 1:
        summary["history_residual_interpretation"] = (
            "fixed_wake_single_solve: reported residual is retained as diagnostic only; "
            "do not treat the single history row as a nonlinear wake-convergence history"
        )
    else:
        summary["history_residual_interpretation"] = "standard_iteration_history_or_unknown"
    return {"summary": summary, "final": final}

def _build_extreme_cells(rows: list[dict], extreme_angle_deg: float, extreme_area_ratio: float, extreme_edge_ratio: float) -> list[dict]:
    out = []
    for row in rows:
        angle = row.get("child_min_angle_deg", row.get("angle_min_deg", math.nan))
        area_ratio = row.get("child_max_neighbor_area_ratio", row.get("neighbor_area_ratio", math.nan))
        edge_ratio = row.get("child_max_edge_ratio", row.get("edge_ratio", math.nan))
        angle_flag = math.isfinite(angle) and angle <= extreme_angle_deg
        area_flag = math.isfinite(area_ratio) and area_ratio >= extreme_area_ratio
        edge_flag = math.isfinite(edge_ratio) and edge_ratio >= extreme_edge_ratio
        cp_flag = bool(row.get("is_negative_cp_spike")) or bool(row.get("is_cp_floor_plateau"))
        if angle_flag or area_flag or edge_flag:
            out.append(
                {
                    **row,
                    "extreme_angle_warning": angle_flag,
                    "extreme_area_transition_warning": area_flag,
                    "extreme_edge_ratio_warning": edge_flag,
                    "coincident_cp_anomaly": cp_flag,
                    "strong_mesh_advisory": cp_flag and (angle_flag or area_flag or edge_flag),
                }
            )
    return out

def _json_value(value):
    """Recursively convert analysis results to strict JSON-compatible values.

    NumPy scalar classes are not all subclasses of the corresponding Python
    scalar classes (notably ``np.bool_``).  Converting only floating/integer
    values therefore leaves nested NumPy booleans in summaries and causes
    ``json.dump`` to fail.  Normalize every ``np.generic`` scalar first, then
    recurse through containers.  Non-finite floating values are represented as
    JSON null rather than non-standard NaN/Infinity literals.
    """
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, dict):
        return {str(_json_value(k)): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_value(x) for x in value]
    return value

def _write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None):
    if fieldnames is None:
        fieldnames = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8-sig") as fp:
        if not fieldnames:
            return
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _json_value(row.get(key)) for key in fieldnames})

def _autodiscover_sidecar(adb_path: Path, explicit: str | Path | None, suffix: str) -> Path | None:
    if explicit is not None:
        p = Path(explicit)
        return p if p.is_file() else None
    candidate = adb_path.with_suffix(suffix)
    return candidate if candidate.is_file() else None

def analyze_vspaero_mesh_quality(
    adb_path: str | Path,
    vspgeom_path: str | Path | None = None,
    output_dir: str | Path | None = None,
    *,
    history_path: str | Path | None = None,
    polar_path: str | Path | None = None,
    lod_path: str | Path | None = None,
    vspaero_path: str | Path | None = None,
    vsp3_path: str | Path | None = None,
    solution_case: int = 1,
    component_ids: Iterable[int] | None = None,
    surface_ids: Iterable[int] | None = None,
    bbox: tuple[float, float, float, float, float, float] | None = None,
    normal_angle_limit_deg: float = 30.0,
    spike_z_threshold: float = -6.0,
    spike_delta_threshold: float = -1.0,
    cp_floor_tolerance: float = 1e-6,
    cp_floor_min_cells: int = 4,
    cp_floor_max_cp: float = -2.0,
    small_angle_advisory_deg: float = 15.0,
    extreme_angle_deg: float = 0.1,
    extreme_area_ratio: float = 1000.0,
    extreme_edge_ratio: float = 100.0,
    kutta_coverage_min: float = 0.98,
    kutta_symmetry_tolerance_fraction: float = 0.02,
    te_param_near_miss_band: float = DEFAULT_TE_PARAM_NEAR_MISS_BAND,
    lod_cl_abs_warning: float = 3.0,
    lod_cmy_abs_warning: float = 10.0,
    lod_symmetry_relative_warning: float = 0.10,
    surface_wake_cl_relative_warning: float = 0.10,
) -> dict:
    """Run v4 mesh/Cp/Kutta/wake/solution verification diagnostics."""
    adb_path = Path(adb_path)
    adb = read_adb_v3(adb_path, solution_case=solution_case)
    surface_meta, surface_metadata_rows = _surface_metadata_from_adb(adb)
    component_metadata_rows = _build_component_metadata(surface_metadata_rows)
    triangle_rows, triangle_neighbors, topology = _build_triangle_rows(adb, surface_meta)
    metadata_check = _validate_surface_metadata(triangle_rows, surface_meta)

    surface_rows = [row for row in triangle_rows if row["surface_id"] > 0]
    surface_triangle_ids = {row["triangle_id"] for row in surface_rows}
    surface_neighbors = {
        triangle_id: {n for n in triangle_neighbors[triangle_id] if n in surface_triangle_ids}
        for triangle_id in surface_triangle_ids
    }
    _add_local_cp_metrics(surface_rows, surface_neighbors, normal_angle_limit_deg, spike_z_threshold, spike_delta_threshold)
    for row in surface_rows:
        row["warn_small_angle"] = math.isfinite(row["angle_min_deg"]) and row["angle_min_deg"] <= small_angle_advisory_deg

    mapping = None
    ngon_rows = []
    vspgeom = None
    kutta_result = None
    if vspgeom_path is not None:
        vspgeom_path = Path(vspgeom_path)
        vspgeom = read_vspgeom_v3(vspgeom_path)
        mapping = _attach_parent_ngons(surface_rows, vspgeom)
        ngon_rows = _build_ngon_rows(surface_rows, vspgeom, normal_angle_limit_deg, spike_z_threshold, spike_delta_threshold)
        for row in ngon_rows:
            row["warn_small_angle"] = math.isfinite(row["child_min_angle_deg"]) and row["child_min_angle_deg"] <= small_angle_advisory_deg
        kutta_result = _analyze_kutta(vspgeom, surface_meta, kutta_coverage_min, kutta_symmetry_tolerance_fraction, te_param_near_miss_band)

    component_ids = set(component_ids) if component_ids is not None else None
    surface_ids = set(surface_ids) if surface_ids is not None else None
    filtered_triangles = [r for r in surface_rows if _matches_filters(r, component_ids, surface_ids, bbox)]
    filtered_ngons = [r for r in ngon_rows if _matches_filters(r, component_ids, surface_ids, bbox)]

    triangle_floor = _mark_cp_floor_plateau(filtered_triangles, cp_floor_tolerance, cp_floor_min_cells, cp_floor_max_cp)
    ngon_floor = _mark_cp_floor_plateau(filtered_ngons, cp_floor_tolerance, cp_floor_min_cells, cp_floor_max_cp)
    anomaly_source = filtered_ngons if filtered_ngons else filtered_triangles
    floor_summary = ngon_floor if filtered_ngons else triangle_floor
    spikes = [row for row in anomaly_source if row["is_negative_cp_spike"]]
    spikes.sort(key=lambda r: (r["cp"], r.get("cp_local_delta", 0.0)))
    plateaus = [row for row in anomaly_source if row["is_cp_floor_plateau"]]

    if kutta_result and vspgeom:
        _attach_nearest_kutta_distance(spikes, kutta_result, vspgeom)
        _attach_nearest_kutta_distance(plateaus, kutta_result, vspgeom)

    extreme_source = filtered_ngons if filtered_ngons else filtered_triangles
    extreme_cells = _build_extreme_cells(extreme_source, extreme_angle_deg, extreme_area_ratio, extreme_edge_ratio)
    junction_rows = _build_junction_quality_rows(topology, surface_rows, adb["geometry"]["nodes"])
    cp_force_total, cp_force_components = _cp_force_integration(filtered_triangles, adb["header"], adb["solution"]["alpha_deg"])

    history_path = _autodiscover_sidecar(adb_path, history_path, ".history")
    polar_path = _autodiscover_sidecar(adb_path, polar_path, ".polar")
    lod_path = _autodiscover_sidecar(adb_path, lod_path, ".lod")
    vspaero_path = _autodiscover_sidecar(adb_path, vspaero_path, ".vspaero")
    vsp3_path = _autodiscover_sidecar(adb_path, vsp3_path, ".vsp3")
    history = read_history(history_path) if history_path else None
    polar = read_polar(polar_path) if polar_path else None
    lod = read_lod(lod_path) if lod_path else None
    vspaero = read_vspaero_config(vspaero_path) if vspaero_path else None
    vsp3_meta = read_vsp3_metadata(vsp3_path) if vsp3_path else None
    solution_diag = _solution_diagnostics(history, polar, vspaero, surface_wake_cl_relative_warning)
    lod_diag = _lod_diagnostics(lod, lod_cl_abs_warning, lod_cmy_abs_warning, lod_symmetry_relative_warning)

    cp_values = [row["cp"] for row in filtered_triangles]
    cp_array = np.asarray(cp_values, dtype=float) if cp_values else np.asarray([], dtype=float)
    structural_failures = []
    if mapping is not None and not mapping["passed"]:
        structural_failures.append("vspgeom_adb_mapping")
    if topology["same_surface_non_manifold_edges"]:
        structural_failures.append("same_surface_non_manifold")
    if not metadata_check["passed"]:
        structural_failures.append("surface_metadata")
    if kutta_result is not None and not kutta_result["passed"]:
        structural_failures.append("kutta_wake_topology")

    summary = {
        "mesh_quality_tool_version": MESH_QUALITY_TOOL_VERSION,
        "openvsp_reference_version": OPENVSP_REFERENCE_VERSION,
        "openvsp_target_version": OPENVSP_REFERENCE_VERSION,  # v3 compatibility alias
        "adb_version": adb["header"]["version"],
        "solution_case": solution_case,
        "mach": adb["solution"]["mach"],
        "alpha_deg": adb["solution"]["alpha_deg"],
        "beta_deg": adb["solution"]["beta_deg"],
        "sref": adb["header"]["sref"],
        "cref": adb["header"]["cref"],
        "bref": adb["header"]["bref"],
        "xcg": adb["header"]["xcg"],
        "ycg": adb["header"]["ycg"],
        "zcg": adb["header"]["zcg"],
        "adb_file_name": adb_path.name,
        "adb_sha256": _sha256(adb_path),
        "vspgeom_file_name": Path(vspgeom_path).name if vspgeom_path else None,
        "vspgeom_sha256": _sha256(vspgeom_path),
        "history_file_name": history_path.name if history_path else None,
        "history_sha256": _sha256(history_path),
        "polar_file_name": polar_path.name if polar_path else None,
        "polar_sha256": _sha256(polar_path),
        "lod_file_name": lod_path.name if lod_path else None,
        "lod_sha256": _sha256(lod_path),
        "vspaero_file_name": vspaero_path.name if vspaero_path else None,
        "vspaero_sha256": _sha256(vspaero_path),
        "vsp3_file_name": vsp3_path.name if vsp3_path else None,
        "vsp3_sha256": _sha256(vsp3_path),
        "n_surface_triangles": len(filtered_triangles),
        "n_ngons": len(filtered_ngons),
        "cp_min_actual": min(cp_values) if cp_values else None,
        "cp_max_actual": max(cp_values) if cp_values else None,
        "cp_p001": float(np.quantile(cp_array, 0.001)) if cp_values else None,
        "cp_p01": float(np.quantile(cp_array, 0.01)) if cp_values else None,
        "n_local_cp_spikes": len(spikes),
        "cp_floor_value": floor_summary["cp_floor_value"],
        "cp_floor_plateau_size": floor_summary["cp_floor_plateau_size"],
        "cp_floor_plateau_fraction": floor_summary["cp_floor_plateau_fraction"],
        "cp_floor_plateau_detected": floor_summary["cp_floor_plateau_detected"],
        "n_cp_limiter_candidates": len(plateaus),
        "cp_floor_cpcrit_equivalent_local_mach": floor_summary["cp_floor_cpcrit_equivalent_local_mach"],
        "n_small_angle_advisories": sum(row["warn_small_angle"] for row in filtered_triangles),
        "n_small_angle_warnings": sum(row["warn_small_angle"] for row in filtered_triangles),  # v3 compatibility alias
        "small_angle_advisory_deg": small_angle_advisory_deg,
        "n_extreme_mesh_cells": len(extreme_cells),
        "n_strong_mesh_advisories": sum(bool(r["strong_mesh_advisory"]) for r in extreme_cells),
        "extreme_angle_deg": extreme_angle_deg,
        "extreme_area_ratio": extreme_area_ratio,
        "extreme_edge_ratio": extreme_edge_ratio,
        "edge_multiplicity_gt2_count": len(topology["edge_multiplicity_gt2"]),
        "same_surface_non_manifold_edge_count": len(topology["same_surface_non_manifold_edges"]),
        "cross_surface_junction_edge_count": len(topology["cross_surface_junction_edges"]),
        "surface_metadata_checks_passed": metadata_check["passed"],
        "surface_metadata_missing_count": len(metadata_check["missing_surface_ids"]),
        "surface_metadata_missing_surface_ids": metadata_check["missing_surface_ids"],
        "surface_component_mismatch_count": len(metadata_check["component_mismatches"]),
        "surface_component_mismatches": metadata_check["component_mismatches"],
        "surface_metadata": surface_metadata_rows,
        "component_metadata": component_metadata_rows,
        "mapping_checks_passed": mapping["passed"] if mapping is not None else None,
        "mapping_mismatch_count": len(mapping["mismatches"]) if mapping is not None else None,
        "adb_trailing_vortex_edge_count": len(adb["solution"]["trailing_vortex_edges"]),
        "vspgeom_kutta_list_count": len(vspgeom["kutta_lists"]) if vspgeom else None,
        "vspgeom_kutta_node_count": sum(len(x["node_ids"]) for x in vspgeom["kutta_lists"]) if vspgeom else None,
        "kutta_checks_passed": kutta_result["passed"] if kutta_result else None,
        "kutta_coverage_failure_count": len(kutta_result["coverage_failures"]) if kutta_result else None,
        "kutta_symmetry_failure_count": len(kutta_result["symmetry_failures"]) if kutta_result else None,
        "kutta_unresolved_lifting_line_count": len(kutta_result["unresolved"]) if kutta_result else None,
        "kutta_te_param_near_miss_edge_count": len(kutta_result["near_miss_edges"]) if kutta_result else None,
        "kutta_te_param_near_miss_on_failed_surface_count": sum(not bool(x.get("surface_kutta_coverage_passed")) for x in kutta_result["near_miss_edges"] if x.get("surface_kutta_coverage_passed") is not None) if kutta_result else None,
        "kutta_surface_summary": kutta_result["surfaces"] if kutta_result else [],
        "cp_force_integral": cp_force_total,
        "structural_checks_passed": not structural_failures,
        "structural_failure_categories": structural_failures,
        "note": (
            "Threshold-based Cp, small-angle, extreme-cell, LOD and surface-vs-wake flags are diagnostic advisories, "
            "not OpenVSP validity limits. Cross-surface junction edges are not a failure by count alone. "
            "Mapping, same-surface non-manifold and incomplete Kutta/wake topology are treated as structural failures."
        ),
        **solution_diag["summary"],
        **lod_diag["summary"],
    }
    if vspaero:
        summary["vspaero_settings"] = vspaero["values"]
    if vsp3_meta:
        summary["vsp3_metadata"] = vsp3_meta

    result = {
        "summary": summary,
        "triangles": filtered_triangles,
        "ngons": filtered_ngons,
        "spikes": spikes,
        "plateaus": plateaus,
        "extreme_cells": extreme_cells,
        "junction_quality": junction_rows,
        "mapping": mapping,
        "topology": topology,
        "kutta": kutta_result,
        "cp_force_components": cp_force_components,
        "history": history,
        "polar": polar,
        "lod": lod,
        "lod_outliers": lod_diag["outliers"],
        "lod_symmetry": lod_diag["symmetry"],
        "solution": solution_diag,
    }

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        _write_csv(output_dir / "mesh_quality_triangles.csv", filtered_triangles)
        _write_csv(output_dir / "mesh_quality_ngons.csv", filtered_ngons)
        _write_csv(output_dir / "cp_spikes.csv", spikes, list(spikes[0].keys()) if spikes else ["ngon_id", "triangle_id", "cp", "nearest_kutta_distance", "nearest_kutta_endpoint_distance", "nearest_kutta_list_id"])
        _write_csv(output_dir / "cp_plateaus.csv", plateaus, list(plateaus[0].keys()) if plateaus else ["ngon_id", "triangle_id", "cp", "cp_floor_plateau_size", "nearest_kutta_distance", "nearest_kutta_endpoint_distance", "nearest_kutta_list_id"])
        _write_csv(output_dir / "mesh_quality_surfaces.csv", surface_metadata_rows, ["surface_id", "surface_name", "component_id"])
        _write_csv(output_dir / "mesh_quality_components.csv", component_metadata_rows, ["component_id", "n_surfaces", "surface_ids", "surface_names"])
        _write_csv(output_dir / "mesh_extreme_cells.csv", extreme_cells)
        _write_csv(output_dir / "junction_quality.csv", junction_rows)
        _write_csv(output_dir / "cp_force_components.csv", cp_force_components + [cp_force_total])
        _write_csv(output_dir / "adb_trailing_vortex_edges.csv", [{k: v for k, v in x.items() if k != "points"} for x in adb["solution"]["trailing_vortex_edges"]])
        if kutta_result:
            _write_csv(output_dir / "kutta_lines.csv", kutta_result["lines"], list(kutta_result["lines"][0].keys()) if kutta_result["lines"] else ["kutta_list_id", "body_wake", "wake_part_num", "n_nodes", "surface_id", "surface_name", "component_id", "outer_coverage_fraction"])
            _write_csv(output_dir / "kutta_nodes.csv", kutta_result["nodes"], list(kutta_result["nodes"][0].keys()) if kutta_result["nodes"] else ["kutta_list_id", "sequence_index", "node_id", "surface_id", "surface_name", "x", "y", "z"])
            _write_csv(output_dir / "kutta_surfaces.csv", kutta_result["surfaces"], list(kutta_result["surfaces"][0].keys()) if kutta_result["surfaces"] else ["surface_id", "surface_name", "component_id", "n_kutta_lists", "n_kutta_nodes", "kutta_outer_coverage_fraction", "kutta_coverage_passed"])
            _write_csv(output_dir / "kutta_symmetry.csv", kutta_result["symmetry"], list(kutta_result["symmetry"][0].keys()) if kutta_result["symmetry"] else ["component_id", "surface_name", "positive_surface_id", "negative_surface_id", "passed"])
            _write_csv(output_dir / "kutta_te_near_miss_edges.csv", kutta_result["near_miss_edges"], list(kutta_result["near_miss_edges"][0].keys()) if kutta_result["near_miss_edges"] else ["surface_id", "surface_name", "node1", "node2", "param_w_offset_max", "surface_kutta_coverage_passed"])
        else:
            for name in ("kutta_lines.csv", "kutta_nodes.csv", "kutta_surfaces.csv", "kutta_symmetry.csv", "kutta_te_near_miss_edges.csv"):
                _write_csv(output_dir / name, [], ["unavailable"])
        _write_csv(output_dir / "solution_history.csv", history["rows"] if history else [], history["header"] if history else ["unavailable"])
        _write_csv(output_dir / "solution_polar.csv", polar["rows"] if polar else [], polar["header"] if polar else ["unavailable"])
        _write_csv(output_dir / "lod_loads.csv", lod["rows"] if lod else [], lod["header"] if lod else ["unavailable"])
        _write_csv(output_dir / "lod_outliers.csv", lod_diag["outliers"], list(lod_diag["outliers"][0].keys()) if lod_diag["outliers"] else ["Iter", "VortexSheet", "TrailVort", "Xavg", "Yavg", "Zavg", "Cl", "Cmy", "warn_abs_cl", "warn_abs_cmy"])
        _write_csv(output_dir / "lod_symmetry.csv", lod_diag["symmetry"], list(lod_diag["symmetry"][0].keys()) if lod_diag["symmetry"] else ["positive_row_index", "negative_row_index", "Cl_relative_difference", "passed"])
        with (output_dir / "mesh_quality_summary.json").open("w", encoding="utf-8") as fp:
            json.dump(_json_value(summary), fp, ensure_ascii=False, indent=2)
    return result

def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="VSPAERO mesh/Cp/Kutta/solution verification diagnostics (v4).")
    p.add_argument("adb", type=Path)
    p.add_argument("--vspgeom", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--history", type=Path)
    p.add_argument("--polar", type=Path)
    p.add_argument("--lod", type=Path)
    p.add_argument("--vspaero", type=Path)
    p.add_argument("--vsp3", type=Path)
    p.add_argument("--solution-case", type=int, default=1)
    return p

def main() -> int:
    args = _build_arg_parser().parse_args()
    result = analyze_vspaero_mesh_quality(
        args.adb,
        args.vspgeom,
        args.output_dir,
        history_path=args.history,
        polar_path=args.polar,
        lod_path=args.lod,
        vspaero_path=args.vspaero,
        vsp3_path=args.vsp3,
        solution_case=args.solution_case,
    )
    print(json.dumps(_json_value(result["summary"]), ensure_ascii=False, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
