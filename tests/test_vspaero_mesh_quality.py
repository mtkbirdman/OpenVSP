import csv
import json
import math
import struct
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import VSPAEROMeshQuality as vmq


def write_synthetic_adb(path: Path):
    endian = "<"

    def w(fmt, value):
        fp.write(struct.pack(endian + fmt, value))

    with path.open("wb") as fp:
        w("i", vmq.ADB_MAGIC_V3)
        w("i", 1)
        w("i", 0)
        w("i", 0)
        w("i", 2)
        w("i", 4)
        w("i", 2)
        w("i", 0)
        for value in (1.0, 1.0, 1.0, 0.0, 0.0, 0.0):
            w("f", value)
        w("i", 0)  # Cart3D surfaces

        for tri in (
            (1, 2, 3, 1, 1, 0, 0.5),
            (1, 3, 4, 1, 1, 0, 0.5),
        ):
            for value in tri[:6]:
                w("i", value)
            w("f", tri[6])
        for xyz in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)):
            for value in xyz:
                w("f", float(value))

        w("i", 0)  # rotors
        w("i", 0)  # nozzles
        w("i", 1)  # mesh levels
        w("i", 0)
        w("i", 0)
        w("i", 0)  # Kutta edges
        w("i", 0)  # Kutta nodes
        w("i", 0)  # controls

        w("f", 0.1)
        w("f", math.radians(2.0))
        w("f", 0.0)
        w("f", -15.0)
        w("f", 1.0)
        for gamma, dcp in ((1.0, 0.0), (2.0, 0.0)):
            w("d", gamma)
            w("d", dcp)
        for uvw in ((1.0, 0.0, 0.0), (1.1, 0.0, 0.0)):
            for value in uvw:
                w("d", value)
        for cp in (-0.8, -4.0):
            w("f", cp)
            w("f", 0.0)
            w("f", 1.0)
        w("i", 0)


def write_synthetic_vspgeom(path: Path):
    path.write_text(
        """# vspgeom v3
1
4 1 0
0 0 0
1 0 0
1 1 0
0 1 0
1
4 1 2 3 4
1 0 0 0 1 0 1 1 0 1
1 1
0
1 2 1 2 3 1 3 4
""",
        encoding="utf-8",
    )


def make_kutta_vspgeom(full_span=True, thick_tip=False):
    if full_span:
        coords = [
            (math.nan, math.nan, math.nan),
            (1.0, -2.0, 0.0),
            (1.0, -1.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, -2.0, 0.0),
            (1.0, 1.0, 0.0),
            (1.0, 2.0, 0.0),
            (0.0, 2.0, 0.0),
        ]
        ngons = [
            {"ngon_id": 1, "node_ids": [1, 2, 3, 4, 5], "vspgeom_surface_id": 1,
             "uvs": [(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1)]},
            {"ngon_id": 2, "node_ids": [3, 6, 7, 8, 4], "vspgeom_surface_id": 2,
             "uvs": [(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1)]},
        ]
        kutta_lists = [{"kutta_list_id": 1, "body_wake": False, "wake_part_num": 1,
                        "node_ids": [1, 2, 3, 6, 7]}]
    else:
        coords = [
            (math.nan, math.nan, math.nan),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (1.0, 2.0, 0.0),
            (0.0, 2.2, 0.0),
            (0.0, 0.0, 0.0),
        ]
        ngons = [
            {"ngon_id": 1, "node_ids": [1, 2, 3, 4, 5], "vspgeom_surface_id": 1,
             "uvs": [(0, 0), (0.5, 0), (1, 0), (1, 1), (0, 1)]},
        ]
        kutta_lists = [{"kutta_list_id": 1, "body_wake": False, "wake_part_num": 1,
                        "node_ids": [1, 2, 3]}]

    node_surface_uv = defaultdict(lambda: defaultdict(list))
    for ngon in ngons:
        for node_id, uv in zip(ngon["node_ids"], ngon["uvs"]):
            node_surface_uv[ngon["vspgeom_surface_id"]][node_id].append(uv)
    return {
        "nodes": np.asarray(coords, dtype=float),
        "ngons": ngons,
        "kutta_lists": kutta_lists,
        "node_surface_uv": node_surface_uv,
    }


def test_triangle_geometry():
    nodes = np.array([
        [np.nan, np.nan, np.nan],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    metrics = vmq._triangle_geometry((1, 2, 3), nodes)
    assert metrics["area"] == pytest.approx(0.5)
    assert metrics["angle_min_deg"] == pytest.approx(45.0)
    assert metrics["angle_max_deg"] == pytest.approx(90.0)


def test_adb_vspgeom_and_report(tmp_path):
    adb = tmp_path / "synthetic.adb"
    vspgeom = tmp_path / "synthetic.vspgeom"
    write_synthetic_adb(adb)
    write_synthetic_vspgeom(vspgeom)

    result = vmq.analyze_vspaero_mesh_quality(
        adb,
        vspgeom,
        tmp_path / "out",
        runtime_openvsp_version="OpenVSP 3.52.1",
    )
    assert result["summary"]["provenance"]["runtime_openvsp_version"] == "OpenVSP 3.52.1"
    assert result["summary"]["topology"]["mapping_checks_passed"] is True
    assert result["summary"]["cp"]["min_actual"] < -3.9
    assert result["summary"]["mesh"]["surface_edge_count"] == 4
    assert (tmp_path / "out" / "mesh_edges.csv").stat().st_size > 0
    assert (tmp_path / "out" / "junction_quality.csv").read_text(encoding="utf-8-sig").startswith("junction_edge_id,")
    json.loads((tmp_path / "out" / "mesh_quality_summary.json").read_text(encoding="utf-8"))


def test_non_manifold_detection():
    nodes = np.array([
        [np.nan, np.nan, np.nan],
        [0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., -1., 0.], [0., 0., 1.],
    ])
    geometry = {
        "nodes": nodes,
        "triangles": [
            {"triangle_id": 1, "node1": 1, "node2": 2, "node3": 3, "component_id": 1, "surface_id": 1, "min_valid_timestep": 0, "stored_area": 0.5},
            {"triangle_id": 2, "node1": 2, "node2": 1, "node3": 4, "component_id": 1, "surface_id": 1, "min_valid_timestep": 0, "stored_area": 0.5},
            {"triangle_id": 3, "node1": 1, "node2": 2, "node3": 5, "component_id": 1, "surface_id": 1, "min_valid_timestep": 0, "stored_area": 0.5},
        ],
    }
    adb = {"geometry": geometry, "solution": {"cp": np.array([-1., -1., -1.]), "cp_unsteady": np.zeros(3), "gamma": np.ones(3)}}
    rows, _, topology = vmq._build_triangle_rows(adb, {1: {"surface_name": "Wing", "component_id": 1}})
    assert len(rows) == 3
    assert len(topology["same_surface_non_manifold_edges"]) == 1


def test_full_span_kutta_line_is_attributed_to_both_mirrored_surfaces():
    vspgeom = make_kutta_vspgeom(full_span=True)
    surface_meta = {
        1: {"surface_id": 1, "surface_name": "Wing", "component_id": 1},
        2: {"surface_id": 2, "surface_name": "Wing", "component_id": 1},
    }
    result = vmq._analyze_kutta(vspgeom, surface_meta, 0.98, 0.02, 1e-3)
    assert result["passed"] is True
    assert result["lines"][0]["surface_ids"] == [1, 2]
    rows = {row["surface_id"]: row for row in result["surfaces"]}
    assert rows[1]["kutta_outer_coverage_fraction"] == pytest.approx(1.0)
    assert rows[2]["kutta_outer_coverage_fraction"] == pytest.approx(1.0)
    assert result["symmetry"][0]["passed"] is True


def test_kutta_coverage_uses_trailing_edge_extent_not_thick_tip_bbox():
    vspgeom = make_kutta_vspgeom(full_span=False, thick_tip=True)
    surface_meta = {1: {"surface_id": 1, "surface_name": "HTail", "component_id": 1}}
    result = vmq._analyze_kutta(vspgeom, surface_meta, 0.98, 0.02, 1e-3)
    row = result["surfaces"][0]
    assert row["span_outer_abs"] == pytest.approx(2.2)
    assert row["te_outer_abs_span"] == pytest.approx(2.0)
    assert row["kutta_outer_coverage_fraction"] == pytest.approx(1.0)
    assert row["kutta_coverage_passed"] is True


def test_surface_mesh_size_reports_physical_edge_lengths_and_parametric_direction():
    vspgeom = make_kutta_vspgeom(full_span=False)
    surface_meta = {1: {"surface_id": 1, "surface_name": "Wing", "component_id": 1}}
    edges, summary = vmq._build_surface_mesh_size_rows(vspgeom, surface_meta, 1e-10)
    assert edges
    assert {row["param_direction"] for row in edges} >= {"u", "w"}
    assert summary[0]["all_edge_p50"] > 0
    assert summary[0]["u_edge_count"] > 0
    assert summary[0]["w_edge_count"] > 0


def test_lod_symmetry_absolute_floor_suppresses_near_zero_false_warning():
    lod = {
        "rows": [
            {"Xavg": 1.0, "Yavg": 2.0, "Zavg": 0.0, "Chord": 1.0, "dSpan": 0.5, "Cl": 0.00315, "Cdi": 0.001, "Cmy": 0.0},
            {"Xavg": 1.0, "Yavg": -2.0, "Zavg": 0.0, "Chord": 1.0, "dSpan": 0.5, "Cl": 0.00254, "Cdi": 0.001, "Cmy": 0.0},
        ]
    }
    result = vmq._lod_diagnostics(lod, 3.0, 10.0, 0.10, 1e-3)
    assert result["symmetry"][0]["Cl_relative_difference"] > 0.10
    assert result["symmetry"][0]["Cl_absolute_difference"] < 1e-3
    assert result["symmetry"][0]["passed"] is True


def test_json_normalization_handles_numpy_scalars_and_nonfinite_values():
    value = {"a": [np.bool_(True)], "b": np.int64(3), "c": np.nan, "d": np.array([True, False])}
    normalized = vmq._json_value(value)
    assert json.dumps(normalized, allow_nan=False) == '{"a": [true], "b": 3, "c": null, "d": [true, false]}'


def test_explicit_missing_sidecar_raises(tmp_path):
    adb = tmp_path / "synthetic.adb"
    write_synthetic_adb(adb)
    with pytest.raises(FileNotFoundError):
        vmq.analyze_vspaero_mesh_quality(adb, history_path=tmp_path / "missing.history")


def test_empty_csv_has_header(tmp_path):
    path = tmp_path / "empty.csv"
    vmq._write_csv(path, [], ["a", "b"])
    with path.open(encoding="utf-8-sig", newline="") as fp:
        assert next(csv.reader(fp)) == ["a", "b"]


def test_vsp3_metadata_reports_mesh_parameters(tmp_path):
    path = tmp_path / "model.vsp3"
    path.write_text(
        """<Vsp_Geometry><Version>4</Version><Vehicle><ParmContainer><Name>Vehicle</Name></ParmContainer>
<Geom><ParmContainer><Name>Wing</Name></ParmContainer><WingGeom><Tess_W Value="33"/><LECluster Value="0.5"/><TECluster Value="1"/><CapUMinOption Value="1"/><CapUMaxOption Value="2"/><CapUMinTess Value="5"/><RotateAirfoilMatchDideralFlag Value="0"/><XSecSurf>
<XSec><SectTess_U Value="6"/><InCluster Value="1"/><OutCluster Value="0.5"/><FwdCluster Value="0.25"/><AftCluster Value="1"/><TE_Close_Type Value="3"/><TE_Close_Thick Value="0"/></XSec>
</XSecSurf></WingGeom></Geom>
<VSPAEROSettings><ParmContainer><VSPAERO><GeomSet Value="4"/><ThinGeomSet Value="3"/><FixedWakeFlag Value="1"/><WakeNumIter Value="12"/><RootWakeNodes Value="64"/></VSPAERO></ParmContainer></VSPAEROSettings>
</Vehicle></Vsp_Geometry>""",
        encoding="utf-8",
    )
    metadata = vmq.read_vsp3_metadata(path)
    wing = metadata["geometry_parameters"][0]
    assert wing["tess_w"] == 33
    assert wing["le_cluster"] == pytest.approx(0.5)
    assert wing["cap_u_min_option"] == 1
    assert wing["cap_u_max_option"] == 2
    assert wing["cap_u_min_tess"] == 5
    assert wing["xsecs"][0]["sect_tess_u"] == 6
    assert wing["xsecs"][0]["fwd_cluster"] == pytest.approx(0.25)
    assert wing["xsecs"][0]["aft_cluster"] == pytest.approx(1.0)
    assert wing["xsecs"][0]["te_close_type"] == 3
    assert metadata["vspaero_settings"] == {
        "geom_set": 4,
        "thin_geom_set": 3,
        "fixed_wake_flag": 1,
        "wake_num_iter": 12,
        "root_wake_nodes": 64,
    }


def test_junction_quality_compares_cut_edge_with_local_physical_edges():
    nodes = np.array([
        [np.nan, np.nan, np.nan],
        [0.0, 0.0, 0.0], [0.01, 0.0, 0.0], [0.0, 1.0, 0.0], [0.01, -1.0, 0.0],
    ])
    triangle_rows = [
        {"triangle_id": 1, "node1": 1, "node2": 2, "node3": 3, "angle_min_deg": 1.0, "edge_ratio": 100.0, "neighbor_area_ratio": 2.0, "cp": -1.0},
        {"triangle_id": 2, "node1": 2, "node2": 1, "node3": 4, "angle_min_deg": 1.0, "edge_ratio": 100.0, "neighbor_area_ratio": 2.0, "cp": -1.1},
    ]
    topology = {
        "cross_surface_junction_edges": [{
            "edge": (1, 2),
            "triangle_ids": [1, 2],
            "component_surface_pairs": [(1, 1), (2, 2)],
        }]
    }
    row = vmq._build_junction_quality_rows(topology, triangle_rows, nodes)[0]
    assert row["edge_length"] == pytest.approx(0.01)
    assert row["local_regular_edge_p50"] > 0.9
    assert row["junction_to_local_p50_ratio"] < 0.02
