import math
import struct
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import VSPAEROMeshQuality as vmq


def write_synthetic_adb(path: Path):
    endian = "<"

    def w(fmt, value):
        fp.write(struct.pack(endian + fmt, value))

    with path.open("wb") as fp:
        w("i", vmq.ADB_MAGIC_V3)
        w("i", 1)  # model type
        w("i", 0)  # symmetry
        w("i", 0)  # time accurate
        w("i", 2)  # computational vortex loops
        w("i", 4)  # nodes
        w("i", 2)  # tris
        w("i", 0)  # surface vortex edges

        for value in (1.0, 1.0, 1.0, 0.0, 0.0, 0.0):
            w("f", value)

        w("i", 0)  # Cart3d surfaces

        # Geometry record
        triangles = [
            (1, 2, 3, 1, 1, 0, 0.5),
            (1, 3, 4, 1, 1, 0, 0.5),
        ]
        for tri in triangles:
            for value in tri[:6]:
                w("i", value)
            w("f", tri[6])

        for xyz in ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)):
            for value in xyz:
                w("f", float(value))

        w("i", 0)  # rotors
        w("i", 0)  # nozzles
        w("i", 1)  # mesh levels
        w("i", 0)  # coarse nodes
        w("i", 0)  # coarse edges
        w("i", 0)  # Kutta edges
        w("i", 0)  # Kutta nodes
        w("i", 0)  # controls

        # Solution
        w("f", 0.1)  # Mach
        w("f", math.radians(2.0))
        w("f", 0.0)
        w("f", -15.0)  # solver clip/min field, not actual min
        w("f", 1.0)

        for gamma, dcp in ((1.0, 0.0), (2.0, 0.0)):
            w("d", gamma)
            w("d", dcp)

        # no edge forces
        for uvw in ((1.0, 0.0, 0.0), (1.1, 0.0, 0.0)):
            for value in uvw:
                w("d", value)

        for cp in (-0.8, -4.0):
            w("f", cp)
            w("f", 0.0)
            w("f", 1.0)

        w("i", 0)  # trailing vortex edges
        # no control deflections


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


def test_triangle_geometry():
    nodes = np.array([
        [np.nan, np.nan, np.nan],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    m = vmq._triangle_geometry((1, 2, 3), nodes)
    assert abs(m["area"] - 0.5) < 1e-12
    assert abs(m["angle_min_deg"] - 45.0) < 1e-12
    assert abs(m["angle_max_deg"] - 90.0) < 1e-12


def test_adb_and_vspgeom(tmp_path):
    adb = tmp_path / "synthetic.adb"
    vspgeom = tmp_path / "synthetic.vspgeom"
    write_synthetic_adb(adb)
    write_synthetic_vspgeom(vspgeom)

    parsed = vmq.read_adb_v3(adb)
    assert parsed["header"]["version"] == 3
    assert parsed["header"]["number_of_tris"] == 2
    assert np.allclose(parsed["solution"]["cp"], [-0.8, -4.0])

    vg = vmq.read_vspgeom_v3(vspgeom)
    assert vg["num_loops"] == 1
    assert vg["child_to_parent"] == {1: 1, 2: 1}

    out = tmp_path / "out"
    result = vmq.analyze_vspaero_mesh_quality(adb, vspgeom, out)
    assert result["summary"]["mapping_checks_passed"] is True
    assert result["summary"]["cp_min_actual"] < -3.9
    assert (out / "mesh_quality_triangles.csv").exists()
    assert (out / "mesh_quality_ngons.csv").exists()
    assert (out / "mesh_quality_summary.json").exists()


def test_non_manifold_detection(tmp_path):
    # Unit-level topology test using three triangles sharing one edge.
    nodes = np.array([
        [np.nan, np.nan, np.nan],
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., -1., 0.],
        [0., 0., 1.],
    ])
    geometry = {
        "nodes": nodes,
        "triangles": [
            {"triangle_id": 1, "node1": 1, "node2": 2, "node3": 3, "component_id": 1, "surface_id": 1, "min_valid_timestep": 0, "stored_area": 0.5},
            {"triangle_id": 2, "node1": 2, "node2": 1, "node3": 4, "component_id": 1, "surface_id": 1, "min_valid_timestep": 0, "stored_area": 0.5},
            {"triangle_id": 3, "node1": 1, "node2": 2, "node3": 5, "component_id": 1, "surface_id": 1, "min_valid_timestep": 0, "stored_area": 0.5},
        ],
    }
    adb = {
        "geometry": geometry,
        "solution": {
            "cp": np.array([-1., -1., -1.]),
            "cp_unsteady": np.zeros(3),
            "gamma": np.ones(3),
        },
    }
    rows, neighbors, topology = vmq._build_triangle_rows(adb)
    assert len(topology["non_manifold_edges"]) == 1
