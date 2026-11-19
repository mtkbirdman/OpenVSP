"""Deterministic tests of the OpenVSP 3.50.5 ADB reader and Cp rasterizer.

The synthetic fixture is written in the exact VSP_Solver.C version-3 order,
including repeated per-case geometry, coarse edges, control surfaces and wake.
"""
from __future__ import annotations

import math
import struct
import sys
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from openvsppy.visualization import VSPAEROCpViewCapture as capture


def _make_adb(tmp_path, *, endian="<", n_cases=2, time_accurate=0, model_type=2):
    """Write small exact-layout ADB v3, two independently positioned cases."""
    path = tmp_path / "wing.adb"
    data = bytearray()
    def write(fmt, *values):
        data.extend(struct.pack(endian + fmt, *values))
    # Header: magic/model/symmetry/time; loop/node/tri/edge; six ref floats.
    write("8i", capture.MAGIC_V3, model_type, 0, time_accurate, 2, 4, 2, 1)
    write("6f", 10., 1., 10., 0., 0., 0.)
    write("i", 1)  # surface records
    write("i", 1)
    data.extend(b"Wing" + b"\0" * 96)
    write("i", 7)
    for k in range(n_cases):
        # Two triangles sharing an edge. Second case has shifted geometry.
        for nodes, area in [((1, 2, 3), .5), ((1, 3, 4), .5)]:
            write("6if", *nodes, 7, 1, 0, area)
        zshift = .25 * k
        for xyz in ((0.,0.,zshift), (1.,0.,zshift),
                    (1.,1.,zshift), (0.,1.,zshift)):
            write("3f", *xyz)
        write("2i", 0, 0)  # no rotor/nozzle
        write("i", 1)  # one coarse level
        write("2i", 4, 1)
        for xyz in ((0.,0.,zshift), (1.,0.,zshift),
                    (1.,1.,zshift), (0.,1.,zshift)):
            write("3f", *xyz)
        write("4i", 1, 0, 1, 2)  # coarse edge
        write("2i", 1, 1)  # one Kutta edge, index 1
        write("2i", 1, 2)  # one Kutta node, node 2
        write("i", 1)     # one control surface
        write("i", 0)     # polygon nodes=0
        write("9f", 0,0,0, 1,0,0, 1,0,0)
        write("2i", 1, 1)  # affected loop 1
        alpha = 2. if k == 0 else 5.
        write("5f", .5, math.radians(alpha), 0., -1.2, 1.3)
        for _ in range(2):
            write("2d", .1, .2)  # loop Gamma, dCp
        write("3d", 1., 2., 3.)  # one edge force vector
        for _ in range(2):
            write("3d", 1., 2., 3.)  # loop U/V/W
        write("3f", -.2-k, .1, 4.)
        write("3f", +.7+k, .0, 5.)
        write("i", 1)  # one wake vortex
        write("idi", 2, .25, 2)  # node, span, two sub vortex nodes
        write("6d", 1,0,0, 3,0,0)
        write("f", .03)  # control deflection
    path.write_bytes(data)
    (tmp_path / "wing.adb.cases").write_text(
        "".join(f"0.5000000 {2. if i == 0 else 5.:.7f} 0.0000000 Case {i+1}\n"
                for i in range(n_cases)), encoding="utf-8")
    return path


@pytest.mark.parametrize("endian", ["<", ">"])
def test_reader_multiple_cases_and_geometry_updates(tmp_path, endian):
    path = _make_adb(tmp_path, endian=endian)
    assert [info.index for info in capture.list_adb_cases(path)] == [1, 2]
    first = capture.read_adb_case(path, case=1)
    second = capture.read_adb_case(path, case=2)
    assert first.model_type == 2 and first.label == "Cp"
    assert np.array_equal(first.triangles, [[0, 1, 2], [0, 2, 3]])
    assert np.allclose(first.cp, [-.2, .7], atol=1e-6)
    assert np.allclose(second.cp, [-1.2, 1.7], atol=1e-6)
    assert second.case.alpha_deg == 5.
    assert np.allclose(second.vertices[:, 2] - first.vertices[:, 2], .25)
    assert np.array_equal(first.component_id, [7, 7])
    assert np.array_equal(first.surface_id, [1, 1])


def test_viewer_quantity_label_for_vlm(tmp_path):
    path = _make_adb(tmp_path, model_type=1)
    assert capture.read_adb_case(path).label == "Delta Cp"


def test_reader_detects_truncation_and_cases_mismatch(tmp_path):
    path = _make_adb(tmp_path)
    bad = path.read_bytes()[:-10]
    path.write_bytes(bad)
    with pytest.raises(capture.ADBFormatError, match="truncated"):
        capture.read_adb_case(path)
    path = _make_adb(tmp_path)
    (tmp_path / "wing.adb.cases").write_text("0.5 2 0 a\n", encoding="utf-8")
    with pytest.raises(capture.ADBFormatError, match="trailing bytes"):
        capture.read_adb_case(path)


def test_reader_rejects_unsteady_and_wrong_layout(tmp_path):
    path = _make_adb(tmp_path, time_accurate=1)
    with pytest.raises(capture.ADBFormatError, match="unsteady"):
        capture.read_adb_case(path)
    path = _make_adb(tmp_path)
    bits = bytearray(path.read_bytes())
    # Corrupt first triangle node 1, matching writer's six-int record layout.
    first_tri_offset = 8*4 + 6*4 + 4 + 4 + 100 + 4
    bits[first_tri_offset:first_tri_offset+4] = struct.pack("<i", 9999)
    path.write_bytes(bits)
    with pytest.raises(capture.ADBFormatError, match="invalid triangle node"):
        capture.read_adb_case(path)


def test_zbuffer_nearest_face_independent_of_draw_order():
    vertices = np.array([
        [0,0,0], [1,0,0], [0,1,0],
        [0,0,1], [1,0,1], [0,1,1],
    ], dtype=float)
    faces = np.array([[0,1,2],[3,4,5]])
    camera = capture._camera(capture.ViewSpec("top"))
    img, ids = capture._rasterize(vertices, faces, np.array([-1.,1.]), camera,
                                  np.array([.5,.5,0.]), 35., 60, 60, (-1.,1.))
    assert ids[30, 24] == 1
    assert tuple(img[30,24,:3]) == (255,0,255)
    img_rev, ids_rev = capture._rasterize(vertices, faces[::-1], np.array([1.,-1.]),
                                         camera, np.array([.5,.5,0.]), 35., 60, 60, (-1.,1.))
    assert ids_rev[30,24] == 0
    assert np.array_equal(img, img_rev)


def test_zbuffer_crossing_faces_changes_front_with_position():
    # Triangle A is tilted in depth; B lies at z=0, overlapping same XY extent.
    verts = np.array([
        [0,0,-1],[2,0,1],[0,2,-1],
        [0,0,0],[2,0,0],[0,2,0],
    ], float)
    tri = np.array([[0,1,2],[3,4,5]])
    camera = capture._camera(capture.ViewSpec("top"))
    _, ids = capture._rasterize(verts,tri,np.array([0.,1.]),camera,
                                np.array([1.,1.,0.]),35,100,100, (0.,1.))
    assert 0 in ids and 1 in ids


def test_mirror_symmetry_duplicates_only_half_geometry(tmp_path):
    path = _make_adb(tmp_path)
    data = capture.read_adb_case(path)
    data.symmetry_flag = 1
    # Original square has y in [0,1], so mirror it once.
    verts, tri, cp = capture._mirror_if_needed(data, "auto")
    assert len(tri) == 4 and len(cp) == 4
    assert verts[:,1].min() == -1
    assert verts[:,1].max() == 1


def test_camera_is_orthonormal_and_rotation_is_preserved():
    for view in ("top", "bottom", "left", "right", "front", "rear", "left_iso", "right_iso"):
        camera = capture._camera(capture.ViewSpec(view,90))
        assert np.allclose(camera @ camera.T, np.eye(3), atol=1e-12)
    top = capture._camera(capture.ViewSpec("top",90))
    # With x pointing downstream, negative-X nose must point down in image.
    assert (np.array([-1.,0.,0.]) @ top[1]) < 0


def test_capture_png_and_input_guards(tmp_path):
    path = _make_adb(tmp_path)
    output = tmp_path / "cp_views.png"
    final = capture.capture_adb_cp_views(path, output, size=(500,400),
                                         cp_limits=(-1.5, 1.5))
    assert final == output.resolve()
    with Image.open(output) as im:
        assert im.size == (500,400)
        assert im.format == "PNG"
        rgb = np.asarray(im.convert("RGB"))
        assert (rgb != 255).any()
    with pytest.raises(FileExistsError):
        capture.capture_adb_cp_views(path, output)
    with pytest.raises(ValueError, match="at least two"):
        capture.capture_adb_cp_views(path, tmp_path/"other.png", view_grid=(("top",),))
    with pytest.raises(ValueError, match="strictly increasing"):
        capture.capture_adb_cp_views(path, tmp_path/"other.png", cp_limits=(1.,-1.))
    with pytest.raises(ValueError, match="quantity"):
        capture.capture_adb_cp_views(path, tmp_path/"other.png", quantity="unsteady")


def _make_cube_adb(tmp_path):
    """Closed 3D Panel mesh for end-to-end visibility in every view."""
    data = bytearray()
    def write(fmt, *args):
        data.extend(struct.pack("<" + fmt, *args))
    verts = np.array([
        [0,-.5,-.5], [1,-.5,-.5], [1,+.5,-.5], [0,+.5,-.5],
        [0,-.5,+.5], [1,-.5,+.5], [1,+.5,+.5], [0,+.5,+.5],
    ])
    quads = [(0,1,2,3),(4,5,6,7),(0,1,5,4),(3,2,6,7),
             (0,3,7,4),(1,2,6,5)]
    faces = [(q[0],q[1],q[2]) for q in quads] + [(q[0],q[2],q[3]) for q in quads]
    ntri = len(faces)
    write("8i", capture.MAGIC_V3, 2, 0, 0, ntri, 8, ntri, 0)
    write("6f", 1,1,1,0,0,0)
    write("i", 0)  # no surface names
    for face in faces:
        p=verts[list(face)]
        area=.5*np.linalg.norm(np.cross(p[1]-p[0],p[2]-p[0]))
        write("6if", *(int(j)+1 for j in face), 1, 1, 0, float(area))
    for xyz in verts:
        write("3f", *map(float,xyz))
    write("2i", 0, 0)
    write("i", 1) # one coarse mesh level
    write("2i", 0, 0)
    write("3i", 0, 0, 0)  # kutta edges, kutta nodes, control surfaces
    write("5f", .4, math.radians(3.), 0, -1., 1.)
    for _ in range(ntri):
        write("2d", 0., 0.)
    for _ in range(ntri):
        write("3d", 0., 0., 0.)
    for i in range(ntri):
        write("3f", float(i)/ntri * 2 - 1, 0., 0.)
    write("i", 0)  # no wake
    p=tmp_path/"cube.adb"
    p.write_bytes(data)
    (tmp_path/"cube.adb.cases").write_text("0.4 3.0 0.0 No Comment\n", encoding="utf-8")
    return p


def test_cube_all_four_views_show_geometry(tmp_path):
    adb=_make_cube_adb(tmp_path)
    png=capture.capture_adb_cp_views(adb, tmp_path/"cube.png", size=(800,600),
                                     cp_limits=(-1,1))
    with Image.open(png) as im:
        img=np.asarray(im.convert("RGB"))
    # Non-white pixels in the central part of each cell mean the front/side
    # faces are visible, not just the planar top-view face in the simpler test.
    assert (img[110:225, 70:300] != 255).any()    # top
    assert (img[110:225, 500:735] != 255).any()   # left
    assert (img[350:470, 70:300] != 255).any()    # front
    assert (img[350:470, 500:735] != 255).any()   # isometric
