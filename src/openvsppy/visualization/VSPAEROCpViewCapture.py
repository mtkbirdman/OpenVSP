"""Render OpenVSP 3.50.5 VSPAERO ADB pressure solutions without the Viewer.

The binary layout follows VSP_Solver.C's WriteOutAerothermalDatabase* and
Viewer/glviewer.C's LoadMeshData / UpdateMeshData / LoadExistingSolutionData.
Only version-3 *steady* ADBs without propulsion records are supported.
Surface pressure is constant on each triangle, as in the Viewer per-triangle
mode. NumPy performs orthographic projection and pixelwise Z-buffering;
Matplotlib's Agg backend composes the final PNG.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .OpenVSPViewCapture import (
    ORTHOGRAPHIC_FOUR_VIEW_LAYOUT,
    ViewSpec,
    _normalize_view_grid,
    _split_extent,
    _validate_size,
)

MAGIC_V3 = -123789453
# VSPAERO Viewer glviewer.C: percent_to_rgb(per, ..., 0).
_CP_RGB = np.array([
    [0., 0., 1.], [0., 1., 1.], [0., 1., 0.],
    [1., 1., 0.], [1., 0., 0.], [1., 0., 1.],
], dtype=np.float64)
_CP_STOPS = np.linspace(0., 1., 6)
_ORTHOGRAPHIC = {"top", "bottom", "front", "rear", "left", "right"}
_SUPPORTED = _ORTHOGRAPHIC | {"left_iso", "right_iso"}


class ADBFormatError(ValueError):
    """Unsupported, truncated, or internally inconsistent ADB binary data."""


@dataclass(frozen=True)
class ADBCaseInfo:
    index: int
    mach: float
    alpha_deg: float
    beta_deg: float
    comment: str


@dataclass
class ADBPressureCase:
    case: ADBCaseInfo
    vertices: np.ndarray
    triangles: np.ndarray
    cp: np.ndarray
    cp_unsteady: np.ndarray
    component_id: np.ndarray
    surface_id: np.ndarray
    model_type: int
    symmetry_flag: int
    source_cp_limits: tuple[float, float]

    @property
    def label(self) -> str:
        return "Delta Cp" if self.model_type == 1 else "Cp"


class _ADBReader:
    """A bounded, sequential binary reader with contextual error messages."""

    def __init__(self, path: Path):
        self.path = path
        self.data = path.read_bytes()
        self.offset = 0
        if len(self.data) < 4:
            raise ADBFormatError(f"{path}: empty/truncated ADB header")
        little = int.from_bytes(self.data[:4], "little", signed=True)
        big = int.from_bytes(self.data[:4], "big", signed=True)
        if little == MAGIC_V3:
            self.endian = "<"
        elif big == MAGIC_V3:
            self.endian = ">"
        else:
            raise ADBFormatError(
                f"{path}: unsupported ADB magic {little} (expected {MAGIC_V3}); "
                "only OpenVSP 3.50.5 version-3 layout is supported"
            )

    def array(self, kind: str, count: int, name: str) -> np.ndarray:
        if count < 0 or count > 100_000_000:
            raise ADBFormatError(f"{self.path}: {name} invalid element count {count} at byte {self.offset}")
        dt = np.dtype(self.endian + kind)
        length = dt.itemsize * count
        end = self.offset + length
        if end > len(self.data):
            raise ADBFormatError(
                f"{self.path}: {name} truncated at byte {self.offset}: "
                f"expected {length} bytes, {len(self.data)-self.offset} available"
            )
        values = np.frombuffer(self.data, dtype=dt, count=count, offset=self.offset)
        self.offset = end
        return values

    def integer(self, name: str) -> int:
        return int(self.array("i4", 1, name)[0])

    def count(self, name: str, *, minimum: int = 0) -> int:
        n = self.integer(name)
        if n < minimum or n > 20_000_000:
            raise ADBFormatError(f"{self.path}: {name} invalid count {n} at byte {self.offset-4}")
        return n

    def scalar(self, name: str) -> float:
        return float(self.array("f4", 1, name)[0])

    def skip(self, kind: str, count: int, name: str) -> None:
        self.array(kind, count, name)

    def raw(self, n: int, name: str) -> bytes:
        if self.offset + n > len(self.data):
            raise ADBFormatError(f"{self.path}: {name} truncated at byte {self.offset}")
        result = self.data[self.offset:self.offset+n]
        self.offset += n
        return result


def _case_file(adb_path: Path) -> list[ADBCaseInfo]:
    filename = Path(str(adb_path) + ".cases")
    if not filename.is_file():
        raise FileNotFoundError(f"Missing corresponding ADB cases file: {filename}")
    result = []
    for line_no, line in enumerate(filename.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if not line.strip():
            raise ADBFormatError(f"{filename}:{line_no}: empty case record")
        fields = line.split(maxsplit=3)
        if len(fields) < 3:
            raise ADBFormatError(f"{filename}:{line_no}: expected Mach Alpha Beta")
        try:
            mach, alpha, beta = [float(v) for v in fields[:3]]
        except ValueError as exc:
            raise ADBFormatError(f"{filename}:{line_no}: invalid case conditions") from exc
        if not all(map(math.isfinite, (mach, alpha, beta))):
            raise ADBFormatError(f"{filename}:{line_no}: non-finite case conditions")
        result.append(ADBCaseInfo(len(result)+1, mach, alpha, beta, fields[3].strip() if len(fields) > 3 else ""))
    if not result:
        raise ADBFormatError(f"{filename}: no solution cases")
    return result


def _read_header(reader: _ADBReader) -> dict:
    if reader.integer("magic") != MAGIC_V3:
        raise ADBFormatError(f"{reader.path}: magic changed unexpectedly")
    model_type = reader.integer("ModelType")
    symmetry = reader.integer("SymmetryFlag")
    time_accurate = reader.integer("TimeAccurate")
    n_loop = reader.count("NumberOfVortexLoops")
    n_node = reader.count("NumberOfNodes", minimum=3)
    n_tri = reader.count("NumberOfTris", minimum=1)
    n_edge = reader.count("NumberOfSurfaceVortexEdges")
    reader.skip("f4", 6, "Sref Cref Bref Xcg Ycg Zcg")
    if model_type not in (1, 2):
        raise ADBFormatError(f"{reader.path}: unsupported model type {model_type} (VLM=1, PANEL=2)")
    if time_accurate:
        raise ADBFormatError(f"{reader.path}: unsteady ADB not supported (TimeAccurate={time_accurate})")
    n_surface = reader.count("NumberOfSurfaces")
    for i in range(n_surface):
        reader.integer(f"surface[{i}].ID")
        reader.raw(100, f"surface[{i}].Name")
        reader.integer(f"surface[{i}].ComponentID")
    return dict(model_type=model_type, symmetry=symmetry,
                n_loop=n_loop, n_node=n_node, n_tri=n_tri, n_edge=n_edge)


def _read_geometry(reader: _ADBReader, hdr: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    n_tri, n_node = hdr["n_tri"], hdr["n_node"]
    # Six 32-bit integer fields and one float area, interleaved per triangle.
    records = reader.array("i4", n_tri * 7, "triangles/area").reshape(-1, 7)
    tri = records[:, :3].astype(np.int64) - 1
    if (tri < 0).any() or (tri >= n_node).any():
        raise ADBFormatError(f"{reader.path}: invalid triangle node ID at byte {reader.offset}")
    component = records[:, 3].astype(np.int32).copy()
    surface = records[:, 4].astype(np.int32).copy()
    # Source has MinValidTimeStep at column 5; the final column is float32 bits.
    areas = records[:, 6].view(np.dtype(reader.endian + "f4")).astype(np.float64)
    vertices = reader.array("f4", n_node * 3, "node coordinates").reshape(-1, 3).astype(np.float64)
    if not np.isfinite(vertices).all() or not np.isfinite(areas).all() or (areas < 0).any():
        raise ADBFormatError(f"{reader.path}: non-finite coordinates or invalid triangle areas")
    calculated = 0.5 * np.linalg.norm(
        np.cross(vertices[tri[:, 1]] - vertices[tri[:, 0]],
                 vertices[tri[:, 2]] - vertices[tri[:, 0]]), axis=1)
    if not np.allclose(calculated, areas, rtol=5e-3, atol=max(1e-8, calculated.max(initial=0) * 1e-7)):
        raise ADBFormatError(f"{reader.path}: triangle area does not match node geometry (possible incompatible ADB layout)")

    rotors = reader.count("NumberOfRotors")
    nozzles = reader.count("NumberOfNozzles")
    if rotors or nozzles:
        raise ADBFormatError(
            f"{reader.path}: rotor/nozzle propulsion records not yet supported ({rotors} rotors, {nozzles} nozzles)"
        )
    n_levels = reader.count("NumberOfMeshLevels")
    coarse_edges = []
    for level in range(n_levels):
        n_cn = reader.count(f"level[{level}].NumberOfNodes")
        n_ce = reader.count(f"level[{level}].NumberOfEdges")
        reader.skip("f4", 3 * n_cn, f"level[{level}].nodes")
        reader.skip("i4", 4 * n_ce, f"level[{level}].edges")
        coarse_edges.append(n_ce)
    n_kutta_edges = reader.count("NumberOfKuttaEdges")
    kutta_edges = reader.array("i4", n_kutta_edges, "KuttaEdges")
    if n_kutta_edges and (not coarse_edges or (kutta_edges < 1).any() or (kutta_edges > coarse_edges[0]).any()):
        raise ADBFormatError(f"{reader.path}: invalid Kutta edge IDs")
    reader.skip("i4", reader.count("NumberOfKuttaNodes"), "KuttaNodes")
    n_control = reader.count("NumberOfControlSurfaces")
    for i in range(n_control):
        n = reader.count(f"control[{i}].NumberOfNodes")
        reader.skip("f4", 3 * n + 9, f"control[{i}].nodes/hinge")
        reader.skip("i4", reader.count(f"control[{i}].NumberOfLoops"), f"control[{i}].loops")
    return vertices, tri, component, surface, n_control


def _read_solution(reader: _ADBReader, hdr: dict, n_control: int, info: ADBCaseInfo) -> tuple[np.ndarray, np.ndarray, tuple[float, float]]:
    vals = reader.array("f4", 5, f"case {info.index} Mach Alpha Beta CpMin CpMax").astype(np.float64)
    mach, alpha_rad, beta_rad, cp_min, cp_max = vals
    if not np.isfinite(vals).all():
        raise ADBFormatError(f"{reader.path}: case {info.index} invalid solution header")
    if not np.allclose((mach, np.degrees(alpha_rad), np.degrees(beta_rad)),
                       (info.mach, info.alpha_deg, info.beta_deg), rtol=1e-5, atol=3e-4):
        raise ADBFormatError(
            f"{reader.path}: case {info.index} solution conditions do not match .adb.cases: "
            f"ADB={(mach, np.degrees(alpha_rad), np.degrees(beta_rad))}, "
            f"cases={(info.mach, info.alpha_deg, info.beta_deg)}"
        )
    reader.skip("f8", 2 * hdr["n_loop"], f"case {info.index} loop Gamma/dCp")
    reader.skip("f8", 3 * hdr["n_edge"], f"case {info.index} edge forces")
    reader.skip("f8", 3 * hdr["n_loop"], f"case {info.index} loop velocities")
    pressure = reader.array("f4", 3 * hdr["n_tri"], f"case {info.index} Cp CpUnsteady Gamma").reshape(-1, 3)
    if not np.isfinite(pressure).all():
        raise ADBFormatError(f"{reader.path}: case {info.index} non-finite pressure values")
    n_wake = reader.count(f"case {info.index} NumberOfTrailingVortexEdges")
    for i in range(n_wake):
        te_node = reader.integer(f"case {info.index} wake[{i}].TE_node")
        if te_node < 1 or te_node > hdr["n_node"]:
            raise ADBFormatError(f"{reader.path}: invalid wake trailing-edge node {te_node}")
        reader.skip("f8", 1, f"case {info.index} wake[{i}].span")
        reader.skip("f8", 3 * reader.count(f"case {info.index} wake[{i}].NumberOfSubVortices"),
                    f"case {info.index} wake[{i}].points")
    reader.skip("f4", n_control, f"case {info.index} control deflection")
    return pressure[:, 0].astype(np.float64), pressure[:, 1].astype(np.float64), (float(cp_min), float(cp_max))


def list_adb_cases(adb_path: str | Path) -> list[ADBCaseInfo]:
    """List case conditions from an ADB's companion ``.adb.cases`` file."""
    path = Path(adb_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    return _case_file(path)


def read_adb_case(adb_path: str | Path, case: int = 1) -> ADBPressureCase:
    """Read a steady OpenVSP 3.50.5 version-3 ADB solution (case is 1-based)."""
    path = Path(adb_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    cases = _case_file(path)
    if isinstance(case, bool) or not isinstance(case, int) or not 1 <= case <= len(cases):
        raise ValueError(f"case must be an integer in 1..{len(cases)}")
    reader = _ADBReader(path)
    hdr = _read_header(reader)
    chosen = None
    for info in cases:
        verts, tris, comps, surfaces, n_control = _read_geometry(reader, hdr)
        cp, cp_unsteady, limits = _read_solution(reader, hdr, n_control, info)
        if info.index == case:
            chosen = ADBPressureCase(info, verts, tris, cp, cp_unsteady, comps, surfaces,
                                     hdr["model_type"], hdr["symmetry"], limits)
    if reader.offset != len(reader.data):
        raise ADBFormatError(
            f"{path}: {len(reader.data)-reader.offset} unexpected trailing bytes after "
            f"{len(cases)} cases at byte {reader.offset}; file format/case count may differ"
        )
    assert chosen is not None
    return chosen


def _mirror_if_needed(data: ADBPressureCase, mode: bool | str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if mode not in (True, False, "auto"):
        raise ValueError("reflect_symmetry must be True, False, or 'auto'")
    visible = data.surface_id > 0
    if not np.any(visible):
        raise ValueError("ADB contains no positive-SurfaceID body/wing triangles to display")
    triangles = data.triangles[visible]
    cp = data.cp[visible]
    vertices = data.vertices
    do_mirror = mode is True or (mode == "auto" and data.symmetry_flag != 0)
    if do_mirror:
        scale = max(1., np.ptp(vertices[:, 1]))
        tolerance = max(1e-8, scale * 1e-6)
        y = vertices[:, 1]
        if np.any(y > tolerance) and np.any(y < -tolerance):
            if mode is True:
                raise ValueError("Cannot mirror geometry already spanning both sides of Y=0")
            do_mirror = False
        elif np.all(np.abs(y) <= tolerance):
            do_mirror = False
    if do_mirror:
        # Do not duplicate triangles lying entirely in the symmetry plane.
        mirrorable = np.any(np.abs(vertices[triangles, 1]) > tolerance, axis=1)
        reflected = vertices.copy()
        reflected[:, 1] *= -1
        offset = len(vertices)
        triangles = np.concatenate((triangles, triangles[mirrorable][:, ::-1] + offset), axis=0)
        cp = np.concatenate((cp, cp[mirrorable]))
        vertices = np.concatenate((vertices, reflected), axis=0)
    return vertices, triangles, cp


def _camera(view: ViewSpec) -> np.ndarray:
    """Rows are camera-right, camera-up, and camera-toward-viewer unit axes."""
    if view.view not in _SUPPORTED:
        raise ValueError(f"view {view.view!r} not supported for ADB capture")
    axes = {
        "top": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        "bottom": ((1, 0, 0), (0, -1, 0), (0, 0, -1)),
        "front": ((0, -1, 0), (0, 0, 1), (-1, 0, 0)),
        "rear": ((0, 1, 0), (0, 0, 1), (1, 0, 0)),
        "left": ((1, 0, 0), (0, 0, 1), (0, -1, 0)),
        "right": ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
    }
    if view.view in axes:
        camera = np.asarray(axes[view.view], dtype=np.float64)
    else:
        eye = np.array([-1., -1. if view.view == "left_iso" else 1., 0.85])
        eye /= np.linalg.norm(eye)
        right = np.cross(np.array([0., 0., 1.]), eye)
        right /= np.linalg.norm(right)
        up = np.cross(eye, right)
        camera = np.stack((right, up, eye))
    # Existing ViewSpec rotation is counter-clockwise as seen on the image.
    radians = np.deg2rad(view.rotate_deg)
    c, s = np.cos(radians), np.sin(radians)
    right, up = camera[:2].copy()
    camera[0] = c * right - s * up
    camera[1] = s * right + c * up
    return camera


def _rasterize(
    vertices: np.ndarray, triangles: np.ndarray, cp: np.ndarray,
    camera: np.ndarray, center: np.ndarray, scale: float,
    width: int, height: int, limits: tuple[float, float],
    *, transparent: bool = False, show_mesh_edges: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Pixelwise orthographic Z-buffer; triangle colors remain constant per face."""
    projected = (vertices - center) @ camera.T
    projected[:, 0] = (width - 1) / 2. + scale * projected[:, 0]
    projected[:, 1] = (height - 1) / 2. - scale * projected[:, 1]
    depth = np.full((height, width), -np.inf, dtype=np.float64)
    tri_id = np.full((height, width), -1, dtype=np.int32)
    # Bound each vectorized operation to a tile; avoid per-pixel Python loops.
    for index, face in enumerate(triangles):
        a, b, c = projected[face]
        den = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
        if abs(den) < 1e-14:
            continue
        xmin = max(0, int(math.floor(min(a[0], b[0], c[0]) - .5)))
        xmax = min(width, int(math.ceil(max(a[0], b[0], c[0]) + .5)))
        ymin = max(0, int(math.floor(min(a[1], b[1], c[1]) - .5)))
        ymax = min(height, int(math.ceil(max(a[1], b[1], c[1]) + .5)))
        if xmax <= xmin or ymax <= ymin:
            continue
        for y0 in range(ymin, ymax, 128):
            y1 = min(y0 + 128, ymax)
            for x0 in range(xmin, xmax, 128):
                x1 = min(x0 + 128, xmax)
                yy, xx = np.ogrid[y0:y1, x0:x1]
                xx = xx + .5
                yy = yy + .5
                w0 = ((b[1]-c[1])*(xx-c[0])+(c[0]-b[0])*(yy-c[1])) / den
                w1 = ((c[1]-a[1])*(xx-c[0])+(a[0]-c[0])*(yy-c[1])) / den
                w2 = 1. - w0 - w1
                inside = (w0 >= -1e-12) & (w1 >= -1e-12) & (w2 >= -1e-12)
                if not inside.any():
                    continue
                candidate = w0*a[2] + w1*b[2] + w2*c[2]
                current = depth[y0:y1, x0:x1]
                old_id = tri_id[y0:y1, x0:x1]
                tolerance = 1e-10 * np.maximum(1., np.abs(candidate))
                update = inside & ((candidate > current + tolerance) |
                                   ((np.abs(candidate - current) <= tolerance) & (index < old_id)))
                current[update] = candidate[update]
                old_id[update] = index
    visible = tri_id >= 0
    rgba = np.full((height, width, 4), 255, dtype=np.uint8)
    if transparent:
        rgba[:, :, 3] = 0
        rgba[:, :, 3][visible] = 255
    if visible.any():
        normal = np.clip((cp[tri_id[visible]] - limits[0]) / (limits[1] - limits[0]), 0., 1.)
        values = np.stack([np.interp(normal, _CP_STOPS, _CP_RGB[:, i]) for i in range(3)], axis=1)
        rgba[visible, :3] = np.rint(255 * values).astype(np.uint8)
    if show_mesh_edges and visible.any():
        boundary = np.zeros_like(visible)
        boundary[1:] |= visible[1:] & (tri_id[1:] != tri_id[:-1])
        boundary[:, 1:] |= visible[:, 1:] & (tri_id[:, 1:] != tri_id[:, :-1])
        rgba[boundary, :3] = 32
    return rgba, tri_id


def capture_adb_cp_views(
    adb_path: str | Path,
    output_path: str | Path,
    *,
    case: int = 1,
    quantity: str = "total",
    size: tuple[int, int] = (1600, 1200),
    view_grid: Sequence[Sequence[ViewSpec | str]] = ORTHOGRAPHIC_FOUR_VIEW_LAYOUT,
    cp_limits: tuple[float, float] | None = None,
    align_orthographic_views: bool = True,
    content_padding_px: int = 12,
    reflect_symmetry: bool | str = "auto",
    show_colorbar: bool = True,
    show_mesh_edges: bool = False,
    transparent_background: bool = False,
    overwrite: bool = False,
) -> Path:
    """Save Viewer-style per-triangle Cp/Delta-Cp as an orthographic PNG view grid.

    No OpenVSP GUI, VSPAERO solver, or native 3D library is invoked. Viewer
    lighting/antialiasing are not replicated. ModelType=2 labels the global
    result Cp; an unidentifiable thin/thick mixed model must be inspected
    separately because ADB does not expose per-face pressure semantics.
    """
    if quantity != "total":
        raise ValueError("Only quantity='total' is supported for steady ADB files")
    width, height = _validate_size(size)
    views = _normalize_view_grid(view_grid)
    if any(item.view not in _SUPPORTED for row in views for item in row):
        raise ValueError("ADB supports top/bottom/front/rear/left/right/left_iso/right_iso views")
    n_ortho = sum(item.view in _ORTHOGRAPHIC for row in views for item in row)
    if align_orthographic_views and n_ortho < 2:
        raise ValueError("align_orthographic_views requires at least two orthographic views")
    if isinstance(content_padding_px, bool) or not isinstance(content_padding_px, int) or content_padding_px < 0:
        raise ValueError("content_padding_px must be a nonnegative integer")
    if not isinstance(align_orthographic_views, bool):
        raise TypeError("align_orthographic_views must be a bool")
    if not isinstance(transparent_background, bool):
        raise TypeError("transparent_background must be a bool")
    destination = Path(output_path).expanduser().resolve()
    if destination.suffix.lower() != ".png":
        raise ValueError("output_path must end in .png")
    if destination.exists() and not overwrite:
        raise FileExistsError(destination)

    result = read_adb_case(adb_path, case)
    vertices, tris, cp = _mirror_if_needed(result, reflect_symmetry)
    if not np.isfinite(cp).all():
        raise ValueError("Non-finite pressures")
    if cp_limits is None:
        low, high = float(cp.min()), float(cp.max())
        if low == high:
            low, high = low - .5, high + .5
    else:
        if len(cp_limits) != 2:
            raise ValueError("cp_limits must be a two-element sequence")
        low, high = map(float, cp_limits)
        if not math.isfinite(low) or not math.isfinite(high) or low >= high:
            raise ValueError("cp_limits must be finite and strictly increasing")
    limits = (low, high)

    # Exact pixel layout: reserved title, labels, and optional shared colorbar.
    left_margin, right_margin, gutter, top_margin, bottom_margin = 16, 16, 8, 50, 12
    colorbar_height = 78 if show_colorbar else 0
    inner_width = width - left_margin - right_margin - (len(views[0]) - 1) * gutter
    inner_height = height - top_margin - bottom_margin - colorbar_height - (len(views) - 1) * gutter
    if inner_width < len(views[0]) * 20 or inner_height < len(views) * 40:
        raise ValueError("Output size is too small for labels, view grid, and colorbar")
    columns = _split_extent(inner_width, len(views[0]))
    rows = _split_extent(inner_height, len(views))
    panel_dims = {(r, c): (columns[c], rows[r]-22) for r in range(len(rows)) for c in range(len(columns))}
    if any(min(w,h) <= 2 * content_padding_px for w,h in panel_dims.values()):
        raise ValueError("content_padding_px leaves no room for a model view")

    center = (vertices.min(axis=0) + vertices.max(axis=0)) / 2
    offsets = vertices - center
    cameras = {position: _camera(views[position[0]][position[1]]) for position in panel_dims}
    extents = {}
    for key, cam in cameras.items():
        xy = offsets @ cam[:2].T
        extents[key] = np.ptp(xy, axis=0)
    def best_scale(key):
        w, h = panel_dims[key]
        a, b = extents[key]
        options = [(w-2*content_padding_px)/a if a > 0 else np.inf,
                   (h-2*content_padding_px)/b if b > 0 else np.inf]
        return min(options)
    ortho_keys = [key for key in panel_dims if views[key[0]][key[1]].view in _ORTHOGRAPHIC]
    common_scale = min(best_scale(key) for key in ortho_keys) if ortho_keys else None

    # Do not choose pyplot's display backend; rendering must also work in CI.
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.figure import Figure
    from matplotlib.cm import ScalarMappable

    dpi = 100
    fig = Figure(figsize=(width / dpi, height / dpi), dpi=dpi,
                 facecolor=(1,1,1,0) if transparent_background else "white")
    FigureCanvasAgg(fig)
    title = (f"{Path(adb_path).name} | Case {result.case.index} | "
             f"M={result.case.mach:g}, alpha={result.case.alpha_deg:g} deg, "
             f"beta={result.case.beta_deg:g} deg | {result.label}")
    fig.text(0.5, 1 - 24/height, title, ha="center", va="center", fontsize=11)
    cmap = LinearSegmentedColormap.from_list("vspaero-3505", list(zip(_CP_STOPS, _CP_RGB)), N=256)
    y = top_margin
    for r, row in enumerate(views):
        x = left_margin
        for c, spec in enumerate(row):
            w, h = panel_dims[(r,c)]
            scale = common_scale if align_orthographic_views and spec.view in _ORTHOGRAPHIC else best_scale((r,c))
            image, _ = _rasterize(vertices, tris, cp, cameras[(r,c)], center, scale,
                                  w, h, limits, transparent=transparent_background,
                                  show_mesh_edges=show_mesh_edges)
            # Pixel-precise axes; no auto layout or rescaling of the output cells.
            axis = fig.add_axes([x/width, (height-y-22-h)/height, w/width, h/height])
            axis.imshow(image, origin="upper", interpolation="nearest", aspect="auto")
            axis.axis("off")
            fig.text((x+w/2)/width, 1-(y+9)/height, spec.view.upper().replace("_", " "),
                     ha="center", va="center", fontsize=9)
            x += w + gutter
        y += rows[r] + gutter
    if show_colorbar:
        cb_w = max(60, min(width-120, int(width * .55)))
        cb_x = (width-cb_w)/2
        cb_y = height - bottom_margin - 45
        cb_ax = fig.add_axes([cb_x/width, (height-cb_y-15)/height, cb_w/width, 15/height])
        mappable = ScalarMappable(norm=Normalize(vmin=low, vmax=high), cmap=cmap)
        cb = fig.colorbar(mappable, cax=cb_ax, orientation="horizontal")
        cb.set_label(result.label, fontsize=9)
        cb.ax.tick_params(labelsize=8)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Save to a temporary sibling; never leave a partially written final PNG.
    temporary = destination.with_name(destination.name + ".tmp")
    try:
        fig.savefig(temporary, format="png", dpi=dpi, transparent=transparent_background)
        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return destination
