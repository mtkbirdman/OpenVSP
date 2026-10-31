"""Capture OpenVSP views and turn an ordered model list into an animation.

OpenVSP's graphics switches must be set before ``openvsp`` is imported. The
public functions therefore perform all OpenVSP work in a short-lived child
process. A single screenshot starts one GUI for one model; an animation starts
one GUI and reloads every model in the supplied order.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

try:
    from PIL import Image
except ImportError as exc:  # pragma: no cover - exercised only without Pillow
    raise ImportError(
        "Pillow is required for OpenVSP screenshot composition. "
        "Install it with: python -m pip install pillow"
    ) from exc

_VIEW_ENUM_NAMES = {
    "top": "CAM_TOP",
    "front": "CAM_FRONT",
    "front_y_up": "CAM_FRONT_YUP",
    "left": "CAM_LEFT",
    "left_iso": "CAM_LEFT_ISO",
    "bottom": "CAM_BOTTOM",
    "rear": "CAM_REAR",
    "right": "CAM_RIGHT",
    "right_iso": "CAM_RIGHT_ISO",
    "center": "CAM_CENTER",
}

_RENDER_ENUM_NAMES = {
    "wire": "GEOM_DRAW_WIRE",
    "hidden": "GEOM_DRAW_HIDDEN",
    "shade": "GEOM_DRAW_SHADE",
    "texture": "GEOM_DRAW_TEXTURE",
}

@dataclass(frozen=True)
class ViewSpec:
    """One output panel.

    ``rotate_deg`` is a counter-clockwise image rotation applied after OpenVSP
    captures the standard camera view. Multiples of 90 degrees are accepted.
    """

    view: str
    rotate_deg: int = 0

    def __post_init__(self) -> None:
        if not isinstance(self.view, str):
            raise TypeError("view must be a string")
        view = self.view.strip().lower()
        if view not in _VIEW_ENUM_NAMES:
            supported = ", ".join(_VIEW_ENUM_NAMES)
            raise ValueError(f"Unsupported view {self.view!r}. Supported: {supported}")
        if isinstance(self.rotate_deg, bool) or not isinstance(self.rotate_deg, int):
            raise TypeError("rotate_deg must be an integer multiple of 90")
        if self.rotate_deg % 90:
            raise ValueError("rotate_deg must be an integer multiple of 90")
        object.__setattr__(self, "view", view)
        object.__setattr__(self, "rotate_deg", self.rotate_deg % 360)

# Projection-style arrangement requested for this project:
#   top view (nose down)       left-side view (nose down)
#   front view                 left isometric view
ORTHOGRAPHIC_FOUR_VIEW_LAYOUT = (
    (ViewSpec("top", 90), ViewSpec("left", 90)),
    (ViewSpec("front"), ViewSpec("left_iso")),
)

def capture_vsp3_views(
    vsp3_path: str | Path,
    output_path: str | Path,
    *,
    size: tuple[int, int] = (1600, 1200),
    view_grid: Sequence[Sequence[ViewSpec | str]] = ORTHOGRAPHIC_FOUR_VIEW_LAYOUT,
    render_mode: str = "preserve",
    fit: bool = True,
    fit_to_content: bool = False,
    content_padding_px: int = 12,
    transparent_background: bool = False,
    show_axis: bool = False,
    show_borders: bool = False,
    overwrite: bool = False,
    python_executable: str | Path | None = None,
    timeout_s: float | None = 120.0,
) -> Path:
    """Capture one ``.vsp3`` model as a rectangular PNG view grid.

    ``render_mode`` accepts ``preserve``, ``wire``, ``hidden``, ``shade``, or
    ``texture``. ``preserve`` keeps every Geom's saved draw type. Other modes
    change only Geoms that are already visible; Geoms saved as
    ``GEOM_DRAW_NONE`` remain hidden.

    When ``fit_to_content`` is true, every rendered panel is cropped to its
    non-transparent pixels after rotation, resized without distortion, and
    centered in its output cell with ``content_padding_px`` pixels of padding.
    """

    model_path = _resolve_vsp3_paths([vsp3_path])[0]
    destination = Path(output_path).expanduser().resolve()
    if destination.suffix.lower() != ".png":
        raise ValueError("output_path must use the .png extension")
    if destination.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {destination}. Pass overwrite=True to replace it."
        )

    width, height = _validate_size(size)
    normalized_grid = _normalize_view_grid(view_grid)
    render_mode = _validate_render_mode(render_mode)
    _validate_content_fit(
        fit_to_content,
        content_padding_px,
        show_borders,
        size=(width, height),
        view_grid=normalized_grid,
    )
    if timeout_s is not None and timeout_s <= 0:
        raise ValueError("timeout_s must be greater than zero or None")

    destination.parent.mkdir(parents=True, exist_ok=True)
    _capture_vsp3_files(
        [model_path],
        [destination],
        size=(width, height),
        view_grid=normalized_grid,
        render_mode=render_mode,
        fit=fit,
        fit_to_content=fit_to_content,
        content_padding_px=content_padding_px,
        transparent_background=transparent_background,
        show_axis=show_axis,
        show_borders=show_borders,
        python_executable=python_executable,
        timeout_s=timeout_s,
    )
    return destination

def create_vsp3_animation(
    vsp3_paths: Sequence[str | Path],
    output_path: str | Path,
    *,
    size: tuple[int, int] = (1600, 1200),
    view_grid: Sequence[Sequence[ViewSpec | str]] = ORTHOGRAPHIC_FOUR_VIEW_LAYOUT,
    render_mode: str = "preserve",
    fps: float = 5.0,
    frames_dir: str | Path | None = None,
    keep_frames: bool = True,
    fit: bool = True,
    fit_to_content: bool = False,
    content_padding_px: int = 12,
    transparent_background: bool = False,
    show_axis: bool = False,
    show_borders: bool = False,
    overwrite: bool = False,
    python_executable: str | Path | None = None,
    timeout_s: float | None = None,
) -> dict[str, object]:
    """Capture an ordered list of ``.vsp3`` files and create GIF or MP4 output.

    Input order and duplicate paths are preserved. The output format is chosen
    from ``output_path`` (``.gif`` or ``.mp4``). When ``keep_frames`` is true,
    numbered PNG files and a capture manifest are retained. Repeating the same
    request reuses valid frames and captures only missing or damaged ones.
    ``fit_to_content`` is applied independently to every view in every frame.
    """

    model_paths = _resolve_vsp3_paths(vsp3_paths)
    destination = Path(output_path).expanduser().resolve()
    output_format = destination.suffix.lower()
    if output_format not in {".gif", ".mp4"}:
        raise ValueError("output_path must use the .gif or .mp4 extension")
    if output_format == ".mp4" and transparent_background:
        raise ValueError("H.264 MP4 output does not support a transparent background")
    if destination.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {destination}. Pass overwrite=True to replace it."
        )
    if fps <= 0:
        raise ValueError("fps must be greater than zero")
    if frames_dir is not None and not keep_frames:
        raise ValueError("frames_dir cannot be used when keep_frames=False")
    if timeout_s is not None and timeout_s <= 0:
        raise ValueError("timeout_s must be greater than zero or None")

    width, height = _validate_size(size)
    normalized_grid = _normalize_view_grid(view_grid)
    render_mode = _validate_render_mode(render_mode)
    _validate_content_fit(
        fit_to_content,
        content_padding_px,
        show_borders,
        size=(width, height),
        view_grid=normalized_grid,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)

    try:
        import imageio_ffmpeg
    except ImportError as exc:
        raise ImportError(
            "imageio-ffmpeg is required for GIF/MP4 encoding. "
            "Install it with: python -m pip install imageio-ffmpeg"
        ) from exc
    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()

    temporary_frames = None
    if keep_frames:
        frame_root = (
            Path(frames_dir).expanduser().resolve()
            if frames_dir is not None
            else destination.with_name(f"{destination.stem}_frames")
        )
        frame_root.mkdir(parents=True, exist_ok=True)
    else:
        temporary_frames = tempfile.TemporaryDirectory(prefix="openvsp-animation-")
        frame_root = Path(temporary_frames.name)

    try:
        frame_paths = [
            frame_root / f"frame_{index:06d}.png"
            for index in range(1, len(model_paths) + 1)
        ]
        capture_manifest = {
            "schema_version": 3,
            "models": [
                {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "modified_ns": path.stat().st_mtime_ns,
                }
                for path in model_paths
            ],
            "size": [width, height],
            "view_grid": [
                [
                    {"view": spec.view, "rotate_deg": spec.rotate_deg}
                    for spec in row
                ]
                for row in normalized_grid
            ],
            "render_mode": render_mode,
            "fit": fit,
            "fit_to_content": fit_to_content,
            "content_padding_px": content_padding_px,
            "transparent_background": transparent_background,
            "show_axis": show_axis,
            "show_borders": show_borders,
        }
        manifest_path = frame_root / "capture_manifest.json"

        if keep_frames and manifest_path.exists():
            previous_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if previous_manifest != capture_manifest:
                if not overwrite:
                    raise ValueError(
                        f"Existing frames use different capture settings: {frame_root}. "
                        "Use another frames_dir or pass overwrite=True."
                    )
                for old_frame in frame_root.glob("frame_*.png"):
                    old_frame.unlink()
        elif keep_frames:
            old_frames = list(frame_root.glob("frame_*.png"))
            if old_frames and not overwrite:
                raise ValueError(
                    f"Frames exist without a matching manifest: {frame_root}. "
                    "Use another frames_dir or pass overwrite=True."
                )
            for old_frame in old_frames:
                old_frame.unlink()

        temporary_manifest = manifest_path.with_suffix(".json.tmp")
        temporary_manifest.write_text(
            json.dumps(capture_manifest, indent=2), encoding="utf-8"
        )
        temporary_manifest.replace(manifest_path)

        missing_models = []
        missing_frames = []
        for model_path, frame_path in zip(model_paths, frame_paths):
            reusable = False
            if frame_path.is_file() and frame_path.stat().st_size > 0:
                try:
                    with Image.open(frame_path) as frame:
                        reusable = frame.size == (width, height) and frame.format == "PNG"
                except OSError:
                    reusable = False
            if not reusable:
                missing_models.append(model_path)
                missing_frames.append(frame_path)

        if missing_models:
            _capture_vsp3_files(
                missing_models,
                missing_frames,
                size=(width, height),
                view_grid=normalized_grid,
                render_mode=render_mode,
                fit=fit,
                fit_to_content=fit_to_content,
                content_padding_px=content_padding_px,
                transparent_background=transparent_background,
                show_axis=show_axis,
                show_borders=show_borders,
                python_executable=python_executable,
                timeout_s=timeout_s,
            )

        input_pattern = frame_root / "frame_%06d.png"
        fps_text = f"{float(fps):g}"
        commands: list[list[str]] = []

        if output_format == ".mp4":
            commands.append(
                [
                    ffmpeg,
                    "-y",
                    "-framerate",
                    fps_text,
                    "-start_number",
                    "1",
                    "-i",
                    str(input_pattern),
                    "-frames:v",
                    str(len(frame_paths)),
                    "-vf",
                    "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                    "-c:v",
                    "libx264",
                    "-pix_fmt",
                    "yuv420p",
                    "-movflags",
                    "+faststart",
                    str(destination),
                ]
            )
        else:
            palette_path = frame_root / "gif_palette.png"
            commands.extend(
                [
                    [
                        ffmpeg,
                        "-y",
                        "-framerate",
                        fps_text,
                        "-start_number",
                        "1",
                        "-i",
                        str(input_pattern),
                        "-frames:v",
                        str(len(frame_paths)),
                        "-vf",
                        f"fps={fps_text},palettegen=stats_mode=diff",
                        str(palette_path),
                    ],
                    [
                        ffmpeg,
                        "-y",
                        "-framerate",
                        fps_text,
                        "-start_number",
                        "1",
                        "-i",
                        str(input_pattern),
                        "-i",
                        str(palette_path),
                        "-frames:v",
                        str(len(frame_paths)),
                        "-lavfi",
                        f"fps={fps_text}[frame];[frame][1:v]paletteuse=dither=sierra2_4a",
                        "-loop",
                        "0",
                        str(destination),
                    ],
                ]
            )

        for command in commands:
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                details = result.stderr.strip() or result.stdout.strip()
                raise RuntimeError(f"FFmpeg animation encoding failed.\n{details}")

        if not destination.is_file() or destination.stat().st_size == 0:
            raise RuntimeError(f"Animation was not created: {destination}")

        return {
            "animation_path": destination,
            "frames_dir": frame_root if keep_frames else None,
            "frame_paths": frame_paths if keep_frames else [],
            "frame_count": len(frame_paths),
            "captured_frame_count": len(missing_frames),
            "reused_frame_count": len(frame_paths) - len(missing_frames),
        }
    finally:
        if temporary_frames is not None:
            temporary_frames.cleanup()

def _resolve_vsp3_paths(paths: Sequence[str | Path]) -> list[Path]:
    if isinstance(paths, (str, bytes, Path)):
        raise TypeError("vsp3_paths must be a non-empty sequence of paths")
    resolved = [Path(path).expanduser().resolve() for path in paths]
    if not resolved:
        raise ValueError("vsp3_paths must contain at least one model")
    for path in resolved:
        if not path.is_file():
            raise FileNotFoundError(f"OpenVSP model does not exist: {path}")
        if path.suffix.lower() != ".vsp3":
            raise ValueError(f"OpenVSP model must have a .vsp3 extension: {path}")
    return resolved

def _validate_size(size: tuple[int, int]) -> tuple[int, int]:
    if (
        not isinstance(size, Sequence)
        or isinstance(size, (str, bytes))
        or len(size) != 2
    ):
        raise TypeError("size must be a two-item (width, height) sequence")
    width, height = size
    if not isinstance(width, int) or isinstance(width, bool):
        raise TypeError("size width must be an integer")
    if not isinstance(height, int) or isinstance(height, bool):
        raise TypeError("size height must be an integer")
    if width <= 0 or height <= 0:
        raise ValueError("size width and height must be greater than zero")
    return width, height

def _normalize_view_grid(
    view_grid: Sequence[Sequence[ViewSpec | str]],
) -> tuple[tuple[ViewSpec, ...], ...]:
    if isinstance(view_grid, (str, bytes)) or not isinstance(view_grid, Sequence):
        raise TypeError("view_grid must be a non-empty rectangular sequence of rows")
    if not view_grid:
        raise ValueError("view_grid must contain at least one row")

    normalized_rows = []
    expected_columns = None
    for row in view_grid:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence):
            raise TypeError("every view_grid row must be a sequence")
        if not row:
            raise ValueError("view_grid rows must not be empty")
        normalized_row = tuple(
            item if isinstance(item, ViewSpec) else ViewSpec(item) for item in row
        )
        if expected_columns is None:
            expected_columns = len(normalized_row)
        elif len(normalized_row) != expected_columns:
            raise ValueError("view_grid must be rectangular")
        normalized_rows.append(normalized_row)
    return tuple(normalized_rows)

def _validate_render_mode(render_mode: str) -> str:
    if not isinstance(render_mode, str):
        raise TypeError("render_mode must be a string")
    render_mode = render_mode.strip().lower()
    if render_mode != "preserve" and render_mode not in _RENDER_ENUM_NAMES:
        supported = ", ".join(["preserve", *_RENDER_ENUM_NAMES])
        raise ValueError(f"Unsupported render_mode {render_mode!r}. Supported: {supported}")
    return render_mode

def _validate_content_fit(
    fit_to_content: bool,
    content_padding_px: int,
    show_borders: bool,
    *,
    size: tuple[int, int],
    view_grid: Sequence[Sequence[ViewSpec]],
) -> None:
    if not isinstance(fit_to_content, bool):
        raise TypeError("fit_to_content must be a bool")
    if isinstance(content_padding_px, bool) or not isinstance(
        content_padding_px, int
    ):
        raise TypeError("content_padding_px must be an integer")
    if content_padding_px < 0:
        raise ValueError("content_padding_px must be nonnegative")
    if fit_to_content and show_borders:
        raise ValueError(
            "fit_to_content=True cannot be combined with show_borders=True"
        )
    if fit_to_content:
        cell_width = min(_split_extent(size[0], len(view_grid[0])))
        cell_height = min(_split_extent(size[1], len(view_grid)))
        if 2 * content_padding_px >= min(cell_width, cell_height):
            raise ValueError(
                f"content_padding_px={content_padding_px} leaves no room in "
                f"the smallest {cell_width}x{cell_height} output cell"
            )

def _split_extent(total: int, count: int) -> list[int]:
    quotient, remainder = divmod(total, count)
    if quotient == 0:
        raise ValueError(f"Output extent {total} px is too small for {count} grid cells")
    return [quotient + (1 if index < remainder else 0) for index in range(count)]

def _capture_vsp3_files(
    model_paths: Sequence[Path],
    output_paths: Sequence[Path],
    *,
    size: tuple[int, int],
    view_grid: Sequence[Sequence[ViewSpec]],
    render_mode: str,
    fit: bool,
    fit_to_content: bool,
    content_padding_px: int,
    transparent_background: bool,
    show_axis: bool,
    show_borders: bool,
    python_executable: str | Path | None,
    timeout_s: float | None,
) -> None:
    if len(model_paths) != len(output_paths):
        raise ValueError("model_paths and output_paths must have the same length")
    for output_path in output_paths:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    request = {
        "frames": [
            {"model_path": str(model), "output_path": str(output)}
            for model, output in zip(model_paths, output_paths)
        ],
        "size": list(size),
        "view_grid": [
            [
                {"view": spec.view, "rotate_deg": spec.rotate_deg}
                for spec in row
            ]
            for row in view_grid
        ],
        "render_mode": render_mode,
        "fit": fit,
        "fit_to_content": fit_to_content,
        "content_padding_px": content_padding_px,
        "transparent_background": transparent_background,
        "show_axis": show_axis,
        "show_borders": show_borders,
    }

    with tempfile.TemporaryDirectory(prefix="openvsp-view-request-") as temp_dir:
        request_path = Path(temp_dir) / "request.json"
        request_path.write_text(json.dumps(request, indent=2), encoding="utf-8")
        interpreter = Path(python_executable or sys.executable).expanduser().resolve()
        if not interpreter.is_file():
            raise FileNotFoundError(f"Python executable does not exist: {interpreter}")
        command = [
            str(interpreter),
            str(Path(__file__).resolve()),
            "--worker",
            str(request_path),
        ]
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout_s,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            duration = f"{timeout_s:g} seconds" if timeout_s is not None else "timeout"
            raise TimeoutError(f"OpenVSP screenshot worker exceeded {duration}") from exc
        if result.returncode != 0:
            details = result.stderr.strip() or result.stdout.strip() or "No worker output."
            raise RuntimeError(
                "OpenVSP screenshot worker failed. Ensure that this Python environment "
                "contains a graphics-capable OpenVSP build and can open the GUI.\n"
                f"{details}"
            )

    for output_path in output_paths:
        if not output_path.is_file() or output_path.stat().st_size == 0:
            raise RuntimeError(f"Screenshot was not created: {output_path}")

def _compose_panels(
    panel_paths: Sequence[Sequence[Path]],
    view_grid: Sequence[Sequence[ViewSpec]],
    destination: Path,
    *,
    size: tuple[int, int],
    fit_to_content: bool,
    content_padding_px: int,
    transparent_background: bool,
) -> None:
    width, height = size
    column_widths = _split_extent(width, len(view_grid[0]))
    row_heights = _split_extent(height, len(view_grid))
    canvas_color = (0, 0, 0, 0) if transparent_background else (255, 255, 255, 255)
    canvas = Image.new("RGBA", (width, height), canvas_color)

    y = 0
    for row_index, (row_paths, row_specs) in enumerate(zip(panel_paths, view_grid)):
        x = 0
        for column_index, (panel_path, spec) in enumerate(zip(row_paths, row_specs)):
            if not panel_path.is_file() or panel_path.stat().st_size == 0:
                raise RuntimeError(f"OpenVSP did not create panel: {panel_path}")
            with Image.open(panel_path) as opened_image:
                panel = opened_image.convert("RGBA")
            if spec.rotate_deg:
                panel = panel.rotate(spec.rotate_deg, expand=True)

            panel_width = column_widths[column_index]
            panel_height = row_heights[row_index]
            if fit_to_content:
                panel = _fit_panel_to_cell(
                    panel,
                    (panel_width, panel_height),
                    content_padding_px,
                    panel_label=(
                        f"view {spec.view!r} at row {row_index}, "
                        f"column {column_index}"
                    ),
                )
            else:
                if panel.width > panel_width or panel.height > panel_height:
                    raise RuntimeError(
                        "OpenVSP panel is larger than its output cell after rotation: "
                        f"{panel.size} does not fit {(panel_width, panel_height)}"
                    )
                cell = Image.new("RGBA", (panel_width, panel_height), (0, 0, 0, 0))
                cell.alpha_composite(
                    panel,
                    (
                        (panel_width - panel.width) // 2,
                        (panel_height - panel.height) // 2,
                    ),
                )
                panel = cell

            canvas.alpha_composite(panel, (x, y))
            x += panel_width
        y += row_heights[row_index]

    if transparent_background:
        canvas.save(destination, format="PNG")
    else:
        canvas.convert("RGB").save(destination, format="PNG")

def _fit_panel_to_cell(
    panel: Image.Image,
    cell_size: tuple[int, int],
    padding_px: int,
    *,
    panel_label: str = "panel",
) -> Image.Image:
    cell_width, cell_height = cell_size
    available_width = cell_width - 2 * padding_px
    available_height = cell_height - 2 * padding_px
    if available_width <= 0 or available_height <= 0:
        raise ValueError(
            f"content_padding_px={padding_px} leaves no room in "
            f"the {cell_width}x{cell_height} output cell"
        )

    rgba_panel = panel.convert("RGBA")
    content_bbox = rgba_panel.getchannel("A").getbbox()
    if content_bbox is None:
        raise RuntimeError(f"OpenVSP rendered no visible content for {panel_label}")

    content = rgba_panel.crop(content_bbox)
    scale = min(
        available_width / content.width,
        available_height / content.height,
    )
    resized_size = (
        max(1, round(content.width * scale)),
        max(1, round(content.height * scale)),
    )
    if content.size != resized_size:
        content = content.resize(resized_size, Image.Resampling.LANCZOS)

    cell = Image.new("RGBA", (cell_width, cell_height), (0, 0, 0, 0))
    cell.alpha_composite(
        content,
        (
            (cell_width - content.width) // 2,
            (cell_height - content.height) // 2,
        ),
    )
    return cell

def _capture_worker(request_path: Path) -> None:
    request = json.loads(request_path.read_text(encoding="utf-8"))

    import openvsp_config

    openvsp_config.LOAD_GRAPHICS = True
    openvsp_config.LOAD_FACADE = True

    # OpenVSP 3.50.4 imports facade_server while loading the facade client, and
    # facade_server reads sys.argv[1] as a port number. Hide this worker's
    # ``--worker`` argument during that import so it is not parsed as the port.
    worker_argv = sys.argv
    try:
        sys.argv = [sys.argv[0], "-1"]
        import openvsp as vsp
    finally:
        sys.argv = worker_argv

    if not vsp.IsGUIBuild():
        raise RuntimeError("The installed OpenVSP Python API is not a GUI build.")

    width, height = request["size"]
    view_grid = tuple(
        tuple(ViewSpec(**spec) for spec in row) for row in request["view_grid"]
    )
    column_widths = _split_extent(width, len(view_grid[0]))
    row_heights = _split_extent(height, len(view_grid))

    gui_started = False
    try:
        vsp.InitGUI()
        vsp.StartGUI()
        gui_started = True
        deadline = time.monotonic() + 10.0
        while not vsp.IsEventLoopRunning():
            if time.monotonic() >= deadline:
                raise RuntimeError("OpenVSP GUI event loop did not start within 10 seconds.")
            time.sleep(0.05)

        vsp.SetWindowLayout(1, 1)
        vsp.SetViewAxis(request["show_axis"])
        vsp.SetShowBorders(request["show_borders"])

        with tempfile.TemporaryDirectory(prefix="openvsp-view-panels-") as temp_dir:
            panel_root = Path(temp_dir)
            for frame_index, frame in enumerate(request["frames"], start=1):
                vsp.ClearVSPModel()
                vsp.ReadVSPFile(frame["model_path"])
                vsp.Update()
                vehicle_id = vsp.GetVehicleID()

                if request["render_mode"] != "preserve":
                    draw_type = getattr(
                        vsp, _RENDER_ENUM_NAMES[request["render_mode"]]
                    )
                    for geom_id in vsp.GetGeomSetAtIndex(vsp.SET_SHOWN):
                        vsp.SetGeomDrawType(geom_id, draw_type)
                    vsp.Update()

                panel_paths = []
                for row_index, row in enumerate(view_grid):
                    path_row = []
                    for column_index, spec in enumerate(row):
                        panel_path = panel_root / (
                            f"frame_{frame_index:06d}_"
                            f"panel_{row_index}_{column_index}.png"
                        )
                        path_row.append(panel_path)

                        panel_width = column_widths[column_index]
                        panel_height = row_heights[row_index]
                        if spec.rotate_deg in {90, 270}:
                            capture_width, capture_height = panel_height, panel_width
                        else:
                            capture_width, capture_height = panel_width, panel_height

                        # OpenVSP 3.50 renders ScreenGrab with the current GL
                        # viewport projection. Match both sizes before fitting so
                        # the capture is not stretched, then rotate without resize.
                        vsp.SetParmVal(
                            vehicle_id,
                            "ViewportX",
                            "AdjustView",
                            capture_width,
                        )
                        vsp.SetParmVal(
                            vehicle_id,
                            "ViewportY",
                            "AdjustView",
                            capture_height,
                        )
                        vsp.UpdateGUI()

                        vsp.SetView(0, getattr(vsp, _VIEW_ENUM_NAMES[spec.view]))
                        if request["fit"]:
                            vsp.FitAllViews()
                        vsp.UpdateGUI()
                        capture_transparent_background = (
                            request["transparent_background"]
                            or request["fit_to_content"]
                        )
                        vsp.ScreenGrab(
                            str(panel_path),
                            capture_width,
                            capture_height,
                            capture_transparent_background,
                            False,
                        )
                    panel_paths.append(path_row)

                _compose_panels(
                    panel_paths,
                    view_grid,
                    Path(frame["output_path"]),
                    size=(width, height),
                    fit_to_content=request["fit_to_content"],
                    content_padding_px=request["content_padding_px"],
                    transparent_background=request["transparent_background"],
                )
    finally:
        if gui_started:
            vsp.StopGUI()

def _main() -> int:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--worker", type=Path)
    arguments = parser.parse_args()
    if arguments.worker is None:
        parser.error("This module is an internal worker; use the public functions.")
    _capture_worker(arguments.worker)
    return 0

if __name__ == "__main__":
    raise SystemExit(_main())
