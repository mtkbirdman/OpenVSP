import builtins
import json
import os
import sys
import types
from pathlib import Path

import pytest
from PIL import Image

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src import OpenVSPViewCapture as capture


def test_view_spec_normalizes_name_and_rotation():
    spec = capture.ViewSpec(" LEFT ", -90)

    assert spec.view == "left"
    assert spec.rotate_deg == 270


@pytest.mark.parametrize(
    ("args", "exception_type"),
    [
        (("unknown", 0), ValueError),
        (("top", 45), ValueError),
        (("top", 90.0), TypeError),
    ],
)
def test_view_spec_rejects_invalid_values(args, exception_type):
    with pytest.raises(exception_type):
        capture.ViewSpec(*args)


def test_normalize_view_grid_accepts_names_and_rejects_ragged_rows():
    normalized = capture._normalize_view_grid(
        (("top", "left"), ("front", "left_iso"))
    )

    assert normalized[0][0] == capture.ViewSpec("top")
    assert normalized[1][1] == capture.ViewSpec("left_iso")

    with pytest.raises(ValueError, match="rectangular"):
        capture._normalize_view_grid((("top",), ("front", "left")))


@pytest.mark.parametrize("mode", ["preserve", "wire", "hidden", "shade", "texture"])
def test_render_modes(mode):
    assert capture._validate_render_mode(f" {mode.upper()} ") == mode


def test_split_extent_preserves_exact_total():
    assert capture._split_extent(10, 3) == [4, 3, 3]
    assert sum(capture._split_extent(1601, 2)) == 1601


def _make_models(tmp_path, names=("first", "second")):
    models = []
    for name in names:
        path = tmp_path / f"{name}.vsp3"
        path.write_text(name, encoding="utf-8")
        models.append(path)
    return models


def _fake_capture(calls):
    def capture_files(model_paths, output_paths, **settings):
        calls.append((list(model_paths), list(output_paths), settings))
        for index, output_path in enumerate(output_paths):
            color = ((index + 1) * 70 % 255, 80, 160)
            Image.new("RGB", settings["size"], color).save(output_path)

    return capture_files


def test_capture_uses_shared_batch_path_and_render_mode(tmp_path, monkeypatch):
    model_path = _make_models(tmp_path, ("aircraft",))[0]
    output_path = tmp_path / "views.png"
    calls = []
    monkeypatch.setattr(capture, "_capture_vsp3_files", _fake_capture(calls))

    result = capture.capture_vsp3_views(
        model_path,
        output_path,
        size=(101, 79),
        render_mode="hidden",
    )

    assert result == output_path.resolve()
    assert calls[0][0] == [model_path.resolve()]
    assert calls[0][2]["render_mode"] == "hidden"
    with Image.open(result) as image:
        assert image.size == (101, 79)
        assert image.format == "PNG"


@pytest.mark.parametrize("extension", [".gif", ".mp4"])
def test_animation_preserves_order_and_encodes_requested_format(
    tmp_path, monkeypatch, extension
):
    first, second = _make_models(tmp_path)
    output_path = tmp_path / f"animation{extension}"
    calls = []
    monkeypatch.setattr(capture, "_capture_vsp3_files", _fake_capture(calls))

    result = capture.create_vsp3_animation(
        [second, first, second],
        output_path,
        size=(64, 48),
        render_mode="shade",
        fps=2,
        keep_frames=False,
    )

    assert result["animation_path"] == output_path.resolve()
    assert result["frame_count"] == 3
    assert result["captured_frame_count"] == 3
    assert result["frames_dir"] is None
    assert output_path.stat().st_size > 0
    assert calls[0][0] == [second.resolve(), first.resolve(), second.resolve()]
    assert calls[0][2]["render_mode"] == "shade"


def test_animation_reuses_matching_frames(tmp_path, monkeypatch):
    models = _make_models(tmp_path)
    frames_dir = tmp_path / "frames"
    calls = []
    monkeypatch.setattr(capture, "_capture_vsp3_files", _fake_capture(calls))

    first = capture.create_vsp3_animation(
        models,
        tmp_path / "first.gif",
        size=(64, 48),
        frames_dir=frames_dir,
        fps=2,
    )
    second = capture.create_vsp3_animation(
        models,
        tmp_path / "second.gif",
        size=(64, 48),
        frames_dir=frames_dir,
        fps=4,
    )

    assert first["captured_frame_count"] == 2
    assert second["captured_frame_count"] == 0
    assert second["reused_frame_count"] == 2
    assert len(calls) == 1


def test_animation_rejects_frames_from_different_capture_settings(
    tmp_path, monkeypatch
):
    models = _make_models(tmp_path)
    frames_dir = tmp_path / "frames"
    monkeypatch.setattr(capture, "_capture_vsp3_files", _fake_capture([]))
    capture.create_vsp3_animation(
        models,
        tmp_path / "wire.gif",
        size=(64, 48),
        frames_dir=frames_dir,
        render_mode="wire",
    )

    with pytest.raises(ValueError, match="different capture settings"):
        capture.create_vsp3_animation(
            models,
            tmp_path / "shade.gif",
            size=(64, 48),
            frames_dir=frames_dir,
            render_mode="shade",
        )


def test_capture_validates_paths_before_starting_worker(tmp_path, monkeypatch):
    worker_called = False

    def fake_capture(*args, **kwargs):
        nonlocal worker_called
        worker_called = True

    monkeypatch.setattr(capture, "_capture_vsp3_files", fake_capture)

    with pytest.raises(FileNotFoundError):
        capture.capture_vsp3_views(tmp_path / "missing.vsp3", tmp_path / "out.png")
    assert worker_called is False


def test_worker_hides_cli_arguments_from_openvsp_3504_facade(tmp_path, monkeypatch):
    output_path = tmp_path / "captured.png"
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "frames": [
                    {"model_path": "model.vsp3", "output_path": str(output_path)}
                ],
                "size": [64, 48],
                "view_grid": [[{"view": "top", "rotate_deg": 0}]],
                "render_mode": "preserve",
                "fit": True,
                "transparent_background": False,
                "autocrop": False,
                "show_axis": False,
                "show_borders": False,
            }
        ),
        encoding="utf-8",
    )

    calls = []
    fake_vsp = types.SimpleNamespace(
        CAM_TOP=0,
        IsGUIBuild=lambda: True,
        InitGUI=lambda: calls.append("init"),
        StartGUI=lambda: calls.append("start"),
        IsEventLoopRunning=lambda: True,
        SetWindowLayout=lambda *_: None,
        SetViewAxis=lambda *_: None,
        SetShowBorders=lambda *_: None,
        ClearVSPModel=lambda: None,
        ReadVSPFile=lambda *_: None,
        Update=lambda: None,
        SetView=lambda *_: None,
        FitAllViews=lambda: None,
        UpdateGUI=lambda: None,
        ScreenGrab=lambda path, width, height, *_: Image.new(
            "RGB", (width, height), "white"
        ).save(path),
        StopGUI=lambda: calls.append("stop"),
    )
    fake_config = types.SimpleNamespace()
    monkeypatch.setitem(sys.modules, "openvsp_config", fake_config)
    original_import = builtins.__import__
    imported_argv = []

    def import_openvsp(name, *args, **kwargs):
        if name == "openvsp":
            imported_argv.append(sys.argv.copy())
            return fake_vsp
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_openvsp)
    monkeypatch.setattr(sys, "argv", ["OpenVSPViewCapture.py", "--worker", "request"])

    capture._capture_worker(request_path)

    assert imported_argv == [["OpenVSPViewCapture.py", "-1"]]
    assert sys.argv == ["OpenVSPViewCapture.py", "--worker", "request"]
    assert calls == ["init", "start", "stop"]
    assert output_path.is_file()


@pytest.mark.slow
@pytest.mark.skipif(
    os.environ.get("RUN_OPENVSP_GUI_TESTS") != "1",
    reason="Set RUN_OPENVSP_GUI_TESTS=1 to run the interactive OpenVSP GUI test.",
)
def test_capture_g103a_with_openvsp_gui(tmp_path):
    model_path = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"
    output_path = tmp_path / "G103A_four_views.png"

    result = capture.capture_vsp3_views(
        model_path,
        output_path,
        size=(800, 600),
        render_mode="hidden",
    )

    with Image.open(result) as image:
        assert image.size == (800, 600)
        assert image.format == "PNG"
