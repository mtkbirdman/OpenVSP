from pathlib import Path
import sys

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.VSPAEROMesh import equalize_vspaero_tessellation

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


def test_equalize_rejects_missing_input(tmp_path):
    with pytest.raises(FileNotFoundError):
        equalize_vspaero_tessellation(
            input_vsp3_path=tmp_path / "missing.vsp3",
            output_vsp3_path=tmp_path / "out.vsp3",
            reference_wing_name="WingGeom",
            reference_tess_w=65,
        )


def test_equalize_rejects_invalid_settings(tmp_path):
    output_path = tmp_path / "out.vsp3"

    with pytest.raises(ValueError, match="reference_tess_w"):
        equalize_vspaero_tessellation(
            input_vsp3_path=MODEL_PATH,
            output_vsp3_path=output_path,
            reference_wing_name="WingGeom",
            reference_tess_w=0,
        )

    with pytest.raises(ValueError, match="tolerance"):
        equalize_vspaero_tessellation(
            input_vsp3_path=MODEL_PATH,
            output_vsp3_path=output_path,
            reference_wing_name="WingGeom",
            reference_tess_w=65,
            tolerance=1.0,
        )

    with pytest.raises(ValueError, match="max_iterations"):
        equalize_vspaero_tessellation(
            input_vsp3_path=MODEL_PATH,
            output_vsp3_path=output_path,
            reference_wing_name="WingGeom",
            reference_tess_w=65,
            max_iterations=0,
        )


def test_equalize_saves_selected_wing_and_active_end_cap(tmp_path):
    vsp = pytest.importorskip("openvsp")
    from src.util import find_container_parm, find_one_geom

    output_path = tmp_path / "G103A_equalized.vsp3"

    report = equalize_vspaero_tessellation(
        input_vsp3_path=MODEL_PATH,
        output_vsp3_path=output_path,
        reference_wing_name="WingGeom",
        reference_tess_w=65,
        geom_names=["WingGeom"],
        max_iterations=2,
    )

    assert output_path.is_file()
    assert set(report["geom_name"]) == {"WingGeom"}

    wing_row = report.loc[report["geom_name"] == "WingGeom"].iloc[0]
    assert wing_row["status"] == "adjusted"
    assert wing_row["before_cap_tess"] is not None
    assert wing_row["after_cap_tess"] is not None
    assert (
        wing_row["after_cap_min_u_edge_median"] is not None
        or wing_row["after_cap_max_u_edge_median"] is not None
    )
    assert wing_row["cap_target_ratio"] is not None
    assert wing_row["u_target_ratio"] > 0.0
    assert wing_row["w_target_ratio"] > 0.0

    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(output_path))
    vsp.Update()

    wing_id = find_one_geom(vsp, "WingGeom")
    cap_tess_id = find_container_parm(vsp, wing_id, "CapUMinTess")
    assert cap_tess_id
    assert int(round(vsp.GetParmVal(cap_tess_id))) == int(wing_row["after_cap_tess"])
