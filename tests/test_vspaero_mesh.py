from pathlib import Path
import sys

import pytest

vsp = pytest.importorskip("openvsp")

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.VSPAEROMesh import equalize_vspaero_tessellation
from src.util import find_container_parm, find_one_geom

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


def test_equalize_includes_active_wing_end_cap_tessellation(tmp_path):
    output_path = tmp_path / "G103A_equalized.vsp3"

    report = equalize_vspaero_tessellation(
        input_vsp3_path=MODEL_PATH,
        output_vsp3_path=output_path,
        reference_wing_name="WingGeom",
        reference_tess_w=65,
        geom_names=["WingGeom"],
        max_iterations=2,
    )

    wing_row = report.loc[report["geom_name"] == "WingGeom"].iloc[0]
    assert wing_row["before_cap_tess"] is not None
    assert wing_row["after_cap_tess"] is not None
    assert (
        wing_row["after_cap_min_u_edge_median"] is not None
        or wing_row["after_cap_max_u_edge_median"] is not None
    )
    assert wing_row["cap_target_ratio"] is not None

    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(output_path))
    vsp.Update()

    wing_id = find_one_geom(vsp, "WingGeom")
    cap_tess_id = find_container_parm(vsp, wing_id, "CapUMinTess")
    assert cap_tess_id
    assert int(round(vsp.GetParmVal(cap_tess_id))) == int(wing_row["after_cap_tess"])
