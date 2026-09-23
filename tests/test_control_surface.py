from pathlib import Path
import sys

import numpy as np
import pytest

vsp = pytest.importorskip("openvsp")

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import vsp_sweep
from src.util import find_container_parm, results_dataframe, set_control_surface

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


def test_aileron_deflection_reverses_rolling_moment():
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    # Keep the other controls explicitly neutral so this test isolates the
    # aileron group rather than depending on whatever deflections were last
    # saved in the model.
    set_control_surface(
        vsp,
        geom_name="HTailGeom",
        deflection=0.0,
        cs_group_name="ELEVATOR_GROUP",
        gains=(1, -1),
        verbose=0,
    )
    set_control_surface(
        vsp,
        geom_name="VTailGeom",
        deflection=0.0,
        cs_group_name="RUDDER_GROUP",
        gains=(-1,),
        verbose=0,
    )

    def run_aileron_case(deflection_deg):
        deflection_parm_id = set_control_surface(
            vsp,
            geom_name="WingGeom",
            deflection=deflection_deg,
            cs_group_name="AILERON_GROUP",
            gains=(1, 1),
            verbose=0,
        )
        assert vsp.GetParmVal(deflection_parm_id) == pytest.approx(deflection_deg)

        result_id = vsp_sweep(
            vsp=vsp,
            alpha=[2.0],
            mach=[0.1],
            reynolds=[4.4e6],
            verbose=0,
        )
        assert result_id

        polar = results_dataframe(vsp, result_id, "VSPAERO_Polar")
        if polar.empty:
            polar = results_dataframe(vsp, result_id, "VSPAERO Polar")

        assert len(polar) == 1
        assert {"CMxtot", "CLtot"} <= set(polar.columns)
        assert np.isfinite(polar[["CMxtot", "CLtot"]].to_numpy()).all()
        return float(polar["CMxtot"].iloc[0]), float(polar["CLtot"].iloc[0])

    neutral_roll, _ = run_aileron_case(0.0)

    # Verify that the existing OpenVSP control group still contains the two
    # mirrored aileron surfaces and that set_control_surface kept both gains at
    # the intended value.
    group_names = [
        vsp.GetVSPAEROControlGroupName(i)
        for i in range(vsp.GetNumControlSurfaceGroups())
    ]
    assert group_names.count("AILERON_GROUP") == 1
    group_index = group_names.index("AILERON_GROUP")
    assert len(vsp.GetActiveCSNameVec(group_index)) == 2

    wing_id = vsp.FindGeomsWithName("WingGeom")[0]
    aileron_subsurf_id = vsp.GetSubSurf(wing_id, 0)
    settings_id = vsp.FindContainer("VSPAEROSettings", 0)
    for reflected_index in (0, 1):
        gain_parm_id = find_container_parm(
            vsp,
            settings_id,
            f"Surf_{aileron_subsurf_id}_{reflected_index}_Gain",
        )
        assert gain_parm_id
        assert vsp.GetParmVal(gain_parm_id) == pytest.approx(1.0)

    positive_roll, positive_cl = run_aileron_case(10.0)
    negative_roll, negative_cl = run_aileron_case(-10.0)

    # Reversing an antisymmetric aileron input must reverse the rolling moment.
    assert positive_roll * negative_roll < 0.0

    roll_effect = max(abs(positive_roll), abs(negative_roll))
    assert roll_effect > 1e-4

    # The neutral case should be small compared with the controlled response,
    # and the two opposite deflections should be approximately antisymmetric.
    assert abs(neutral_roll) <= 0.10 * roll_effect + 1e-4
    assert abs(positive_roll + negative_roll) <= 0.20 * roll_effect + 1e-4

    # Opposite aileron deflections primarily redistribute lift between the two
    # wings; they should not create a large change in total aircraft lift.
    lift_scale = max(abs(positive_cl), abs(negative_cl), 1.0)
    assert abs(positive_cl - negative_cl) <= 0.02 * lift_scale

    set_control_surface(
        vsp,
        geom_name="WingGeom",
        deflection=0.0,
        cs_group_name="AILERON_GROUP",
        gains=(1, 1),
        verbose=0,
    )
