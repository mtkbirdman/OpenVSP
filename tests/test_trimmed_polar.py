from pathlib import Path
import sys

import numpy as np
import pytest

vsp = pytest.importorskip("openvsp")

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import vsp_trimmed_sweep

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


def test_trimmed_polar_converges_pitch_moment_and_restores_elevator():
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    group_names = [
        vsp.GetVSPAEROControlGroupName(i)
        for i in range(vsp.GetNumControlSurfaceGroups())
    ]
    assert group_names.count("ELEVATOR_GROUP") == 1
    elevator_group_index = group_names.index("ELEVATOR_GROUP")
    assert len(vsp.GetActiveCSNameVec(elevator_group_index)) == 2

    settings_id = vsp.FindContainer("VSPAEROSettings", 0)
    elevator_parm_id = vsp.FindParm(
        settings_id,
        "DeflectionAngle",
        f"ControlSurfaceGroup_{elevator_group_index}",
    )
    assert elevator_parm_id and str(elevator_parm_id).upper() != "NONE"
    original_elevator_deg = float(vsp.GetParmVal(elevator_parm_id))

    alpha = [0.0, 4.0, 8.0]
    cmy_tolerance = 5e-4
    elevator_bounds = (-20.0, 20.0)
    max_trim_evaluations = 30

    polar = vsp_trimmed_sweep(
        vsp=vsp,
        alpha_list=alpha,
        Weight=580,
        elevator_bounds=elevator_bounds,
        cmy_tolerance=cmy_tolerance,
        max_trim_evaluations=max_trim_evaluations,
        verbose=0,
    )

    required_columns = {
        "InputAlpha_deg",
        "InputMach",
        "InputReCref",
        "Elevator_deg",
        "CMyResidual",
        "CMyElevatorSlope_per_deg",
        "TrimFunctionEvaluations",
        "CL",
        "CMytot",
        "CMy",
        "Velocity",
    }
    assert required_columns <= set(polar.columns)
    assert len(polar) == len(alpha)
    np.testing.assert_allclose(polar["InputAlpha_deg"], alpha, rtol=0.0, atol=1e-10)

    numeric_columns = [
        "InputMach",
        "InputReCref",
        "Elevator_deg",
        "CMyResidual",
        "CMyElevatorSlope_per_deg",
        "CL",
        "CMytot",
        "CMy",
        "Velocity",
    ]
    assert np.isfinite(polar[numeric_columns].to_numpy()).all()
    assert np.all(polar["InputMach"].to_numpy() > 0.0)
    assert np.all(polar["InputReCref"].to_numpy() > 0.0)

    # The defining contract of a pitch-trimmed polar is zero pitching moment
    # within the requested numerical tolerance.
    assert np.all(np.abs(polar["CMyResidual"].to_numpy()) <= cmy_tolerance)
    assert np.all(np.abs(polar["CMytot"].to_numpy()) <= cmy_tolerance)
    np.testing.assert_allclose(
        polar["CMyResidual"], polar["CMytot"], rtol=0.0, atol=1e-12
    )
    np.testing.assert_allclose(polar["CMy"], polar["CMytot"], rtol=0.0, atol=1e-12)

    assert np.all(polar["Elevator_deg"].to_numpy() >= elevator_bounds[0])
    assert np.all(polar["Elevator_deg"].to_numpy() <= elevator_bounds[1])
    assert np.all(polar["TrimFunctionEvaluations"].to_numpy() >= 1)
    assert np.all(
        polar["TrimFunctionEvaluations"].to_numpy() <= max_trim_evaluations
    )

    # With this G103A model and its saved elevator gain convention, increasing
    # elevator deflection produces a negative pitching-moment derivative. A
    # sign reversal here usually means the control group/gain convention or
    # result interpretation has been broken.
    assert np.all(polar["CMyElevatorSlope_per_deg"].to_numpy() < 0.0)

    # These cases are deliberately within the attached-flow range. Lift should
    # rise with alpha, and the speed required to support fixed aircraft weight
    # should therefore fall.
    assert np.all(np.diff(polar["CL"].to_numpy()) > 0.0)
    assert np.all(np.diff(polar["Velocity"].to_numpy()) < 0.0)

    # vsp_trimmed_sweep is responsible for restoring the model state after the
    # continuation sweep; later tests must not inherit the final trim angle.
    assert vsp.GetParmVal(elevator_parm_id) == pytest.approx(original_elevator_deg)
