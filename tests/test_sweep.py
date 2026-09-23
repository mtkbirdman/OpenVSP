from pathlib import Path
import sys

import numpy as np
import pytest

vsp = pytest.importorskip("openvsp")

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import vsp_sweep
from src.util import results_dataframe

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


def test_sweep_returns_requested_polar_with_physical_lift_trend():
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    alpha = [-2.0, 2.0, 6.0]
    mach = [0.1]
    reynolds = [4.4e6]

    result_id = vsp_sweep(
        vsp=vsp,
        alpha=alpha,
        mach=mach,
        reynolds=reynolds,
        verbose=0,
    )

    assert result_id

    polar = results_dataframe(vsp, result_id, "VSPAERO_Polar")
    if polar.empty:
        polar = results_dataframe(vsp, result_id, "VSPAERO Polar")

    assert not polar.empty
    assert len(polar) == len(alpha)

    alpha_column = "Alpha" if "Alpha" in polar.columns else "AoA"
    required_columns = {alpha_column, "Mach", "Re_1e6", "CLtot", "CDi", "CMytot"}
    assert required_columns <= set(polar.columns)

    np.testing.assert_allclose(polar[alpha_column], alpha, atol=1e-10, rtol=0.0)
    np.testing.assert_allclose(polar["Mach"], mach[0], atol=1e-10, rtol=0.0)
    np.testing.assert_allclose(
        polar["Re_1e6"], reynolds[0] / 1e6, atol=1e-10, rtol=0.0
    )

    aerodynamic_columns = ["CLtot", "CDi", "CMytot"]
    assert np.isfinite(polar[aerodynamic_columns].to_numpy()).all()

    # In the attached-flow range used here, a conventional glider must have a
    # positive lift-curve slope. This catches broken alpha setup, result
    # extraction, or a severely corrupted aerodynamic model without pinning the
    # test to one OpenVSP/VSPAERO version's exact coefficient values.
    assert np.all(np.diff(polar["CLtot"].to_numpy()) > 0.0)

    # Induced drag should not become materially negative in this positive-lift
    # operating range. Keep a tiny numerical allowance rather than testing an
    # exact solver value.
    assert np.all(polar["CDi"].to_numpy() >= -1e-8)


def test_sweep_applies_explicit_thick_thin_set_override():
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    thin_geom_set = int(vsp.GetSetIndex("ThinGeom"))
    assert thin_geom_set != int(vsp.SET_NONE)
    assert list(vsp.GetGeomSetAtIndex(thin_geom_set))

    result_id = vsp_sweep(
        vsp=vsp,
        alpha=[2.0],
        mach=[0.1],
        reynolds=[4.4e6],
        verbose=0,
        thick_geom_set=vsp.SET_NONE,
        thin_geom_set=thin_geom_set,
    )

    assert result_id

    for analysis_name in ("VSPAEROComputeGeometry", "VSPAEROSweep"):
        assert int(vsp.GetIntAnalysisInput(analysis_name, "GeomSet")[0]) == int(
            vsp.SET_NONE
        )
        assert int(vsp.GetIntAnalysisInput(analysis_name, "ThinGeomSet")[0]) == thin_geom_set
