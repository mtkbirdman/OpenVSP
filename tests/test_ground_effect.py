from pathlib import Path
import sys

import numpy as np
import pytest

vsp = pytest.importorskip("openvsp")

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import make_CDo_correction, vsp_sweep_wig
from src.util import get_container_parm_value

MODEL_PATH = REPO_ROOT / "examples" / "models" / "SampleGlider" / "SampleGlider.vsp3"


def test_ground_effect_reduces_induced_drag_and_preserves_output_contract():
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    settings_id = vsp.FindContainer("VSPAEROSettings", 0)
    bref, bref_parm_id = get_container_parm_value(vsp, settings_id, "bref")
    assert bref_parm_id
    assert bref is not None and bref > 0.0

    height_ratios = np.array([0.10, 0.25, 1.00])
    heights = (height_ratios * bref).tolist()

    polar = vsp_sweep_wig(
        vsp=vsp,
        alpha_deg=[0.0],
        mach=[0.1],
        reynolds=[1e6],
        height=heights,
        verbose=0,
    )

    # vsp_sweep_wig adds exactly one out-of-ground-effect reference followed by
    # the requested finite heights in caller order. No artificial height such
    # as 999 m is needed to emulate the reference case.
    assert len(polar) == 1 + len(heights)
    required_columns = {
        "GroundEffectEnabled",
        "CGHeight",
        "CGHeight_bref",
        "InputAlpha_deg",
        "InputMach",
        "InputReCref",
        "CL",
        "CDi",
    }
    assert required_columns <= set(polar.columns)

    oge = polar.loc[~polar["GroundEffectEnabled"].astype(bool)]
    finite_height = polar.loc[polar["GroundEffectEnabled"].astype(bool)]

    assert len(oge) == 1
    assert len(finite_height) == len(heights)
    assert np.isnan(float(oge["CGHeight"].iloc[0]))
    assert np.isnan(float(oge["CGHeight_bref"].iloc[0]))

    np.testing.assert_allclose(finite_height["CGHeight"], heights, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(
        finite_height["CGHeight_bref"], height_ratios, rtol=0.0, atol=1e-10
    )
    np.testing.assert_allclose(polar["InputAlpha_deg"], 0.0, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(polar["InputMach"], 0.1, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(polar["InputReCref"], 1e6, rtol=0.0, atol=1e-6)
    assert np.isfinite(polar[["CL", "CDi"]].to_numpy()).all()

    finite_cdi = finite_height["CDi"].to_numpy(dtype=float)
    oge_cdi = float(oge["CDi"].iloc[0])

    # As the wing approaches the ground, downwash and induced drag decrease.
    # Therefore CDi must rise monotonically as h/b increases back toward the
    # out-of-ground-effect condition.
    assert np.all(np.diff(finite_cdi) > 0.0)
    assert finite_cdi[-1] < oge_cdi

    # At h/b = 1 the ground-effect correction should already be small. The
    # tolerance is intentionally loose so this remains a physics sanity check,
    # not a VSPAERO-version-specific golden-value test.
    assert abs(finite_cdi[-1] - oge_cdi) / oge_cdi < 0.05

    corrected = make_CDo_correction(
        vsp=vsp,
        trimmed_polar=polar.copy(),
        Weight=580,
        xTr=(0.5, 0.7),
        CDpCL=0.0016,
        thickness=0.19,
        interference_factor=1.14,
        altitude=0,
        dT=0,
    )

    # These are identities defined by make_CDo_correction itself, so they can
    # be tested tightly rather than with a loose aerodynamic tolerance.
    np.testing.assert_allclose(
        corrected["CDtot_corr"],
        corrected["CDo_corr"] + corrected["CDi"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        corrected["L_D_corr"],
        corrected["CL"] / corrected["CDtot_corr"],
        rtol=1e-12,
        atol=1e-12,
    )
    assert np.isfinite(
        corrected[["CDo_corr", "CDtot_corr", "L_D_corr", "Velocity"]].to_numpy()
    ).all()
