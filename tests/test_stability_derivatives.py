from pathlib import Path
import sys

import numpy as np
import pytest

vsp = pytest.importorskip("openvsp")

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import (
    validate_vsp3_for_stability_derivatives,
    vsp_stability_derivatives,
)

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"
WAKE_NUM_NODES = 32


def test_stability_derivative_preflight_accepts_saved_g103a_configuration():
    report = validate_vsp3_for_stability_derivatives(MODEL_PATH, verbose=0)

    assert report["passed"], report["errors"]
    assert report["errors"] == []

    selected_geom_count = sum(
        len(report["set_summary"][name]["geom_ids"])
        for name in ("GeomSet", "ThinGeomSet")
    )
    assert selected_geom_count > 0

    settings = report["vspaero_settings_summary"]
    for name in ("Sref", "bref", "cref"):
        assert np.isfinite(settings[name])
        assert settings[name] > 0.0
    for name in ("Xcg", "Ycg", "Zcg"):
        assert np.isfinite(settings[name])

    # Full stability derivatives require the complete aircraft rather than an
    # XZ-symmetry half model.  The preflight validator explicitly guarantees
    # this saved-model contract.
    assert report["symmetry_summary"]["passed"] is True
    assert int(round(float(settings["Symmetry"]))) == 0


@pytest.mark.parametrize(
    ("stability_type", "expected_name", "expected_suffix"),
    [
        (vsp.STABILITY_DEFAULT, "STABILITY_DEFAULT", ".stab"),
        pytest.param(
            vsp.STABILITY_ADJOINT,
            "STABILITY_ADJOINT",
            ".adjoint.stab",
            marks=pytest.mark.slow,
            id="adjoint",
        ),
    ],
    ids=["default", "adjoint"],
)
def test_stability_derivatives_return_finite_data_and_preserve_requested_settings(
    stability_type,
    expected_name,
    expected_suffix,
):
    report = vsp_stability_derivatives(
        MODEL_PATH,
        alpha=2.0,
        mach=0.1,
        reynolds=4.4e6,
        stability_type=stability_type,
        wake_num_nodes=WAKE_NUM_NODES,
        verbose=0,
        vspaero_verbose=0,
    )

    assert report["passed"], report["errors"]
    assert report["errors"] == []
    assert report["compute_geometry_result_id"]
    assert report["wrapper_result_id"]
    assert report["stab_result_id"]
    assert "VSPAERO_Stab" in report["result_names"]

    settings = report["vspaero_settings"]
    assert settings["effective_stability_type"] == stability_type
    assert settings["effective_stability_name"] == expected_name
    assert settings["effective_stab_file_suffix"] == expected_suffix
    assert settings["effective_wake_num_nodes"] == WAKE_NUM_NODES

    derivatives = report["derivatives"]
    assert not derivatives.empty
    assert derivatives.shape[1] > 0
    assert np.isfinite(derivatives.to_numpy(dtype=float)).all()
