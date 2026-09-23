from pathlib import Path
import math
import sys

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.TrimTurnSolver import (
    solve_rudder_limit_turn,
    solve_steady_gliding_turn,
    solve_steady_level_turn,
)

STAB_PATH = REPO_ROOT / "examples" / "models" / "SampleGlider" / "SampleGlider.stab"
MASS_KG = 100.0
RHO_KG_M3 = 1.225
SPEED_M_S = 10.0
TURN_RATE_RAD_S = 0.03
RESIDUAL_TOL = 1e-6


def test_steady_gliding_turn_satisfies_force_and_moment_balance_without_thrust():
    result = solve_steady_gliding_turn(
        fixed={
            "V": SPEED_M_S,
            "Omega": TURN_RATE_RAD_S,
            "beta": 0.0,
        },
        stab_path=STAB_PATH,
        mass=MASS_KG,
        rho=RHO_KG_M3,
        initial_guess={
            "alpha": math.radians(0.0),
            "phi": math.radians(2.0),
            "theta": math.radians(-1.0),
        },
        residual_tol=RESIDUAL_TOL,
    )

    assert result["passed"], result["message"]
    assert result["max_abs_residual"] <= RESIDUAL_TOL
    assert all(abs(value) <= RESIDUAL_TOL for value in result["residuals"].values())

    solution = result["solution"]
    assert solution["V"] == pytest.approx(SPEED_M_S)
    assert solution["Omega"] == pytest.approx(TURN_RATE_RAD_S)
    assert solution["beta"] == pytest.approx(0.0)
    assert solution["T"] == pytest.approx(0.0)

    # A steady unpowered glide must lose altitude while retaining positive
    # horizontal speed.  These checks catch sign errors in the kinematics or a
    # regression that accidentally reintroduces thrust into the gliding model.
    assert result["derived"]["horizontal_speed"] > 0.0
    assert result["derived"]["sink_rate"] > 0.0
    assert result["derived"]["h_dot"] < 0.0

    finite_values = [
        solution[name]
        for name in ("alpha", "phi", "theta", "delta_a", "delta_e", "delta_r")
    ]
    assert np.isfinite(finite_values).all()


def test_steady_level_rudder_only_turn_holds_altitude_with_aileron_fixed_neutral():
    result = solve_steady_level_turn(
        fixed={
            "V": SPEED_M_S,
            "Omega": TURN_RATE_RAD_S,
            "delta_a": 0.0,
        },
        stab_path=STAB_PATH,
        mass=MASS_KG,
        rho=RHO_KG_M3,
        initial_guess={
            "alpha": math.radians(0.0),
            "beta": math.radians(2.0),
            "phi": math.radians(2.0),
            "theta": 0.0,
        },
        residual_tol=RESIDUAL_TOL,
    )

    assert result["passed"], result["message"]
    assert result["max_abs_residual"] <= RESIDUAL_TOL
    assert all(abs(value) <= RESIDUAL_TOL for value in result["residuals"].values())

    solution = result["solution"]
    assert solution["V"] == pytest.approx(SPEED_M_S)
    assert solution["Omega"] == pytest.approx(TURN_RATE_RAD_S)
    assert solution["delta_a"] == pytest.approx(0.0)

    # solve_steady_level_turn includes the height constraint.  A converged
    # result must therefore have essentially zero vertical velocity.  Positive
    # thrust is expected because this is a level, not gliding, equilibrium.
    assert result["derived"]["h_dot"] == pytest.approx(0.0, abs=RESIDUAL_TOL)
    assert result["derived"]["sink_rate"] == pytest.approx(0.0, abs=RESIDUAL_TOL)
    assert solution["T"] > 0.0

    # With the aileron fixed at neutral, this SampleGlider model must use
    # sideslip and rudder to sustain the requested non-zero turn rate.
    assert abs(solution["beta"]) > 1e-4
    assert abs(solution["delta_r"]) > 1e-4


def test_rudder_limit_turn_solves_both_endpoints_and_selects_larger_bank_angle():
    delta_r_max = math.radians(10.0)

    result = solve_rudder_limit_turn(
        stab_path=STAB_PATH,
        mass=MASS_KG,
        delta_r_max=delta_r_max,
        mode="gliding",
        V=SPEED_M_S,
        delta_a=0.0,
        rho=RHO_KG_M3,
        residual_tol=RESIDUAL_TOL,
    )

    assert result["passed"], result["message"]
    assert result["rudder_limit_complete"] is True
    assert result["negative_trim"]["passed"] is True
    assert result["positive_trim"]["passed"] is True
    assert result["selected_side"] in {"negative", "positive"}

    assert result["fixed_V"] == pytest.approx(SPEED_M_S)
    assert result["fixed_delta_a"] == pytest.approx(0.0)
    assert result["delta_r_max"] == pytest.approx(delta_r_max)
    assert abs(result["limiting_delta_r"]) == pytest.approx(delta_r_max)
    assert result["solution"]["delta_a"] == pytest.approx(0.0)
    assert result["solution"]["T"] == pytest.approx(0.0)
    assert result["max_abs_residual"] <= RESIDUAL_TOL

    negative_phi = abs(result["negative_trim"]["solution"]["phi"])
    positive_phi = abs(result["positive_trim"]["solution"]["phi"])
    assert result["max_abs_phi"] == pytest.approx(max(negative_phi, positive_phi))
    assert abs(result["solution"]["phi"]) == pytest.approx(result["max_abs_phi"])

    selected_trim = result[f'{result["selected_side"]}_trim']
    assert result["limiting_delta_r"] == pytest.approx(selected_trim["solution"]["delta_r"])
