from pathlib import Path
import math
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.TrimTurnSolver import (
    solve_rudder_limit_turn,
    solve_steady_gliding_turn,
    solve_steady_level_turn,
)

STAB_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.stab"
MASS_KG = 580.0
RHO_KG_M3 = 1.225
SPEED_M_S = 30.0
TURN_RATE_RAD_S = 0.06


if __name__ == "__main__":
    gliding_turn = solve_steady_gliding_turn(
        fixed={
            "V": SPEED_M_S,
            "Omega": TURN_RATE_RAD_S,
            "beta": 0.0,
        },
        stab_path=STAB_PATH,
        mass=MASS_KG,
        rho=RHO_KG_M3,
        initial_guess={
            "alpha": math.radians(2.0),
            "phi": math.radians(10.0),
            "theta": 0.0,
        },
    )

    if not gliding_turn["passed"]:
        print("gliding coordinated turn failed")
        print("max_abs_residual:", gliding_turn["max_abs_residual"])
        sys.exit(1)

    rudder_only_level_turn = solve_steady_level_turn(
        fixed={
            "V": SPEED_M_S,
            "Omega": TURN_RATE_RAD_S,
            "delta_a": 0.0,
        },
        stab_path=STAB_PATH,
        mass=MASS_KG,
        rho=RHO_KG_M3,
        initial_guess={
            "alpha": math.radians(2.0),
            "beta": math.radians(2.0),
            "phi": math.radians(10.0),
            "theta": math.radians(2.0),
        },
    )

    rudder_limit_turn = solve_rudder_limit_turn(
        stab_path=STAB_PATH,
        mass=MASS_KG,
        delta_r_max=math.radians(10.0),
        mode="gliding",
        V=SPEED_M_S,
        delta_a=0.0,
        rho=RHO_KG_M3,
        initial_guess={
            "alpha": gliding_turn["solution"]["alpha"],
            "beta": 0.0,
            "phi": gliding_turn["solution"]["phi"],
            "theta": gliding_turn["solution"]["theta"],
            "Omega": gliding_turn["solution"]["Omega"],
            "delta_e": gliding_turn["solution"]["delta_e"],
        },
    )

    cases = (
        ("gliding coordinated turn", gliding_turn),
        ("level rudder-only turn", rudder_only_level_turn),
        ("gliding rudder-limit turn", rudder_limit_turn),
    )

    for case_name, result in cases:
        print(f"\n{case_name}")
        print("-" * len(case_name))
        print("passed:", result["passed"])
        print("max_abs_residual:", result["max_abs_residual"])

        if result.get("solution"):
            solution = result["solution"]
            print("phi_deg:", math.degrees(solution["phi"]))
            print("beta_deg:", math.degrees(solution["beta"]))
            print("delta_a_deg:", math.degrees(solution["delta_a"]))
            print("delta_r_deg:", math.degrees(solution["delta_r"]))

        if not result["passed"]:
            sys.exit(1)
