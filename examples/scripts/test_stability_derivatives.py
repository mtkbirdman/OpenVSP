from pathlib import Path
import sys

import openvsp as vsp

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import (
    validate_vsp3_for_stability_derivatives,
    vsp_stability_derivatives,
)

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


if __name__ == "__main__":
    validation_report = validate_vsp3_for_stability_derivatives(
        MODEL_PATH,
        verbose=2,
    )

    print("\nValidation result")
    print("-----------------")
    print("passed:", validation_report["passed"])

    if validation_report["errors"]:
        print("\nErrors")
        for error in validation_report["errors"]:
            print("-", error["code"], error["message"])

    if validation_report["warnings"]:
        print("\nWarnings")
        for warning in validation_report["warnings"]:
            print("-", warning["code"], warning["message"])

    if not validation_report["passed"]:
        sys.exit(1)

    for stability_type, stability_name in (
        (vsp.STABILITY_DEFAULT, "STABILITY_DEFAULT"),
        (vsp.STABILITY_ADJOINT, "STABILITY_ADJOINT"),
    ):
        print(f"\nRunning {stability_name}")
        print("-" * (8 + len(stability_name)))

        stability_report = vsp_stability_derivatives(
            MODEL_PATH,
            alpha=2.0,
            mach=0.1,
            reynolds=4.4e6,
            stability_type=stability_type,
            wake_num_nodes=64,
            verbose=2,
        )

        print("\nStability derivative result")
        print("---------------------------")
        print("passed:", stability_report["passed"])

        if stability_report["errors"]:
            print("\nErrors")
            for error in stability_report["errors"]:
                print("-", error["code"], error["message"])

        if stability_report["warnings"]:
            print("\nWarnings")
            for warning in stability_report["warnings"]:
                print("-", warning["code"], warning["message"])

        if not stability_report["passed"]:
            sys.exit(1)

        effective_settings = stability_report["vspaero_settings"]
        if effective_settings.get("effective_stability_type") != stability_type:
            print(
                "Unexpected stability type:",
                effective_settings.get("effective_stability_type"),
            )
            sys.exit(1)

        if effective_settings.get("effective_wake_num_nodes") != 64:
            print(
                "Unexpected NumWakeNodes:",
                effective_settings.get("effective_wake_num_nodes"),
            )
            sys.exit(1)

        print("stab file suffix:", effective_settings.get("effective_stab_file_suffix"))
        print("\nStability derivatives")
        print(stability_report["derivatives"].T)
