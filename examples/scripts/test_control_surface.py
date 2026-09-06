from pathlib import Path
import sys

import numpy as np
import openvsp as vsp

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import vsp_sweep
from src.util import set_control_surface

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


if __name__ == "__main__":
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    set_control_surface(
        vsp,
        geom_name="WingGeom",
        deflection=10,
        cs_group_name="AILERON_GROUP",
        gains=(1, 1),
    )
    set_control_surface(
        vsp,
        geom_name="HTailGeom",
        deflection=0,
        cs_group_name="ELEVATOR_GROUP",
        gains=(1, -1),
    )
    set_control_surface(
        vsp,
        geom_name="VTailGeom",
        deflection=0,
        cs_group_name="RUDDER_GROUP",
    )
    vsp.Update()

    alpha = np.linspace(-4, 12, 9)
    mach = [0.1]
    vsp_sweep(vsp=vsp, alpha=alpha, mach=mach)
