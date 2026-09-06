from pathlib import Path
import sys

import openvsp as vsp

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import vsp_sweep

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


if __name__ == "__main__":
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    alpha = list(range(-4, 13, 2))
    mach = [0.1]
    vsp_sweep(vsp=vsp, alpha=alpha, mach=mach)
