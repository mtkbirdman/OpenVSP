from pathlib import Path
import sys

import openvsp as vsp

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import make_CDo_correction, vsp_trimmed_sweep

MODEL_PATH = REPO_ROOT / "examples" / "models" / "G103A" / "G103A.vsp3"


if __name__ == "__main__":
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    alpha_list = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12]
    weight_kg = 580

    trimmed_polar = vsp_trimmed_sweep(
        vsp=vsp,
        alpha_list=alpha_list,
        Weight=weight_kg,
    )
    trimmed_polar = make_CDo_correction(
        vsp,
        trimmed_polar,
        Weight=weight_kg,
        xTr=(0.5, 0.7),
        CDpCL=0.0016,
        thickness=0.19,
        interference_factor=1.14,
        altitude=0,
        dT=0,
    )
    trimmed_polar.to_csv("G103A_DegenGeom_trimmed.polar", index=False)
