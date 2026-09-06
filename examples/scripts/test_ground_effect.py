from pathlib import Path
import sys

import numpy as np
import openvsp as vsp

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.AnalysisVSPAERO import make_CDo_correction, vsp_sweep_wig

MODEL_PATH = REPO_ROOT / "examples" / "models" / "SampleGlider" / "SampleGlider.vsp3"


if __name__ == "__main__":
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile(str(MODEL_PATH))
    vsp.Update()

    alpha = [0]
    mach = [0.1]
    reynolds = [1e6]

    bref = 27
    height = list((1 - np.cos(np.linspace(0, 1, 12) * np.pi / 2)) * bref)[1:] + [999]

    polar = vsp_sweep_wig(
        vsp,
        alpha,
        mach,
        reynolds,
        height,
        AnalysisMethod=0,
        verbose=1,
    )
    polar = make_CDo_correction(
        vsp,
        polar,
        Weight=580,
        xTr=(0.5, 0.7),
        CDpCL=0.0016,
        thickness=0.19,
        interference_factor=1.14,
        altitude=0,
        dT=0,
    )
    polar.to_csv("SampleGlider_DegenGeom.polar", sep="\t")
