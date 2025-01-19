import sys
import os

import numpy as np
import pandas as pd

# ../bin/AnalysisVSPAERO.py をモジュールとしてインポート
sys.path.append(os.path.join('..')) # 親ディレクトリをモジュール探索パスに追加
from bin.AnalysisVSPAERO import *
from bin.AnalysisVSPAERO import *

import openvsp as vsp

if __name__=='__main__':
    
    # Clear the current VSP model and read the new VSP file
    vsp.ClearVSPModel()
    vsp.Update()
    vsp.ReadVSPFile('G103A.vsp3')
    vsp.Update()

    # Define the list of alpha angles, the Mach number and Reynolds number for the sweep
    alpha = [2]
    mach = [0.1]
    reynolds = [1e6]

    # Define the list of height
    bref = 17.536552870371967
    height = list((1-np.cos(np.linspace(0,1,12)*np.pi/2))*bref)[1:] + [999]
    
    df = vsp_sweep_wig(vsp, alpha, mach, reynolds, height, AnalysisMethod=0, verbose=1)
    df = make_CDo_correction(vsp, df, Weight=580, xTr=(0.5, 0.7), CDpCL=0.0016, thickness=0.19, interference_factor=1.14, altitude=0, dT=0)
    
    df.to_csv('G103A_DegenGeom.polar', sep='\t')